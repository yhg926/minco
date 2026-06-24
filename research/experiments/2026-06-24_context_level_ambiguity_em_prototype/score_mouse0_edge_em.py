#!/usr/bin/env python3
"""Prototype context-level ambiguity EM on one Toy Mouse Gut sample."""

from __future__ import annotations

import importlib.util
import math
import mmap
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
RUN = Path("/tmp/minco_context_em_20260624")
REF = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup")
PROFILE = RUN / "mouse_s0_profile.tsv"
EDGES = RUN / "mouse_s0_edges.tsv"
STAT = REF / "minco.stat"
INTEGRATED_SCORER = ROOT / "research/experiments/2026-06-24_integrated_dual_evidence_benchmark/score_integrated_dual_evidence.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


integrated = load_module("context_em_integrated_scorer", INTEGRATED_SCORER)


def refname_for_gid(mm: mmap.mmap, gid: int) -> str:
    header = struct.Struct("I??2x6i")
    infile_num = header.unpack_from(mm, 0)[-1]
    if gid < 0 or gid >= infile_num:
        return ""
    begin = header.size + gid * 256
    raw = mm[begin : begin + 256]
    return raw.split(b"\0", 1)[0].decode("utf-8", errors="replace")


def gid_species_map(gids: np.ndarray, by_accession, by_core) -> dict[int, str]:
    out: dict[int, str] = {}
    with open(STAT, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        for gid in sorted(int(x) for x in np.unique(gids)):
            ref = refname_for_gid(mm, gid)
            acc = integrated.truth.extract_accession(ref)
            rec, _, _ = integrated.truth.lookup_accession(acc, by_accession, by_core)
            out[gid] = rec["gtdb_species"] if rec else ""
        mm.close()
    return out


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0 and math.isfinite(v))
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0 and math.isfinite(v)}


def score_prediction(truth: dict[str, float], pred: dict[str, float], method: str) -> dict[str, object]:
    truth_set = set(truth)
    pred_set = set(pred)
    tp = truth_set & pred_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set
    precision = len(tp) / len(pred_set) if pred_set else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    keys = sorted(truth)
    y_true = [truth[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    pearson = float(pd.Series(y_pred, dtype=float).corr(pd.Series(y_true, dtype=float), method="pearson"))
    return {
        "method": method,
        "truth_taxa": len(truth_set),
        "pred_taxa": len(pred_set),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) * 100.0,
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) / max(len(keys), 1) * 100.0,
        "pearson": pearson,
        "missing_truth_mass": sum(truth[k] for k in fn),
        "fp_pred_mass": sum(pred.get(k, 0.0) for k in fp),
    }


def load_marker_values():
    by_accession, by_core = integrated.truth.load_gtdb_metadata(integrated.truth.GTDB_METADATA)
    rows = pd.read_csv(PROFILE, sep="\t")
    marker = integrated.add_mouse_metadata_and_gate(integrated.marker_view(rows), by_accession, by_core)
    full = integrated.add_mouse_metadata_and_gate(integrated.full_view(rows), by_accession, by_core)
    marker_prepared = integrated.abn.prev.prepare_rows(
        marker,
        dataset="mouse_gtdb",
        sample_id=0,
        target_col="gtdb_species",
    )
    full_prepared = integrated.abn.prev.prepare_rows(
        full,
        dataset="mouse_gtdb",
        sample_id=0,
        target_col="gtdb_species",
    )
    truth = integrated.mouse_ev.load_truth_profile(0)
    sample = integrated.abn.prev.PreparedSample(
        dataset="mouse_gtdb",
        sample_id=0,
        truth=truth,
        marker=marker_prepared,
        full=full_prepared,
    )
    marker_raw = integrated.abn.raw_values_for_sample(sample, integrated.abn.MARKER_L1)
    marker_norm = normalize(marker_raw)
    return marker_raw, marker_norm, truth, by_accession, by_core


def build_edge_species(marker_species: set[str], by_accession, by_core) -> tuple[pd.DataFrame, dict[str, object]]:
    usecols = [
        "unit_id",
        "qctx",
        "gid",
        "diff",
        "best_diff",
        "candidate_refs",
        "selected_refs",
        "selected_by_mode",
    ]
    edges = pd.read_csv(
        EDGES,
        sep="\t",
        usecols=usecols,
        dtype={
            "unit_id": "uint64",
            "qctx": "uint64",
            "gid": "uint32",
            "diff": "uint16",
            "best_diff": "uint16",
            "candidate_refs": "uint16",
            "selected_refs": "uint16",
            "selected_by_mode": "uint8",
        },
    )
    gid_map = gid_species_map(edges["gid"].to_numpy(), by_accession, by_core)
    edges["target_id"] = edges["gid"].map(gid_map).fillna("")
    mapped = edges.loc[edges["target_id"].astype(bool)].copy()
    filtered = mapped.loc[mapped["target_id"].isin(marker_species)].copy()
    grouped = (
        filtered.groupby(["unit_id", "qctx", "target_id"], as_index=False)
        .agg(
            diff=("diff", "min"),
            best_diff=("best_diff", "min"),
            selected_by_mode=("selected_by_mode", "max"),
            candidate_refs=("candidate_refs", "max"),
            selected_refs=("selected_refs", "max"),
            edge_refs=("gid", "nunique"),
        )
    )
    stats = {
        "edge_rows": len(edges),
        "edge_groups": edges.groupby(["unit_id", "qctx"]).ngroups,
        "mapped_edge_rows": len(mapped),
        "marker_filtered_edge_rows": len(filtered),
        "marker_filtered_groups": filtered.groupby(["unit_id", "qctx"]).ngroups if len(filtered) else 0,
        "species_edge_rows": len(grouped),
        "species_groups": grouped.groupby(["unit_id", "qctx"]).ngroups if len(grouped) else 0,
        "unique_edge_gids": int(edges["gid"].nunique()),
        "mapped_species": int(mapped["target_id"].nunique()),
        "marker_edge_species": int(grouped["target_id"].nunique()) if len(grouped) else 0,
    }
    return grouped, stats


def em_counts(
    edge_species: pd.DataFrame,
    marker_norm: dict[str, float],
    *,
    alpha: float,
    prior_power: float,
    selected_only: bool,
    iterations: int = 30,
) -> dict[str, float]:
    rows = edge_species
    if selected_only:
        rows = rows.loc[rows["selected_by_mode"] > 0].copy()
    if rows.empty:
        return {}
    marker_species = sorted(marker_norm)
    target_to_code = {target: i for i, target in enumerate(marker_species)}
    rows = rows.loc[rows["target_id"].isin(target_to_code)].copy()
    if rows.empty:
        return {}
    target_codes = rows["target_id"].map(target_to_code).to_numpy(dtype=np.int64)
    group_codes, _ = pd.factorize(pd.MultiIndex.from_frame(rows[["unit_id", "qctx"]]), sort=False)
    n_groups = int(group_codes.max()) + 1
    diff_delta = rows["diff"].to_numpy(dtype=float) - rows["best_diff"].to_numpy(dtype=float)
    diff_delta = np.maximum(diff_delta, 0.0)
    likelihood = np.exp(-alpha * diff_delta)
    prior = np.array([max(marker_norm.get(target, 0.0), 1e-12) for target in marker_species], dtype=float)
    prior /= prior.sum()
    abundance = prior.copy()
    for _ in range(iterations):
        weights = (abundance[target_codes] ** prior_power) * likelihood
        denom = np.bincount(group_codes, weights=weights, minlength=n_groups)
        valid = denom[group_codes] > 0.0
        post = np.zeros_like(weights)
        post[valid] = weights[valid] / denom[group_codes][valid]
        new_abundance = np.bincount(target_codes, weights=post, minlength=len(marker_species))
        total = new_abundance.sum()
        if total <= 0.0:
            break
        new_abundance /= total
        if np.max(np.abs(new_abundance - abundance)) < 1e-10:
            abundance = new_abundance
            break
        abundance = new_abundance
    return {target: float(abundance[i]) for i, target in enumerate(marker_species) if abundance[i] > 0.0}


def blend(marker_norm: dict[str, float], edge_norm: dict[str, float], beta: float) -> dict[str, float]:
    keys = set(marker_norm) | set(edge_norm)
    vals = {
        k: (1.0 - beta) * marker_norm.get(k, 0.0) + beta * edge_norm.get(k, 0.0)
        for k in keys
    }
    # Keep the marker callset fixed for this abundance-only prototype.
    for k in marker_norm:
        vals[k] = vals.get(k, 0.0) + 1e-15
    return normalize(vals)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    marker_raw, marker_norm, truth, by_accession, by_core = load_marker_values()
    edge_species_path = RESULTS / "mouse0_marker_candidate_edge_species.tsv"
    edge_stats_path = RESULTS / "mouse0_edge_stats.tsv"
    if edge_species_path.exists() and edge_stats_path.exists():
        edge_species = pd.read_csv(edge_species_path, sep="\t")
    else:
        edge_species, edge_stats = build_edge_species(set(marker_raw), by_accession, by_core)
        pd.DataFrame([edge_stats]).to_csv(edge_stats_path, sep="\t", index=False)
        edge_species.to_csv(edge_species_path, sep="\t", index=False)

    rows = [score_prediction(truth, marker_norm, "marker_l1_per_read_profile")]
    beta_values = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.1]
    for selected_only in [True, False]:
        for alpha in [0.0, 0.5, 1.0, 2.0, 4.0]:
            edge_norm = em_counts(
                edge_species,
                marker_norm,
                alpha=alpha,
                prior_power=1.0,
                selected_only=selected_only,
                iterations=20,
            )
            mode = "selected" if selected_only else "allcand"
            rows.append(score_prediction(
                truth,
                edge_norm,
                f"edge_only_{mode}_alpha{alpha:g}",
            ))
            for beta in beta_values:
                pred = blend(marker_norm, edge_norm, beta)
                rows.append(score_prediction(
                    truth,
                    pred,
                    f"blend_{mode}_alpha{alpha:g}_beta{beta:g}",
                ))
    out = pd.DataFrame(rows).sort_values(["l1_pct_points", "F1"], ascending=[True, False])
    out.to_csv(RESULTS / "mouse0_edge_em_focused_grid.tsv", sep="\t", index=False)
    out.head(30).to_csv(EXP / "summary.tsv", sep="\t", index=False)
    print(out.head(30).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
