#!/usr/bin/env python3
"""Score the context-level edge-EM prototype on external spot samples."""

from __future__ import annotations

import importlib.util
import math
import mmap
import struct
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
RUN = Path("/tmp/minco_context_em_external_20260624")
REF = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup")
STAT = REF / "minco.stat"

PROTO = ROOT / "research/experiments/2026-06-24_context_level_ambiguity_em_prototype"
SEARCH = PROTO / "search_mouse0_filtered_edge_em.py"
CAMI3_BASE = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"

BETA_VALUES = [
    0.0,
    0.0001,
    0.0003,
    0.0005,
    0.0007,
    0.001,
    0.0015,
    0.002,
    0.0025,
    0.003,
    0.0035,
    0.004,
    0.0045,
    0.005,
    0.006,
    0.008,
    0.01,
    0.02,
]

DROP_TRUTH_NAMES = {"unidentified", "unidentified plasmid", "unidentified virus"}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


search = load_module("external_edge_em_search", SEARCH)
edge_base = search.base
integrated = edge_base.integrated
cami3_base = integrated.cami3_base


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0 and math.isfinite(v))
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0 and math.isfinite(v)}


def parse_profile_truth(path: Path, sample_id: str, *, scope: str) -> dict[str, float]:
    rows: list[dict[str, object]] = []
    active = sample_id == ""
    seen_sample = False
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("@SampleID:"):
                sid = line.split(":", 1)[1].strip()
                if seen_sample and sid != sample_id:
                    break
                active = sid == sample_id
                seen_sample = seen_sample or active
                continue
            if not active or not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] == "TAXID":
                continue
            taxid, rank, taxpath, taxpathsn, pct = fields[:5]
            if rank != "species":
                continue
            try:
                pct_f = float(pct)
            except ValueError:
                continue
            if pct_f <= 0.0:
                continue
            name = taxpathsn.split("|")[-1].strip()
            if name.lower() in DROP_TRUTH_NAMES:
                continue
            if scope == "bacteria" and not (
                taxpath.startswith("2") or taxpathsn.startswith("Bacteria")
            ):
                continue
            rows.append({"taxid": str(taxid), "pct": pct_f})
    if not rows:
        raise RuntimeError(f"no truth species rows found in {path} for sample {sample_id!r}")
    df = pd.DataFrame(rows)
    grouped = df.groupby("taxid", as_index=False)["pct"].sum()
    total = float(grouped["pct"].sum())
    return {
        str(row.taxid): float(row.pct) / total
        for row in grouped.itertuples(index=False)
        if total > 0.0
    }


def load_cami3_truth() -> dict[str, float]:
    truth_rows = cami3_base.load_truth(
        Path("/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt")
    )
    return dict(
        zip(
            truth_rows["ncbi_species_taxid"].astype(str),
            truth_rows["truth_abundance_bacterial_norm"].astype(float),
        )
    )


def refname_for_gid(mm: mmap.mmap, gid: int) -> str:
    header = struct.Struct("I??2x6i")
    infile_num = header.unpack_from(mm, 0)[-1]
    if gid < 0 or gid >= infile_num:
        return ""
    begin = header.size + gid * 256
    raw = mm[begin : begin + 256]
    return raw.split(b"\0", 1)[0].decode("utf-8", errors="replace")


def gid_taxid_map(gids: np.ndarray, by_accession, by_core) -> dict[int, str]:
    out: dict[int, str] = {}
    with STAT.open("rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        for gid in sorted(int(x) for x in np.unique(gids)):
            ref = refname_for_gid(mm, gid)
            acc = integrated.truth.extract_accession(ref)
            rec, _, _ = integrated.truth.lookup_accession(acc, by_accession, by_core)
            out[gid] = (
                str(rec.get("ncbi_species_taxid") or rec.get("ncbi_taxid") or "")
                if rec
                else ""
            )
        mm.close()
    return out


def load_marker_values(
    profile_path: Path,
    dataset: str,
    sample_id: int,
    truth: dict[str, float],
    by_accession,
    by_core,
) -> tuple[dict[str, float], dict[str, float]]:
    rows = pd.read_csv(profile_path, sep="\t")
    marker = integrated.add_cami3_metadata_and_gate(
        integrated.marker_view(rows),
        by_accession,
        by_core,
    )
    prepared = integrated.abn.prev.prepare_rows(
        marker,
        dataset=dataset,
        sample_id=sample_id,
        target_col="ncbi_species_taxid",
    )
    sample = integrated.abn.prev.PreparedSample(
        dataset=dataset,
        sample_id=sample_id,
        truth=truth,
        marker=prepared,
        full=pd.DataFrame(),
    )
    marker_raw = integrated.abn.raw_values_for_sample(sample, integrated.abn.MARKER_L1)
    return marker_raw, normalize(marker_raw)


def build_edge_taxids(
    edge_path: Path,
    marker_taxids: set[str],
    by_accession,
    by_core,
) -> tuple[pd.DataFrame, dict[str, object]]:
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
        edge_path,
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
    gid_map = gid_taxid_map(edges["gid"].to_numpy(), by_accession, by_core)
    edges["target_id"] = edges["gid"].map(gid_map).fillna("")
    mapped = edges.loc[edges["target_id"].astype(bool)].copy()
    filtered = mapped.loc[mapped["target_id"].isin(marker_taxids)].copy()
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
        "edge_groups": edges.groupby(["unit_id", "qctx"]).ngroups if len(edges) else 0,
        "mapped_edge_rows": len(mapped),
        "marker_filtered_edge_rows": len(filtered),
        "marker_filtered_groups": filtered.groupby(["unit_id", "qctx"]).ngroups
        if len(filtered)
        else 0,
        "taxid_edge_rows": len(grouped),
        "taxid_groups": grouped.groupby(["unit_id", "qctx"]).ngroups if len(grouped) else 0,
        "unique_edge_gids": int(edges["gid"].nunique()) if len(edges) else 0,
        "mapped_taxids": int(mapped["target_id"].nunique()) if len(mapped) else 0,
        "marker_edge_taxids": int(grouped["target_id"].nunique()) if len(grouped) else 0,
    }
    return grouped, stats


def score_prediction(
    truth: dict[str, float],
    pred: dict[str, float],
    *,
    dataset: str,
    sample_id: int,
    method: str,
    beta: float,
) -> dict[str, object]:
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
    tp_keys = sorted(tp)
    tp_true = [truth[k] for k in tp_keys]
    tp_pred = [pred.get(k, 0.0) for k in tp_keys]
    pearson = float(pd.Series(y_pred, dtype=float).corr(pd.Series(y_true, dtype=float), method="pearson"))
    spearman = float(pd.Series(y_pred, dtype=float).corr(pd.Series(y_true, dtype=float), method="spearman"))
    tp_pearson = (
        float(pd.Series(tp_pred, dtype=float).corr(pd.Series(tp_true, dtype=float), method="pearson"))
        if len(tp_keys) >= 2
        else float("nan")
    )
    tp_mae = (
        sum(abs(a - b) for a, b in zip(tp_pred, tp_true)) / len(tp_keys)
        if tp_keys
        else float("nan")
    )
    return {
        "dataset": dataset,
        "sample_id": sample_id,
        "method": method,
        "beta": beta,
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
        "spearman": spearman,
        "tp_abundance_pearson": tp_pearson,
        "tp_abundance_mae": tp_mae,
        "missing_truth_mass": sum(truth[k] for k in fn),
        "fp_pred_mass": sum(pred.get(k, 0.0) for k in fp),
    }


def sample_specs() -> list[dict[str, object]]:
    return [
        {
            "dataset": "cami3_toygut",
            "sample_id": 0,
            "profile": RUN / "cami3_s0_profile.tsv",
            "edges": RUN / "cami3_s0_edges.tsv",
            "truth": load_cami3_truth(),
            "truth_scope": "NCBI bacterial species taxid",
        },
        {
            "dataset": "marine",
            "sample_id": 0,
            "profile": RUN / "marine_s0_profile.tsv",
            "edges": RUN / "marine_s0_edges.tsv",
            "truth": parse_profile_truth(
                Path("/tmp/gs_marine_short.profile"),
                "marmgCAMI2_short_read_sample_0",
                scope="all",
            ),
            "truth_scope": "NCBI species taxid excluding unidentified/unidentified plasmid",
        },
        {
            "dataset": "strain_madness",
            "sample_id": 0,
            "profile": RUN / "strain_s0_profile.tsv",
            "edges": RUN / "strain_s0_edges.tsv",
            "truth": parse_profile_truth(
                Path("/mnt/new3T/minco_cami2_strain_20260621/short_read/taxonomic_profile_0.txt"),
                "",
                scope="bacteria",
            ),
            "truth_scope": "NCBI bacterial species taxid",
        },
    ]


def score_sample(
    spec: dict[str, object],
    by_accession,
    by_core,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    profile = Path(spec["profile"])
    edges = Path(spec["edges"])
    if not profile.exists() or not edges.exists():
        raise FileNotFoundError(f"missing profile or edges for {spec['dataset']}: {profile}, {edges}")
    dataset = str(spec["dataset"])
    sample_id = int(spec["sample_id"])
    truth = dict(spec["truth"])

    marker_raw, marker_norm = load_marker_values(
        profile,
        dataset,
        sample_id,
        truth,
        by_accession,
        by_core,
    )
    edge_taxids, stats = build_edge_taxids(edges, set(marker_raw), by_accession, by_core)
    edge_taxids = search.enrich_edges(edge_taxids)
    cfg = next(c for c in search.make_configs() if c.name == "selected_group_species2")
    filtered = search.apply_filter(edge_taxids, cfg)
    edge_norm = search.em_counts_filtered(
        filtered,
        marker_norm,
        alpha=0.0,
        weight_mode="event",
    )
    filtered_groups = filtered.groupby(["unit_id", "qctx"]).ngroups if len(filtered) else 0
    filtered_taxids = int(filtered["target_id"].nunique()) if len(filtered) else 0

    rows: list[dict[str, object]] = []
    for beta in BETA_VALUES:
        pred = edge_base.blend(marker_norm, edge_norm, beta)
        method = (
            "marker_l1_per_read_profile"
            if beta == 0.0
            else f"edge_em_selected_group_species2_beta{beta:g}"
        )
        row = score_prediction(
            truth,
            pred,
            dataset=dataset,
            sample_id=sample_id,
            method=method,
            beta=beta,
        )
        row.update(
            {
                "truth_scope": spec["truth_scope"],
                "marker_raw_targets": len(marker_raw),
                "edge_norm_targets": len(edge_norm),
                "filtered_rows": len(filtered) if beta > 0.0 else 0,
                "filtered_groups": filtered_groups if beta > 0.0 else 0,
                "filtered_taxids": filtered_taxids if beta > 0.0 else 0,
            }
        )
        rows.append(row)

    stats.update(
        {
            "dataset": dataset,
            "sample_id": sample_id,
            "truth_taxa": len(truth),
            "marker_raw_targets": len(marker_raw),
            "marker_norm_targets": len(marker_norm),
            "edge_norm_targets": len(edge_norm),
            "selected_group_species2_rows": len(filtered),
            "selected_group_species2_groups": filtered_groups,
            "selected_group_species2_taxids": filtered_taxids,
        }
    )
    return rows, stats


def summarize(sample_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, group in sample_metrics.groupby("dataset", sort=False):
        baseline = group.loc[group["beta"] == 0.0].iloc[0].to_dict()
        baseline["summary_role"] = "marker baseline from same edge run"
        rows.append(baseline)

        positive = group.loc[group["beta"] > 0.0].copy()
        if not positive.empty:
            best = positive.sort_values(["l1_pct_points", "F1"], ascending=[True, False]).iloc[0].to_dict()
            best["summary_role"] = "best positive beta diagnostic"
            rows.append(best)

            fixed = positive.loc[np.isclose(positive["beta"].astype(float), 0.0001)]
            if not fixed.empty:
                fixed_row = fixed.iloc[0].to_dict()
                fixed_row["summary_role"] = "fixed smallest positive beta"
                rows.append(fixed_row)
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core = integrated.truth.load_gtdb_metadata(integrated.truth.GTDB_METADATA)
    metric_rows: list[dict[str, object]] = []
    stat_rows: list[dict[str, object]] = []
    for spec in sample_specs():
        rows, stats = score_sample(spec, by_accession, by_core)
        metric_rows.extend(rows)
        stat_rows.append(stats)
        print(f"scored {spec['dataset']}", file=sys.stderr)

    metrics = pd.DataFrame(metric_rows)
    diagnostics = pd.DataFrame(stat_rows)
    summary = summarize(metrics)
    metrics.to_csv(RESULTS / "external_edge_em_sample_metrics.tsv", sep="\t", index=False)
    diagnostics.to_csv(RESULTS / "external_edge_em_diagnostics.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "external_edge_em_summary.tsv", sep="\t", index=False)
    print(summary.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
