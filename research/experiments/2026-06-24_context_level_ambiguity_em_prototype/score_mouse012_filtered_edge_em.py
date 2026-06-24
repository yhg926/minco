#!/usr/bin/env python3
"""Validate the filtered edge-EM abundance blend on Toy Mouse samples 0-2."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
RUN = Path("/tmp/minco_context_em_20260624")
SEARCH = EXP / "search_mouse0_filtered_edge_em.py"

PROFILE = {
    0: RUN / "mouse_s0_profile.tsv",
    1: RUN / "mouse_s1_profile.tsv",
    2: RUN / "mouse_s2_profile.tsv",
}
EDGES = {
    0: RUN / "mouse_s0_edges.tsv",
    1: RUN / "mouse_s1_edges.tsv",
    2: RUN / "mouse_s2_edges.tsv",
}

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


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


search = load_module("mouse012_filtered_edge_search", SEARCH)
base = search.base
integrated = base.integrated


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0 and math.isfinite(v))
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0 and math.isfinite(v)}


def sample_cache(sample_id: int) -> tuple[Path, Path]:
    return (
        RESULTS / f"mouse{sample_id}_marker_candidate_edge_species.tsv",
        RESULTS / f"mouse{sample_id}_edge_stats.tsv",
    )


def load_marker_values(sample_id: int, by_accession, by_core):
    rows = pd.read_csv(PROFILE[sample_id], sep="\t")
    marker = integrated.add_mouse_metadata_and_gate(integrated.marker_view(rows), by_accession, by_core)
    full = integrated.add_mouse_metadata_and_gate(integrated.full_view(rows), by_accession, by_core)
    marker_prepared = integrated.abn.prev.prepare_rows(
        marker,
        dataset="mouse_gtdb",
        sample_id=sample_id,
        target_col="gtdb_species",
    )
    full_prepared = integrated.abn.prev.prepare_rows(
        full,
        dataset="mouse_gtdb",
        sample_id=sample_id,
        target_col="gtdb_species",
    )
    truth = integrated.mouse_ev.load_truth_profile(sample_id)
    sample = integrated.abn.prev.PreparedSample(
        dataset="mouse_gtdb",
        sample_id=sample_id,
        truth=truth,
        marker=marker_prepared,
        full=full_prepared,
    )
    marker_raw = integrated.abn.raw_values_for_sample(sample, integrated.abn.MARKER_L1)
    marker_norm = normalize(marker_raw)
    return marker_raw, marker_norm, truth


def build_edge_species(sample_id: int, marker_species: set[str], by_accession, by_core) -> tuple[pd.DataFrame, dict[str, object]]:
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
        EDGES[sample_id],
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
    gid_map = base.gid_species_map(edges["gid"].to_numpy(), by_accession, by_core)
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
        "sample_id": sample_id,
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


def load_or_build_edge_species(sample_id: int, marker_species: set[str], by_accession, by_core) -> tuple[pd.DataFrame, dict[str, object]]:
    edge_path, stats_path = sample_cache(sample_id)
    if edge_path.exists() and stats_path.exists():
        edge = pd.read_csv(edge_path, sep="\t")
        stats = pd.read_csv(stats_path, sep="\t").iloc[0].to_dict()
        return edge, stats
    edge, stats = build_edge_species(sample_id, marker_species, by_accession, by_core)
    edge.to_csv(edge_path, sep="\t", index=False)
    pd.DataFrame([stats]).to_csv(stats_path, sep="\t", index=False)
    return edge, stats


def score_sample(sample_id: int, by_accession, by_core) -> tuple[list[dict[str, object]], dict[str, object]]:
    marker_raw, marker_norm, truth = load_marker_values(sample_id, by_accession, by_core)
    edge, stats = load_or_build_edge_species(sample_id, set(marker_raw), by_accession, by_core)
    edge = search.enrich_edges(edge)
    cfg = next(c for c in search.make_configs() if c.name == "selected_group_species2")
    filtered = search.apply_filter(edge, cfg)
    filtered_groups = filtered.groupby(["unit_id", "qctx"]).ngroups if len(filtered) else 0
    filtered_species = int(filtered["target_id"].nunique()) if len(filtered) else 0
    edge_norm = search.em_counts_filtered(filtered, marker_norm, alpha=0.0, weight_mode="event")

    rows: list[dict[str, object]] = []
    for beta in BETA_VALUES:
        pred = base.blend(marker_norm, edge_norm, beta)
        method = "marker_l1_per_read_profile" if beta == 0.0 else f"filtered_group_species2_beta{beta:g}"
        score = base.score_prediction(truth, pred, method)
        score.update(
            {
                "sample_id": sample_id,
                "filter_id": cfg.name if beta > 0.0 else "baseline",
                "alpha": 0.0 if beta > 0.0 else math.nan,
                "weight_mode": "event" if beta > 0.0 else "none",
                "beta": beta,
                "filtered_rows": len(filtered) if beta > 0.0 else 0,
                "filtered_groups": filtered_groups if beta > 0.0 else 0,
                "filtered_species": filtered_species if beta > 0.0 else 0,
            }
        )
        rows.append(score)

    diag = dict(stats)
    diag.update(
        {
            "sample_id": sample_id,
            "marker_species": len(marker_norm),
            "truth_species": len(truth),
            "selected_group_species2_rows": len(filtered),
            "selected_group_species2_groups": filtered_groups,
            "selected_group_species2_species": filtered_species,
        }
    )
    return rows, diag


def summarize(sample_rows: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [
        "truth_taxa",
        "pred_taxa",
        "TP",
        "FP",
        "FN",
        "precision",
        "recall",
        "F1",
        "l1_pct_points",
        "mae_pct_points",
        "pearson",
        "missing_truth_mass",
        "fp_pred_mass",
        "beta",
    ]
    rows = []
    for method, group in sample_rows.groupby("method", sort=False):
        rec = {"method": method, "samples": len(group)}
        for col in metric_cols:
            rec[f"mean_{col}"] = float(group[col].mean())
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["mean_l1_pct_points", "mean_F1"], ascending=[True, False])


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core = integrated.truth.load_gtdb_metadata(integrated.truth.GTDB_METADATA)
    sample_rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []
    for sample_id in [0, 1, 2]:
        rows, diag = score_sample(sample_id, by_accession, by_core)
        sample_rows.extend(rows)
        diagnostics.append(diag)
        print(f"scored sample {sample_id}", file=sys.stderr)

    sample_df = pd.DataFrame(sample_rows)
    summary = summarize(sample_df)
    diag_df = pd.DataFrame(diagnostics)
    sample_df.to_csv(RESULTS / "mouse012_filtered_edge_em_sample_metrics.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "mouse012_filtered_edge_em_summary.tsv", sep="\t", index=False)
    diag_df.to_csv(RESULTS / "mouse012_filtered_edge_em_diagnostics.tsv", sep="\t", index=False)
    summary.head(20).to_csv(EXP / "mouse012_filtered_summary.tsv", sep="\t", index=False)
    print(summary.head(20).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
