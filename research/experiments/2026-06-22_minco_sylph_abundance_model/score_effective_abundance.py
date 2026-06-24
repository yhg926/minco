#!/usr/bin/env python3
"""Focused abundance comparison for MinCO robust effective depth vs Sylph."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(EXP_DIR))
sys.path.insert(0, str(MORE_EXP))
sys.path.insert(0, str(OLD_EXP))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import evaluate_abundance_estimators as ev  # noqa: E402
import score_more_toymouse_gtdb as more  # noqa: E402


SAMPLES = [0, 1, 2]
MEDIAN_DEPTH_CUTOFF = 20.0
AF_EXPONENT = 1.05


def ctx_median_path(sample: int) -> Path:
    return EXP_DIR / f"toymouse_sample{sample}_ctxmarker_median_col.tsv"


def ctxobj_median_path(sample: int) -> Path:
    return EXP_DIR / f"toymouse_sample{sample}_ctxobjmarker_median_col.tsv"


def robust_effective_depth(row: pd.Series) -> float:
    med = ev.safe_num(row.get("Reliable_Ref_hit_median_depth"))
    if med >= MEDIAN_DEPTH_CUTOFF:
        return med
    mean = ev.safe_num(row.get("Reliable_Ref_mean_depth"))
    af = max(ev.safe_num(row.get("Reliable_Ref_zip_af")), ev.EPS)
    if mean <= 0.0:
        return 0.0
    return mean / (af**AF_EXPONENT)


def collapse_selected(
    rows: pd.DataFrame,
    value_fn: Callable[[pd.Series], float],
) -> dict[str, float]:
    selected = rows.loc[rows["active_gate_pass"] & rows["gtdb_species"].astype(bool)].copy()
    values: dict[str, float] = {}
    for _, row in selected.iterrows():
        species = str(row["gtdb_species"])
        value = value_fn(row)
        if math.isfinite(value) and value > 0.0:
            values[species] = max(values.get(species, 0.0), value)
    return values


def sylph_values(rows: pd.DataFrame) -> dict[str, float]:
    values: dict[str, float] = {}
    for _, row in rows.loc[rows["gtdb_species"].astype(bool)].iterrows():
        species = str(row["gtdb_species"])
        value = ev.safe_num(row.get("Taxonomic_abundance")) / 100.0
        if value > 0.0:
            values[species] = max(values.get(species, 0.0), value)
    return values


def score_values(
    sample: int,
    method: str,
    estimator: str,
    values: dict[str, float],
    gold: dict[str, float],
    renorm: bool,
) -> dict[str, object]:
    return ev.estimate_metrics(sample, method, estimator, values, gold, renorm)


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    rows_out: list[dict[str, object]] = []

    for sample in SAMPLES:
        gold = ev.load_truth_profile(sample)
        method_specs = [
            (
                "minco_ctx_marker",
                "current_normalized_depth",
                more.minco_output_path(sample, "ctx_only_current_best"),
                lambda r: ev.safe_num(r.get("Normalized_abundance_depth")),
                "minco",
            ),
            (
                "minco_ctx_marker",
                f"robust_effective_depth_T{MEDIAN_DEPTH_CUTOFF:g}_AFexp{AF_EXPONENT:g}",
                ctx_median_path(sample),
                robust_effective_depth,
                "minco",
            ),
            (
                "minco_ctxobj_marker",
                "current_normalized_depth",
                more.minco_output_path(sample, "ctxobj_markerdb_product0_active"),
                lambda r: ev.safe_num(r.get("Normalized_abundance_depth")),
                "minco",
            ),
            (
                "minco_ctxobj_marker",
                f"robust_effective_depth_T{MEDIAN_DEPTH_CUTOFF:g}_AFexp{AF_EXPONENT:g}",
                ctxobj_median_path(sample),
                robust_effective_depth,
                "minco",
            ),
            (
                "sylph_gtdb_profile",
                "reported_taxonomic_abundance",
                more.sylph_profile_path(sample),
                None,
                "sylph",
            ),
        ]
        for method, estimator, path, value_fn, kind in method_specs:
            if kind == "minco":
                rows = more.load_minco_active_rows(path, by_accession, by_core)
                values = collapse_selected(rows, value_fn)
            else:
                rows = more.load_sylph_rows(path, by_accession, by_core)
                values = sylph_values(rows)
            for renorm in [False, True]:
                rows_out.append(score_values(sample, method, estimator, values, gold, renorm))

    metrics = pd.DataFrame(rows_out)
    metrics.to_csv(EXP_DIR / "effective_abundance_sample_metrics.tsv", sep="\t", index=False)
    mean_cols = [
        "TP",
        "FP",
        "FN",
        "missing_truth_mass",
        "fp_pred_mass",
        "pred_sum_on_truth",
        "pred_sum_all",
        "pearson",
        "spearman",
        "mae_pct_points",
        "l1_pct_points",
    ]
    mean = (
        metrics.groupby(["method", "estimator", "renorm_pred"], as_index=False)[mean_cols]
        .mean()
        .sort_values(["renorm_pred", "l1_pct_points", "method"])
    )
    mean.to_csv(EXP_DIR / "effective_abundance_mean_metrics.tsv", sep="\t", index=False)
    print(mean.loc[mean["renorm_pred"]].to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
