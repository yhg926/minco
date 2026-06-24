#!/usr/bin/env python3
"""Diagnostic coden15 gate sweep for abundance after mechanism analysis."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
OLD_ABUND_EXP = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"
MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
OLD_TRUTH_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"

sys.path.insert(0, str(OLD_ABUND_EXP))
sys.path.insert(0, str(MORE_EXP))
sys.path.insert(0, str(OLD_TRUTH_EXP))
sys.path.insert(0, str(EXP_DIR))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import evaluate_abundance_estimators as ev  # noqa: E402
import score_coden15_toymouse_three_samples as c15  # noqa: E402
import score_intragenus_winner_rescue as rescue  # noqa: E402
import score_more_toymouse_gtdb as more  # noqa: E402
from analyze_coden15_mechanism import values_with_rescue  # noqa: E402


SAMPLES = [0, 1, 2]


def gate_mask(
    rows: pd.DataFrame,
    xny_min: float,
    af_floor: float,
    vmr_min: float,
    delta_max: float,
    ani_min: float = 0.95,
    hit_mean_min: float = 3.0,
) -> pd.Series:
    delta_trigger = (
        (rows["Reliable_Ref_hit_mean_depth"] > hit_mean_min)
        & (rows["Reliable_depth_vmr"] > vmr_min)
    )
    delta_pass = ~delta_trigger | (rows["ANI_AF_delta"] < delta_max)
    return (
        (rows["XnY_ctx"] >= xny_min)
        & (rows["ANI_naive_calc"] > ani_min)
        & (rows["Reliable_ztp_af"] >= af_floor)
        & delta_pass
        & rows["gtdb_species"].astype(bool)
    )


def score_one(sample: int, rows: pd.DataFrame, gold: dict[str, float]) -> dict[str, object]:
    values, rescued = values_with_rescue(
        rows,
        formula="pinned",
        exponent=1.0,
        median_cutoff=20.0,
        rescue_enabled=True,
    )
    return rescue.score_values(sample, "grid", values, gold, rescued)


def mean_score(sample_rows: list[dict[str, object]]) -> dict[str, float]:
    df = pd.DataFrame(sample_rows)
    cols = [
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
    return df[cols].mean().to_dict()


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    rows_by_sample = {
        sample: more.load_minco_active_rows(c15.coden15_output_path(sample), by_accession, by_core)
        for sample in SAMPLES
    }
    gold_by_sample = {sample: ev.load_truth_profile(sample) for sample in SAMPLES}

    grid_rows = []
    best_sample_rows: list[dict[str, object]] = []
    best_key = None
    best_l1 = math.inf

    for xny_min in [15.0, 20.0]:
        for af_floor in [0.35, 0.37, 0.40]:
            for vmr_min in [10.0, 15.0, 20.0, 50.0]:
                for delta_max in [0.03]:
                    sample_scores = []
                    for sample in SAMPLES:
                        rows = rows_by_sample[sample].copy()
                        rows["active_gate_pass"] = gate_mask(
                            rows,
                            xny_min=xny_min,
                            af_floor=af_floor,
                            vmr_min=vmr_min,
                            delta_max=delta_max,
                        )
                        sample_scores.append(score_one(sample, rows, gold_by_sample[sample]))
                    mean = mean_score(sample_scores)
                    rec = {
                        "xny_min": xny_min,
                        "af_floor": af_floor,
                        "vmr_min": vmr_min,
                        "delta_max": delta_max,
                        **mean,
                    }
                    grid_rows.append(rec)
                    if mean["l1_pct_points"] < best_l1:
                        best_l1 = mean["l1_pct_points"]
                        best_key = (xny_min, af_floor, vmr_min, delta_max)
                        best_sample_rows = [
                            {
                                **row,
                                "xny_min": xny_min,
                                "af_floor": af_floor,
                                "vmr_min": vmr_min,
                                "delta_max": delta_max,
                            }
                            for row in sample_scores
                        ]

    grid_df = pd.DataFrame(grid_rows).sort_values("l1_pct_points")
    grid_df.to_csv(EXP_DIR / "coden15_gate_grid.tsv", sep="\t", index=False)
    pd.DataFrame(best_sample_rows).to_csv(
        EXP_DIR / "coden15_gate_grid_best_sample_metrics.tsv", sep="\t", index=False
    )
    print("best", best_key, "l1", best_l1)
    print(grid_df.head(30).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
