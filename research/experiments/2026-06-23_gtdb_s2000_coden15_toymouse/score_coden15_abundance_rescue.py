#!/usr/bin/env python3
"""Score coden15 Toy Mouse abundance with the pinned robust-depth rescue rule."""

from __future__ import annotations

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


SAMPLES = [0, 1, 2]


def score_minco_path(
    sample: int,
    method: str,
    path: Path,
    by_accession,
    by_core,
    relaxed: bool = False,
) -> dict[str, object]:
    gold = ev.load_truth_profile(sample)
    rows = more.load_minco_active_rows(path, by_accession, by_core)
    if relaxed:
        rows = rows.copy()
        rows["active_gate_pass"] = c15.coden15_relaxed_mask(rows)
    values, rescued = rescue.minco_values_with_rescue(rows)
    return rescue.score_values(sample, method, values, gold, rescued)


def score_sylph(sample: int, by_accession, by_core) -> dict[str, object]:
    gold = ev.load_truth_profile(sample)
    rows = more.load_sylph_rows(more.sylph_profile_path(sample), by_accession, by_core)
    return rescue.score_values(
        sample,
        "sylph_gtdb_profile_reported_taxonomic_abundance",
        rescue.sylph_values(rows),
        gold,
        [],
    )


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    sample_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        sample_rows.append(
            score_minco_path(
                sample,
                "old_ctx_marker_robust_intragenus_winner_rescue",
                rescue.ctx_marker_path(sample),
                by_accession,
                by_core,
            )
        )
        sample_rows.append(
            score_minco_path(
                sample,
                "coden15_ctxmarker_active_robust_intragenus_winner_rescue",
                c15.coden15_output_path(sample),
                by_accession,
                by_core,
            )
        )
        sample_rows.append(
            score_minco_path(
                sample,
                "coden15_ctxmarker_relaxed_robust_intragenus_winner_rescue",
                c15.coden15_output_path(sample),
                by_accession,
                by_core,
                relaxed=True,
            )
        )
        sample_rows.append(score_sylph(sample, by_accession, by_core))

    sample_df = pd.DataFrame(sample_rows)
    sample_path = EXP_DIR / "coden15_abundance_rescue_sample_metrics.tsv"
    sample_df.to_csv(sample_path, sep="\t", index=False)

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
    mean_df = sample_df.groupby("method", as_index=False)[mean_cols].mean()
    mean_df = mean_df.sort_values("l1_pct_points")
    mean_path = EXP_DIR / "coden15_abundance_rescue_mean_metrics.tsv"
    mean_df.to_csv(mean_path, sep="\t", index=False)
    print(mean_df.to_csv(sep="\t", index=False), end="")
    print(f"wrote {sample_path}")
    print(f"wrote {mean_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
