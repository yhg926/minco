#!/usr/bin/env python3
"""Score final baseline-plus full-sketch fallback strategies on GTDB toy mouse."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXP_DIR))

import score_hybrid_mouse as hybrid  # noqa: E402


STRATEGIES = [
    ("hybrid_baseline_plus_full_xny500_scale04", 500, 0.4),
    ("hybrid_baseline_plus_full_xny700_scale04", 700, 0.4),
    ("hybrid_baseline_plus_full_xny1000_scale04", 1000, 0.4),
]


def values_scaled(rows: pd.DataFrame, scale: float) -> tuple[dict[str, float], list[str]]:
    rows = rows.loc[rows["gtdb_species"].astype(bool)].copy()
    rows["abund_depth"] = rows["effective_depth"]
    rows.loc[rows["source_mode"] == "full_fallback", "abund_depth"] *= scale

    active = rows.loc[rows["active_gate_pass"]].copy()
    marker_species = set(active.loc[active["source_mode"] == "ctx_marker", "gtdb_species"].astype(str))
    if marker_species:
        active = active.loc[
            ~(
                (active["source_mode"] == "full_fallback")
                & active["gtdb_species"].astype(str).isin(marker_species)
            )
        ].copy()

    values = active.loc[active["abund_depth"] > 0.0].groupby("gtdb_species")["abund_depth"].max().to_dict()
    active_genus_max = active.loc[active["abund_depth"] > 0.0].groupby("species_genus")["abund_depth"].max()

    candidates = rows.loc[
        ~rows["active_gate_pass"]
        & (rows["ANI_naive_calc"] >= hybrid.rescue.RESCUE_ANI_MIN)
        & (rows["XnY_ctx"] >= hybrid.rescue.RESCUE_XNY_MIN)
        & (rows["abund_depth"] >= hybrid.rescue.RESCUE_EFFECTIVE_MIN)
        & (rows["Reliable_Ref_zip_af"] <= hybrid.rescue.RESCUE_ZIP_AF_MAX)
    ].copy()

    rescued: list[str] = []
    if not candidates.empty and not active_genus_max.empty:
        candidates["active_genus_value"] = candidates["species_genus"].map(active_genus_max).fillna(0.0)
        candidates = candidates.loc[
            (candidates["active_genus_value"] > 0.0)
            & (candidates["abund_depth"] >= hybrid.rescue.RESCUE_ACTIVE_GENUS_RATIO * candidates["active_genus_value"])
        ].copy()
        if not candidates.empty:
            winners = (
                candidates.sort_values(["species_genus", "abund_depth", "XnY_ctx"])
                .groupby("species_genus", as_index=False)
                .tail(1)
            )
            for _, row in winners.iterrows():
                species = str(row["gtdb_species"])
                if row.get("source_mode", "") == "full_fallback" and species in values:
                    continue
                value = float(row["abund_depth"])
                values[species] = max(float(values.get(species, 0.0)), value)
                rescued.append(species)
    return values, rescued


def score_values(sample: int, method: str, values: dict[str, float], rescued: list[str]) -> dict[str, object]:
    return hybrid.rescue.score_values(sample, method, values, hybrid.ev.load_truth_profile(sample), rescued)


def main() -> int:
    ctx_marker_sizes = hybrid.load_psmp_sizes(hybrid.CTX_MARKER_PSM_PATH)
    by_accession, by_core = hybrid.truth.load_gtdb_metadata(hybrid.truth.GTDB_METADATA)
    views = {
        sample: hybrid.load_sample_views(sample, by_accession, by_core, ctx_marker_sizes)
        for sample in hybrid.SAMPLES
    }

    sample_rows = []
    for sample, (marker, full) in views.items():
        old_values, old_rescued = values_scaled(marker, 1.0)
        sample_rows.append(score_values(sample, "old_ctxmarker_robust_rescue_recalc", old_values, old_rescued))

        for method, fallback_xny_min, scale in STRATEGIES:
            rows = pd.concat(
                [
                    marker.copy(),
                    full.loc[
                        (full["ctx_marker_size"] < 200)
                        & (full["XnY_ctx"] >= fallback_xny_min)
                    ].copy(),
                ],
                ignore_index=True,
                sort=False,
            ).fillna(0)
            values, rescued = values_scaled(rows, scale)
            sample_rows.append(score_values(sample, method, values, rescued))

    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(EXP_DIR / "hybrid_mouse_baseline_plus_sample_metrics.tsv", sep="\t", index=False)

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
    mean_df["precision"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FP"])
    mean_df["recall"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FN"])
    mean_df["F1"] = 2.0 * mean_df["precision"] * mean_df["recall"] / (mean_df["precision"] + mean_df["recall"])
    mean_df = mean_df.sort_values(["F1", "l1_pct_points"], ascending=[False, True])
    mean_df.to_csv(EXP_DIR / "hybrid_mouse_baseline_plus_mean_metrics.tsv", sep="\t", index=False)
    print(mean_df.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
