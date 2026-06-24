#!/usr/bin/env python3
"""Mechanism diagnostics for coden15 ctx-marker abundance regression."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Callable

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
EPS = 1e-12


def species_genus(species: str) -> str:
    return rescue.species_genus(species)


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0}


def safe_num(value: object) -> float:
    return ev.safe_num(value)


def row_value(row: pd.Series, formula: str, exponent: float, median_cutoff: float) -> float:
    median = safe_num(row.get("Reliable_Ref_hit_median_depth"))
    mean = safe_num(row.get("Reliable_Ref_mean_depth"))
    hit_mean = safe_num(row.get("Reliable_Ref_hit_mean_depth"))
    af = max(safe_num(row.get("Reliable_Ref_zip_af")), EPS)
    ztp_af = max(safe_num(row.get("Reliable_ztp_af")), EPS)
    if formula == "pinned":
        if median >= median_cutoff:
            return median
        return mean / (af**exponent) if mean > 0.0 else 0.0
    if formula == "mean_over_zipaf":
        return mean / (af**exponent) if mean > 0.0 else 0.0
    if formula == "mean_over_ztpaf":
        return mean / (ztp_af**exponent) if mean > 0.0 else 0.0
    if formula == "hit_mean":
        return hit_mean
    if formula == "mean":
        return mean
    if formula == "median":
        return median
    raise ValueError(formula)


def values_with_rescue(
    rows: pd.DataFrame,
    formula: str,
    exponent: float,
    median_cutoff: float,
    rescue_enabled: bool,
    rescue_ratio: float = rescue.RESCUE_ACTIVE_GENUS_RATIO,
    rescue_zip_af_max: float = rescue.RESCUE_ZIP_AF_MAX,
) -> tuple[dict[str, float], list[str]]:
    rows = rows.loc[rows["gtdb_species"].astype(bool)].copy()
    rows["mechanism_value"] = [
        row_value(row, formula, exponent, median_cutoff) for _, row in rows.iterrows()
    ]
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(species_genus)

    values: dict[str, float] = {}
    active_genus_max: dict[str, float] = {}
    active = rows.loc[rows["active_gate_pass"]].copy()
    for _, row in active.iterrows():
        species = str(row["gtdb_species"])
        genus = str(row["species_genus"])
        value = safe_num(row["mechanism_value"])
        if value <= 0.0:
            continue
        values[species] = max(values.get(species, 0.0), value)
        active_genus_max[genus] = max(active_genus_max.get(genus, 0.0), value)

    if not rescue_enabled:
        return values, []

    rescue_best: dict[str, pd.Series] = {}
    for _, row in rows.loc[~rows["active_gate_pass"]].iterrows():
        species = str(row["gtdb_species"])
        genus = str(row["species_genus"])
        active_value = active_genus_max.get(genus, 0.0)
        value = safe_num(row["mechanism_value"])
        if active_value <= 0.0 or value <= 0.0:
            continue
        if safe_num(row.get("ANI_naive_calc")) < rescue.RESCUE_ANI_MIN:
            continue
        if safe_num(row.get("XnY_ctx")) < rescue.RESCUE_XNY_MIN:
            continue
        if value < rescue.RESCUE_EFFECTIVE_MIN:
            continue
        if safe_num(row.get("Reliable_Ref_zip_af")) > rescue_zip_af_max:
            continue
        if value < rescue_ratio * active_value:
            continue
        old = rescue_best.get(genus)
        if old is None or (
            value,
            safe_num(row.get("XnY_ctx")),
        ) > (
            safe_num(old.get("mechanism_value")),
            safe_num(old.get("XnY_ctx")),
        ):
            rescue_best[genus] = row

    rescued: list[str] = []
    for row in rescue_best.values():
        species = str(row["gtdb_species"])
        value = safe_num(row["mechanism_value"])
        values[species] = max(values.get(species, 0.0), value)
        rescued.append(species)
    return values, rescued


def score_pred(
    sample: int,
    method: str,
    values: dict[str, float],
    gold: dict[str, float],
    rescued: list[str],
) -> dict[str, object]:
    return rescue.score_values(sample, method, values, gold, rescued)


def species_value_table(
    sample: int,
    rows: pd.DataFrame,
    values: dict[str, float],
    gold: dict[str, float],
    method: str,
) -> pd.DataFrame:
    pred = normalize(values)
    selected_species = {s for s, v in pred.items() if v > 0.0}
    active_rows = rows.loc[rows["gtdb_species"].astype(str).isin(selected_species)].copy()
    summaries = []
    for species, sub in active_rows.groupby("gtdb_species"):
        best = sub.sort_values(
            ["active_gate_pass", "XnY_ctx", "Reliable_Ref_hit_mean_depth"],
            ascending=[False, False, False],
        ).iloc[0]
        summaries.append(
            {
                "sample": sample,
                "method": method,
                "species": species,
                "truth_abundance": gold.get(str(species), 0.0),
                "pred_abundance": pred.get(str(species), 0.0),
                "is_truth": str(species) in gold,
                "active_gate_pass": bool(best.get("active_gate_pass")),
                "XnY_ctx": safe_num(best.get("XnY_ctx")),
                "ANI_naive_calc": safe_num(best.get("ANI_naive_calc")),
                "Reliable_ztp_af": safe_num(best.get("Reliable_ztp_af")),
                "Reliable_Ref_breadth": safe_num(best.get("Reliable_Ref_breadth")),
                "Reliable_Ref_mean_depth": safe_num(best.get("Reliable_Ref_mean_depth")),
                "Reliable_Ref_hit_mean_depth": safe_num(best.get("Reliable_Ref_hit_mean_depth")),
                "Reliable_Ref_hit_median_depth": safe_num(best.get("Reliable_Ref_hit_median_depth")),
                "Reliable_Ref_zip_af": safe_num(best.get("Reliable_Ref_zip_af")),
                "Ref": str(best.get("Ref", "")),
            }
        )
    return pd.DataFrame(summaries)


def aggregate_score(sample_rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(sample_rows)
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
    return df.groupby("method", as_index=False)[mean_cols].mean().sort_values("l1_pct_points")


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)

    old_rows_by_sample = {
        sample: more.load_minco_active_rows(rescue.ctx_marker_path(sample), by_accession, by_core)
        for sample in SAMPLES
    }
    c15_rows_by_sample = {
        sample: more.load_minco_active_rows(c15.coden15_output_path(sample), by_accession, by_core)
        for sample in SAMPLES
    }
    gold_by_sample = {sample: ev.load_truth_profile(sample) for sample in SAMPLES}

    row_count_rows = []
    for sample in SAMPLES:
        for method, rows in [
            ("old_ctx_marker", old_rows_by_sample[sample]),
            ("coden15_ctxmarker", c15_rows_by_sample[sample]),
        ]:
            row_count_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "rows": len(rows),
                    "mapped_rows": int(rows["gtdb_species"].astype(bool).sum()),
                    "active_rows": int(rows["active_gate_pass"].sum()),
                    "weak_rows": int((~rows["active_gate_pass"]).sum()),
                }
            )
    pd.DataFrame(row_count_rows).to_csv(
        EXP_DIR / "coden15_mechanism_row_counts.tsv", sep="\t", index=False
    )

    grid_rows: list[dict[str, object]] = []
    formulas = ["pinned", "mean_over_zipaf", "mean_over_ztpaf", "hit_mean", "median"]
    exponents = [0.0, 0.75, 1.0, 1.05, 1.25]
    median_cutoffs = [0.0, 20.0, 1e9]
    rescue_options = [False, True]
    gate_options: list[tuple[str, Callable[[pd.DataFrame], pd.Series] | None]] = [
        ("active", None),
        ("relaxed", c15.coden15_relaxed_mask),
    ]
    for gate_name, mask_fn in gate_options:
        for formula in formulas:
            for exponent in exponents:
                if formula not in {"pinned", "mean_over_zipaf", "mean_over_ztpaf"} and exponent != 0.0:
                    continue
                for median_cutoff in median_cutoffs:
                    if formula != "pinned" and median_cutoff != 20.0:
                        continue
                    if formula == "pinned" and median_cutoff == 20.0 and exponent not in {1.0, 1.05, 1.25}:
                        continue
                    for rescue_enabled in rescue_options:
                        sample_scores = []
                        for sample in SAMPLES:
                            rows = c15_rows_by_sample[sample].copy()
                            if mask_fn is not None:
                                rows["active_gate_pass"] = mask_fn(rows)
                            values, rescued = values_with_rescue(
                                rows,
                                formula=formula,
                                exponent=exponent,
                                median_cutoff=median_cutoff,
                                rescue_enabled=rescue_enabled,
                            )
                            sample_scores.append(
                                score_pred(
                                    sample,
                                    "grid",
                                    values,
                                    gold_by_sample[sample],
                                    rescued,
                                )
                            )
                        mean = aggregate_score(sample_scores).iloc[0].to_dict()
                        grid_rows.append(
                            {
                                "gate": gate_name,
                                "formula": formula,
                                "exponent": exponent,
                                "median_cutoff": median_cutoff,
                                "rescue_enabled": rescue_enabled,
                                "TP": mean["TP"],
                                "FP": mean["FP"],
                                "FN": mean["FN"],
                                "pearson": mean["pearson"],
                                "spearman": mean["spearman"],
                                "mae_pct_points": mean["mae_pct_points"],
                                "l1_pct_points": mean["l1_pct_points"],
                            }
                        )
    grid_df = pd.DataFrame(grid_rows).sort_values("l1_pct_points")
    grid_df.to_csv(EXP_DIR / "coden15_abundance_estimator_grid.tsv", sep="\t", index=False)

    # Detailed comparison for pinned old and pinned coden15 active.
    detail_tables = []
    sample_scores = []
    for sample in SAMPLES:
        gold = gold_by_sample[sample]
        old_values, old_rescued = values_with_rescue(
            old_rows_by_sample[sample], "pinned", 1.05, 20.0, True
        )
        c15_values, c15_rescued = values_with_rescue(
            c15_rows_by_sample[sample], "pinned", 1.05, 20.0, True
        )
        sample_scores.append(score_pred(sample, "old_pinned", old_values, gold, old_rescued))
        sample_scores.append(score_pred(sample, "coden15_pinned", c15_values, gold, c15_rescued))
        detail_tables.append(
            species_value_table(sample, old_rows_by_sample[sample], old_values, gold, "old_pinned")
        )
        detail_tables.append(
            species_value_table(sample, c15_rows_by_sample[sample], c15_values, gold, "coden15_pinned")
        )

    detail_df = pd.concat(detail_tables, ignore_index=True)
    detail_df.to_csv(EXP_DIR / "coden15_old_vs_coden15_species_values.tsv", sep="\t", index=False)

    deltas = []
    for sample in SAMPLES:
        gold = gold_by_sample[sample]
        old_sub = detail_df[(detail_df["sample"] == sample) & (detail_df["method"] == "old_pinned")]
        c15_sub = detail_df[(detail_df["sample"] == sample) & (detail_df["method"] == "coden15_pinned")]
        old_pred = dict(zip(old_sub["species"], old_sub["pred_abundance"]))
        c15_pred = dict(zip(c15_sub["species"], c15_sub["pred_abundance"]))
        for species in sorted(set(gold) | set(old_pred) | set(c15_pred)):
            old_err = abs(old_pred.get(species, 0.0) - gold.get(species, 0.0))
            c15_err = abs(c15_pred.get(species, 0.0) - gold.get(species, 0.0))
            deltas.append(
                {
                    "sample": sample,
                    "species": species,
                    "truth_abundance": gold.get(species, 0.0),
                    "old_pred": old_pred.get(species, 0.0),
                    "coden15_pred": c15_pred.get(species, 0.0),
                    "old_abs_error": old_err,
                    "coden15_abs_error": c15_err,
                    "error_delta_coden15_minus_old": c15_err - old_err,
                }
            )
    delta_df = pd.DataFrame(deltas).sort_values("error_delta_coden15_minus_old", ascending=False)
    delta_df.to_csv(EXP_DIR / "coden15_error_delta_species.tsv", sep="\t", index=False)
    delta_df.head(40).to_csv(EXP_DIR / "coden15_error_delta_top40.tsv", sep="\t", index=False)

    print("row counts")
    print(pd.DataFrame(row_count_rows).to_csv(sep="\t", index=False), end="")
    print("\nbest coden15 estimator grid")
    print(grid_df.head(20).to_csv(sep="\t", index=False), end="")
    print("\ntop coden15 error increases")
    print(delta_df.head(20).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
