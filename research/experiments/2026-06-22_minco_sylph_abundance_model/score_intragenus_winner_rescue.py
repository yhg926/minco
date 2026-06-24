#!/usr/bin/env python3
"""Score MinCO robust abundance with an intra-genus winner rescue rule."""

from __future__ import annotations

import math
import sys
from pathlib import Path

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
RESCUE_ANI_MIN = 0.999
RESCUE_XNY_MIN = 100.0
RESCUE_EFFECTIVE_MIN = 5.0
RESCUE_ZIP_AF_MAX = 0.25
RESCUE_ACTIVE_GENUS_RATIO = 3.0


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def ctx_marker_path(sample: int) -> Path:
    return EXP_DIR / f"toymouse_sample{sample}_ctxmarker_median_col.tsv"


def robust_effective_depth(row: pd.Series) -> float:
    median = ev.safe_num(row.get("Reliable_Ref_hit_median_depth"))
    if median >= MEDIAN_DEPTH_CUTOFF:
        return median
    mean = ev.safe_num(row.get("Reliable_Ref_mean_depth"))
    af = max(ev.safe_num(row.get("Reliable_Ref_zip_af")), ev.EPS)
    return mean / (af**AF_EXPONENT) if mean > 0.0 else 0.0


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0}


def score_values(
    sample: int,
    method: str,
    values: dict[str, float],
    gold: dict[str, float],
    rescued: list[str],
) -> dict[str, object]:
    pred = normalize(values)
    keys = sorted(gold)
    y_true = [gold[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    pearson, spearman = ev.pearson_spearman(y_true, y_pred)
    pred_species = set(pred)
    gold_species = set(gold)
    tp = pred_species & gold_species
    fp = pred_species - gold_species
    fn = gold_species - pred_species
    return {
        "sample": sample,
        "method": method,
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "missing_truth_mass": sum(gold[s] for s in fn),
        "fp_pred_mass": sum(pred.get(s, 0.0) for s in fp),
        "pred_sum_on_truth": sum(y_pred),
        "pred_sum_all": sum(pred.values()),
        "pearson": pearson,
        "spearman": spearman,
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) / len(keys) * 100.0,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) * 100.0,
        "rescued_species": ",".join(sorted(set(rescued))),
    }


def minco_values_with_rescue(rows: pd.DataFrame) -> tuple[dict[str, float], list[str]]:
    rows = rows.loc[rows["gtdb_species"].astype(bool)].copy()
    rows["effective_depth"] = rows.apply(robust_effective_depth, axis=1)
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(species_genus)

    values: dict[str, float] = {}
    active_genus_max: dict[str, float] = {}
    active = rows.loc[rows["active_gate_pass"]].copy()
    for _, row in active.iterrows():
        species = str(row["gtdb_species"])
        genus = str(row["species_genus"])
        value = ev.safe_num(row["effective_depth"])
        if value <= 0.0:
            continue
        values[species] = max(values.get(species, 0.0), value)
        active_genus_max[genus] = max(active_genus_max.get(genus, 0.0), value)

    rescue_best: dict[str, pd.Series] = {}
    for _, row in rows.loc[~rows["active_gate_pass"]].iterrows():
        species = str(row["gtdb_species"])
        genus = str(row["species_genus"])
        active_value = active_genus_max.get(genus, 0.0)
        value = ev.safe_num(row["effective_depth"])
        if active_value <= 0.0 or value <= 0.0:
            continue
        if ev.safe_num(row.get("ANI_naive_calc")) < RESCUE_ANI_MIN:
            continue
        if ev.safe_num(row.get("XnY_ctx")) < RESCUE_XNY_MIN:
            continue
        if value < RESCUE_EFFECTIVE_MIN:
            continue
        if ev.safe_num(row.get("Reliable_Ref_zip_af")) > RESCUE_ZIP_AF_MAX:
            continue
        if value < RESCUE_ACTIVE_GENUS_RATIO * active_value:
            continue
        old = rescue_best.get(genus)
        if old is None or (
            value,
            ev.safe_num(row.get("XnY_ctx")),
        ) > (
            ev.safe_num(old.get("effective_depth")),
            ev.safe_num(old.get("XnY_ctx")),
        ):
            rescue_best[genus] = row

    rescued: list[str] = []
    for row in rescue_best.values():
        species = str(row["gtdb_species"])
        value = ev.safe_num(row["effective_depth"])
        values[species] = max(values.get(species, 0.0), value)
        rescued.append(species)

    return values, rescued


def sylph_values(rows: pd.DataFrame) -> dict[str, float]:
    values: dict[str, float] = {}
    for _, row in rows.loc[rows["gtdb_species"].astype(bool)].iterrows():
        species = str(row["gtdb_species"])
        value = ev.safe_num(row.get("Taxonomic_abundance")) / 100.0
        if value > 0.0:
            values[species] = max(values.get(species, 0.0), value)
    return values


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    sample_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        gold = ev.load_truth_profile(sample)
        minco_rows = more.load_minco_active_rows(ctx_marker_path(sample), by_accession, by_core)
        minco_values, rescued = minco_values_with_rescue(minco_rows)
        sample_rows.append(
            score_values(
                sample,
                "minco_ctx_marker_robust_intragenus_winner_rescue",
                minco_values,
                gold,
                rescued,
            )
        )

        sylph_rows = more.load_sylph_rows(more.sylph_profile_path(sample), by_accession, by_core)
        sample_rows.append(
            score_values(
                sample,
                "sylph_gtdb_profile_reported_taxonomic_abundance",
                sylph_values(sylph_rows),
                gold,
                [],
            )
        )

    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(EXP_DIR / "intragenus_winner_rescue_sample_metrics.tsv", sep="\t", index=False)

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
    mean_df.to_csv(EXP_DIR / "intragenus_winner_rescue_mean_metrics.tsv", sep="\t", index=False)
    print(mean_df.sort_values("l1_pct_points").to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
