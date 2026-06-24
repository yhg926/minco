#!/usr/bin/env python3
"""Score a hybrid full-sketch / unique-marker fallback strategy on GTDB toy mouse.

Prototype rule:
- use ctx-marker rows when the reference has enough ctx-marker contexts;
- use full-sketch split rows only when the reference is marker-poor;
- score calls with the existing active gate and abundance with robust effective
  depth plus the previous conservative intra-genus rescue.

This tests the idea before changing MinCO's internal data structures.
"""

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

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import evaluate_abundance_estimators as ev  # noqa: E402
import score_intragenus_winner_rescue as rescue  # noqa: E402
import score_more_toymouse_gtdb as more  # noqa: E402


SAMPLES = [0, 1, 2]
CTX_MARKER_PSM_PATH = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv")


def full_split_path(sample: int) -> Path:
    return OLD_ABUND_EXP / f"toymouse_sample{sample}_s2000_dedup_full_split_naive_product_topfrac_median025.tsv"


def load_psmp_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip():
                continue
            size_s, ref = raw.rstrip("\n").split("\t")[:2]
            acc = truth.extract_accession(ref)
            try:
                size = int(float(size_s))
            except ValueError:
                continue
            if acc:
                sizes[acc] = size
    return sizes


def robust_effective_depth(row: pd.Series) -> float:
    median = ev.safe_num(row.get("Reliable_Ref_hit_median_depth"))
    if median >= rescue.MEDIAN_DEPTH_CUTOFF:
        return median
    mean = ev.safe_num(row.get("Reliable_Ref_mean_depth"))
    af = max(ev.safe_num(row.get("Reliable_Ref_zip_af")), ev.EPS)
    return mean / (af**rescue.AF_EXPONENT) if mean > 0.0 else 0.0


def numeric(rows: pd.DataFrame, col: str) -> pd.Series:
    if col not in rows.columns:
        return pd.Series([0.0] * len(rows), index=rows.index)
    return pd.to_numeric(rows[col], errors="coerce").fillna(0.0)


def load_minco_active_rows_fast(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows["accession"] = rows["Ref"].map(truth.extract_accession)
    if "Ref_annotation" in rows.columns:
        missing = rows["accession"] == ""
        rows.loc[missing, "accession"] = rows.loc[missing, "Ref_annotation"].map(truth.extract_accession)

    accession_map = {}
    for acc in rows["accession"].dropna().astype(str).unique():
        rec, method, key = truth.lookup_accession(acc, by_accession, by_core)
        accession_map[acc] = (
            rec["gtdb_species"] if rec else "",
            rec["gtdb_taxonomy"] if rec else "",
            method,
            key,
        )
    rows["gtdb_species"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[0])
    rows["gtdb_taxonomy"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[1])
    rows["gtdb_mapping_method"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[2])
    rows["gtdb_mapping_key"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[3])

    rows = truth.add_naive_ani(rows)
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_median_depth",
        "Reliable_Ref_hit_depth_variance",
        "Normalized_abundance_depth",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
    ]:
        rows[col] = numeric(rows, col)

    rows["Reliable_ztp_af"] = [
        more.ztp_adjusted_af(float(b), float(m))
        for b, m in zip(rows["Reliable_Ref_breadth"], rows["Reliable_Ref_hit_mean_depth"])
    ]
    rows["ANI_from_Reliable_ztp_af"] = [
        1.0 + math.log(max(float(af), 1e-300)) / more.EFFECTIVE_CTX_LENGTH
        for af in rows["Reliable_ztp_af"]
    ]
    rows["ANI_AF_delta"] = rows["ANI_naive_calc"] - rows["ANI_from_Reliable_ztp_af"]
    rows["Reliable_depth_vmr"] = [
        (float(var) / float(mean)) if float(mean) > 0.0 else 0.0
        for mean, var in zip(rows["Reliable_Ref_hit_mean_depth"], rows["Reliable_Ref_hit_depth_variance"])
    ]
    rows["active_delta_trigger"] = (
        (rows["Reliable_Ref_hit_mean_depth"] > more.ACTIVE_MEAN_DEPTH_MIN)
        & (rows["Reliable_depth_vmr"] > more.ACTIVE_VMR_MIN)
    )
    rows["active_delta_pass"] = ~rows["active_delta_trigger"] | (rows["ANI_AF_delta"] < more.ACTIVE_DELTA_MAX)
    rows["active_gate_pass"] = (
        (rows["XnY_ctx"] >= more.ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > more.ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= more.ACTIVE_RELIABLE_ZTP_AF_FLOOR)
        & rows["active_delta_pass"]
        & rows["gtdb_species"].astype(bool)
    )
    return rows


def attach_mode(rows: pd.DataFrame, mode: str, ctx_marker_sizes: dict[str, int]) -> pd.DataFrame:
    rows = rows.copy()
    rows["source_mode"] = mode
    if "accession" not in rows.columns:
        rows["accession"] = rows["Ref"].map(truth.extract_accession)
    rows["ctx_marker_size"] = rows["accession"].map(ctx_marker_sizes).fillna(0).astype(int)
    median = numeric(rows, "Reliable_Ref_hit_median_depth")
    mean = numeric(rows, "Reliable_Ref_mean_depth")
    af = numeric(rows, "Reliable_Ref_zip_af").clip(lower=ev.EPS)
    rows["effective_depth"] = mean.where(median < rescue.MEDIAN_DEPTH_CUTOFF, median)
    rows.loc[median < rescue.MEDIAN_DEPTH_CUTOFF, "effective_depth"] = (
        mean.loc[median < rescue.MEDIAN_DEPTH_CUTOFF] / (af.loc[median < rescue.MEDIAN_DEPTH_CUTOFF] ** rescue.AF_EXPONENT)
    )
    rows.loc[mean <= 0.0, "effective_depth"] = 0.0
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(rescue.species_genus)
    return rows


def load_sample_views(
    sample: int,
    by_accession,
    by_core,
    ctx_marker_sizes: dict[str, int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker = attach_mode(
        load_minco_active_rows_fast(rescue.ctx_marker_path(sample), by_accession, by_core),
        "ctx_marker",
        ctx_marker_sizes,
    )
    full = attach_mode(
        load_minco_active_rows_fast(full_split_path(sample), by_accession, by_core),
        "full_fallback",
        ctx_marker_sizes,
    )
    return marker, full


def hybrid_rows_from_views(
    marker: pd.DataFrame,
    full: pd.DataFrame,
    marker_size_cutoff: int,
) -> pd.DataFrame:
    marker_part = marker.loc[marker["ctx_marker_size"] >= marker_size_cutoff].copy()
    full_part = full.loc[full["ctx_marker_size"] < marker_size_cutoff].copy()
    rows = pd.concat([marker_part, full_part], ignore_index=True, sort=False).fillna(0)
    return rows


def minco_values_with_rescue(rows: pd.DataFrame) -> tuple[dict[str, float], list[str]]:
    rows = rows.loc[rows["gtdb_species"].astype(bool)].copy()
    active = rows.loc[rows["active_gate_pass"]].copy()
    marker_species = set(
        active.loc[active["source_mode"] == "ctx_marker", "gtdb_species"].astype(str)
    )
    if marker_species:
        active = active.loc[
            ~(
                (active["source_mode"] == "full_fallback")
                & active["gtdb_species"].astype(str).isin(marker_species)
            )
        ].copy()
    values = active.loc[active["effective_depth"] > 0.0].groupby("gtdb_species")["effective_depth"].max().to_dict()
    active_genus_max = active.loc[active["effective_depth"] > 0.0].groupby("species_genus")["effective_depth"].max()

    candidates = rows.loc[
        ~rows["active_gate_pass"]
        & (rows["ANI_naive_calc"] >= rescue.RESCUE_ANI_MIN)
        & (rows["XnY_ctx"] >= rescue.RESCUE_XNY_MIN)
        & (rows["effective_depth"] >= rescue.RESCUE_EFFECTIVE_MIN)
        & (rows["Reliable_Ref_zip_af"] <= rescue.RESCUE_ZIP_AF_MAX)
    ].copy()
    rescued: list[str] = []
    if not candidates.empty and not active_genus_max.empty:
        candidates["active_genus_value"] = candidates["species_genus"].map(active_genus_max).fillna(0.0)
        candidates = candidates.loc[
            (candidates["active_genus_value"] > 0.0)
            & (candidates["effective_depth"] >= rescue.RESCUE_ACTIVE_GENUS_RATIO * candidates["active_genus_value"])
        ].copy()
        if not candidates.empty:
            winners = (
                candidates.sort_values(["species_genus", "effective_depth", "XnY_ctx"])
                .groupby("species_genus", as_index=False)
                .tail(1)
            )
            for _, row in winners.iterrows():
                species = str(row["gtdb_species"])
                if row.get("source_mode", "") == "full_fallback" and species in values:
                    continue
                value = float(row["effective_depth"])
                values[species] = max(float(values.get(species, 0.0)), value)
                rescued.append(species)
    return values, rescued


def score_one(sample: int, method: str, rows: pd.DataFrame) -> dict[str, object]:
    gold = ev.load_truth_profile(sample)
    values, rescued = minco_values_with_rescue(rows)
    scored = rescue.score_values(sample, method, values, gold, rescued)
    selected = rows.loc[rows["active_gate_pass"] & rows["gtdb_species"].astype(bool)].copy()
    scored["selected_rows"] = int(len(selected))
    scored["selected_species"] = int(selected["gtdb_species"].astype(str).nunique())
    scored["selected_full_fallback_species"] = int(
        selected.loc[selected["source_mode"] == "full_fallback", "gtdb_species"].astype(str).nunique()
    )
    scored["selected_marker_species"] = int(
        selected.loc[selected["source_mode"] == "ctx_marker", "gtdb_species"].astype(str).nunique()
    )
    return scored


def main() -> int:
    if not CTX_MARKER_PSM_PATH.exists():
        raise SystemExit(f"missing ctx-marker psmp size file: {CTX_MARKER_PSM_PATH}")
    ctx_marker_sizes = load_psmp_sizes(CTX_MARKER_PSM_PATH)
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)

    sample_rows: list[dict[str, object]] = []
    diagnostic_rows: list[dict[str, object]] = []
    cutoffs = [0, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500]
    sample_views: dict[int, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for sample in SAMPLES:
        sample_views[sample] = load_sample_views(sample, by_accession, by_core, ctx_marker_sizes)

    for cutoff in cutoffs:
        for sample in SAMPLES:
            marker, full = sample_views[sample]
            rows = hybrid_rows_from_views(marker, full, cutoff)
            method = f"hybrid_marker_ge{cutoff}_else_full"
            sample_rows.append(score_one(sample, method, rows))
            diagnostic_rows.append(
                {
                    "sample": sample,
                    "marker_size_cutoff": cutoff,
                    "rows_total": int(len(rows)),
                    "rows_marker": int((rows["source_mode"] == "ctx_marker").sum()),
                    "rows_full_fallback": int((rows["source_mode"] == "full_fallback").sum()),
                    "species_total": int(rows["gtdb_species"].astype(str).nunique()),
                    "species_full_fallback": int(
                        rows.loc[rows["source_mode"] == "full_fallback", "gtdb_species"].astype(str).nunique()
                    ),
                }
            )

    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(EXP_DIR / "hybrid_mouse_sample_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(EXP_DIR / "hybrid_mouse_diagnostics.tsv", sep="\t", index=False)

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
        "selected_full_fallback_species",
        "selected_marker_species",
    ]
    mean_df = sample_df.groupby("method", as_index=False)[mean_cols].mean()
    mean_df["precision"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FP"])
    mean_df["recall"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FN"])
    mean_df["F1"] = 2.0 * mean_df["precision"] * mean_df["recall"] / (mean_df["precision"] + mean_df["recall"])
    mean_df = mean_df.sort_values(["F1", "l1_pct_points"], ascending=[False, True])
    mean_df.to_csv(EXP_DIR / "hybrid_mouse_mean_metrics.tsv", sep="\t", index=False)
    print(mean_df.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
