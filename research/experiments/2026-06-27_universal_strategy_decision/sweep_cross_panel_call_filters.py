#!/usr/bin/env python3
"""Sweep output-only call filters across cached GTDB panels.

This is a post-hoc diagnostic. It starts from the current calibrated calls and
tests simple universal filters that only remove called rows using fields already
present in the MinCO profile output. A filter is promotable only if it improves
F1 safely across panels and samples; abundance is secondary.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import score_cami3_gtdb_source_readmap as cami3
import score_cami3_gtdb_taxid_transfer as taxid_score
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp
import sweep_cross_panel_abundance_variants as abundance


warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in divide",
    category=RuntimeWarning,
)

EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Sweep output-only call filters across cached GTDB panels.",
    )
    ap.add_argument(
        "--include-hmp-omitted",
        action="store_true",
        help="Add restored HMP omitted samples 2, 8, and 26 to a separate audit output set.",
    )
    ap.add_argument(
        "--output-prefix",
        default="cross_panel_call_filter",
        help="Output prefix under results/. Use a non-default prefix for opt-in extensions.",
    )
    return ap.parse_args()


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def arr(df: pd.DataFrame, col: str, default: float = 0.0) -> np.ndarray:
    return numeric(df, col, default).to_numpy(dtype=float)


def max_feature(calls: pd.DataFrame, left: str, right: str, default: float = 0.0) -> np.ndarray:
    return np.maximum(arr(calls, left, default), arr(calls, right, default))


def min_positive_feature(calls: pd.DataFrame, left: str, right: str, default: float = 0.0) -> np.ndarray:
    a = arr(calls, left, default)
    b = arr(calls, right, default)
    out = np.maximum(a, b)
    both = (a > 0.0) & (b > 0.0)
    out[both] = np.minimum(a[both], b[both])
    return out


def feature_arrays(calls: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        "prob": arr(calls, "calibrated_probability", 1.0),
        "reported_ani": arr(calls, "reported_ani", 0.0),
        "s_ani": arr(calls, "s_ANI_max", 0.0),
        "max_ani": max_feature(calls, "s_ANI_max", "u_ANI_max"),
        "s_xny": arr(calls, "s_XnY_ctx_max", 0.0),
        "u_xny": arr(calls, "u_XnY_ctx_max", 0.0),
        "max_xny": max_feature(calls, "s_XnY_ctx_max", "u_XnY_ctx_max"),
        "min_xny": min_positive_feature(calls, "s_XnY_ctx_max", "u_XnY_ctx_max"),
        "s_breadth": arr(calls, "s_Ref_breadth_max", 0.0),
        "u_breadth": arr(calls, "u_Ref_breadth_max", 0.0),
        "max_breadth": max_feature(calls, "s_Ref_breadth_max", "u_Ref_breadth_max"),
        "min_breadth": min_positive_feature(calls, "s_Ref_breadth_max", "u_Ref_breadth_max"),
        "s_real_min_af": arr(calls, "s_Real_min_align_fraction_max", 0.0),
        "u_real_min_af": arr(calls, "u_Real_min_align_fraction_max", 0.0),
        "max_real_min_af": max_feature(
            calls,
            "s_Real_min_align_fraction_max",
            "u_Real_min_align_fraction_max",
        ),
        "s_zip_af": arr(calls, "s_Ref_zip_af_max", 0.0),
        "u_zip_af": arr(calls, "u_Ref_zip_af_max", 0.0),
        "max_zip_af": max_feature(calls, "s_Ref_zip_af_max", "u_Ref_zip_af_max"),
        "s_mean_depth": arr(calls, "s_Ref_mean_depth_max", 0.0),
        "u_mean_depth": arr(calls, "u_Ref_mean_depth_max", 0.0),
        "max_mean_depth": max_feature(calls, "s_Ref_mean_depth_max", "u_Ref_mean_depth_max"),
        "s_hit_mean_depth": arr(calls, "s_Ref_hit_mean_depth_max", 0.0),
        "u_hit_mean_depth": arr(calls, "u_Ref_hit_mean_depth_max", 0.0),
        "max_hit_mean_depth": max_feature(
            calls,
            "s_Ref_hit_mean_depth_max",
            "u_Ref_hit_mean_depth_max",
        ),
        "s_depth_cv": arr(calls, "s_Ref_depth_cv_max", 0.0),
        "u_depth_cv": arr(calls, "u_Ref_depth_cv_max", 0.0),
        "max_depth_cv": max_feature(calls, "s_Ref_depth_cv_max", "u_Ref_depth_cv_max"),
        "s_read_frac": arr(calls, "s_Read_match_fraction_max", 0.0),
        "u_read_frac": arr(calls, "u_Read_match_fraction_max", 0.0),
        "max_read_frac": max_feature(
            calls,
            "s_Read_match_fraction_max",
            "u_Read_match_fraction_max",
        ),
        "s_block_frac": arr(calls, "s_Block_match_fraction_max", 0.0),
        "u_block_frac": arr(calls, "u_Block_match_fraction_max", 0.0),
        "max_block_frac": max_feature(
            calls,
            "s_Block_match_fraction_max",
            "u_Block_match_fraction_max",
        ),
        "split_unique_xny_ratio": arr(calls, "split_unique_xny_ratio", 0.0),
        "split_unique_breadth_ratio": arr(calls, "split_unique_breadth_ratio", 0.0),
        "calibrated_abundance": arr(calls, "calibrated_abundance", 0.0),
        "calibrated_abundance_raw": arr(calls, "calibrated_abundance_raw", 0.0),
    }


def threshold_filters(calls: pd.DataFrame) -> dict[str, np.ndarray]:
    features = feature_arrays(calls)
    n = len(calls)
    out: dict[str, np.ndarray] = {"current_calls": np.ones(n, dtype=bool)}

    ge_thresholds: dict[str, list[float]] = {
        "prob": [0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99],
        "reported_ani": [0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99],
        "s_ani": [0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99],
        "max_ani": [0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99],
        "s_xny": [10, 25, 50, 100, 200, 300, 500, 750, 1000],
        "max_xny": [10, 25, 50, 100, 200, 300, 500, 750, 1000],
        "min_xny": [10, 25, 50, 100, 200, 300, 500],
        "s_breadth": [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "max_breadth": [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "min_breadth": [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "s_real_min_af": [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "max_real_min_af": [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "s_zip_af": [0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "max_zip_af": [0.05, 0.10, 0.20, 0.30, 0.40, 0.50],
        "max_mean_depth": [0.25, 0.50, 1.0, 2.0, 3.0, 5.0],
        "max_hit_mean_depth": [0.50, 1.0, 2.0, 3.0, 5.0],
    }
    for name, thresholds in ge_thresholds.items():
        values = features[name]
        for threshold in thresholds:
            out[f"{name}_ge_{threshold:g}"] = values >= threshold

    le_thresholds: dict[str, list[float]] = {
        "max_depth_cv": [5, 10, 20, 50, 100],
        "split_unique_xny_ratio": [1.1, 1.25, 1.5, 2.0, 3.0, 5.0],
        "split_unique_breadth_ratio": [1.1, 1.25, 1.5, 2.0, 3.0],
    }
    for name, thresholds in le_thresholds.items():
        values = features[name]
        for threshold in thresholds:
            out[f"{name}_le_{threshold:g}"] = values <= threshold

    prob = features["prob"]
    max_xny = features["max_xny"]
    max_breadth = features["max_breadth"]
    max_real = features["max_real_min_af"]
    max_zip = features["max_zip_af"]
    max_depth = features["max_mean_depth"]
    reported_ani = features["reported_ani"]
    for p in [0.50, 0.70, 0.90, 0.95]:
        for xny in [25, 50, 100, 200, 500]:
            out[f"drop_prob_lt_{p:g}_and_xny_lt_{xny:g}"] = ~((prob < p) & (max_xny < xny))
        for breadth in [0.05, 0.10, 0.20, 0.30]:
            out[f"drop_prob_lt_{p:g}_and_breadth_lt_{breadth:g}"] = ~(
                (prob < p) & (max_breadth < breadth)
            )
        for real_af in [0.05, 0.10, 0.20, 0.30]:
            out[f"drop_prob_lt_{p:g}_and_realaf_lt_{real_af:g}"] = ~(
                (prob < p) & (max_real < real_af)
            )
        for zip_af in [0.10, 0.20, 0.30, 0.40]:
            out[f"drop_prob_lt_{p:g}_and_zipaf_lt_{zip_af:g}"] = ~(
                (prob < p) & (max_zip < zip_af)
            )

    for xny in [25, 50, 100, 200, 500]:
        for breadth in [0.05, 0.10, 0.20, 0.30]:
            out[f"drop_xny_lt_{xny:g}_and_breadth_lt_{breadth:g}"] = ~(
                (max_xny < xny) & (max_breadth < breadth)
            )
    for depth in [0.5, 1.0, 2.0]:
        for breadth in [0.05, 0.10, 0.20, 0.30]:
            out[f"drop_depth_lt_{depth:g}_and_breadth_lt_{breadth:g}"] = ~(
                (max_depth < depth) & (max_breadth < breadth)
            )
    for ani in [0.95, 0.96, 0.97, 0.98]:
        for breadth in [0.05, 0.10, 0.20]:
            out[f"drop_ani_lt_{ani:g}_and_breadth_lt_{breadth:g}"] = ~(
                (reported_ani < ani) & (max_breadth < breadth)
            )

    return out


def sample_records(include_hmp_omitted: bool = False) -> list[tuple[str, int]]:
    records = [
        ("cami2_toy_mouse_gut", 5),
        ("cami2_toy_mouse_gut", 6),
        ("cami2_toy_mouse_gut", 7),
        ("hmp_airskin_gtdb_source_abundance", 0),
        ("hmp_airskin_gtdb_source_abundance", 1),
        ("hmp_airskin_gtdb_source_abundance", 3),
        ("hmp_airskin_gtdb_source_abundance", 4),
        ("hmp_airskin_gtdb_source_abundance", 5),
        ("hmp_airskin_gtdb_source_abundance", 6),
        ("hmp_airskin_gtdb_source_abundance", 7),
        ("hmp_airskin_gtdb_source_abundance", 9),
        ("hmp_airskin_gtdb_source_abundance", 10),
        ("hmp_airskin_gtdb_source_abundance", 11),
        ("hmp_airskin_gtdb_source_abundance", 13),
        ("hmp_airskin_gtdb_source_abundance", 14),
        ("hmp_airskin_gtdb_source_abundance", 15),
        ("hmp_airskin_gtdb_source_abundance", 16),
        ("hmp_airskin_gtdb_source_abundance", 17),
        ("hmp_airskin_gtdb_source_abundance", 18),
        ("hmp_airskin_gtdb_source_abundance", 19),
        ("hmp_airskin_gtdb_source_abundance", 20),
        ("hmp_airskin_gtdb_source_abundance", 21),
        ("hmp_airskin_gtdb_source_abundance", 22),
        ("hmp_airskin_gtdb_source_abundance", 23),
        ("hmp_airskin_gtdb_source_abundance", 24),
        ("hmp_airskin_gtdb_source_abundance", 25),
        ("hmp_airskin_gtdb_source_abundance", 28),
        ("hmp_gastrooral_gtdb_source_abundance", 0),
        ("hmp_gastrooral_gtdb_source_abundance", 6),
        ("cami3_toy_human_gut_gtdb_source_readmap", 0),
        ("cami3_toy_human_gut_gtdb_source_readmap", 1),
        ("cami3_toy_human_gut_gtdb_source_readmap", 2),
    ]
    if include_hmp_omitted:
        records.extend(
            [
                ("hmp_airskin_gtdb_source_abundance", 2),
                ("hmp_airskin_gtdb_source_abundance", 8),
                ("hmp_airskin_gtdb_source_abundance", 26),
            ]
        )
    return records


def load_panel_inputs(
    panel: str,
    sample: int,
    loaders: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, str, Path]:
    if panel == "cami2_toy_mouse_gut":
        truth = abundance.toy_truth(sample)
        calls, collapse, path = abundance.load_toy_calls(sample)
    elif panel == "hmp_airskin_gtdb_source_abundance":
        truth = abundance.hmp_truth(sample)
        calls, collapse, path = abundance.load_hmp_calls(
            sample,
            loaders["hmp_taxmap"],
            loaders["by_accession"],
            loaders["by_core"],
        )
    elif panel == "hmp_gastrooral_gtdb_source_abundance":
        truth = abundance.hmp_gastro_truth(sample)
        calls, collapse, path = abundance.load_hmp_gastro_calls(
            sample,
            loaders["hmp_gastro_taxmap"],
            loaders["by_accession"],
            loaders["by_core"],
        )
    elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
        truth = abundance.cami3_truth(sample)
        calls, collapse, path = abundance.load_cami3_calls(
            sample,
            loaders["cami3_taxid_to_species"],
            loaders["cami3_name_to_species"],
            loaders["cami3_taxmap"],
            loaders["by_accession"],
            loaders["by_core"],
        )
    else:
        raise AssertionError(panel)
    return truth, calls, collapse, path


def official_cols(panel: str) -> tuple[str, str]:
    if panel == "cami2_toy_mouse_gut":
        return "L1_truth_only_pp", "Pearson_truth_only"
    return "L1_union_pp", "Pearson_union"


def summarize_scores(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    current_sample = scores.loc[scores["method"].eq("current_calls")][
        ["panel", "sample", "F1", "L1_union_pp", "L1_truth_only_pp", "Pearson_union", "Pearson_truth_only"]
    ].rename(
        columns={
            "F1": "current_F1",
            "L1_union_pp": "current_L1_union_pp",
            "L1_truth_only_pp": "current_L1_truth_only_pp",
            "Pearson_union": "current_Pearson_union",
            "Pearson_truth_only": "current_Pearson_truth_only",
        }
    )
    sample_summary = scores.merge(current_sample, on=["panel", "sample"], how="left")
    sample_summary["delta_F1"] = sample_summary["F1"] - sample_summary["current_F1"]
    sample_summary["official_L1_pp"] = np.where(
        sample_summary["panel"].eq("cami2_toy_mouse_gut"),
        sample_summary["L1_truth_only_pp"],
        sample_summary["L1_union_pp"],
    )
    sample_summary["current_official_L1_pp"] = np.where(
        sample_summary["panel"].eq("cami2_toy_mouse_gut"),
        sample_summary["current_L1_truth_only_pp"],
        sample_summary["current_L1_union_pp"],
    )
    sample_summary["delta_official_L1_pp"] = (
        sample_summary["official_L1_pp"] - sample_summary["current_official_L1_pp"]
    )
    sample_summary["official_Pearson"] = np.where(
        sample_summary["panel"].eq("cami2_toy_mouse_gut"),
        sample_summary["Pearson_truth_only"],
        sample_summary["Pearson_union"],
    )
    sample_summary["current_official_Pearson"] = np.where(
        sample_summary["panel"].eq("cami2_toy_mouse_gut"),
        sample_summary["current_Pearson_truth_only"],
        sample_summary["current_Pearson_union"],
    )
    sample_summary["delta_official_Pearson"] = (
        sample_summary["official_Pearson"] - sample_summary["current_official_Pearson"]
    )

    panel_rows: list[dict[str, object]] = []
    for (panel, method), sub in sample_summary.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        panel_rows.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique(), key=str))),
                "mean_F1": float(sub["F1"].mean()),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_delta_F1": float(sub["delta_F1"].mean()),
                "min_delta_F1": float(sub["delta_F1"].min()),
                "worsened_samples": int((sub["delta_F1"] < -1e-12).sum()),
                "improved_samples": int((sub["delta_F1"] > 1e-12).sum()),
                "mean_official_L1_pp": float(sub["official_L1_pp"].mean()),
                "mean_delta_official_L1_pp": float(sub["delta_official_L1_pp"].mean()),
                "max_worse_official_L1_pp": float(sub["delta_official_L1_pp"].max()),
                "mean_official_Pearson": float(sub["official_Pearson"].mean()),
                "mean_delta_official_Pearson": float(sub["delta_official_Pearson"].mean()),
            }
        )
    panel_summary = pd.DataFrame(panel_rows)

    current_panel = panel_summary.loc[panel_summary["method"].eq("current_calls")][
        ["panel", "mean_F1", "pooled_F1", "mean_official_L1_pp", "mean_official_Pearson"]
    ].rename(
        columns={
            "mean_F1": "current_panel_mean_F1",
            "pooled_F1": "current_panel_pooled_F1",
            "mean_official_L1_pp": "current_panel_L1_pp",
            "mean_official_Pearson": "current_panel_Pearson",
        }
    )
    panel_summary = panel_summary.merge(current_panel, on="panel", how="left")
    panel_summary["delta_panel_mean_F1"] = panel_summary["mean_F1"] - panel_summary["current_panel_mean_F1"]
    panel_summary["delta_panel_pooled_F1"] = panel_summary["pooled_F1"] - panel_summary["current_panel_pooled_F1"]
    panel_summary["delta_panel_L1_pp"] = (
        panel_summary["mean_official_L1_pp"] - panel_summary["current_panel_L1_pp"]
    )
    panel_summary["delta_panel_Pearson"] = (
        panel_summary["mean_official_Pearson"] - panel_summary["current_panel_Pearson"]
    )

    overall_rows: list[dict[str, object]] = []
    expected_panels = int(panel_summary["panel"].nunique())
    for method, psub in panel_summary.groupby("method", sort=True):
        ssub = sample_summary.loc[sample_summary["method"].eq(method)]
        panel_count = int(psub["panel"].nunique())
        sample_worse = int((ssub["delta_F1"] < -1e-12).sum())
        sample_improve = int((ssub["delta_F1"] > 1e-12).sum())
        panel_worse = int((psub["delta_panel_mean_F1"] < -1e-12).sum())
        panel_improve = int((psub["delta_panel_mean_F1"] > 1e-12).sum())
        mean_delta_f1 = float(ssub["delta_F1"].mean())
        min_delta_f1 = float(ssub["delta_F1"].min())
        mean_delta_l1 = float(ssub["delta_official_L1_pp"].mean())
        decision = "worse_or_equal"
        if method == "current_calls":
            decision = "current_default"
        elif panel_count != expected_panels:
            decision = "incomplete"
        elif sample_worse == 0 and sample_improve > 0:
            decision = "sample_safe_f1_candidate"
        elif panel_worse == 0 and panel_improve > 0:
            decision = "panel_safe_f1_candidate"
        elif mean_delta_f1 > 0.0 and sample_worse > 0:
            decision = "f1_tradeoff_not_default"
        elif abs(mean_delta_f1) <= 1e-12 and mean_delta_l1 < 0.0:
            decision = "f1_tie_l1_tradeoff"
        overall_rows.append(
            {
                "method": method,
                "panel_count": panel_count,
                "sample_count": int(len(ssub)),
                "mean_delta_F1": mean_delta_f1,
                "min_delta_F1": min_delta_f1,
                "improved_samples": sample_improve,
                "worsened_samples": sample_worse,
                "mean_panel_delta_F1": float(psub["delta_panel_mean_F1"].mean()),
                "min_panel_delta_F1": float(psub["delta_panel_mean_F1"].min()),
                "improved_panels": panel_improve,
                "worsened_panels": panel_worse,
                "mean_delta_L1_pp": mean_delta_l1,
                "max_worse_L1_pp": float(ssub["delta_official_L1_pp"].max()),
                "mean_delta_Pearson": float(ssub["delta_official_Pearson"].mean()),
                "decision": decision,
            }
        )
    overall = pd.DataFrame(overall_rows).sort_values(
        [
            "decision",
            "worsened_samples",
            "worsened_panels",
            "mean_delta_F1",
            "mean_delta_L1_pp",
            "method",
        ],
        ascending=[True, True, True, False, True, True],
    )
    return sample_summary, panel_summary, overall


def validation(panel_summary: pd.DataFrame, skip_official: bool = False) -> pd.DataFrame:
    if skip_official:
        return pd.DataFrame(
            [
                {
                    "panel": "all",
                    "F1_delta": np.nan,
                    "L1_delta_pp": np.nan,
                    "Pearson_delta": np.nan,
                    "source": "official validation skipped for opt-in extended sample set",
                }
            ]
        )
    official = pd.read_csv(RESULTS / "abundance_error_decomposition_summary.tsv", sep="\t")
    rows: list[dict[str, object]] = []
    for panel in sorted(panel_summary["panel"].unique()):
        ours = panel_summary.loc[
            panel_summary["panel"].eq(panel) & panel_summary["method"].eq("current_calls")
        ].iloc[0]
        theirs = official.loc[
            official["panel"].eq(panel) & ~official["method"].astype(str).str.startswith("sylph")
        ].iloc[0]
        official_l1 = (
            float(theirs["mean_truth_only_L1_pp"])
            if panel == "cami2_toy_mouse_gut"
            else float(theirs["mean_union_L1_pp"])
        )
        official_pearson = (
            float(theirs["mean_Pearson_truth_only"])
            if panel == "cami2_toy_mouse_gut"
            else float(theirs["mean_Pearson_union"])
        )
        rows.append(
            {
                "panel": panel,
                "F1_delta": float(ours["mean_F1"]) - float(theirs["mean_F1"]),
                "L1_delta_pp": float(ours["mean_official_L1_pp"]) - official_l1,
                "Pearson_delta": float(ours["mean_official_Pearson"]) - official_pearson,
                "source": "abundance_error_decomposition_summary.tsv",
            }
        )
    return pd.DataFrame(rows)


def safety_audit(overall: pd.DataFrame, evidence_prefix: str = "cross_panel_call_filter") -> pd.DataFrame:
    noncurrent = overall.loc[~overall["method"].eq("current_calls")].copy()
    if noncurrent.empty:
        return pd.DataFrame()
    evidence = f"{evidence_prefix}_overall.tsv"
    ranked = noncurrent.sort_values(
        ["mean_delta_F1", "worsened_samples", "mean_delta_L1_pp"],
        ascending=[False, True, True],
    )
    best_mean = ranked.iloc[0]
    sample_safe = noncurrent.loc[noncurrent["decision"].eq("sample_safe_f1_candidate")]
    panel_safe = noncurrent.loc[noncurrent["decision"].eq("panel_safe_f1_candidate")]
    tradeoff = noncurrent.loc[noncurrent["decision"].eq("f1_tradeoff_not_default")]
    rows = [
        {
            "metric": "tested_filters",
            "value": int(len(noncurrent)),
            "evidence": evidence,
            "decision": "output_only_posthoc",
        },
        {
            "metric": "sample_safe_f1_filters",
            "value": int(len(sample_safe)),
            "evidence": evidence,
            "decision": "must_be_positive_to_consider_default",
        },
        {
            "metric": "panel_safe_f1_filters",
            "value": int(len(panel_safe)),
            "evidence": evidence,
            "decision": "needs_sample_safety_or_external_validation",
        },
        {
            "metric": "f1_tradeoff_filters",
            "value": int(len(tradeoff)),
            "evidence": evidence,
            "decision": "not_default",
        },
        {
            "metric": "best_mean_f1_filter",
            "value": (
                f"{best_mean['method']};mean_delta_F1={best_mean['mean_delta_F1']:.6f};"
                f"worsened_samples={int(best_mean['worsened_samples'])};"
                f"worsened_panels={int(best_mean['worsened_panels'])};"
                f"mean_delta_L1_pp={best_mean['mean_delta_L1_pp']:.6f}"
            ),
            "evidence": evidence,
            "decision": str(best_mean["decision"]),
        },
    ]
    if len(sample_safe):
        best_safe = sample_safe.sort_values(
            ["mean_delta_F1", "mean_delta_L1_pp"],
            ascending=[False, True],
        ).iloc[0]
        promotion = f"candidate_requires_raw_wrapper_validation:{best_safe['method']}"
    else:
        promotion = "reject_output_only_call_filter_replacement"
    rows.append(
        {
            "metric": "promotion_decision",
            "value": promotion,
            "evidence": evidence,
            "decision": "current_call_gate_remains_default",
        }
    )
    return pd.DataFrame(rows)


def compact_methods(overall: pd.DataFrame) -> set[str]:
    """Keep repo artifacts small while preserving the decision-critical rows."""
    methods = {"current_calls"}
    noncurrent = overall.loc[~overall["method"].eq("current_calls")].copy()
    if noncurrent.empty:
        return methods
    best = noncurrent.sort_values(
        ["mean_delta_F1", "worsened_samples", "mean_delta_L1_pp"],
        ascending=[False, True, True],
    ).iloc[0]
    methods.add(str(best["method"]))
    safe = noncurrent.loc[
        noncurrent["decision"].isin(
            ["sample_safe_f1_candidate", "panel_safe_f1_candidate"]
        )
    ]
    methods.update(safe["method"].astype(str))
    return methods


def main() -> int:
    args = parse_args()
    output_prefix = args.output_prefix
    if args.include_hmp_omitted and output_prefix == "cross_panel_call_filter":
        output_prefix = "cross_panel_call_filter_with_hmp_omitted"
    RESULTS.mkdir(parents=True, exist_ok=True)

    toy_mod = abundance.decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    _cami3_wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _map_diag = (
        cami3.build_transfer_maps()
    )
    loaders = {
        "by_accession": by_accession,
        "by_core": by_core,
        "hmp_taxmap": hmp.parse_species_taxmap(hmp.TAXMAP),
        "hmp_gastro_taxmap": hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP),
        "cami3_taxid_to_species": cami3_taxid_to_species,
        "cami3_name_to_species": cami3_name_to_species,
        "cami3_taxmap": cami3.parse_species_taxmap(cami3.TAXMAP),
    }

    score_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    label_rows: list[dict[str, object]] = []
    for panel, sample in sample_records(args.include_hmp_omitted):
        truth, calls, collapse, path = load_panel_inputs(panel, sample, loaders)
        truth_species = set(truth["gtdb_species"].astype(str))
        filters = threshold_filters(calls)
        base_abundance = arr(calls, "calibrated_abundance", 0.0)
        for method, keep in filters.items():
            keep = np.asarray(keep, dtype=bool)
            pred = abundance.collapse_prediction(calls.loc[keep].copy(), base_abundance[keep], collapse)
            row = taxid_score.score_prediction(sample, method, pred, truth, {})
            row["panel"] = panel
            row["collapse_rule"] = collapse
            row["kept_rows"] = int(keep.sum())
            row["dropped_rows"] = int(len(keep) - keep.sum())
            score_rows.append(row)
        map_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "profile": str(path),
                "collapse_rule": collapse,
                "called_rows_mapped": len(calls),
                "called_species_mapped": calls["gtdb_species"].nunique(),
                "truth_species": len(truth_species),
            }
        )
        row_features = feature_arrays(calls)
        for idx, row in calls.reset_index(drop=True).iterrows():
            label_rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "gtdb_species": row["gtdb_species"],
                    "is_truth_species": str(row["gtdb_species"]) in truth_species,
                    "calibrated_probability": row_features["prob"][idx],
                    "calibrated_abundance": row_features["calibrated_abundance"][idx],
                    "reported_ani": row_features["reported_ani"][idx],
                    "max_xny": row_features["max_xny"][idx],
                    "max_breadth": row_features["max_breadth"][idx],
                    "max_real_min_af": row_features["max_real_min_af"][idx],
                    "max_zip_af": row_features["max_zip_af"][idx],
                    "max_mean_depth": row_features["max_mean_depth"][idx],
                    "max_depth_cv": row_features["max_depth_cv"][idx],
                    "profile": str(path),
                }
            )

    scores = pd.DataFrame(score_rows)
    sample_summary, panel_summary, overall = summarize_scores(scores)
    validation_df = validation(panel_summary, skip_official=args.include_hmp_omitted)
    audit_df = safety_audit(overall, output_prefix)
    kept_methods = compact_methods(overall)
    compact_scores = scores.loc[scores["method"].isin(kept_methods)].copy()
    compact_sample_summary = sample_summary.loc[sample_summary["method"].isin(kept_methods)].copy()

    compact_scores.to_csv(RESULTS / f"{output_prefix}_scores.tsv", sep="\t", index=False)
    compact_sample_summary.to_csv(
        RESULTS / f"{output_prefix}_sample_summary.tsv",
        sep="\t",
        index=False,
    )
    panel_summary.to_csv(RESULTS / f"{output_prefix}_panel_summary.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / f"{output_prefix}_overall.tsv", sep="\t", index=False)
    validation_df.to_csv(RESULTS / f"{output_prefix}_validation.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / f"{output_prefix}_safety_audit.tsv", sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(RESULTS / f"{output_prefix}_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(label_rows).to_csv(
        RESULTS / f"{output_prefix}_row_labels.tsv",
        sep="\t",
        index=False,
    )

    print("VALIDATION")
    print(validation_df.to_string(index=False))
    print("\nSAFETY")
    print(audit_df.to_string(index=False))
    print("\nTOP OVERALL")
    print(overall.head(25).to_string(index=False))
    print("\nPANEL TOP")
    print(
        panel_summary.sort_values(
            ["panel", "delta_panel_mean_F1", "mean_delta_official_L1_pp"],
            ascending=[True, False, True],
        )
        .groupby("panel", group_keys=False)
        .head(8)
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
