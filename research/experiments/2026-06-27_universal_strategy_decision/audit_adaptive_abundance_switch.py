#!/usr/bin/env python3
"""Audit whether output-derived switches can safely replace abundance defaults.

The fixed-call abundance sweep found variants that improve mean L1 but regress
some panels. This diagnostic asks whether a simple automatic switch, based only
on MinCO output features and not dataset labels or truth, can choose when to use
one of those variants.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

SCORES = RESULTS / "cross_panel_abundance_variant_scores.tsv"
OVERALL = RESULTS / "cross_panel_abundance_variant_overall.tsv"
MAPPING = RESULTS / "cross_panel_abundance_variant_mapping.tsv"

CURRENT = "current_calibrated_abundance"
TRIVIAL_EQUIVALENTS = {CURRENT, "current_calibrated_raw", "current_raw_power_p1"}
MAX_METHODS = 12


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def official_l1(row: pd.Series) -> float:
    if str(row["panel"]) == "cami2_toy_mouse_gut":
        return finite(row["L1_truth_only_pp"])
    return finite(row["L1_union_pp"])


def entropy(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    arr = np.maximum(arr.to_numpy(dtype=float), 0.0)
    total = float(arr.sum())
    if total <= 0.0:
        return 0.0
    p = arr / total
    p = p[p > 0.0]
    return float(-(p * np.log(p)).sum())


def quantile_thresholds(values: pd.Series) -> list[float]:
    clean = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if clean.empty or clean.nunique() <= 1:
        return []
    thresholds = []
    for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        value = float(clean.quantile(q))
        if math.isfinite(value):
            thresholds.append(value)
    return sorted(set(thresholds))


def sample_features(profile_path: Path) -> dict[str, float]:
    raw = pd.read_csv(profile_path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    calls = raw.loc[call].copy()
    if calls.empty:
        return {
            "n_calls": 0.0,
            "abundance_entropy": 0.0,
            "top1_abundance": 0.0,
            "top3_abundance": 0.0,
            "probability_mean": 0.0,
            "probability_min": 0.0,
        }

    abundance = numeric(calls, "calibrated_abundance")
    abundance_sorted = np.sort(np.maximum(abundance.to_numpy(dtype=float), 0.0))[::-1]
    prob = numeric(calls, "calibrated_probability", 1.0)
    s_xny = numeric(calls, "s_XnY_ctx_max")
    u_xny = numeric(calls, "u_XnY_ctx_max")
    s_zip = numeric(calls, "s_Ref_zip_af_max")
    u_zip = numeric(calls, "u_Ref_zip_af_max")
    s_depth = numeric(calls, "s_Ref_mean_depth_max")
    u_depth = numeric(calls, "u_Ref_mean_depth_max")
    s_breadth = numeric(calls, "s_Ref_breadth_max")
    u_breadth = numeric(calls, "u_Ref_breadth_max")
    s_cv = numeric(calls, "s_Ref_depth_cv_max")
    u_cv = numeric(calls, "u_Ref_depth_cv_max")
    xny_ratio = numeric(calls, "split_unique_xny_ratio", 1.0)
    breadth_ratio = numeric(calls, "split_unique_breadth_ratio", 1.0)
    ani_delta = numeric(calls, "split_unique_ani_delta")
    raw_mass = numeric(calls, "calibrated_abundance_raw")

    return {
        "n_calls": float(len(calls)),
        "abundance_entropy": entropy(abundance),
        "top1_abundance": float(abundance_sorted[0]) if abundance_sorted.size else 0.0,
        "top3_abundance": float(abundance_sorted[:3].sum()) if abundance_sorted.size else 0.0,
        "probability_mean": float(prob.mean()),
        "probability_median": float(prob.median()),
        "probability_min": float(prob.min()),
        "split_xny_median": float(s_xny.median()),
        "unique_xny_median": float(u_xny.median()),
        "xny_ratio_median": float(xny_ratio.median()),
        "breadth_ratio_median": float(breadth_ratio.median()),
        "zip_af_split_median": float(s_zip.median()),
        "zip_af_unique_median": float(u_zip.median()),
        "depth_split_median": float(s_depth.median()),
        "depth_unique_median": float(u_depth.median()),
        "breadth_split_median": float(s_breadth.median()),
        "breadth_unique_median": float(u_breadth.median()),
        "depth_cv_mean": float(pd.concat([s_cv, u_cv], ignore_index=True).mean()),
        "abs_ani_delta_mean": float(np.abs(ani_delta.to_numpy(dtype=float)).mean()),
        "raw_mass_sum": float(raw_mass.sum()),
    }


def load_features() -> pd.DataFrame:
    rows = []
    mapping = pd.read_csv(MAPPING, sep="\t")
    for item in mapping.itertuples(index=False):
        profile = Path(getattr(item, "profile"))
        row = {
            "panel": str(getattr(item, "panel")),
            "sample": str(getattr(item, "sample")),
            "profile": str(profile),
            "collapse_rule": getattr(item, "collapse_rule"),
            "called_rows_mapped": int(getattr(item, "called_rows_mapped")),
            "called_species_mapped": int(getattr(item, "called_species_mapped")),
        }
        row.update(sample_features(profile))
        rows.append(row)
    return pd.DataFrame(rows)


def choose_methods() -> list[str]:
    overall = pd.read_csv(OVERALL, sep="\t")
    candidates = overall.loc[~overall["method"].isin(TRIVIAL_EQUIVALENTS)].copy()
    candidates = candidates.sort_values(
        ["mean_delta_current_L1_pp", "max_worse_current_L1_pp", "method"],
        kind="mergesort",
    )
    return list(candidates["method"].head(MAX_METHODS))


def load_score_matrix(methods: list[str]) -> pd.DataFrame:
    scores = pd.read_csv(SCORES, sep="\t")
    scores["sample"] = scores["sample"].astype(str)
    scores["official_L1_pp"] = scores.apply(official_l1, axis=1)
    keep = scores.loc[scores["method"].isin([CURRENT, *methods])].copy()
    current = keep.loc[
        keep["method"].eq(CURRENT),
        ["panel", "sample", "official_L1_pp", "Pearson_union", "Pearson_truth_only", "F1"],
    ].rename(
        columns={
            "official_L1_pp": "current_official_L1_pp",
            "Pearson_union": "current_Pearson_union",
            "Pearson_truth_only": "current_Pearson_truth_only",
            "F1": "current_F1",
        }
    )
    return keep.merge(current, on=["panel", "sample"], how="left")


def evaluate_rule(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    method: str,
    feature: str,
    op: str,
    threshold: float,
) -> dict[str, object]:
    condition = features[feature].astype(float) >= threshold if op == ">=" else features[feature].astype(float) <= threshold
    switched = set(
        zip(
            features.loc[condition, "panel"].astype(str),
            features.loc[condition, "sample"].astype(str),
        )
    )
    rows = []
    for row in scores.loc[scores["method"].isin([CURRENT, method])].itertuples(index=False):
        key = (str(getattr(row, "panel")), str(getattr(row, "sample")))
        should_use_variant = key in switched
        if should_use_variant != (getattr(row, "method") == method):
            continue
        current_l1 = finite(getattr(row, "current_official_L1_pp"))
        official_l1 = finite(getattr(row, "official_L1_pp"))
        rows.append(
            {
                "panel": key[0],
                "sample": key[1],
                "method_used": getattr(row, "method"),
                "official_L1_pp": official_l1,
                "current_official_L1_pp": current_l1,
                "delta_current_L1_pp": official_l1 - current_l1,
                "switched": should_use_variant,
            }
        )
    chosen = pd.DataFrame(rows)
    sample_count = int(chosen.shape[0])
    if sample_count == 0:
        raise RuntimeError("empty rule evaluation")
    panel_delta = (
        chosen.groupby("panel", as_index=False)["delta_current_L1_pp"]
        .mean()
        .rename(columns={"delta_current_L1_pp": "panel_delta_current_L1_pp"})
    )
    return {
        "variant_method": method,
        "feature": feature,
        "op": op,
        "threshold": threshold,
        "switched_samples": int(chosen["switched"].sum()),
        "sample_count": sample_count,
        "mean_sample_delta_L1_pp": float(chosen["delta_current_L1_pp"].mean()),
        "max_sample_worse_L1_pp": float(chosen["delta_current_L1_pp"].max()),
        "improved_samples": int((chosen["delta_current_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((chosen["delta_current_L1_pp"] > 1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["panel_delta_current_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["panel_delta_current_L1_pp"].max()),
        "improved_panels": int((panel_delta["panel_delta_current_L1_pp"] < -1e-9).sum()),
        "worsened_panels": int((panel_delta["panel_delta_current_L1_pp"] > 1e-9).sum()),
    }


def sample_outcomes_for_rule(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    rule: pd.Series,
) -> pd.DataFrame:
    feature = str(rule["feature"])
    threshold = finite(rule["threshold"])
    op = str(rule["op"])
    method = str(rule["variant_method"])
    condition = features[feature].astype(float) >= threshold if op == ">=" else features[feature].astype(float) <= threshold
    switched = set(
        zip(
            features.loc[condition, "panel"].astype(str),
            features.loc[condition, "sample"].astype(str),
        )
    )
    rows = []
    for row in scores.loc[scores["method"].isin([CURRENT, method])].itertuples(index=False):
        key = (str(getattr(row, "panel")), str(getattr(row, "sample")))
        should_use_variant = key in switched
        if should_use_variant != (getattr(row, "method") == method):
            continue
        current_l1 = finite(getattr(row, "current_official_L1_pp"))
        official_l1 = finite(getattr(row, "official_L1_pp"))
        rows.append(
            {
                "panel": key[0],
                "sample": key[1],
                "method_used": getattr(row, "method"),
                "official_L1_pp": official_l1,
                "current_official_L1_pp": current_l1,
                "delta_current_L1_pp": official_l1 - current_l1,
                "switched": should_use_variant,
            }
        )
    return pd.DataFrame(rows)


def panel_stats(rows: pd.DataFrame) -> dict[str, float | int]:
    panel_delta = (
        rows.groupby("panel", as_index=False)["delta_current_L1_pp"]
        .mean()
        .rename(columns={"delta_current_L1_pp": "panel_delta_current_L1_pp"})
    )
    return {
        "sample_count": int(rows.shape[0]),
        "mean_sample_delta_L1_pp": float(rows["delta_current_L1_pp"].mean()),
        "max_sample_worse_L1_pp": float(rows["delta_current_L1_pp"].max()),
        "improved_samples": int((rows["delta_current_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((rows["delta_current_L1_pp"] > 1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["panel_delta_current_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["panel_delta_current_L1_pp"].max()),
        "improved_panels": int((panel_delta["panel_delta_current_L1_pp"] < -1e-9).sum()),
        "worsened_panels": int((panel_delta["panel_delta_current_L1_pp"] > 1e-9).sum()),
    }


def leave_one_panel_out(
    rules: pd.DataFrame,
    scores: pd.DataFrame,
    features: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    panels = sorted(features["panel"].astype(str).unique())
    for holdout in panels:
        candidates = []
        for rule in rules.itertuples(index=False):
            rule_series = pd.Series(rule._asdict())
            outcomes = sample_outcomes_for_rule(scores, features, rule_series)
            train = outcomes.loc[~outcomes["panel"].eq(holdout)].copy()
            test = outcomes.loc[outcomes["panel"].eq(holdout)].copy()
            if train.empty or test.empty:
                continue
            train_stats = panel_stats(train)
            test_stats = panel_stats(test)
            eligible = (
                int(train_stats["worsened_panels"]) == 0
                and int(train_stats["improved_panels"]) > 0
            )
            candidates.append(
                {
                    "holdout_panel": holdout,
                    "variant_method": rule_series["variant_method"],
                    "feature": rule_series["feature"],
                    "op": rule_series["op"],
                    "threshold": rule_series["threshold"],
                    "train_strict_panel_safe": eligible,
                    "train_mean_panel_delta_L1_pp": train_stats["mean_panel_delta_L1_pp"],
                    "train_max_panel_worse_L1_pp": train_stats["max_panel_worse_L1_pp"],
                    "train_improved_panels": train_stats["improved_panels"],
                    "train_worsened_panels": train_stats["worsened_panels"],
                    "holdout_mean_panel_delta_L1_pp": test_stats["mean_panel_delta_L1_pp"],
                    "holdout_max_sample_worse_L1_pp": test_stats["max_sample_worse_L1_pp"],
                    "holdout_improved_samples": test_stats["improved_samples"],
                    "holdout_worsened_samples": test_stats["worsened_samples"],
                    "holdout_switched_samples": int(test["switched"].sum()),
                    "holdout_sample_count": int(test.shape[0]),
                }
            )
        cand_df = pd.DataFrame(candidates)
        strict = cand_df.loc[cand_df["train_strict_panel_safe"]].copy()
        pool = strict if not strict.empty else cand_df
        selected = pool.sort_values(
            ["train_mean_panel_delta_L1_pp", "train_max_panel_worse_L1_pp"],
            kind="mergesort",
        ).iloc[0].to_dict()
        selected["selection_pool"] = "strict_train_panel_safe" if not strict.empty else "best_train_mean"
        rows.append(selected)
    return pd.DataFrame(rows)


def audit_rows(rules: pd.DataFrame, lopo: pd.DataFrame, methods: list[str]) -> list[dict[str, object]]:
    nontrivial = rules.loc[rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])]
    strict_sample = nontrivial.loc[
        nontrivial["worsened_samples"].eq(0) & nontrivial["improved_samples"].gt(0)
    ]
    strict_panel = nontrivial.loc[
        nontrivial["worsened_panels"].eq(0) & nontrivial["improved_panels"].gt(0)
    ]
    best_sample = nontrivial.sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"],
        kind="mergesort",
    ).iloc[0]
    best_panel = nontrivial.sort_values(
        ["mean_panel_delta_L1_pp", "max_panel_worse_L1_pp"],
        kind="mergesort",
    ).iloc[0]
    lopo_mean = float(lopo["holdout_mean_panel_delta_L1_pp"].mean())
    lopo_worse = int((lopo["holdout_mean_panel_delta_L1_pp"] > 1e-9).sum())
    lopo_improved = int((lopo["holdout_mean_panel_delta_L1_pp"] < -1e-9).sum())
    lopo_max_worse = float(lopo["holdout_mean_panel_delta_L1_pp"].max())
    if int(strict_panel.shape[0]) and lopo_worse == 0 and lopo_improved > 0:
        promotion_value = "candidate_requires_raw_wrapper_validation"
        promotion_decision = "do_not_promote_without_independent_holdout"
    elif int(strict_panel.shape[0]):
        promotion_value = "reject_as_overfit_by_leave_one_panel_out"
        promotion_decision = "current_calibrated_abundance_remains_default"
    else:
        promotion_value = "reject_simple_output_feature_abundance_switch"
        promotion_decision = "current_calibrated_abundance_remains_default"
    return [
        {
            "metric": "candidate_methods",
            "value": len(methods),
            "evidence": ",".join(methods),
            "decision": "top_fixed_call_variants_by_mean_L1",
        },
        {
            "metric": "tested_switch_rules",
            "value": int(rules.shape[0]),
            "evidence": "output_feature_thresholds",
            "decision": "truth_used_only_for_scoring",
        },
        {
            "metric": "nontrivial_switch_rules",
            "value": int(nontrivial.shape[0]),
            "evidence": "0 < switched_samples < sample_count",
            "decision": "eligible_for_default_audit",
        },
        {
            "metric": "strict_sample_safe_switches",
            "value": int(strict_sample.shape[0]),
            "evidence": "adaptive_abundance_switch_rules.tsv",
            "decision": "must_be_zero_to_reject_simple_output_switch",
        },
        {
            "metric": "strict_panel_safe_switches",
            "value": int(strict_panel.shape[0]),
            "evidence": "adaptive_abundance_switch_rules.tsv",
            "decision": "must_be_zero_to_reject_simple_output_switch",
        },
        {
            "metric": "lopo_holdout_panels_improved",
            "value": lopo_improved,
            "evidence": (
                f"mean_delta={lopo_mean:.6f};worsened={lopo_worse};"
                f"max_worse={lopo_max_worse:.6f}"
            ),
            "decision": "leave_one_panel_out_validation",
        },
        {
            "metric": "best_sample_mean_switch",
            "value": (
                f"{best_sample['variant_method']};{best_sample['feature']}{best_sample['op']}"
                f"{best_sample['threshold']:.6g};mean_delta={best_sample['mean_sample_delta_L1_pp']:.6f};"
                f"worsened_samples={int(best_sample['worsened_samples'])};"
                f"max_worse={best_sample['max_sample_worse_L1_pp']:.6f};"
                f"switched={int(best_sample['switched_samples'])}"
            ),
            "evidence": "adaptive_abundance_switch_rules.tsv",
            "decision": "tradeoff_not_default",
        },
        {
            "metric": "best_panel_mean_switch",
            "value": (
                f"{best_panel['variant_method']};{best_panel['feature']}{best_panel['op']}"
                f"{best_panel['threshold']:.6g};mean_delta={best_panel['mean_panel_delta_L1_pp']:.6f};"
                f"worsened_panels={int(best_panel['worsened_panels'])};"
                f"max_worse={best_panel['max_panel_worse_L1_pp']:.6f};"
                f"switched={int(best_panel['switched_samples'])}"
            ),
            "evidence": "adaptive_abundance_switch_rules.tsv",
            "decision": "tradeoff_not_default",
        },
        {
            "metric": "promotion_decision",
            "value": promotion_value,
            "evidence": "adaptive_abundance_switch_rules.tsv;adaptive_abundance_switch_lopo.tsv",
            "decision": promotion_decision,
        },
    ]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    methods = choose_methods()
    features = load_features()
    scores = load_score_matrix(methods)
    feature_cols = [
        col
        for col in features.columns
        if col not in {"panel", "sample", "profile", "collapse_rule"}
        and pd.api.types.is_numeric_dtype(features[col])
        and features[col].nunique(dropna=True) > 1
    ]
    rules = []
    for method in methods:
        for feature in feature_cols:
            for threshold in quantile_thresholds(features[feature]):
                for op in [">=", "<="]:
                    rules.append(evaluate_rule(scores, features, method, feature, op, threshold))
    rule_df = pd.DataFrame(rules).sort_values(
        ["mean_panel_delta_L1_pp", "max_panel_worse_L1_pp", "mean_sample_delta_L1_pp"],
        kind="mergesort",
    )
    feature_df = features.sort_values(["panel", "sample"], kind="mergesort")
    lopo = leave_one_panel_out(rule_df, scores, features)
    audit = audit_rows(rule_df, lopo, methods)

    rule_df.to_csv(RESULTS / "adaptive_abundance_switch_rules.tsv", sep="\t", index=False)
    lopo.to_csv(RESULTS / "adaptive_abundance_switch_lopo.tsv", sep="\t", index=False)
    feature_df.to_csv(RESULTS / "adaptive_abundance_switch_features.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "adaptive_abundance_switch_audit.tsv",
        audit,
        ["metric", "value", "evidence", "decision"],
    )
    print(pd.DataFrame(audit).to_string(index=False))
    print("\nTOP RULES")
    print(rule_df.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
