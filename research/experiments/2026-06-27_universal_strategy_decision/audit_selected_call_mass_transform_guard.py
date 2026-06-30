#!/usr/bin/env python3
"""Audit output-derived guards for selected-call abundance transforms."""

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SCORES_TSV = RESULTS / "selected_call_mass_transform_scores.tsv"
OVERALL_TSV = RESULTS / "selected_call_mass_transform_overall.tsv"
SELECTED_SCORES_TSV = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"

FEATURES_TSV = RESULTS / "selected_call_mass_transform_guard_features.tsv"
RULES_TSV = RESULTS / "selected_call_mass_transform_guard_rules.tsv"
LOPO_TSV = RESULTS / "selected_call_mass_transform_guard_lopo.tsv"
AUDIT_TSV = RESULTS / "selected_call_mass_transform_guard_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_MASS_TRANSFORM_GUARD.md"

MAX_METHODS = 8
MAX_SINGLE_RULES_FOR_PAIRS = 24
MAX_LOPO_RULES = 512
MAX_WRITTEN_RULES = 1000
PROFILE_FEATURE_COLUMNS = {
    "calibrated_call",
    "candidate_rescue_added",
    "candidate_surface_added",
    "calibrated_abundance",
    "calibrated_abundance_raw",
    "s_XnY_ctx_max",
    "u_XnY_ctx_max",
    "s_Ref_breadth_max",
    "u_Ref_breadth_max",
    "s_Ref_mean_depth_max",
    "u_Ref_mean_depth_max",
    "s_Ref_zip_af_max",
    "u_Ref_zip_af_max",
    "calibrated_probability",
    "reported_ani",
    "species_name",
}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def entropy(values: pd.Series) -> float:
    arr = np.maximum(pd.to_numeric(values, errors="coerce").fillna(0.0).to_numpy(dtype=float), 0.0)
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
    thresholds: list[float] = []
    for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        value = float(clean.quantile(q))
        if math.isfinite(value):
            thresholds.append(value)
    return sorted(set(thresholds))


def genus_from_name(value: object) -> str:
    text = str(value or "").strip()
    if text.startswith("s__"):
        text = text[3:]
    return text.split()[0] if text else ""


def candidate_methods() -> list[str]:
    overall = pd.read_csv(OVERALL_TSV, sep="\t")
    candidates = overall.loc[
        overall["mean_sample_L1_delta_pp"].lt(0.0)
        & overall["improved_samples"].gt(overall["worsened_samples"])
    ].copy()
    return candidates.sort_values(
        ["worsened_samples", "mean_sample_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
        kind="mergesort",
    )["method"].astype(str).head(MAX_METHODS).tolist()


def profile_features(path: Path) -> dict[str, float]:
    with path.open() as handle:
        available = handle.readline().rstrip("\n").split("\t")
    usecols = [col for col in available if col in PROFILE_FEATURE_COLUMNS]
    raw = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    called = bool_series(raw, "calibrated_call")
    calls = raw.loc[called].copy()
    if calls.empty:
        return {"n_calls": 0.0}
    rescue = bool_series(calls, "candidate_rescue_added")
    surface = bool_series(calls, "candidate_surface_added")
    candidate = rescue | surface
    base = ~candidate

    abundance = numeric(calls, "calibrated_abundance")
    raw_mass = numeric(calls, "calibrated_abundance_raw")
    base_raw = raw_mass.loc[base]
    candidate_raw = raw_mass.loc[candidate]
    s_xny = numeric(calls, "s_XnY_ctx_max")
    u_xny = numeric(calls, "u_XnY_ctx_max")
    max_xny = pd.concat([s_xny, u_xny], axis=1).max(axis=1)
    s_breadth = numeric(calls, "s_Ref_breadth_max")
    u_breadth = numeric(calls, "u_Ref_breadth_max")
    s_depth = numeric(calls, "s_Ref_mean_depth_max")
    u_depth = numeric(calls, "u_Ref_mean_depth_max")
    s_zip = numeric(calls, "s_Ref_zip_af_max")
    u_zip = numeric(calls, "u_Ref_zip_af_max")
    prob = numeric(calls, "calibrated_probability", 1.0)
    reported_ani = numeric(calls, "reported_ani")
    species = calls.get("species_name", pd.Series("", index=calls.index)).astype(str)
    genus = species.map(genus_from_name)
    genus_counts = genus.loc[base].value_counts()
    multi_genus = set(genus_counts.loc[genus_counts > 1].index.astype(str))
    base_genus = genus.loc[base]
    base_multi = base_genus.astype(str).isin(multi_genus)
    base_mass = float(base_raw.sum())
    candidate_mass = float(candidate_raw.sum())
    raw_total = float(raw_mass.sum())
    abundance_sorted = np.sort(np.maximum(abundance.to_numpy(dtype=float), 0.0))[::-1]

    return {
        "n_calls": float(len(calls)),
        "n_base_calls": float(int(base.sum())),
        "n_candidate_calls": float(int(candidate.sum())),
        "candidate_call_frac": float(candidate.mean()) if len(candidate) else 0.0,
        "candidate_raw_frac": candidate_mass / raw_total if raw_total > 0.0 else 0.0,
        "abundance_entropy": entropy(abundance),
        "top1_abundance": float(abundance_sorted[0]) if abundance_sorted.size else 0.0,
        "top3_abundance": float(abundance_sorted[:3].sum()) if abundance_sorted.size else 0.0,
        "raw_mass_sum": raw_total,
        "base_raw_mass": base_mass,
        "s_xny_median": float(s_xny.median()),
        "s_xny_mean": float(s_xny.mean()),
        "u_xny_median": float(u_xny.median()),
        "xny_max_median": float(max_xny.median()),
        "xny_max_mean": float(max_xny.mean()),
        "s_breadth_median": float(s_breadth.median()),
        "u_breadth_median": float(u_breadth.median()),
        "s_depth_median": float(s_depth.median()),
        "u_depth_median": float(u_depth.median()),
        "s_zip_median": float(s_zip.median()),
        "u_zip_median": float(u_zip.median()),
        "probability_mean": float(prob.mean()),
        "probability_min": float(prob.min()),
        "reported_ani_median": float(reported_ani.median()),
        "base_genus_count": float(genus_counts.shape[0]),
        "multi_species_genus_count": float(len(multi_genus)),
        "max_species_per_genus": float(genus_counts.max()) if not genus_counts.empty else 0.0,
        "base_mass_multi_genus_frac": float(base_raw.loc[base_multi.to_numpy()].sum() / base_mass)
        if base_mass > 0.0 and len(base_raw)
        else 0.0,
    }


def build_features() -> pd.DataFrame:
    selected = pd.read_csv(SELECTED_SCORES_TSV, sep="\t")
    selected = selected.loc[selected["method"].eq("wrapper_candidate_abundance_normalized_depth_alpha2")].copy()
    rows: list[dict[str, object]] = []
    for row in selected.to_dict("records"):
        rec = {
            "panel": str(row["panel"]),
            "sample": str(int(float(row["sample"]))),
            "profile": str(row["profile"]),
            "selected_F1": finite(row.get("F1")),
            "selected_L1_pp": finite(row.get("L1_union_pp"))
            if str(row["panel"]) != "cami2_toy_mouse_gut"
            else finite(row.get("L1_truth_only_pp")),
            "selected_Pearson": finite(row.get("Pearson_union"), float("nan"))
            if str(row["panel"]) != "cami2_toy_mouse_gut"
            else finite(row.get("Pearson_truth_only"), float("nan")),
        }
        rec.update(profile_features(Path(str(row["profile"]))))
        rows.append(rec)
    return pd.DataFrame(rows)


def load_deltas(methods: list[str]) -> pd.DataFrame:
    scores = pd.read_csv(SCORES_TSV, sep="\t")
    scores["sample"] = scores["sample"].astype(str)
    keep = scores.loc[scores["method"].isin(methods)].copy()
    return keep[["panel", "sample", "method", "delta_L1_pp", "delta_Pearson"]].copy()


def panel_stats(rows: pd.DataFrame) -> dict[str, float | int]:
    panel_delta = (
        rows.groupby("panel", as_index=False)["guard_delta_L1_pp"]
        .mean()
        .rename(columns={"guard_delta_L1_pp": "panel_delta_L1_pp"})
    )
    return {
        "sample_count": int(rows.shape[0]),
        "switched_samples": int(rows["switch_transform"].sum()),
        "mean_sample_delta_L1_pp": float(rows["guard_delta_L1_pp"].mean()),
        "max_sample_worse_L1_pp": float(rows["guard_delta_L1_pp"].max()),
        "improved_samples": int((rows["guard_delta_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((rows["guard_delta_L1_pp"] > 1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["panel_delta_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["panel_delta_L1_pp"].max()),
        "improved_panels": int((panel_delta["panel_delta_L1_pp"] < -1e-9).sum()),
        "worsened_panels": int((panel_delta["panel_delta_L1_pp"] > 1e-9).sum()),
    }


def evaluate_mask(features: pd.DataFrame, deltas: pd.DataFrame, method: str, mask: pd.Series) -> pd.DataFrame:
    keys = features[["panel", "sample"]].copy()
    keys["method"] = method
    keys["switch_transform"] = mask.to_numpy(dtype=bool)
    method_deltas = deltas.loc[deltas["method"].eq(method)].copy()
    rows = keys.merge(method_deltas, on=["panel", "sample", "method"], how="left")
    rows["guard_delta_L1_pp"] = np.where(rows["switch_transform"], rows["delta_L1_pp"], 0.0)
    rows["guard_delta_Pearson"] = np.where(rows["switch_transform"], rows["delta_Pearson"], 0.0)
    return rows


def feature_columns(features: pd.DataFrame) -> list[str]:
    return [
        col
        for col in features.columns
        if col
        not in {
            "panel",
            "sample",
            "profile",
            "selected_F1",
            "selected_L1_pp",
            "selected_Pearson",
        }
        and pd.api.types.is_numeric_dtype(features[col])
        and features[col].nunique(dropna=True) > 1
    ]


def single_rule_masks(features: pd.DataFrame, method: str) -> list[dict[str, object]]:
    rules: list[dict[str, object]] = []
    for feature in feature_columns(features):
        for threshold in quantile_thresholds(features[feature]):
            values = features[feature].astype(float)
            for op in [">=", "<="]:
                mask = values >= threshold if op == ">=" else values <= threshold
                rules.append(
                    {
                        "method": method,
                        "rule_type": "single",
                        "feature1": feature,
                        "op1": op,
                        "threshold1": threshold,
                        "feature2": "",
                        "op2": "",
                        "threshold2": "",
                        "mask": mask,
                    }
                )
    return rules


def score_rule(rule: dict[str, object], features: pd.DataFrame, deltas: pd.DataFrame) -> dict[str, object]:
    outcomes = evaluate_mask(features, deltas, str(rule["method"]), rule["mask"])
    stats = panel_stats(outcomes)
    return {
        "method": rule["method"],
        "rule_type": rule["rule_type"],
        "feature1": rule["feature1"],
        "op1": rule["op1"],
        "threshold1": rule["threshold1"],
        "feature2": rule["feature2"],
        "op2": rule["op2"],
        "threshold2": rule["threshold2"],
        **stats,
    }


def build_rules(features: pd.DataFrame, deltas: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    scored: list[dict[str, object]] = []
    pair_rules: list[dict[str, object]] = []
    for method in methods:
        single_rules = single_rule_masks(features, method)
        scored_single = [score_rule(rule, features, deltas) for rule in single_rules]
        scored.extend(scored_single)
        single_df = pd.DataFrame(scored_single).sort_values(
            ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"], kind="mergesort"
        )
        pool = []
        for row in single_df.head(MAX_SINGLE_RULES_FOR_PAIRS).itertuples(index=False):
            feature = str(row.feature1)
            op = str(row.op1)
            threshold = finite(row.threshold1)
            values = features[feature].astype(float)
            mask = values >= threshold if op == ">=" else values <= threshold
            pool.append(
                {
                    "method": method,
                    "rule_type": "single",
                    "feature1": feature,
                    "op1": op,
                    "threshold1": threshold,
                    "feature2": "",
                    "op2": "",
                    "threshold2": "",
                    "mask": mask,
                }
            )
        for left, right in itertools.combinations(pool, 2):
            if left["feature1"] == right["feature1"] and left["op1"] == right["op1"]:
                continue
            pair_rules.append(
                {
                    "method": method,
                    "rule_type": "two_feature_and",
                    "feature1": left["feature1"],
                    "op1": left["op1"],
                    "threshold1": left["threshold1"],
                    "feature2": right["feature1"],
                    "op2": right["op1"],
                    "threshold2": right["threshold1"],
                    "mask": left["mask"] & right["mask"],
                }
            )
    scored.extend([score_rule(rule, features, deltas) for rule in pair_rules])
    return pd.DataFrame(scored).sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp", "method", "rule_type"],
        kind="mergesort",
    )


def mask_from_rule(features: pd.DataFrame, rule: pd.Series) -> pd.Series:
    values1 = features[str(rule["feature1"])].astype(float)
    mask1 = values1 >= finite(rule["threshold1"]) if str(rule["op1"]) == ">=" else values1 <= finite(rule["threshold1"])
    if str(rule.get("rule_type", "")) != "two_feature_and" or not str(rule.get("feature2", "")):
        return mask1
    values2 = features[str(rule["feature2"])].astype(float)
    mask2 = values2 >= finite(rule["threshold2"]) if str(rule["op2"]) == ">=" else values2 <= finite(rule["threshold2"])
    return mask1 & mask2


def leave_one_panel_out(rules: pd.DataFrame, features: pd.DataFrame, deltas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    nontrivial = rules.loc[
        rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])
    ].copy()
    candidate_rules = nontrivial.head(MAX_LOPO_RULES)
    for holdout in sorted(features["panel"].astype(str).unique()):
        candidates: list[dict[str, object]] = []
        for rule in candidate_rules.itertuples(index=False):
            rule_series = pd.Series(rule._asdict())
            outcomes = evaluate_mask(
                features,
                deltas,
                str(rule_series["method"]),
                mask_from_rule(features, rule_series),
            )
            train = outcomes.loc[~outcomes["panel"].eq(holdout)].copy()
            test = outcomes.loc[outcomes["panel"].eq(holdout)].copy()
            train_stats = panel_stats(train)
            test_stats = panel_stats(test)
            train_sample_safe = (
                int(train_stats["worsened_samples"]) == 0
                and int(train_stats["improved_samples"]) > 0
            )
            train_panel_safe = (
                int(train_stats["worsened_panels"]) == 0
                and int(train_stats["improved_panels"]) > 0
            )
            candidates.append(
                {
                    "holdout_panel": holdout,
                    "method": rule_series["method"],
                    "rule_type": rule_series["rule_type"],
                    "feature1": rule_series["feature1"],
                    "op1": rule_series["op1"],
                    "threshold1": rule_series["threshold1"],
                    "feature2": rule_series.get("feature2", ""),
                    "op2": rule_series.get("op2", ""),
                    "threshold2": rule_series.get("threshold2", ""),
                    "train_sample_safe": train_sample_safe,
                    "train_panel_safe": train_panel_safe,
                    "train_mean_sample_delta_L1_pp": train_stats["mean_sample_delta_L1_pp"],
                    "train_max_sample_worse_L1_pp": train_stats["max_sample_worse_L1_pp"],
                    "train_improved_samples": train_stats["improved_samples"],
                    "train_worsened_samples": train_stats["worsened_samples"],
                    "holdout_mean_sample_delta_L1_pp": test_stats["mean_sample_delta_L1_pp"],
                    "holdout_max_sample_worse_L1_pp": test_stats["max_sample_worse_L1_pp"],
                    "holdout_improved_samples": test_stats["improved_samples"],
                    "holdout_worsened_samples": test_stats["worsened_samples"],
                    "holdout_switched_samples": test_stats["switched_samples"],
                    "holdout_sample_count": test_stats["sample_count"],
                }
            )
        cand_df = pd.DataFrame(candidates)
        strict_sample = cand_df.loc[cand_df["train_sample_safe"]].copy()
        strict_panel = cand_df.loc[cand_df["train_panel_safe"]].copy()
        if not strict_sample.empty:
            pool = strict_sample
            selection_pool = "strict_train_sample_safe"
        elif not strict_panel.empty:
            pool = strict_panel
            selection_pool = "strict_train_panel_safe"
        else:
            pool = cand_df
            selection_pool = "best_train_mean"
        selected = pool.sort_values(
            ["train_mean_sample_delta_L1_pp", "train_max_sample_worse_L1_pp"],
            kind="mergesort",
        ).iloc[0].to_dict()
        selected["selection_pool"] = selection_pool
        rows.append(selected)
    return pd.DataFrame(rows)


def rule_text(row: pd.Series) -> str:
    first = f"{row['feature1']}{row['op1']}{finite(row['threshold1']):.6g}"
    if str(row.get("rule_type", "")) == "two_feature_and" and str(row.get("feature2", "")):
        second = f"{row['feature2']}{row['op2']}{finite(row['threshold2']):.6g}"
        return f"{first} AND {second}"
    return first


def audit_rows(rules: pd.DataFrame, lopo: pd.DataFrame, methods: list[str]) -> list[dict[str, object]]:
    nontrivial = rules.loc[rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])]
    sample_safe = nontrivial.loc[
        nontrivial["worsened_samples"].eq(0) & nontrivial["improved_samples"].gt(0)
    ]
    panel_safe = nontrivial.loc[
        nontrivial["worsened_panels"].eq(0) & nontrivial["improved_panels"].gt(0)
    ]
    best_mean = nontrivial.sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"],
        kind="mergesort",
    ).iloc[0]
    best_safe = (
        sample_safe.sort_values(["mean_sample_delta_L1_pp", "switched_samples"], kind="mergesort").iloc[0]
        if not sample_safe.empty
        else pd.Series(dtype=object)
    )
    lopo_mean_worse = int((lopo["holdout_mean_sample_delta_L1_pp"] > 1e-9).sum())
    lopo_sample_regress = int((lopo["holdout_worsened_samples"] > 0).sum())
    lopo_mean = float(lopo["holdout_mean_sample_delta_L1_pp"].mean())
    lopo_max_worse = float(lopo["holdout_max_sample_worse_L1_pp"].max())
    if not sample_safe.empty and lopo_mean_worse == 0 and lopo_sample_regress == 0:
        decision = "candidate_guard_needs_wrapper_and_independent_holdout"
    elif not sample_safe.empty:
        decision = "sample_safe_in_panel_but_fails_lopo"
    elif not panel_safe.empty:
        decision = "panel_safe_only_not_default"
    else:
        decision = "reject_selected_mass_transform_guard"

    return [
        {
            "metric": "candidate_methods",
            "value": ",".join(methods),
            "evidence": "top selected-call transforms by cached overall ranking",
            "decision": "truth_used_only_for_method_scoring",
        },
        {
            "metric": "tested_guard_rules",
            "value": int(nontrivial.shape[0]),
            "evidence": f"method plus output-derived threshold guards;top_saved={MAX_WRITTEN_RULES}",
            "decision": "truth_used_only_for_scoring",
        },
        {
            "metric": "sample_safe_guards",
            "value": int(sample_safe.shape[0]),
            "evidence": f"full in-memory scoring;top_saved={MAX_WRITTEN_RULES}",
            "decision": "strict_sample_gate",
        },
        {
            "metric": "panel_safe_guards",
            "value": int(panel_safe.shape[0]),
            "evidence": f"full in-memory scoring;top_saved={MAX_WRITTEN_RULES}",
            "decision": "panel_gate",
        },
        {
            "metric": "best_mean_guard",
            "value": (
                f"{best_mean['method']};{best_mean['rule_type']};{rule_text(best_mean)};"
                f"mean_delta={best_mean['mean_sample_delta_L1_pp']:.6f};"
                f"worsened_samples={int(best_mean['worsened_samples'])};"
                f"max_worse={best_mean['max_sample_worse_L1_pp']:.6f};"
                f"switched={int(best_mean['switched_samples'])}"
            ),
            "evidence": f"selected_call_mass_transform_guard_rules.tsv;top_saved={MAX_WRITTEN_RULES}",
            "decision": "tradeoff",
        },
        {
            "metric": "best_sample_safe_guard",
            "value": (
                "none"
                if best_safe.empty
                else f"{best_safe['method']};{best_safe['rule_type']};{rule_text(best_safe)};"
                f"mean_delta={best_safe['mean_sample_delta_L1_pp']:.6f};"
                f"switched={int(best_safe['switched_samples'])}"
            ),
            "evidence": f"selected_call_mass_transform_guard_rules.tsv;top_saved={MAX_WRITTEN_RULES}",
            "decision": "candidate" if not best_safe.empty else "none",
        },
        {
            "metric": "lopo_guard_result",
            "value": (
                f"mean_delta={lopo_mean:.6f};holdouts_with_mean_regression={lopo_mean_worse};"
                f"holdouts_with_sample_regression={lopo_sample_regress};"
                f"max_sample_worse={lopo_max_worse:.6f}"
            ),
            "evidence": f"selected_call_mass_transform_guard_lopo.tsv;top_rules={MAX_LOPO_RULES}",
            "decision": "leave_one_panel_out_validation",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": "guard rules plus leave-one-panel-out",
            "decision": decision,
        },
    ]


def write_markdown(audit: list[dict[str, object]], rules: pd.DataFrame, lopo: pd.DataFrame) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Selected Call-Set Mass Transform Guard Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached selected-candidate profiles and",
        "the selected-call mass-transform sweep. It tests whether output-derived",
        "guards can apply a simple transform only on samples where it is stable.",
        "Cached score columns (`selected_F1`, `selected_L1_pp`, and",
        "`selected_Pearson`) are excluded from guard rule features.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "tested_guard_rules",
        "sample_safe_guards",
        "panel_safe_guards",
        "best_mean_guard",
        "best_sample_safe_guard",
        "lopo_guard_result",
        "promotion_decision",
    ]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Top Rules",
            "",
            "| Method | Rule | Mean sample delta pp | Worsened samples | Max worse pp | Switched |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in rules.head(8).to_dict("records"):
        series = pd.Series(row)
        lines.append(
            "| {method} | {rule} | {mean:.6f} | {worse} | {max_worse:.6f} | {switched} |".format(
                method=row["method"],
                rule=rule_text(series),
                mean=float(row["mean_sample_delta_L1_pp"]),
                worse=int(row["worsened_samples"]),
                max_worse=float(row["max_sample_worse_L1_pp"]),
                switched=int(row["switched_samples"]),
            )
        )

    lines.extend(
        [
            "",
            "## LOPO",
            "",
            "| Holdout | Method | Selection | Mean delta pp | Worsened samples | Switched |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for row in lopo.to_dict("records"):
        lines.append(
            "| {holdout} | {method} | {pool} | {mean:.6f} | {worse} | {switched} |".format(
                holdout=row["holdout_panel"],
                method=row["method"],
                pool=row["selection_pool"],
                mean=float(row["holdout_mean_sample_delta_L1_pp"]),
                worse=int(row["holdout_worsened_samples"]),
                switched=int(row["holdout_switched_samples"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- This is a diagnostic guard search, not a promoted allocator.",
            "- A guard is promotable only if it remains sample-safe in",
            "  leave-one-panel-out validation of top-ranked rules and then survives",
            "  independent holdouts.",
            "",
            "## Outputs",
            "",
            "- `results/selected_call_mass_transform_guard_features.tsv`",
            f"- `results/selected_call_mass_transform_guard_rules.tsv` top {MAX_WRITTEN_RULES} scored rules",
            "- `results/selected_call_mass_transform_guard_lopo.tsv`",
            "- `results/selected_call_mass_transform_guard_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    methods = candidate_methods()
    features = build_features()
    deltas = load_deltas(methods)
    rules = build_rules(features, deltas, methods)
    lopo = leave_one_panel_out(rules, features, deltas)
    audit = audit_rows(rules, lopo, methods)

    features.sort_values(["panel", "sample"], kind="mergesort").to_csv(FEATURES_TSV, sep="\t", index=False)
    rules.head(MAX_WRITTEN_RULES).to_csv(RULES_TSV, sep="\t", index=False)
    lopo.to_csv(LOPO_TSV, sep="\t", index=False)
    write_tsv(AUDIT_TSV, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(audit, rules, lopo)
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
