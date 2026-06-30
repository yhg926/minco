#!/usr/bin/env python3
"""Audit output-only guards for selected-default genus-XnY blend.

The unguarded alpha 0.25 abundance blend improves panel means but regresses
some individual samples. This script asks whether a simple guard based only on
MinCO output and the deterministic blend perturbation can keep the gains
without sample regressions.
"""

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

DELTA_TSV = RESULTS / "candidate_preset_genus_xny_blend_vs_baseline.tsv"
METADATA_TSV = RESULTS / "candidate_preset_genus_xny_blend_metadata.tsv"

FEATURES_TSV = RESULTS / "candidate_preset_genus_xny_blend_guard_features.tsv"
RULES_TSV = RESULTS / "candidate_preset_genus_xny_blend_guard_rules.tsv"
LOPO_TSV = RESULTS / "candidate_preset_genus_xny_blend_guard_lopo.tsv"
AUDIT_TSV = RESULTS / "candidate_preset_genus_xny_blend_guard_audit.tsv"

MAX_SINGLE_RULES_FOR_PAIRS = 80


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


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


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
    thresholds: list[float] = []
    for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        value = float(clean.quantile(q))
        if math.isfinite(value):
            thresholds.append(value)
    return sorted(set(thresholds))


def genus_from_name(name: object) -> str:
    text = str(name or "").strip()
    if text.startswith("s__"):
        text = text[3:]
    parts = text.split()
    return parts[0] if parts else text


def profile_features(src: Path, adjusted: Path) -> dict[str, float]:
    raw = pd.read_csv(src, sep="\t", low_memory=False)
    adj = pd.read_csv(
        adjusted,
        sep="\t",
        usecols=lambda col: col in {"calibrated_abundance_raw"},
        low_memory=False,
    )
    called = bool_series(raw.get("calibrated_call", pd.Series(False, index=raw.index)))
    rescue = bool_series(raw.get("candidate_rescue_added", pd.Series(False, index=raw.index)))
    surface = bool_series(raw.get("candidate_surface_added", pd.Series(False, index=raw.index)))
    candidate = rescue | surface
    base = called & ~candidate
    calls = raw.loc[called].copy()
    base_calls = raw.loc[base].copy()
    if calls.empty:
        return {"n_calls": 0.0}

    abundance = numeric(calls, "calibrated_abundance")
    abundance_sorted = np.sort(np.maximum(abundance.to_numpy(dtype=float), 0.0))[::-1]
    raw_mass = numeric(calls, "calibrated_abundance_raw")
    base_raw = numeric(raw, "calibrated_abundance_raw").loc[base]
    adj_raw = pd.to_numeric(adj["calibrated_abundance_raw"], errors="coerce").fillna(0.0)
    raw_delta = (adj_raw - numeric(raw, "calibrated_abundance_raw")).abs()
    perturb = raw_delta.loc[base]
    s_xny = numeric(calls, "s_XnY_ctx_max")
    u_xny = numeric(calls, "u_XnY_ctx_max")
    s_breadth = numeric(calls, "s_Ref_breadth_max")
    u_breadth = numeric(calls, "u_Ref_breadth_max")
    s_depth = numeric(calls, "s_Ref_mean_depth_max")
    u_depth = numeric(calls, "u_Ref_mean_depth_max")
    s_zip = numeric(calls, "s_Ref_zip_af_max")
    u_zip = numeric(calls, "u_Ref_zip_af_max")
    prob = numeric(calls, "calibrated_probability", 1.0)

    genus_series = base_calls.get("species_name", pd.Series("", index=base_calls.index)).map(genus_from_name)
    genus_counts = genus_series.value_counts()
    multispecies_genus = set(genus_counts.loc[genus_counts > 1].index.astype(str))
    base_mass_by_row = numeric(base_calls, "calibrated_abundance_raw")
    if len(base_calls):
        in_multi = genus_series.astype(str).isin(multispecies_genus).to_numpy()
        multi_mass = float(base_mass_by_row.loc[in_multi].sum())
    else:
        multi_mass = 0.0
    base_mass = float(base_raw.sum())
    perturb_sum = float(perturb.sum())

    return {
        "n_calls": float(len(calls)),
        "n_base_calls": float(int(base.sum())),
        "n_candidate_calls": float(int(candidate.loc[called].sum())),
        "candidate_call_frac": float(candidate.loc[called].mean()) if called.any() else 0.0,
        "abundance_entropy": entropy(abundance),
        "top1_abundance": float(abundance_sorted[0]) if abundance_sorted.size else 0.0,
        "top3_abundance": float(abundance_sorted[:3].sum()) if abundance_sorted.size else 0.0,
        "raw_mass_sum": float(raw_mass.sum()),
        "base_raw_mass": base_mass,
        "candidate_raw_mass": float(numeric(raw, "calibrated_abundance_raw").loc[called & candidate].sum()),
        "candidate_raw_frac": float(numeric(raw, "calibrated_abundance_raw").loc[called & candidate].sum() / raw_mass.sum())
        if float(raw_mass.sum()) > 0.0
        else 0.0,
        "s_xny_median": float(s_xny.median()),
        "s_xny_mean": float(s_xny.mean()),
        "u_xny_median": float(u_xny.median()),
        "xny_max_median": float(pd.concat([s_xny, u_xny], axis=1).max(axis=1).median()),
        "s_breadth_median": float(s_breadth.median()),
        "u_breadth_median": float(u_breadth.median()),
        "s_depth_median": float(s_depth.median()),
        "u_depth_median": float(u_depth.median()),
        "s_zip_median": float(s_zip.median()),
        "u_zip_median": float(u_zip.median()),
        "probability_mean": float(prob.mean()),
        "probability_min": float(prob.min()),
        "base_genus_count": float(genus_counts.shape[0]),
        "multi_species_genus_count": float(len(multispecies_genus)),
        "max_species_per_genus": float(genus_counts.max()) if not genus_counts.empty else 0.0,
        "base_mass_multi_genus_frac": multi_mass / base_mass if base_mass > 0.0 else 0.0,
        "blend_perturb_raw_sum": perturb_sum,
        "blend_perturb_base_frac": perturb_sum / base_mass if base_mass > 0.0 else 0.0,
        "blend_perturb_row_frac": float((perturb > 1e-12).mean()) if len(perturb) else 0.0,
        "blend_perturb_max_frac": float(perturb.max() / base_mass) if len(perturb) and base_mass > 0.0 else 0.0,
    }


def build_features() -> pd.DataFrame:
    metadata = pd.read_csv(METADATA_TSV, sep="\t")
    rows: list[dict[str, object]] = []
    for row in metadata.itertuples(index=False):
        panel = str(row.panel)
        sample = str(row.sample)
        source = Path(str(row.source_profile))
        adjusted = Path(str(row.adjusted_profile))
        rec = {
            "panel": panel,
            "sample": sample,
            "source_profile": str(source),
            "adjusted_profile": str(adjusted),
            "called_rows": finite(getattr(row, "called_rows", 0)),
            "candidate_added_rows": finite(getattr(row, "candidate_added_rows", 0)),
            "eligible_base_rows": finite(getattr(row, "eligible_base_rows", 0)),
            "adjusted_base_rows": finite(getattr(row, "adjusted_base_rows", 0)),
            "adjusted_base_row_frac": finite(getattr(row, "adjusted_base_rows", 0))
            / max(1.0, finite(getattr(row, "eligible_base_rows", 0))),
            "metadata_candidate_raw_mass_before": finite(
                getattr(row, "candidate_raw_mass_before", 0.0)
            ),
        }
        rec.update(profile_features(source, adjusted))
        rows.append(rec)
    return pd.DataFrame(rows)


def load_deltas() -> pd.DataFrame:
    delta = pd.read_csv(DELTA_TSV, sep="\t")
    delta["sample"] = delta["sample"].astype(str)
    return delta[["panel", "sample", "official_L1_delta_pp", "official_Pearson_delta"]].copy()


def panel_stats(rows: pd.DataFrame) -> dict[str, float | int]:
    panel_delta = (
        rows.groupby("panel", as_index=False)["guard_delta_L1_pp"]
        .mean()
        .rename(columns={"guard_delta_L1_pp": "panel_delta_L1_pp"})
    )
    return {
        "sample_count": int(rows.shape[0]),
        "switched_samples": int(rows["switch_blend"].sum()),
        "mean_sample_delta_L1_pp": float(rows["guard_delta_L1_pp"].mean()),
        "max_sample_worse_L1_pp": float(rows["guard_delta_L1_pp"].max()),
        "improved_samples": int((rows["guard_delta_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((rows["guard_delta_L1_pp"] > 1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["panel_delta_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["panel_delta_L1_pp"].max()),
        "improved_panels": int((panel_delta["panel_delta_L1_pp"] < -1e-9).sum()),
        "worsened_panels": int((panel_delta["panel_delta_L1_pp"] > 1e-9).sum()),
    }


def evaluate_mask(features: pd.DataFrame, deltas: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    keys = features[["panel", "sample"]].copy()
    keys["switch_blend"] = mask.to_numpy(dtype=bool)
    rows = keys.merge(deltas, on=["panel", "sample"], how="left")
    rows["guard_delta_L1_pp"] = np.where(rows["switch_blend"], rows["official_L1_delta_pp"], 0.0)
    rows["guard_delta_Pearson"] = np.where(rows["switch_blend"], rows["official_Pearson_delta"], 0.0)
    return rows


def single_rule_masks(features: pd.DataFrame) -> list[dict[str, object]]:
    feature_cols = [
        col
        for col in features.columns
        if col not in {"panel", "sample", "source_profile", "adjusted_profile"}
        and pd.api.types.is_numeric_dtype(features[col])
        and features[col].nunique(dropna=True) > 1
    ]
    rules: list[dict[str, object]] = []
    for feature in feature_cols:
        for threshold in quantile_thresholds(features[feature]):
            for op in [">=", "<="]:
                mask = features[feature].astype(float) >= threshold if op == ">=" else features[feature].astype(float) <= threshold
                rules.append(
                    {
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
    outcomes = evaluate_mask(features, deltas, rule["mask"])
    stats = panel_stats(outcomes)
    return {
        "rule_type": rule["rule_type"],
        "feature1": rule["feature1"],
        "op1": rule["op1"],
        "threshold1": rule["threshold1"],
        "feature2": rule["feature2"],
        "op2": rule["op2"],
        "threshold2": rule["threshold2"],
        **stats,
    }


def build_rules(features: pd.DataFrame, deltas: pd.DataFrame) -> pd.DataFrame:
    single_rules = single_rule_masks(features)
    scored_single = [score_rule(rule, features, deltas) for rule in single_rules]
    single_df = pd.DataFrame(scored_single).sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"],
        kind="mergesort",
    )
    pool = []
    for row in single_df.head(MAX_SINGLE_RULES_FOR_PAIRS).itertuples(index=False):
        feature = str(row.feature1)
        op = str(row.op1)
        threshold = finite(row.threshold1)
        mask = features[feature].astype(float) >= threshold if op == ">=" else features[feature].astype(float) <= threshold
        pool.append(
            {
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
    pair_rules: list[dict[str, object]] = []
    for left, right in itertools.combinations(pool, 2):
        if (
            left["feature1"] == right["feature1"]
            and left["op1"] == right["op1"]
            and left["threshold1"] == right["threshold1"]
        ):
            continue
        pair_rules.append(
            {
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
    scored_pairs = [score_rule(rule, features, deltas) for rule in pair_rules]
    return pd.concat([single_df, pd.DataFrame(scored_pairs)], ignore_index=True, sort=False).sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp", "rule_type"],
        kind="mergesort",
    )


def mask_from_rule(features: pd.DataFrame, rule: pd.Series) -> pd.Series:
    mask1 = (
        features[str(rule["feature1"])].astype(float) >= finite(rule["threshold1"])
        if str(rule["op1"]) == ">="
        else features[str(rule["feature1"])].astype(float) <= finite(rule["threshold1"])
    )
    if str(rule.get("rule_type", "")) != "two_feature_and" or not str(rule.get("feature2", "")):
        return mask1
    mask2 = (
        features[str(rule["feature2"])].astype(float) >= finite(rule["threshold2"])
        if str(rule["op2"]) == ">="
        else features[str(rule["feature2"])].astype(float) <= finite(rule["threshold2"])
    )
    return mask1 & mask2


def leave_one_panel_out(rules: pd.DataFrame, features: pd.DataFrame, deltas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    panels = sorted(features["panel"].astype(str).unique())
    for holdout in panels:
        candidates: list[dict[str, object]] = []
        for rule in rules.itertuples(index=False):
            rule_series = pd.Series(rule._asdict())
            outcomes = evaluate_mask(features, deltas, mask_from_rule(features, rule_series))
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


def audit_rows(rules: pd.DataFrame, lopo: pd.DataFrame) -> list[dict[str, object]]:
    nontrivial = rules.loc[rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])]
    sample_safe = nontrivial.loc[
        nontrivial["worsened_samples"].eq(0) & nontrivial["improved_samples"].gt(0)
    ]
    panel_safe = nontrivial.loc[
        nontrivial["worsened_panels"].eq(0) & nontrivial["improved_panels"].gt(0)
    ]
    best_sample = nontrivial.sort_values(
        ["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"],
        kind="mergesort",
    ).iloc[0]
    best_safe = (
        sample_safe.sort_values(["mean_sample_delta_L1_pp", "switched_samples"], kind="mergesort").iloc[0]
        if not sample_safe.empty
        else pd.Series(dtype=object)
    )
    lopo_worse_holdouts = int((lopo["holdout_mean_sample_delta_L1_pp"] > 1e-9).sum())
    lopo_sample_regress_holdouts = int((lopo["holdout_worsened_samples"] > 0).sum())
    lopo_mean = float(lopo["holdout_mean_sample_delta_L1_pp"].mean())
    lopo_max_worse = float(lopo["holdout_max_sample_worse_L1_pp"].max())
    if not sample_safe.empty and lopo_worse_holdouts == 0 and lopo_sample_regress_holdouts == 0:
        decision = "candidate_guard_needs_wrapper_and_independent_holdout"
    elif not sample_safe.empty:
        decision = "sample_safe_in_panel_but_fails_lopo"
    elif not panel_safe.empty:
        decision = "panel_safe_only_not_default"
    else:
        decision = "reject_simple_blend_guard"

    rows = [
        {
            "metric": "tested_guard_rules",
            "value": int(nontrivial.shape[0]),
            "evidence": "single-feature plus top-pair output-derived threshold guards",
            "decision": "truth_used_only_for_scoring",
        },
        {
            "metric": "sample_safe_guards",
            "value": int(sample_safe.shape[0]),
            "evidence": "candidate_preset_genus_xny_blend_guard_rules.tsv",
            "decision": "strict_sample_gate",
        },
        {
            "metric": "panel_safe_guards",
            "value": int(panel_safe.shape[0]),
            "evidence": "candidate_preset_genus_xny_blend_guard_rules.tsv",
            "decision": "panel_gate",
        },
        {
            "metric": "best_mean_guard",
            "value": (
                f"{best_sample['rule_type']};{rule_text(best_sample)};"
                f"mean_delta={best_sample['mean_sample_delta_L1_pp']:.6f};"
                f"worsened_samples={int(best_sample['worsened_samples'])};"
                f"max_worse={best_sample['max_sample_worse_L1_pp']:.6f};"
                f"switched={int(best_sample['switched_samples'])}"
            ),
            "evidence": "candidate_preset_genus_xny_blend_guard_rules.tsv",
            "decision": "tradeoff",
        },
        {
            "metric": "best_sample_safe_guard",
            "value": (
                "none"
                if best_safe.empty
                else f"{best_safe['rule_type']};{rule_text(best_safe)};"
                f"mean_delta={best_safe['mean_sample_delta_L1_pp']:.6f};"
                f"switched={int(best_safe['switched_samples'])}"
            ),
            "evidence": "candidate_preset_genus_xny_blend_guard_rules.tsv",
            "decision": "candidate" if not best_safe.empty else "none",
        },
        {
            "metric": "lopo_guard_result",
            "value": (
                f"mean_delta={lopo_mean:.6f};holdouts_with_mean_regression={lopo_worse_holdouts};"
                f"holdouts_with_sample_regression={lopo_sample_regress_holdouts};"
                f"max_sample_worse={lopo_max_worse:.6f}"
            ),
            "evidence": "candidate_preset_genus_xny_blend_guard_lopo.tsv",
            "decision": "leave_one_panel_out_validation",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": "guard rules plus leave-one-panel-out",
            "decision": decision,
        },
    ]
    return rows


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    features = build_features()
    deltas = load_deltas()
    rules = build_rules(features, deltas)
    lopo = leave_one_panel_out(rules, features, deltas)
    audit = audit_rows(rules, lopo)

    features.sort_values(["panel", "sample"], kind="mergesort").to_csv(
        FEATURES_TSV, sep="\t", index=False
    )
    rules.to_csv(RULES_TSV, sep="\t", index=False)
    lopo.to_csv(LOPO_TSV, sep="\t", index=False)
    write_tsv(AUDIT_TSV, audit, ["metric", "value", "evidence", "decision"])

    print(pd.DataFrame(audit).to_string(index=False))
    print("\nTOP RULES")
    print(rules.head(20).to_string(index=False))
    print("\nLOPO")
    print(lopo.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
