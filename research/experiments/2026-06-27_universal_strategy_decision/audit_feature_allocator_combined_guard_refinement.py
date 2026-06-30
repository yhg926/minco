#!/usr/bin/env python3
"""Search output-only guard refinements for the guarded feature allocator.

The current `guarded-genus-hit-breadth-a002` allocator is useful but not
promoted because two external exact-split diagnostic samples have small L1
regressions. This cached audit asks whether a simple output-derived guard can
avoid those regressions while preserving most selected-panel gains.
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import audit_selected_call_mass_transform_guard as guard


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SELECTED_FEATURES = RESULTS / "selected_call_feature_allocator_guard_features.tsv"
SELECTED_SCORES = RESULTS / "feature_allocator_wrapper_parity_scores.tsv"
EXTERNAL_SCORES = RESULTS / "feature_allocator_external_exactsplit_scores.tsv"

OUT_FEATURES = RESULTS / "feature_allocator_combined_guard_features.tsv"
OUT_RULES = RESULTS / "feature_allocator_combined_guard_rules.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_combined_guard_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_COMBINED_GUARD_REFINEMENT.md"

METHOD = "guarded_feature_allocator_a002_refined_guard"
SELECTED_METHOD = "wrapper_guarded_genus_hit_breadth_a002"
EXTERNAL_CURRENT = "current_exactsplit_abundance"
EXTERNAL_SWITCH = "guarded_feature_allocator_exactsplit"

MAX_SINGLE_RULES_FOR_PAIRS = 40
MAX_WRITTEN_RULES = 1000
CACHED_GAIN_FRACTION_MIN = 0.80


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    return guard.finite(value, default)


def selected_panel(panel: object) -> bool:
    return str(panel) in {
        "cami2_toy_mouse_gut",
        "cami3_toy_human_gut_gtdb_source_readmap",
        "hmp_airskin_gtdb_source_abundance",
        "hmp_gastrooral_gtdb_source_abundance",
    }


def build_features() -> pd.DataFrame:
    selected = pd.read_csv(SELECTED_FEATURES, sep="\t")
    selected = selected.copy()
    selected["sample"] = selected["sample"].astype(str)
    selected["evidence_group"] = "selected_cached_wrapper"

    external_score_rows = pd.DataFrame(read_tsv(EXTERNAL_SCORES))
    current = external_score_rows.loc[external_score_rows["method"].eq(EXTERNAL_CURRENT)].copy()
    external_rows: list[dict[str, object]] = []
    for row in current.to_dict("records"):
        rec = {
            "panel": str(row["dataset"]),
            "sample": str(int(float(row["sample"]))),
            "profile": str(row["profile_path"]),
            "selected_F1": finite(row.get("F1")),
            "selected_L1_pp": finite(row.get("bacteria_scope_l1")),
            "selected_Pearson": finite(row.get("tp_abundance_pearson"), float("nan")),
            "evidence_group": "external_exactsplit_diagnostic",
        }
        rec.update(guard.profile_features(Path(str(row["profile_path"]))))
        external_rows.append(rec)

    combined = pd.concat([selected, pd.DataFrame(external_rows)], ignore_index=True, sort=False)
    combined["sample"] = combined["sample"].astype(str)
    return combined


def build_deltas() -> pd.DataFrame:
    selected = pd.read_csv(SELECTED_SCORES, sep="\t")
    selected = selected.loc[selected["method"].eq(SELECTED_METHOD)].copy()
    selected_rows = selected[["panel", "sample", "delta_L1_pp", "delta_Pearson"]].copy()
    selected_rows["panel"] = selected_rows["panel"].astype(str)
    selected_rows["sample"] = selected_rows["sample"].astype(str)
    selected_rows["method"] = METHOD

    external = pd.read_csv(EXTERNAL_SCORES, sep="\t")
    current = external.loc[external["method"].eq(EXTERNAL_CURRENT)].copy()
    switched = external.loc[external["method"].eq(EXTERNAL_SWITCH)].copy()
    current = current.set_index(["dataset", "sample"])
    rows: list[dict[str, object]] = []
    for row in switched.to_dict("records"):
        key = (row["dataset"], row["sample"])
        base = current.loc[key]
        rows.append(
            {
                "panel": str(row["dataset"]),
                "sample": str(int(float(row["sample"]))),
                "method": METHOD,
                "delta_L1_pp": finite(row.get("bacteria_scope_l1"))
                - finite(base.get("bacteria_scope_l1")),
                "delta_Pearson": finite(row.get("tp_abundance_pearson"), float("nan"))
                - finite(base.get("tp_abundance_pearson"), float("nan")),
            }
        )
    external_rows = pd.DataFrame(rows)
    return pd.concat([selected_rows, external_rows], ignore_index=True, sort=False)


def evaluate_mask(features: pd.DataFrame, deltas: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    keys = features[["panel", "sample", "evidence_group"]].copy()
    keys["rule_selected"] = mask.to_numpy(dtype=bool)
    keys["method"] = METHOD
    rows = keys.merge(deltas, on=["panel", "sample", "method"], how="left")
    rows["delta_L1_pp"] = pd.to_numeric(rows["delta_L1_pp"], errors="coerce").fillna(0.0)
    rows["delta_Pearson"] = pd.to_numeric(rows["delta_Pearson"], errors="coerce").fillna(0.0)
    rows["guard_delta_L1_pp"] = np.where(rows["rule_selected"], rows["delta_L1_pp"], 0.0)
    rows["guard_delta_Pearson"] = np.where(rows["rule_selected"], rows["delta_Pearson"], 0.0)
    return rows


def panel_delta_stats(outcomes: pd.DataFrame) -> dict[str, float | int]:
    panel_delta = outcomes.groupby("panel", as_index=False)["guard_delta_L1_pp"].mean()
    selected = outcomes.loc[outcomes["panel"].map(selected_panel)].copy()
    external = outcomes.loc[~outcomes["panel"].map(selected_panel)].copy()
    return {
        "sample_count": int(outcomes.shape[0]),
        "selected_by_rule": int(outcomes["rule_selected"].sum()),
        "mean_sample_delta_L1_pp": float(outcomes["guard_delta_L1_pp"].mean()),
        "max_sample_worse_L1_pp": float(outcomes["guard_delta_L1_pp"].max()),
        "improved_samples": int((outcomes["guard_delta_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((outcomes["guard_delta_L1_pp"] > 1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["guard_delta_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["guard_delta_L1_pp"].max()),
        "improved_panels": int((panel_delta["guard_delta_L1_pp"] < -1e-9).sum()),
        "worsened_panels": int((panel_delta["guard_delta_L1_pp"] > 1e-9).sum()),
        "selected_cached_mean_delta_L1_pp": float(selected["guard_delta_L1_pp"].mean()),
        "selected_cached_improved_samples": int((selected["guard_delta_L1_pp"] < -1e-9).sum()),
        "selected_cached_worsened_samples": int((selected["guard_delta_L1_pp"] > 1e-9).sum()),
        "external_mean_delta_L1_pp": float(external["guard_delta_L1_pp"].mean()),
        "external_improved_samples": int((external["guard_delta_L1_pp"] < -1e-9).sum()),
        "external_worsened_samples": int((external["guard_delta_L1_pp"] > 1e-9).sum()),
    }


def feature_columns(features: pd.DataFrame) -> list[str]:
    return guard.feature_columns(features)


def base_rule_masks(features: pd.DataFrame) -> list[dict[str, object]]:
    rules: list[dict[str, object]] = []
    for feature in feature_columns(features):
        for threshold in guard.quantile_thresholds(features[feature]):
            values = features[feature].astype(float)
            for op in [">=", "<="]:
                rules.append(
                    {
                        "rule_type": "single",
                        "feature1": feature,
                        "op1": op,
                        "threshold1": threshold,
                        "feature2": "",
                        "op2": "",
                        "threshold2": "",
                        "mask": values >= threshold if op == ">=" else values <= threshold,
                    }
                )
    return rules


def score_rule(
    rule: dict[str, object],
    features: pd.DataFrame,
    deltas: pd.DataFrame,
    current_selected_mean: float,
) -> dict[str, object]:
    outcomes = evaluate_mask(features, deltas, rule["mask"])
    stats = panel_delta_stats(outcomes)
    cached_mean = float(stats["selected_cached_mean_delta_L1_pp"])
    preserved = cached_mean / current_selected_mean if current_selected_mean < 0.0 else 0.0
    return {
        "method": METHOD,
        "rule_type": rule["rule_type"],
        "feature1": rule["feature1"],
        "op1": rule["op1"],
        "threshold1": rule["threshold1"],
        "feature2": rule.get("feature2", ""),
        "op2": rule.get("op2", ""),
        "threshold2": rule.get("threshold2", ""),
        **stats,
        "cached_gain_preserved_frac": preserved,
        "strict_combined_sample_safe": (
            int(stats["worsened_samples"]) == 0
            and int(stats["improved_samples"]) > 0
            and preserved >= CACHED_GAIN_FRACTION_MIN
        ),
    }


def build_rules(features: pd.DataFrame, deltas: pd.DataFrame) -> pd.DataFrame:
    current = evaluate_mask(features, deltas, pd.Series(True, index=features.index))
    current_stats = panel_delta_stats(current)
    current_selected_mean = float(current_stats["selected_cached_mean_delta_L1_pp"])

    singles = base_rule_masks(features)
    scored_single = [score_rule(rule, features, deltas, current_selected_mean) for rule in singles]
    scored = list(scored_single)
    single_df = pd.DataFrame(scored_single).sort_values(
        [
            "external_worsened_samples",
            "worsened_samples",
            "mean_sample_delta_L1_pp",
            "max_sample_worse_L1_pp",
        ],
        kind="mergesort",
    )
    pool = []
    for row in single_df.head(MAX_SINGLE_RULES_FOR_PAIRS).to_dict("records"):
        values = features[str(row["feature1"])].astype(float)
        threshold = finite(row["threshold1"])
        pool.append(
            {
                "rule_type": "single",
                "feature1": row["feature1"],
                "op1": row["op1"],
                "threshold1": threshold,
                "feature2": "",
                "op2": "",
                "threshold2": "",
                "mask": values >= threshold if row["op1"] == ">=" else values <= threshold,
            }
        )
    pairs: list[dict[str, object]] = []
    for left, right in itertools.combinations(pool, 2):
        if left["feature1"] == right["feature1"] and left["op1"] == right["op1"]:
            continue
        pairs.append(
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
    scored.extend([score_rule(rule, features, deltas, current_selected_mean) for rule in pairs])
    return pd.DataFrame(scored).sort_values(
        [
            "strict_combined_sample_safe",
            "mean_sample_delta_L1_pp",
            "max_sample_worse_L1_pp",
            "cached_gain_preserved_frac",
        ],
        ascending=[False, True, True, False],
        kind="mergesort",
    )


def rule_text(row: pd.Series | dict[str, object]) -> str:
    feature1 = str(row["feature1"])
    text = f"{feature1}{row['op1']}{float(row['threshold1']):.6g}"
    if str(row.get("rule_type", "")) == "two_feature_and" and str(row.get("feature2", "")):
        text += f" AND {row['feature2']}{row['op2']}{float(row['threshold2']):.6g}"
    return text


def build_audit(features: pd.DataFrame, deltas: pd.DataFrame, rules: pd.DataFrame) -> list[dict[str, object]]:
    current = evaluate_mask(features, deltas, pd.Series(True, index=features.index))
    current_stats = panel_delta_stats(current)
    strict = rules.loc[rules["strict_combined_sample_safe"].astype(bool)].copy()
    best_strict = strict.iloc[0].to_dict() if not strict.empty else {}
    best_overall = rules.iloc[0].to_dict()
    if best_strict:
        decision = "refined_guard_candidate_needs_wrapper_validation"
    else:
        decision = "no_simple_output_guard_meets_strict_combined_gate"
    return [
        {
            "metric": "combined_feature_rows",
            "value": f"samples={features.shape[0]};features={len(feature_columns(features))}",
            "evidence": str(OUT_FEATURES.relative_to(EXP)),
            "decision": "cached_profiles_only",
        },
        {
            "metric": "current_guard_effect",
            "value": (
                f"samples={current_stats['sample_count']};"
                f"improved={current_stats['improved_samples']};"
                f"worsened={current_stats['worsened_samples']};"
                f"mean_L1_delta_pp={current_stats['mean_sample_delta_L1_pp']:.6f};"
                f"selected_mean_L1_delta_pp={current_stats['selected_cached_mean_delta_L1_pp']:.6f};"
                f"external_worsened={current_stats['external_worsened_samples']}"
            ),
            "evidence": f"{SELECTED_SCORES.relative_to(EXP)};{EXTERNAL_SCORES.relative_to(EXP)}",
            "decision": "baseline_candidate",
        },
        {
            "metric": "tested_rules",
            "value": int(rules.shape[0]),
            "evidence": str(OUT_RULES.relative_to(EXP)),
            "decision": "single_and_two_feature_output_guards",
        },
        {
            "metric": "strict_safe_rules",
            "value": int(strict.shape[0]),
            "evidence": str(OUT_RULES.relative_to(EXP)),
            "decision": "cached_gain_fraction_min_0.80",
        },
        {
            "metric": "best_overall_rule",
            "value": (
                f"{rule_text(best_overall)};"
                f"mean_L1_delta_pp={float(best_overall['mean_sample_delta_L1_pp']):.6f};"
                f"worsened={int(best_overall['worsened_samples'])};"
                f"external_worsened={int(best_overall['external_worsened_samples'])};"
                f"cached_gain_frac={float(best_overall['cached_gain_preserved_frac']):.6f}"
            ),
            "evidence": str(OUT_RULES.relative_to(EXP)),
            "decision": "best_ranked_rule",
        },
        {
            "metric": "best_strict_rule",
            "value": (
                "none"
                if not best_strict
                else f"{rule_text(best_strict)};"
                f"mean_L1_delta_pp={float(best_strict['mean_sample_delta_L1_pp']):.6f};"
                f"selected_mean_L1_delta_pp={float(best_strict['selected_cached_mean_delta_L1_pp']):.6f};"
                f"external_mean_L1_delta_pp={float(best_strict['external_mean_delta_L1_pp']):.6f};"
                f"cached_gain_frac={float(best_strict['cached_gain_preserved_frac']):.6f}"
            ),
            "evidence": str(OUT_RULES.relative_to(EXP)),
            "decision": "strict_candidate" if best_strict else "none",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": decision,
        },
    ]


def write_markdown(audit: list[dict[str, object]], rules: pd.DataFrame) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Feature Allocator Combined Guard Refinement",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note searches simple output-only guard refinements for",
        "`guarded-genus-hit-breadth-a002` using cached selected-default profiles",
        "plus external exact-split diagnostic profiles. It does not rerun raw",
        "profiling jobs.",
        "",
        "## Decision",
        "",
        f"- Current guard: `{audit_by_metric['current_guard_effect']['value']}`.",
        f"- Tested rules: `{audit_by_metric['tested_rules']['value']}`.",
        f"- Strict safe rules: `{audit_by_metric['strict_safe_rules']['value']}`.",
        f"- Best strict rule: `{audit_by_metric['best_strict_rule']['value']}`.",
        f"- Promotion decision: `{audit_by_metric['promotion_decision']['decision']}`.",
        "",
        "## Top Rules",
        "",
        "| Rule | Mean L1 delta pp | Worsened | External worsened | Cached gain frac | Strict |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in rules.head(12).to_dict("records"):
        lines.append(
            "| {rule} | {mean:.6f} | {worse} | {external_worse} | {gain:.6f} | {strict} |".format(
                rule=rule_text(row),
                mean=float(row["mean_sample_delta_L1_pp"]),
                worse=int(row["worsened_samples"]),
                external_worse=int(row["external_worsened_samples"]),
                gain=float(row["cached_gain_preserved_frac"]),
                strict=str(bool(row["strict_combined_sample_safe"])),
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_FEATURES.relative_to(EXP)}`",
            f"- `{OUT_RULES.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    features = build_features()
    deltas = build_deltas()
    rules = build_rules(features, deltas)
    audit = build_audit(features, deltas, rules)

    features.sort_values(["panel", "sample"], kind="mergesort").to_csv(
        OUT_FEATURES, sep="\t", index=False
    )
    rules.head(MAX_WRITTEN_RULES).to_csv(OUT_RULES, sep="\t", index=False)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(audit, rules)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
