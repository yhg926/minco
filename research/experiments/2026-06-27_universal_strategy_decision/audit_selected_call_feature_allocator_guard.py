#!/usr/bin/env python3
"""Audit output-derived guards for selected-call feature allocators."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import audit_selected_call_mass_transform_guard as guard


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SCORES_TSV = RESULTS / "selected_call_feature_allocator_scores.tsv"
OVERALL_TSV = RESULTS / "selected_call_feature_allocator_overall.tsv"
SELECTED_SCORES_TSV = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"

FEATURES_TSV = RESULTS / "selected_call_feature_allocator_guard_features.tsv"
RULES_TSV = RESULTS / "selected_call_feature_allocator_guard_rules.tsv"
LOPO_TSV = RESULTS / "selected_call_feature_allocator_guard_lopo.tsv"
AUDIT_TSV = RESULTS / "selected_call_feature_allocator_guard_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_GUARD.md"


def configure_guard_module() -> None:
    guard.SCORES_TSV = SCORES_TSV
    guard.OVERALL_TSV = OVERALL_TSV
    guard.SELECTED_SCORES_TSV = SELECTED_SCORES_TSV


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def rule_text(row: pd.Series) -> str:
    return guard.rule_text(row)


def audit_rows(rules: pd.DataFrame, lopo: pd.DataFrame, methods: list[str]) -> list[dict[str, object]]:
    nontrivial = rules.loc[
        rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])
    ]
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
        decision = "reject_selected_feature_allocator_guard"

    return [
        {
            "metric": "candidate_methods",
            "value": ",".join(methods),
            "evidence": "top selected-call feature allocators by cached overall ranking",
            "decision": "truth_used_only_for_method_scoring",
        },
        {
            "metric": "tested_guard_rules",
            "value": int(nontrivial.shape[0]),
            "evidence": f"method plus output-derived threshold guards;top_saved={guard.MAX_WRITTEN_RULES}",
            "decision": "truth_used_only_for_scoring",
        },
        {
            "metric": "sample_safe_guards",
            "value": int(sample_safe.shape[0]),
            "evidence": f"full in-memory scoring;top_saved={guard.MAX_WRITTEN_RULES}",
            "decision": "strict_sample_gate",
        },
        {
            "metric": "panel_safe_guards",
            "value": int(panel_safe.shape[0]),
            "evidence": f"full in-memory scoring;top_saved={guard.MAX_WRITTEN_RULES}",
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
            "evidence": f"selected_call_feature_allocator_guard_rules.tsv;top_saved={guard.MAX_WRITTEN_RULES}",
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
            "evidence": f"selected_call_feature_allocator_guard_rules.tsv;top_saved={guard.MAX_WRITTEN_RULES}",
            "decision": "candidate" if not best_safe.empty else "none",
        },
        {
            "metric": "lopo_guard_result",
            "value": (
                f"mean_delta={lopo_mean:.6f};holdouts_with_mean_regression={lopo_mean_worse};"
                f"holdouts_with_sample_regression={lopo_sample_regress};"
                f"max_sample_worse={lopo_max_worse:.6f}"
            ),
            "evidence": f"selected_call_feature_allocator_guard_lopo.tsv;top_rules={guard.MAX_LOPO_RULES}",
            "decision": "leave_one_panel_out_validation",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": "feature allocator guard rules plus leave-one-panel-out",
            "decision": decision,
        },
    ]


def write_markdown(audit: list[dict[str, object]], rules: pd.DataFrame, lopo: pd.DataFrame) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Selected Call-Set Feature Allocator Guard Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached selected-candidate profiles",
        "and the selected-call feature allocator sweep. It tests whether",
        "output-derived guards can apply a feature allocator only on samples",
        "where it is stable.",
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
            "- `results/selected_call_feature_allocator_guard_features.tsv`",
            f"- `results/selected_call_feature_allocator_guard_rules.tsv` top {guard.MAX_WRITTEN_RULES} scored rules",
            "- `results/selected_call_feature_allocator_guard_lopo.tsv`",
            "- `results/selected_call_feature_allocator_guard_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    configure_guard_module()
    methods = guard.candidate_methods()
    features = guard.build_features()
    deltas = guard.load_deltas(methods)
    rules = guard.build_rules(features, deltas, methods)
    lopo = guard.leave_one_panel_out(rules, features, deltas)
    audit = audit_rows(rules, lopo, methods)

    features.sort_values(["panel", "sample"], kind="mergesort").to_csv(FEATURES_TSV, sep="\t", index=False)
    rules.head(guard.MAX_WRITTEN_RULES).to_csv(RULES_TSV, sep="\t", index=False)
    lopo.to_csv(LOPO_TSV, sep="\t", index=False)
    guard.write_tsv(AUDIT_TSV, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(audit, rules, lopo)
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
