#!/usr/bin/env python3
"""Validate the abundance release blocker for the selected MinCO default.

This uses cached score/audit TSVs only. It separates three claims:

1. The selected candidate preset improves over the previous MinCO default.
2. The selected preset still should not claim broad abundance superiority.
3. The next abundance target is sample-safe allocation plus one call-recovery
   case, not another panel-mean-only rescaling.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

CANDIDATE_VS_CURRENT = RESULTS / "candidate_default_vs_current.tsv"
CANDIDATE_VS_SYLPH = RESULTS / "candidate_default_vs_sylph.tsv"
BLEND_AUDIT = RESULTS / "candidate_preset_genus_xny_blend_audit.tsv"
ORACLE_AUDIT = RESULTS / "candidate_callset_oracle_feasibility_audit.tsv"
ORACLE_SUMMARY = RESULTS / "candidate_callset_oracle_feasibility_summary.tsv"
FEATURE_ALLOCATOR_CANDIDATE = RESULTS / "feature_allocator_release_candidate_audit.tsv"
FEATURE_ALLOCATOR_REFINED_GUARD = RESULTS / "feature_allocator_combined_guard_audit.tsv"
FEATURE_ALLOCATOR_REFINED_PARITY = RESULTS / "feature_allocator_refined_guard_parity_audit.tsv"
FEATURE_ALLOCATOR_REFINED_SCORE = RESULTS / "feature_allocator_refined_score_replay_audit.tsv"
FEATURE_ALLOCATOR_REFINED_MARINE = RESULTS / "feature_allocator_refined_marine_diagnostic_audit.tsv"
FEATURE_ALLOCATOR_REFINED_HOLDOUT_INVENTORY = (
    RESULTS / "feature_allocator_refined_independent_holdout_inventory_audit.tsv"
)
FEATURE_ALLOCATOR_REFINED_RECOVERY = RESULTS / "feature_allocator_refined_holdout_recovery_audit.tsv"
CAMI3_SOURCE_PROFILE_EXTENSION = RESULTS / "cami3_source_profile_binomial_extension_audit.tsv"
CAMI3_SOURCE_READMAP_RECOVERY_INPUTS = RESULTS / "cami3_source_readmap_recovery_inputs_audit.tsv"

OUT_PANEL = RESULTS / "abundance_release_blocker_panel.tsv"
OUT_AUDIT = RESULTS / "abundance_release_blocker_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_RELEASE_BLOCKER.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    return {row[key]: row for row in read_tsv(path)}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return 0.0
    return float(value)


def build_panel_rows() -> list[dict[str, object]]:
    current_by_panel = {row["panel"]: row for row in read_tsv(CANDIDATE_VS_CURRENT)}
    rows = []
    for row in read_tsv(CANDIDATE_VS_SYLPH):
        panel = row["panel"]
        current = current_by_panel.get(panel, {})
        rows.append(
            {
                "panel": panel,
                "panel_label": row.get("panel_label", panel),
                "samples": row.get("samples", ""),
                "candidate_minus_previous_F1": as_float(current, "delta_pooled_F1"),
                "candidate_minus_previous_L1_pp": as_float(current, "delta_official_L1_pp"),
                "candidate_minus_previous_Pearson": as_float(current, "delta_official_Pearson"),
                "candidate_minus_sylph_F1": as_float(row, "delta_pooled_F1"),
                "candidate_minus_sylph_L1_pp": as_float(row, "delta_official_L1_pp"),
                "candidate_minus_sylph_Pearson": as_float(row, "delta_official_Pearson"),
                "abundance_release_status": (
                    "abundance_gap"
                    if as_float(row, "delta_official_L1_pp") > 0.0
                    else "abundance_not_worse"
                ),
            }
        )
    return rows


def build_audit(panel_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    blend = table_by(BLEND_AUDIT, "metric")
    oracle = table_by(ORACLE_AUDIT, "metric")
    oracle_summary = {row["panel"]: row for row in read_tsv(ORACLE_SUMMARY)}
    feature_allocator = table_by(FEATURE_ALLOCATOR_CANDIDATE, "metric")
    refined_guard = table_by(FEATURE_ALLOCATOR_REFINED_GUARD, "metric")
    refined_parity = table_by(FEATURE_ALLOCATOR_REFINED_PARITY, "metric")
    refined_score = table_by(FEATURE_ALLOCATOR_REFINED_SCORE, "metric")
    refined_marine = table_by(FEATURE_ALLOCATOR_REFINED_MARINE, "metric")
    refined_holdout_inventory = table_by(FEATURE_ALLOCATOR_REFINED_HOLDOUT_INVENTORY, "metric")
    refined_recovery = table_by(FEATURE_ALLOCATOR_REFINED_RECOVERY, "metric")
    source_profile_extension = table_by(CAMI3_SOURCE_PROFILE_EXTENSION, "metric")
    source_readmap_recovery = table_by(CAMI3_SOURCE_READMAP_RECOVERY_INPUTS, "metric")
    all_oracle = oracle_summary.get("all", {})

    n = len(panel_rows)
    prev_f1_regressions = sum(1 for row in panel_rows if float(row["candidate_minus_previous_F1"]) < -1e-12)
    prev_l1_regressions = sum(1 for row in panel_rows if float(row["candidate_minus_previous_L1_pp"]) > 1e-12)
    prev_pearson_regressions = sum(
        1 for row in panel_rows if float(row["candidate_minus_previous_Pearson"]) < -1e-12
    )
    sylph_l1_losses = sum(1 for row in panel_rows if float(row["candidate_minus_sylph_L1_pp"]) > 0.0)
    sylph_pearson_losses = sum(1 for row in panel_rows if float(row["candidate_minus_sylph_Pearson"]) < 0.0)
    sylph_f1_wins = sum(1 for row in panel_rows if float(row["candidate_minus_sylph_F1"]) > 0.0)

    blend_decision = blend.get("promotion_decision", {}).get("value", "")
    blend_direction = blend.get("sample_L1_direction", {}).get("value", "")
    blend_max_worse = blend.get("max_worse_L1_delta_vs_candidate_preset_pp", {}).get("value", "")

    allocation_only = oracle.get("allocation_only_can_close_samples", {}).get("value", "")
    call_recovery = oracle.get("call_recovery_required_samples", {}).get("value", "")
    residual = oracle.get("residual_required_samples", {}).get("value", "")
    oracle_decision = oracle.get("promotion_decision", {}).get("decision", "")
    oracle_headroom = all_oracle.get("mean_allocation_headroom_L1_pp", "")
    truth_detected = all_oracle.get("mean_truth_mass_detected_pct", "")
    feature_allocator_effect = feature_allocator.get("combined_candidate_effect", {}).get("value", "")
    feature_allocator_decision = feature_allocator.get("promotion_decision", {}).get(
        "decision", "not_evaluated"
    )
    refined_guard_effect = refined_guard.get("best_strict_rule", {}).get("value", "")
    refined_guard_decision = refined_guard.get("promotion_decision", {}).get(
        "decision", "not_evaluated"
    )
    refined_parity_decision = refined_parity.get("promotion_decision", {}).get(
        "decision", "not_evaluated"
    )
    refined_score_effect = refined_score.get("score_replay_effect", {}).get("value", "")
    refined_score_match = refined_score.get("score_replay_matches_guard_estimate", {}).get(
        "decision", "not_evaluated"
    )
    refined_score_decision = refined_score.get("promotion_decision", {}).get(
        "value", "not_evaluated"
    )
    refined_marine_effect = refined_marine.get("refined_allocator_marine_effect", {}).get("value", "")
    refined_marine_decision = refined_marine.get("promotion_decision", {}).get(
        "value", "not_evaluated"
    )
    refined_holdout_overlap = refined_holdout_inventory.get(
        "release_grade_manifest_overlap", {}
    ).get("value", "")
    refined_holdout_missing = refined_holdout_inventory.get(
        "hmp_omitted_release_candidate_inputs", {}
    ).get("value", "")
    refined_holdout_decision = refined_holdout_inventory.get("promotion_decision", {}).get(
        "value", "not_evaluated"
    )
    refined_recovery_top = refined_recovery.get("top_recovery_route", {}).get(
        "value", "not_evaluated"
    )
    refined_recovery_step = refined_recovery.get("completed_local_step", {}).get("value", "")
    refined_recovery_blocker = refined_recovery.get("remaining_blocker", {}).get("value", "")
    refined_recovery_runner = refined_recovery.get("post_recovery_runner", {}).get("value", "")
    refined_recovery_decision = refined_recovery.get("default_decision", {}).get(
        "decision", "not_evaluated"
    )
    source_profile_scope = source_profile_extension.get("source_profile_truth_scope", {}).get(
        "value", "not_evaluated"
    )
    source_profile_effect = source_profile_extension.get("refined_allocator_effect", {}).get(
        "value", "not_evaluated"
    )
    source_profile_decision = source_profile_extension.get("promotion_decision", {}).get(
        "decision", "not_evaluated"
    )
    source_readmap_inputs = source_readmap_recovery.get("recovery_blockers", {}).get(
        "value", "not_evaluated"
    )
    source_readmap_artifacts = source_readmap_recovery.get("rescoring_artifacts", {}).get(
        "value", ""
    )
    source_readmap_decision = source_readmap_recovery.get("promotion_decision", {}).get(
        "decision", "not_evaluated"
    )

    return [
        {
            "metric": "candidate_vs_previous_minco_regressions",
            "value": (
                f"F1={prev_f1_regressions}/{n};L1={prev_l1_regressions}/{n};"
                f"Pearson={prev_pearson_regressions}/{n}"
            ),
            "evidence": str(CANDIDATE_VS_CURRENT.relative_to(EXP)),
            "decision": "selected_candidate_dominates_previous_minco"
            if prev_f1_regressions == prev_l1_regressions == prev_pearson_regressions == 0
            else "selected_candidate_regression_present",
        },
        {
            "metric": "candidate_vs_sylph_abundance_losses",
            "value": f"L1_losses={sylph_l1_losses}/{n};Pearson_losses={sylph_pearson_losses}/{n}",
            "evidence": str(CANDIDATE_VS_SYLPH.relative_to(EXP)),
            "decision": "abundance_release_gap" if sylph_l1_losses else "abundance_not_worse",
        },
        {
            "metric": "candidate_vs_sylph_F1_wins",
            "value": f"F1_wins={sylph_f1_wins}/{n}",
            "evidence": str(CANDIDATE_VS_SYLPH.relative_to(EXP)),
            "decision": "F1_mixed_not_broad_sylph_beating",
        },
        {
            "metric": "posthoc_blend_promotion",
            "value": f"{blend_decision};{blend_direction};max_worse={blend_max_worse}",
            "evidence": str(BLEND_AUDIT.relative_to(EXP)),
            "decision": "do_not_promote_allocator_with_sample_regressions"
            if blend_decision == "do_not_promote_abundance_blend"
            else "allocator_status_changed",
        },
        {
            "metric": "oracle_next_target",
            "value": (
                f"allocation_only={allocation_only}/32;call_recovery={call_recovery}/32;"
                f"residual={residual}/32;headroom_pp={oracle_headroom};"
                f"detected_truth_pct={truth_detected}"
            ),
            "evidence": f"{ORACLE_AUDIT.relative_to(EXP)};{ORACLE_SUMMARY.relative_to(EXP)}",
            "decision": oracle_decision or "diagnostic_only_not_default",
        },
        {
            "metric": "feature_allocator_release_candidate",
            "value": feature_allocator_effect,
            "evidence": str(FEATURE_ALLOCATOR_CANDIDATE.relative_to(EXP)),
            "decision": feature_allocator_decision,
        },
        {
            "metric": "feature_allocator_refined_guard_candidate",
            "value": refined_guard_effect,
            "evidence": (
                f"{FEATURE_ALLOCATOR_REFINED_GUARD.relative_to(EXP)};"
                f"{FEATURE_ALLOCATOR_REFINED_PARITY.relative_to(EXP)};"
                f"{FEATURE_ALLOCATOR_REFINED_SCORE.relative_to(EXP)}"
            ),
            "decision": refined_score_decision,
        },
        {
            "metric": "feature_allocator_refined_score_replay",
            "value": refined_score_effect,
            "evidence": str(FEATURE_ALLOCATOR_REFINED_SCORE.relative_to(EXP)),
            "decision": f"{refined_score_match};{refined_score_decision}",
        },
        {
            "metric": "feature_allocator_refined_marine_diagnostic",
            "value": refined_marine_effect,
            "evidence": str(FEATURE_ALLOCATOR_REFINED_MARINE.relative_to(EXP)),
            "decision": refined_marine_decision,
        },
        {
            "metric": "feature_allocator_refined_holdout_inventory",
            "value": f"{refined_holdout_overlap};{refined_holdout_missing}",
            "evidence": str(FEATURE_ALLOCATOR_REFINED_HOLDOUT_INVENTORY.relative_to(EXP)),
            "decision": refined_holdout_decision,
        },
        {
            "metric": "feature_allocator_refined_recovery_plan",
            "value": (
                f"top_route={refined_recovery_top};"
                f"completed={refined_recovery_step};"
                f"remaining={refined_recovery_blocker};"
                f"post_recovery_runner={refined_recovery_runner}"
            ),
            "evidence": str(FEATURE_ALLOCATOR_REFINED_RECOVERY.relative_to(EXP)),
            "decision": refined_recovery_decision,
        },
        {
            "metric": "cami3_source_profile_extension_diagnostic",
            "value": f"{source_profile_scope};{source_profile_effect}",
            "evidence": str(CAMI3_SOURCE_PROFILE_EXTENSION.relative_to(EXP)),
            "decision": source_profile_decision,
        },
        {
            "metric": "cami3_source_readmap_recovery_inputs",
            "value": f"{source_readmap_inputs};{source_readmap_artifacts}",
            "evidence": str(CAMI3_SOURCE_READMAP_RECOVERY_INPUTS.relative_to(EXP)),
            "decision": source_readmap_decision,
        },
        {
            "metric": "abundance_release_decision",
            "value": (
                f"selected_default_L1_losses={sylph_l1_losses}/{n};"
                f"posthoc_blend={blend_decision};oracle={oracle_decision};"
                f"feature_allocator={feature_allocator_decision};"
                f"refined_score_replay={refined_score_decision};"
                f"marine_diagnostic={refined_marine_decision};"
                f"holdout_inventory={refined_holdout_decision};"
                f"recovery_plan={refined_recovery_decision};"
                f"source_profile_extension={source_profile_decision};"
                f"source_readmap_inputs={source_readmap_decision}"
            ),
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "expected_gap_keep_default_no_abundance_claim",
        },
    ]


def write_markdown(panel_rows: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Abundance Release Blocker",
        "",
        "Date: 2026-06-29",
        "",
        "This generated audit validates the abundance claim boundary for the",
        "selected default. It uses cached score and allocator-audit TSVs only.",
        "",
        "## Decision",
        "",
        f"- Candidate versus previous MinCO: `{audit_by_metric['candidate_vs_previous_minco_regressions']['value']}`.",
        f"- Candidate versus external baseline abundance: `{audit_by_metric['candidate_vs_sylph_abundance_losses']['value']}`.",
        f"- Posthoc allocator status: `{audit_by_metric['posthoc_blend_promotion']['value']}`.",
        f"- Oracle next target: `{audit_by_metric['oracle_next_target']['value']}`.",
        f"- Feature allocator candidate: `{audit_by_metric['feature_allocator_release_candidate']['value']}` "
        f"({audit_by_metric['feature_allocator_release_candidate']['decision']}).",
        f"- Refined feature allocator guard: `{audit_by_metric['feature_allocator_refined_guard_candidate']['value']}` "
        f"({audit_by_metric['feature_allocator_refined_guard_candidate']['decision']}).",
        f"- Refined feature allocator score replay: `{audit_by_metric['feature_allocator_refined_score_replay']['value']}` "
        f"({audit_by_metric['feature_allocator_refined_score_replay']['decision']}).",
        f"- Refined feature allocator marine diagnostic: `{audit_by_metric['feature_allocator_refined_marine_diagnostic']['value']}` "
        f"({audit_by_metric['feature_allocator_refined_marine_diagnostic']['decision']}).",
        f"- Refined feature allocator holdout inventory: `{audit_by_metric['feature_allocator_refined_holdout_inventory']['value']}` "
        f"({audit_by_metric['feature_allocator_refined_holdout_inventory']['decision']}).",
        f"- Refined feature allocator recovery plan: `{audit_by_metric['feature_allocator_refined_recovery_plan']['value']}` "
        f"({audit_by_metric['feature_allocator_refined_recovery_plan']['decision']}).",
        f"- CAMI3 source-profile extension diagnostic: `{audit_by_metric['cami3_source_profile_extension_diagnostic']['value']}` "
        f"({audit_by_metric['cami3_source_profile_extension_diagnostic']['decision']}).",
        f"- CAMI3 source-readmap recovery inputs: `{audit_by_metric['cami3_source_readmap_recovery_inputs']['value']}` "
        f"({audit_by_metric['cami3_source_readmap_recovery_inputs']['decision']}).",
        f"- Release decision: `{audit_by_metric['abundance_release_decision']['decision']}`.",
        "",
        "## Panel Deltas",
        "",
        "| Panel | Candidate - Previous L1 pp | Candidate - External L1 pp | Candidate - External F1 | Status |",
        "|---|---:|---:|---:|---|",
    ]
    for row in panel_rows:
        lines.append(
            "| {panel} | {prev_l1:.6f} | {base_l1:.6f} | {base_f1:.6f} | {status} |".format(
                panel=row["panel"],
                prev_l1=float(row["candidate_minus_previous_L1_pp"]),
                base_l1=float(row["candidate_minus_sylph_L1_pp"]),
                base_f1=float(row["candidate_minus_sylph_F1"]),
                status=row["abundance_release_status"],
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_PANEL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    panel_rows = build_panel_rows()
    audit = build_audit(panel_rows)
    write_tsv(
        OUT_PANEL,
        panel_rows,
        [
            "panel",
            "panel_label",
            "samples",
            "candidate_minus_previous_F1",
            "candidate_minus_previous_L1_pp",
            "candidate_minus_previous_Pearson",
            "candidate_minus_sylph_F1",
            "candidate_minus_sylph_L1_pp",
            "candidate_minus_sylph_Pearson",
            "abundance_release_status",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(panel_rows, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
