#!/usr/bin/env python3
"""Rank remaining local holdout routes after the marine truth-mapping audit.

This script consumes cached manifests and transfer summaries only. It does not
run profilers, download data, or inspect large read files. Its purpose is to
make the next evidence path explicit once the obvious local routes have either
been scored negative or proven nonrelease because truth mapping is insufficient.
"""

from __future__ import annotations

import csv
from pathlib import Path
from statistics import mean
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

MANIFEST = EXP / "holdout_bundle_manifest.tsv"
HOLDOUT_PLAN = RESULTS / "holdout_gap_action_plan.tsv"
NEXT_ROUTES = RESULTS / "universal_strategy_next_evidence_routes.tsv"
MARINE_SOURCE_INVENTORY = RESULTS / "marine_source_mapping_local_inventory_audit.tsv"
MARINE_TRANSFER = RESULTS / "marine_truth_upgrade_candidate_audit.tsv"
MARINE_SETUP_METADATA = RESULTS / "marine_setup_metadata_truth_upgrade_audit.tsv"
MARINE_PROFILE_RESCORE = RESULTS / "marine_setup_truth_profile_rescore_audit.tsv"
PLANT_TRANSFER = RESULTS / "plant_gtdb_transfer_feasibility_summary.tsv"
STRAIN_TRANSFER = RESULTS / "strain_gtdb_transfer_feasibility_summary.tsv"
CAMI3_EXTENSION = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
RELEASE_GATE = RESULTS / "universal_strategy_release_gate.tsv"

OUT_TSV = RESULTS / "remaining_holdout_route_options.tsv"
OUT_AUDIT = RESULTS / "remaining_holdout_route_options_audit.tsv"
OUT_MD = EXP / "REMAINING_HOLDOUT_ROUTE_OPTIONS.md"

ROUTE_FIELDS = [
    "rank",
    "route_id",
    "route_family",
    "local_status",
    "local_ready",
    "best_available_truth_pct",
    "profile_status",
    "current_evidence",
    "missing_evidence",
    "decision",
]

AUDIT_FIELDS = ["metric", "value", "evidence", "decision"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def optional_read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    return read_tsv(path)


def keyed(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in optional_read_tsv(path)}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fnum(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def min_pct(path: Path) -> tuple[float, str]:
    rows = optional_read_tsv(path)
    if not rows:
        return 0.0, "summary_missing"
    pcts = [fnum(row.get("mapped_abundance_pct", "0")) for row in rows]
    ready = sum(str(row.get("release_grade_ready", "")).lower() == "true" for row in rows)
    return min(pcts), f"samples={len(rows)};release_ready={ready}/{len(rows)};mean={mean(pcts):.6f}"


def route(
    rank: int,
    route_id: str,
    route_family: str,
    local_status: str,
    local_ready: bool,
    best_available_truth_pct: str,
    profile_status: str,
    current_evidence: str,
    missing_evidence: str,
    decision: str,
) -> dict[str, object]:
    return {
        "rank": rank,
        "route_id": route_id,
        "route_family": route_family,
        "local_status": local_status,
        "local_ready": str(local_ready).lower(),
        "best_available_truth_pct": best_available_truth_pct,
        "profile_status": profile_status,
        "current_evidence": current_evidence,
        "missing_evidence": missing_evidence,
        "decision": decision,
    }


def manifest_release_count() -> str:
    rows = optional_read_tsv(MANIFEST)
    if not rows:
        return "manifest_missing"
    release = sum(str(row.get("release_grade_ready", "")).lower() == "true" for row in rows)
    return f"release={release};nonrelease={len(rows) - release};total={len(rows)}"


def build_routes() -> list[dict[str, object]]:
    next_routes = keyed(NEXT_ROUTES, "route_id")
    holdout_plan = keyed(HOLDOUT_PLAN, "panel_id")
    marine_inventory = keyed(MARINE_SOURCE_INVENTORY, "metric")
    marine_transfer = keyed(MARINE_TRANSFER, "metric")
    marine_setup = keyed(MARINE_SETUP_METADATA, "metric")
    marine_profile = keyed(MARINE_PROFILE_RESCORE, "metric")
    cami3_extension = keyed(CAMI3_EXTENSION, "metric")

    plant_min, plant_detail = min_pct(PLANT_TRANSFER)
    strain_min, strain_detail = min_pct(STRAIN_TRANSFER)

    marine = next_routes.get("marine_truth_upgrade", {})
    refined = next_routes.get("refined_allocator_candidate", {})
    hmp = next_routes.get("hmp_airskin_omitted_samples_2_8_12_26_27", {})
    cami3 = next_routes.get("cami3_source_readmap_samples3_5_completed", {})
    marine_plan = holdout_plan.get("cami2_marine_0_3_5_gtdb_taxid_transfer", {})
    mixed_plan = holdout_plan.get("mixed_readiness_26", {})
    cami3_plan = holdout_plan.get("cami3_toy_human_gut_0_5", {})

    marine_setup_decision = marine_setup.get("promotion_decision", {}).get("value", "")
    marine_setup_ready = marine_setup_decision in {
        "strict_setup_source_rule_clears_threshold_profiles_needed",
        "strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed",
    }
    marine_route_id = (
        "marine_selected_default_profile_rescore"
        if marine_setup_ready
        else "marine_source_mapping_or_new_clean_holdout"
    )
    marine_truth_pct = (
        marine_setup.get("strict_source_name_or_assembly_sources_unique_summary", {}).get(
            "value", "NA"
        )
        if marine_setup_ready
        else marine_transfer.get("diagnostic_unique_epithet_rescues", {}).get("value", "NA")
    )
    marine_profile_status = (
        "truth_mapping_ready_selected_default_profiles_missing"
        if marine_setup_ready
        else marine.get("profile_status", "selected_default_same_namespace_incomplete")
    )
    marine_missing = (
        "selected-default same-namespace profiles and scoring are needed for this route"
        if marine_setup_ready
        else marine.get(
            "missing_evidence",
            "source mapping or different clean holdout required",
        )
    )
    marine_decision = (
        "top_local_profile_rescore_route"
        if marine_setup_ready
        else "top_external_or_new_holdout_route"
    )
    marine_profile_scored = (
        marine_profile.get("promotion_decision", {}).get("value")
        == "selected_default_profile_scored_review_metrics"
    )
    if marine_profile_scored:
        marine_route_id = "marine_selected_default_profile_scored"
        marine_profile_status = "selected_default_profiles_scored"
        marine_missing = "none"
        marine_decision = (
            "completed_negative_do_not_repeat_without_strategy_change"
            if marine_profile.get("minco_vs_sylph_summary", {}).get("decision")
            != "minco_higher_F1"
            else "scored_supportive_review_before_default_claim"
        )
    marine_rank = 35 if marine_decision.startswith("completed_negative") else 10

    rows = [
        route(
            marine_rank,
            marine_route_id,
            "new_metagenome_type_holdout",
            marine.get("local_status", "unknown"),
            marine_setup_ready or marine_profile_scored,
            marine_truth_pct,
            marine_profile_status,
            ";".join(
                part
                for part in [
                    marine.get("current_evidence", ""),
                    "local_inventory="
                    + marine_inventory.get("promotion_decision", {}).get("decision", "NA"),
                    "setup_metadata="
                    + marine_setup.get("promotion_decision", {}).get("decision", "NA"),
                    "profile_rescore="
                    + marine_profile.get("minco_vs_sylph_summary", {}).get("value", "NA"),
                    "holdout_action="
                    + marine_plan.get("action_type", "NA"),
                ]
                if part
            ),
            marine_missing,
            marine_decision,
        ),
        route(
            20,
            "refined_allocator_independent_release_holdout",
            "abundance_L1_Pearson",
            refined.get("local_status", "unknown"),
            False,
            "NA",
            "opt_in_validated_cached_but_not_independent_release",
            refined.get("current_evidence", ""),
            refined.get(
                "missing_evidence",
                "independent release-grade holdout with unchanged F1 and nonworse abundance",
            ),
            "not_promotable_without_independent_release_holdout",
        ),
        route(
            30,
            "hmp_omitted_samples_completed",
            "abundance_independent_holdout",
            hmp.get("local_status", "unknown"),
            True,
            "high_mapped_source_abundance_same_release",
            "same_release_profiles_scored",
            hmp.get("current_evidence", ""),
            "none",
            "completed_negative_do_not_repeat_without_strategy_change",
        ),
        route(
            40,
            "cami3_samples3_5_completed",
            "completed_independent_holdout",
            cami3.get("local_status", "unknown"),
            True,
            "source_readmap_binomial_fallback_completed",
            "same_namespace_profiles_scored",
            cami3_extension.get("candidate_minco_vs_sylph", {}).get("value", ""),
            "none",
            "completed_negative_do_not_repeat_without_strategy_change",
        ),
        route(
            50,
            "plant_local_transfer",
            "existing_local_diagnostic",
            "truth_transfer_far_below_release_threshold",
            False,
            f"min={plant_min:.6f};{plant_detail}",
            "cached_exactsplit_profiles_exist_but_nonrelease_truth",
            "plant_gtdb_transfer_feasibility_summary.tsv",
            "GTDB mapped truth mass must reach >=95%; current local transfer is too low",
            "not_viable_release_route_locally",
        ),
        route(
            60,
            "strain_local_transfer",
            "existing_local_diagnostic",
            "truth_transfer_far_below_release_threshold",
            False,
            f"min={strain_min:.6f};{strain_detail}",
            "cached_exactsplit_profiles_exist_but_nonrelease_truth",
            "strain_gtdb_transfer_feasibility_summary.tsv",
            "GTDB mapped truth mass must reach >=95%; current local transfer is too low",
            "not_viable_release_route_locally",
        ),
        route(
            70,
            "mixed_readiness_panel",
            "existing_local_diagnostic",
            "mixed_truth_namespaces",
            False,
            "NA",
            "not_single_profile_set",
            mixed_plan.get("blocking_reason", "truth_not_clean_gtdb_species"),
            "clean GTDB species truth and same-release profiles needed",
            "not_viable_release_route_locally",
        ),
        route(
            80,
            "cami3_older_taxid_transfer",
            "existing_local_diagnostic",
            "completed_or_superseded_by_source_readmap",
            False,
            "samples0-2 promoted; samples3-5 completed negative",
            "not_needed_without_strategy_change",
            cami3_plan.get("next_action", ""),
            "none for current route; use only if strategy changes",
            "superseded_by_source_readmap_evidence",
        ),
    ]
    return sorted(rows, key=lambda row: int(row["rank"]))


def build_audit(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    local_viable = [
        row["route_id"]
        for row in rows
        if row["local_ready"] == "true" and not str(row["decision"]).startswith("completed_negative")
    ]
    exhausted = [
        row["route_id"]
        for row in rows
        if "not_viable" in str(row["decision"]) or "completed_negative" in str(row["decision"])
    ]
    top = rows[0]
    local_decision = (
        "local_profile_rescore_route_available"
        if local_viable
        else "no_local_clean_holdout_ready"
    )
    if str(top["decision"]) == "top_local_profile_rescore_route":
        promotion_decision = "profile_rescore_needed_before_closing_gap"
    elif str(top["route_id"]) == "refined_allocator_independent_release_holdout":
        promotion_decision = "new_strategy_or_new_holdout_needed"
    else:
        promotion_decision = "external_source_mapping_or_new_clean_holdout_needed"
    return [
        {
            "metric": "manifest_release_status",
            "value": manifest_release_count(),
            "evidence": str(MANIFEST.relative_to(EXP)),
            "decision": "release_gap_remains",
        },
        {
            "metric": "local_viable_uncompleted_routes",
            "value": ",".join(local_viable) or "none",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": local_decision,
        },
        {
            "metric": "locally_exhausted_routes",
            "value": ",".join(exhausted),
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "do_not_repeat_without_new_truth_or_strategy_change",
        },
        {
            "metric": "top_remaining_route",
            "value": str(top["route_id"]),
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": str(top["decision"]),
        },
        {
            "metric": "promotion_decision",
            "value": "keep_current_default_with_expected_gaps",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": promotion_decision,
        },
    ]


def write_md(rows: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    lines = [
        "# Remaining Holdout Route Options",
        "",
        "Date: 2026-06-29",
        "",
        "This audit ranks remaining local evidence routes after the HMP omitted",
        "route and CAMI3 samples3-5 route were scored negative. It incorporates",
        "the current setup-metadata truth mapping state for the marine route.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Routes",
            "",
            "| Rank | Route | Local Status | Best Truth % | Decision |",
            "|---:|---|---|---:|---|",
        ]
    )
    for row in rows:
        lines.append(
            "| {rank} | `{route_id}` | {local_status} | {best_available_truth_pct} | {decision} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "No threshold-only sweep can close the release-holdout gap. If the",
            "top route is `marine_selected_default_profile_rescore`, local truth",
            "mapping has cleared the threshold and the next evidence step is",
            "same-namespace selected-default profile scoring. Otherwise the next",
            "step remains stronger source mapping or a different clean holdout.",
            "The current default therefore stays selected, with the refined",
            "abundance allocator opt-in.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_routes()
    audit = build_audit(rows)
    write_tsv(OUT_TSV, rows, ROUTE_FIELDS)
    write_tsv(OUT_AUDIT, audit, AUDIT_FIELDS)
    write_md(rows, audit)
    for row in audit:
        print(row["metric"], row["value"], row["decision"], sep="\t")


if __name__ == "__main__":
    main()
