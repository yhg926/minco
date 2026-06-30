#!/usr/bin/env python3
"""Rank the remaining evidence routes for the universal strategy goal.

This script does not download data or rerun profilers. It turns the current
expected gaps into concrete next-evidence routes so the default strategy cannot
be promoted on a vague or circular benchmark claim.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

READMAP_INPUTS = RESULTS / "cami3_source_readmap_recovery_inputs_audit.tsv"
READMAP_AFTER_RECOVERY = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
HOLDOUT_INVENTORY = RESULTS / "feature_allocator_refined_independent_holdout_inventory_audit.tsv"
HOLDOUT_INVENTORY_DETAIL = RESULTS / "feature_allocator_refined_independent_holdout_inventory.tsv"
HOLDOUT_PLAN = RESULTS / "holdout_gap_action_plan_audit.tsv"
MARINE_TRANSFER = RESULTS / "marine_binomial_transfer_audit.tsv"
MARINE_UPGRADE_CANDIDATES = RESULTS / "marine_truth_upgrade_candidate_audit.tsv"
MARINE_SOURCE_READMAP = RESULTS / "marine_source_readmap_feasibility_audit.tsv"
MARINE_SOURCE_MAPPING_INVENTORY = RESULTS / "marine_source_mapping_local_inventory_audit.tsv"
MARINE_SETUP_METADATA = RESULTS / "marine_setup_metadata_truth_upgrade_audit.tsv"
MARINE_PROFILE_RESCORE = RESULTS / "marine_setup_truth_profile_rescore_audit.tsv"
ABUNDANCE = RESULTS / "abundance_release_blocker_audit.tsv"
RELEASE_GATE = RESULTS / "universal_strategy_release_gate.tsv"
TMP_HEADROOM = RESULTS / "tmp_headroom_audit.tsv"
NEXT_ACTIONS = RESULTS / "next_release_grade_actions.tsv"
HMP_OMITTED_SAMPLES = {2, 8, 12, 26, 27}

OUT_TSV = RESULTS / "universal_strategy_next_evidence_routes.tsv"
OUT_AUDIT = RESULTS / "universal_strategy_next_evidence_routes_audit.tsv"
OUT_MD = EXP / "UNIVERSAL_STRATEGY_NEXT_EVIDENCE_ROUTES.md"


ROUTE_FIELDS = [
    "priority",
    "route_id",
    "gap_target",
    "local_status",
    "ready_without_external_restore",
    "current_evidence",
    "missing_evidence",
    "success_criterion",
    "promotion_if_success",
    "decision",
]

AUDIT_FIELDS = ["metric", "value", "evidence", "decision"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_key(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def optional_by_key(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    return by_key(path, key)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def route(
    priority: int,
    route_id: str,
    gap_target: str,
    local_status: str,
    ready_without_external_restore: bool,
    current_evidence: str,
    missing_evidence: str,
    success_criterion: str,
    promotion_if_success: str,
    decision: str,
) -> dict[str, object]:
    return {
        "priority": priority,
        "route_id": route_id,
        "gap_target": gap_target,
        "local_status": local_status,
        "ready_without_external_restore": str(bool(ready_without_external_restore)).lower(),
        "current_evidence": current_evidence,
        "missing_evidence": missing_evidence,
        "success_criterion": success_criterion,
        "promotion_if_success": promotion_if_success,
        "decision": decision,
    }


def hmp_omitted_counts() -> dict[str, int]:
    if not HOLDOUT_INVENTORY_DETAIL.exists():
        return {"total": 5, "truth": 5, "minco": 0, "sylph": 0, "read_input": 0, "ready": 0}
    omitted = [
        row
        for row in read_tsv(HOLDOUT_INVENTORY_DETAIL)
        if str(row.get("candidate_set", "")).startswith("hmp_airskin_omitted_sample")
    ]
    if not omitted:
        return {"total": 5, "truth": 5, "minco": 0, "sylph": 0, "read_input": 0, "ready": 0}
    return {
        "total": len(omitted),
        "truth": sum(row.get("truth_status") == "truth_file_present" for row in omitted),
        "minco": sum("minco=missing" not in str(row.get("evidence", "")) for row in omitted),
        "sylph": sum("sylph=missing" not in str(row.get("evidence", "")) for row in omitted),
        "read_input": sum("read_input=missing" not in str(row.get("evidence", "")) for row in omitted),
        "ready": sum(int(row.get("independent_ready_n", 0)) > 0 for row in omitted),
    }


def parse_samples(value: str) -> set[int]:
    samples: set[int] = set()
    for item in str(value).split(","):
        item = item.strip()
        if item.isdigit():
            samples.add(int(item))
    return samples


def hmp_completed_score_signal() -> dict[str, str]:
    candidates: list[tuple[int, str, Path, list[dict[str, str]]]] = []
    for path in RESULTS.glob("hmp_current_refresh_r232_source_abundance_sample*_summary.tsv"):
        rows = read_tsv(path)
        if not rows:
            continue
        samples = parse_samples(rows[0].get("samples", ""))
        if HMP_OMITTED_SAMPLES.issubset(samples):
            candidates.append((len(samples), path.name, path, rows))
    if not candidates:
        return {"status": "missing", "value": "hmp_omitted_score_summary=missing"}

    _, _, path, rows = sorted(candidates, key=lambda item: (item[0], item[1]))[-1]
    baseline = next(
        (row for row in rows if row.get("method") == "sylph_gtdb_source_abundance"),
        None,
    )
    minco = next(
        (row for row in rows if row.get("method", "").startswith("minco_current_default")),
        None,
    )
    if baseline is None or minco is None:
        return {
            "status": "missing",
            "value": f"hmp_omitted_score_summary={path.name};baseline_or_minco_row=missing",
        }

    baseline_f1 = float(baseline["mean_F1"])
    minco_f1 = float(minco["mean_F1"])
    baseline_l1 = float(baseline["mean_L1_union_pp"])
    minco_l1 = float(minco["mean_L1_union_pp"])
    status = (
        "completed_supportive"
        if minco_f1 >= baseline_f1 and minco_l1 <= baseline_l1
        else "completed_negative"
    )
    return {
        "status": status,
        "value": (
            f"hmp_omitted_score_summary={path.name};"
            f"samples={baseline.get('samples', 'NA')};"
            f"minco_mean_F1={minco_f1:.6f};baseline_mean_F1={baseline_f1:.6f};"
            f"minco_L1={minco_l1:.6f};baseline_L1={baseline_l1:.6f}"
        ),
    }


def build_routes() -> list[dict[str, object]]:
    readmap = by_key(READMAP_INPUTS, "metric")
    after = by_key(READMAP_AFTER_RECOVERY, "metric")
    inventory = by_key(HOLDOUT_INVENTORY, "metric")
    holdout = by_key(HOLDOUT_PLAN, "metric")
    marine = by_key(MARINE_TRANSFER, "metric")
    marine_upgrade = optional_by_key(MARINE_UPGRADE_CANDIDATES, "metric")
    marine_source_readmap = optional_by_key(MARINE_SOURCE_READMAP, "metric")
    marine_source_inventory = optional_by_key(MARINE_SOURCE_MAPPING_INVENTORY, "metric")
    marine_setup = optional_by_key(MARINE_SETUP_METADATA, "metric")
    marine_profile = optional_by_key(MARINE_PROFILE_RESCORE, "metric")
    abundance = by_key(ABUNDANCE, "metric")
    release = by_key(RELEASE_GATE, "gate_id")
    headroom = optional_by_key(TMP_HEADROOM, "metric")
    next_actions = read_tsv(NEXT_ACTIONS) if NEXT_ACTIONS.exists() else []

    readmap_blocker = readmap["recovery_blockers"]["decision"]
    after_status = after["post_recovery_input_status"]["value"]
    after_scored = after.get("promotion_decision", {}).get("value") == "post_recovery_extension_scored"
    cami3_route = (
        route(
            60,
            "cami3_source_readmap_samples3_5_completed",
            "completed_independent_holdout",
            "scored_negative",
            True,
            after.get("candidate_minco_vs_sylph", {}).get("value", "NA")
            + ";"
            + after.get("refined_allocator_effect", {}).get("value", "NA"),
            "none",
            "Route already scored; result does not support refined-allocator default promotion",
            "Counts as independent holdout evidence against promoting the refined allocator from this route.",
            "completed_negative_holdout_result",
        )
        if after_scored
        else route(
            10,
            "cami3_source_readmap_samples3_5",
            "release_holdout_gap;abundance_independent_holdout",
            "profiles_ready_truth_missing",
            False,
            readmap["rescoring_artifacts"]["value"] + ";" + after_status,
            readmap_blocker,
            "readmaps=3/3 and post-recovery scorer produces selected-default, refined-allocator, and external-baseline scores in one GTDB source-readmap namespace",
            "Can close the nearest independent holdout test for the opt-in refined allocator; default promotion still requires no F1/sample regressions and improved or nonworse abundance.",
            "top_route_requires_readmap_restore_only",
        )
    )
    hmp_headroom = ""
    hmp_counts = hmp_omitted_counts()
    hmp_scores = hmp_completed_score_signal()
    hmp_missing = (
        f"minco_profile={hmp_counts['minco']}/{hmp_counts['total']};"
        f"sylph_profile={hmp_counts['sylph']}/{hmp_counts['total']};"
        f"read_input={hmp_counts['read_input']}/{hmp_counts['total']}"
    )
    hmp_local_status = "truth_present_profiles_missing"
    hmp_decision = "not_local_without_reads_or_profiles"
    if headroom:
        hmp_headroom = (
            "next_sample="
            + headroom.get("next_sample", {}).get("value", "NA")
            + ";tmp_free_GiB="
            + headroom.get("tmp_free_GiB", {}).get("value", "NA")
            + ";estimated_required_GiB="
            + headroom.get("estimated_required_free_GiB", {}).get("value", "NA")
            + ";headroom_pass="
            + headroom.get("headroom_pass", {}).get("value", "NA")
            + ";alternate_free_GiB="
            + headroom.get("alternate_free_GiB", {}).get("value", "NA")
            + ";alternate_headroom_pass="
            + headroom.get("alternate_headroom_pass", {}).get("value", "NA")
            + ";recommended_action="
            + headroom.get("recommended_action", {}).get("value", "NA")
        )
        hmp_missing += (
            ";scratch_headroom="
            + headroom.get("recommended_action", {}).get("value", "NA")
        )
        if headroom.get("headroom_pass", {}).get("value") != "True":
            hmp_local_status = "truth_present_profiles_missing_headroom_blocked"
            hmp_decision = "not_local_without_reads_or_headroom"
    action_status = ""
    for action in next_actions:
        if action.get("action_id", "").startswith("restore_profile_hmp_airskin_sample"):
            action_status = (
                "action_id="
                + action.get("action_id", "NA")
                + ";"
                "action_status="
                + action.get("status", "NA")
                + ";action_artifact="
                + action.get("artifact", "NA")
            )
            break
    hmp_current_evidence = inventory["hmp_omitted_release_candidate_inputs"]["value"]
    if hmp_headroom:
        hmp_current_evidence += ";" + hmp_headroom
    if action_status:
        hmp_current_evidence += ";" + action_status
    hmp_complete = (
        hmp_counts["total"] > 0
        and hmp_counts["ready"] == hmp_counts["total"]
        and hmp_counts["minco"] == hmp_counts["total"]
        and hmp_counts["sylph"] == hmp_counts["total"]
        and hmp_counts["read_input"] == hmp_counts["total"]
    )
    if hmp_complete:
        hmp_current_evidence += ";" + hmp_scores["value"]
        hmp_local_status = (
            "scored_negative"
            if hmp_scores["status"] == "completed_negative"
            else "profiles_ready_scores_missing"
        )
        hmp_missing = "none" if hmp_scores["status"] != "missing" else hmp_scores["value"]
        hmp_decision = (
            "completed_negative_holdout_result"
            if hmp_scores["status"] == "completed_negative"
            else "score_ready_profiles_before_claim"
        )
    marine_current = marine["all_samples_min_mapped_pct_bacteria_archaea"]["value"]
    marine_missing = (
        "mapped Bacteria/Archaea truth mass below 95%; selected-default same-namespace profile set incomplete"
    )
    marine_status = "truth_transfer_partial"
    marine_decision = "truth_upgrade_needed_before_release_use"
    if (
        marine_upgrade.get("promotion_decision", {}).get("value")
        == "do_not_promote_marine_truth_upgrade"
    ):
        marine_current += (
            ";"
            + marine_upgrade.get("diagnostic_unique_epithet_rescues", {}).get("value", "NA")
        )
        marine_missing = (
            "accepted and diagnostic candidate rules remain below 95% mapped Bacteria/Archaea truth mass; "
            "requires accession-level/source-specific truth mapping or a different clean holdout; "
            "selected-default same-namespace profile set incomplete"
        )
        marine_status = "truth_transfer_candidate_rules_negative"
        marine_decision = "truth_upgrade_needs_external_truth_source"
    if (
        marine_source_readmap.get("promotion_decision", {}).get("value")
        == "do_not_promote_marine_source_readmap_from_local_ids"
    ):
        marine_current += (
            ";source_readmap_rows="
            + marine_source_readmap.get("sampled_readmap_rows", {}).get("value", "NA")
            + ";"
            + marine_source_readmap.get("direct_sequence_to_gtdb_matches", {}).get(
                "value", "NA"
            )
            + ";unresolved_rows_sampled="
            + marine_source_readmap.get("unresolved_transfer_rows_sampled", {}).get(
                "value", "NA"
            )
        )
        marine_missing = (
            "accepted and diagnostic candidate rules remain below 95% mapped Bacteria/Archaea truth mass; "
            "local readmaps cover unresolved rows but source IDs do not directly map to GTDB assemblies; "
            "requires contig/OTU-to-assembly metadata, source FASTA provenance, accession-level/source-specific truth mapping, "
            "or a different clean holdout; selected-default same-namespace profile set incomplete"
        )
        marine_status = "truth_transfer_candidate_rules_and_readmaps_negative"
        marine_decision = "truth_upgrade_needs_source_mapping"
    if (
        marine_source_inventory.get("promotion_decision", {}).get("value")
        == "do_not_promote_marine_from_local_inventory"
    ):
        marine_current += (
            ";local_source_mapping="
            + marine_source_inventory.get("local_marine_crosswalk_files_missing", {}).get(
                "decision", "NA"
            )
            + ";"
            + marine_source_inventory.get("archive_relevant_member_status", {}).get(
                "value", "NA"
            )
        )
        marine_missing = (
            "accepted and diagnostic candidate rules remain below 95% mapped Bacteria/Archaea truth mass; "
            "local readmaps cover unresolved rows but source IDs do not directly map to GTDB assemblies; "
            "local inventory found no marine genome_to_id.tsv, metadata.tsv, or pooled mapping crosswalk and read archives are readmap-only; "
            "requires contig/OTU-to-assembly metadata, source FASTA provenance, accession-level/source-specific truth mapping, "
            "or a different clean holdout; selected-default same-namespace profile set incomplete"
        )
        marine_status = "truth_transfer_source_mapping_absent_locally"
        marine_decision = "truth_upgrade_needs_source_mapping"
    if marine_setup.get("promotion_decision", {}).get("value"):
        marine_current += (
            ";setup_metadata="
            + marine_setup.get("setup_metadata_inputs", {}).get("value", "NA")
            + ";"
            + marine_setup.get("exact_source_name_match_summary", {}).get("value", "NA")
            + ";"
            + marine_setup.get("assembly_summary_match_summary", {}).get("value", "NA")
            + ";"
            + marine_setup.get("strict_source_name_or_assembly_sources_unique_summary", {}).get(
                "value", "NA"
            )
            + ";"
            + marine_setup.get("diagnostic_partial_setup_unique_source_name_summary", {}).get(
                "value", "NA"
            )
        )
        if (
            marine_setup.get("promotion_decision", {}).get("value")
            in {
                "strict_setup_source_rule_clears_threshold_profiles_needed",
                "strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed",
            }
        ):
            marine_missing = (
                "strict setup source-name or local assembly-summary rule reaches >=95% mapped Bacteria/Archaea truth mass; "
                "selected-default same-namespace profiles still need scoring before this route can close the holdout gap"
            )
            marine_status = "truth_transfer_setup_metadata_ready_profiles_needed"
            marine_decision = "truth_upgrade_profiles_needed"
        else:
            marine_missing = (
                "local setup metadata is present and exact source-name matching partially resolves GTDB species, "
                "but the release-grade strict setup-source rule remains below 95% mapped Bacteria/Archaea truth mass; "
                "the partial source-name rule is diagnostic only because unmatched setup sources remain; "
                "requires stronger source-specific mapping or a different clean holdout; "
                "selected-default same-namespace profile set incomplete"
            )
            marine_status = "truth_transfer_setup_metadata_incomplete"
            marine_decision = "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout"

    marine_ready = False
    if (
        marine_profile.get("promotion_decision", {}).get("value")
        == "selected_default_profile_scored_review_metrics"
    ):
        marine_current += (
            ";profile_rescore="
            + marine_profile.get("minco_vs_sylph_summary", {}).get("value", "NA")
        )
        marine_missing = "none"
        marine_ready = True
        if marine_profile.get("minco_vs_sylph_summary", {}).get("decision") == "minco_higher_F1":
            marine_status = "scored_supportive_review"
            marine_decision = "scored_supportive_review_before_default_claim"
        else:
            marine_status = "scored_negative"
            marine_decision = "completed_negative_holdout_result"

    routes = [
        cami3_route,
        route(
            50 if hmp_decision == "completed_negative_holdout_result" else 20,
            "hmp_airskin_omitted_samples_2_8_12_26_27",
            "abundance_independent_holdout",
            hmp_local_status,
            hmp_decision == "completed_negative_holdout_result",
            hmp_current_evidence,
            hmp_missing,
            "same-release MinCO and external-baseline profiles exist for all five omitted samples and score under GTDB source-abundance truth",
            "Counts as independent HMP same-domain abundance stress evidence; a negative result keeps the refined allocator opt-in and blocks broad default-promotion claims.",
            hmp_decision,
        ),
        route(
            30,
            "marine_truth_upgrade",
            "new_metagenome_type_holdout",
            marine_status,
            marine_ready,
            marine_current,
            marine_missing,
            "clean GTDB species truth reaches >=95% mapped in-scope mass and selected-default profiles are scored in the same namespace",
            "Would strengthen broad dataset-type coverage; not sufficient by itself to promote the abundance allocator.",
            marine_decision,
        ),
        route(
            40,
            "refined_allocator_candidate",
            "abundance_L1_Pearson",
            "opt_in_validated_cached_but_not_independent_release",
            False,
            abundance["feature_allocator_refined_score_replay"]["value"]
            + ";"
            + abundance["feature_allocator_refined_marine_diagnostic"]["value"],
            inventory["promotion_decision"]["value"],
            "independent release-grade holdout shows unchanged F1 and no sample-level L1/Pearson regression under the guard",
            "If successful on the next independent release holdout, it becomes the first viable default-promotion candidate for abundance.",
            "keep_opt_in_until_independent_release_holdout_passes",
        ),
        route(
            90,
            "current_default_candidate",
            "user_default_now",
            "supported_current_default_with_expected_gaps",
            True,
            release["summary"]["observed"] + ";goal_not_complete_expected_gaps_remain",
            "final universal claim still lacks release-holdout and abundance evidence",
            "Continue passing release gate while expected gaps are explicitly recorded and not overclaimed",
            "Current user-facing default remains scripts/minco_profile_default.py --profile-preset candidate.",
            "keep_current_default_do_not_overclaim",
        ),
    ]
    return sorted(routes, key=lambda row: int(row["priority"]))


def build_audit(routes: list[dict[str, object]]) -> list[dict[str, object]]:
    action_routes = [row for row in routes if row["decision"] != "completed_negative_holdout_result"]
    top = action_routes[0] if action_routes else routes[0]
    local_ready = [row for row in routes if row["ready_without_external_restore"] == "true"]
    external_needed = [
        row
        for row in routes
        if row["ready_without_external_restore"] == "false"
        and row["decision"] != "completed_negative_holdout_result"
    ]
    completed = [row for row in routes if row["decision"] == "completed_negative_holdout_result"]
    completed_ids = [str(row["route_id"]) for row in completed]
    rows: list[dict[str, object]] = [
        {
            "metric": "completed_evidence_routes",
            "value": ",".join(completed_ids) or "none",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": (
                "completed_routes_scored_negative"
                if completed_ids
                else "no_completed_gap_closing_route_yet"
            ),
        },
        {
            "metric": "top_next_evidence_route",
            "value": top["route_id"],
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": top["decision"],
        },
        {
            "metric": "routes_ready_without_external_restore",
            "value": ",".join(str(row["route_id"]) for row in local_ready) or "none",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "only_current_default_audit_ready_locally"
            if len(local_ready) == 1 and local_ready[0]["route_id"] == "current_default_candidate"
            else "review_local_routes",
        },
        {
            "metric": "routes_requiring_external_restore",
            "value": ",".join(str(row["route_id"]) for row in external_needed),
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "expected_gap_requires_input_recovery_not_threshold_tuning",
        },
        {
            "metric": "default_decision",
            "value": "current_default_candidate_with_refined_allocator_opt_in",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "keep_goal_active",
        },
    ]
    rows.insert(
        2,
        {
            "metric": "top_route_blocker_detail",
            "value": str(top["missing_evidence"]),
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": str(top["decision"]),
        },
    )
    return rows


def write_md(routes: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Universal Strategy Next Evidence Routes",
        "",
        "Generated by `audit_universal_strategy_next_evidence_routes.py` from cached decision TSVs. No profiling or downloads are run.",
        "",
        f"Top route: `{audit_by_metric['top_next_evidence_route']['value']}`.",
        "",
        f"Default decision: `{audit_by_metric['default_decision']['decision']}`.",
        "",
        "| Priority | Route | Local Status | Missing Evidence | Decision |",
        "| --- | --- | --- | --- | --- |",
    ]
    for item in routes:
        lines.append(
            "| {priority} | {route_id} | {local_status} | {missing_evidence} | {decision} |".format(
                **item
            )
        )
    lines.extend(
        [
            "",
            "The machine-readable route table is `results/universal_strategy_next_evidence_routes.tsv`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    routes = build_routes()
    audit = build_audit(routes)
    write_tsv(OUT_TSV, routes, ROUTE_FIELDS)
    write_tsv(OUT_AUDIT, audit, AUDIT_FIELDS)
    write_md(routes, audit)
    print(f"wrote {OUT_TSV}")
    print(f"wrote {OUT_AUDIT}")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
