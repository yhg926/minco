#!/usr/bin/env python3
"""Write a concrete action plan for remaining universal-strategy evidence gaps."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

MANIFEST = EXP / "holdout_bundle_manifest.tsv"
AUDIT = RESULTS / "holdout_bundle_audit.tsv"
SUMMARY = RESULTS / "holdout_bundle_summary.tsv"
RELEASE_GATE = RESULTS / "universal_strategy_release_gate.tsv"
DEFAULT_SNAPSHOT = RESULTS / "default_release_snapshot_audit.tsv"
CAMI3_SCOPE_AUDIT = RESULTS / "cami3_source_readmap_scope_audit.tsv"
CAMI3_RESOLVER_SUMMARY = RESULTS / "cami3_source_readmap_resolver_candidate_summary.tsv"
CAMI3_RESOLVER_AUDIT = RESULTS / "cami3_source_readmap_resolver_candidate_audit.tsv"
CAMI3_BINOMIAL_RESCORING_AUDIT = RESULTS / "cami3_source_readmap_binomial_fallback_audit.tsv"
CAMI3_TRUTH_POLICY = RESULTS / "cami3_binomial_fallback_truth_policy_audit.tsv"
CAMI3_EXTENSION_CACHE_AUDIT = RESULTS / "cami3_source_readmap_extension_cache_audit.tsv"
CAMI3_EXTENSION_AFTER_RECOVERY_AUDIT = (
    RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
)
MARINE_BINOMIAL_AUDIT = RESULTS / "marine_binomial_transfer_audit.tsv"
MARINE_CANDIDATE_AUDIT = RESULTS / "marine_truth_upgrade_candidate_audit.tsv"
MARINE_SOURCE_READMAP_AUDIT = RESULTS / "marine_source_readmap_feasibility_audit.tsv"
MARINE_SOURCE_MAPPING_INVENTORY_AUDIT = RESULTS / "marine_source_mapping_local_inventory_audit.tsv"
MARINE_SETUP_METADATA_AUDIT = RESULTS / "marine_setup_metadata_truth_upgrade_audit.tsv"
MARINE_PROFILE_RESCORE_AUDIT = RESULTS / "marine_setup_truth_profile_rescore_audit.tsv"

OUT_PLAN = RESULTS / "holdout_gap_action_plan.tsv"
OUT_AUDIT = RESULTS / "holdout_gap_action_plan_audit.tsv"
OUT_MD = EXP / "HOLDOUT_GAP_ACTION_PLAN.md"


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


def table_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def sample_count(samples: str) -> int:
    value = samples.strip()
    if not value:
        return 0
    if value.isdigit():
        return int(value)
    return len([part for part in value.split(",") if part.strip()])


def action_for_panel(
    manifest: dict[str, str],
    audit: dict[str, str],
    cami3_scope: dict[str, dict[str, str]],
    cami3_resolver: dict[str, str],
    cami3_extension: dict[str, str],
    marine_binomial: dict[str, str],
    marine_candidate: dict[str, str],
    marine_source_readmap: dict[str, str],
    marine_source_inventory: dict[str, str],
    marine_setup_metadata: dict[str, str],
    marine_profile_rescore: dict[str, str],
) -> dict[str, object]:
    panel_id = manifest["panel_id"]
    grade = audit.get("audit_grade", "nonrelease")
    truth_status = manifest.get("truth_status", "")
    family = manifest.get("dataset_family", "")
    samples = manifest.get("samples", "")
    blocker = audit.get("blocking_reason", "")

    priority = 90
    action_type = "maintain_release_panel"
    next_action = "Keep as release-grade evidence and rerun only after default-strategy changes."
    release_upgrade_path = "already_release_grade"
    if grade == "diagnostic":
        priority = 20
        action_type = "upgrade_truth_mapping"
        next_action = (
            "Improve source-to-GTDB species mapping or add missing source genome mappings; "
            "target >=95% mapped truth mass before treating as release-grade."
        )
        release_upgrade_path = "diagnostic_to_release_if_truth_mapping_reaches_threshold"
    elif grade == "nonrelease":
        priority = 50
        action_type = "replace_or_rebuild_panel"
        next_action = (
            "Do not use for release claims. Rebuild this dataset with clean GTDB species truth "
            "and same-release MinCO/Sylph profiles, or keep it diagnostic only."
        )
        release_upgrade_path = "requires_clean_gtdb_truth_and_comparable_profiles"

    if "CAMI3" in family and grade != "release":
        priority = 35
        if cami3_resolver.get("truth_policy_accepted") == "true":
            extension_note = ""
            if cami3_extension:
                if cami3_extension.get("after_recovery_scored") == "true":
                    priority = 60
                    action_type = "completed_negative_holdout_review"
                    release_upgrade_path = "completed_extension_does_not_support_default_promotion"
                    extension_note = (
                        " Post-recovery extension is scored for samples3-5: "
                        f"{cami3_extension.get('candidate_vs_sylph', 'NA')}; "
                        f"allocator effect={cami3_extension.get('allocator_effect', 'NA')}."
                    )
                    next_action = (
                        "Do not repeat this completed CAMI3 extension unless inputs or strategy "
                        "change. Keep the samples0-2 accepted truth policy as release evidence, "
                        "treat samples3-5 as a negative holdout for allocator promotion, and use "
                        f"another independent route for the next release-gap closure.{extension_note}"
                    )
                else:
                    extension_note = (
                        " Local extension-cache audit: accepted subset ready="
                        f"{cami3_extension.get('accepted_ready', 'NA')}; extension ready="
                        f"{cami3_extension.get('extension_ready', 'NA')}; blockers="
                        f"{cami3_extension.get('extension_blockers', 'NA')}."
                    )
                    next_action = (
                        "Do not use this older CAMI3 panel for release claims. CAMI3 "
                        "source-readmap samples0-2 are already promoted under the accepted "
                        "exact-binomial fallback truth policy. Extending CAMI3 beyond samples0-2 "
                        "requires source-readmap truth for samples3-5 and selected-default profiles "
                        f"in the same namespace.{extension_note}"
                    )
        else:
            priority = 10 if grade == "diagnostic" else 25
            scope_status = cami3_scope.get("profile_scope_mapping_status", {}).get("value", "NA")
            rows_needed = cami3_scope.get("additional_rows_needed_for_95pct", {}).get("value", "NA")
            top_sources = cami3_scope.get("top_in_scope_sources_to_review", {}).get("value", "NA")
            next_action = (
                "CAMI3 is the highest-value next release candidate: extend the source-readmap "
                "truth transfer, account for unmapped read rows, and rerun same-release MinCO/Sylph "
                "only if mapped truth mass can reach release-grade coverage. Cached scope audit: "
                f"{scope_status}; additional rows needed for 95%={rows_needed}; "
                f"top in-scope sources={top_sources}."
            )
    elif "marine" in family and grade != "release":
        priority = 30
        marine_note = ""
        if marine_binomial:
            marine_note = (
                " Exact-binomial fallback audit rescues truth mass but remains below "
                f"release threshold: all-sample min mapped B/A %={marine_binomial.get('min_all', 'NA')}, "
                f"scored-sample min={marine_binomial.get('min_scored', 'NA')}; "
                f"decision={marine_binomial.get('decision', 'NA')}."
            )
        if marine_candidate.get("promotion_decision") == "do_not_promote_marine_truth_upgrade":
            action_type = "external_truth_source_needed"
            release_upgrade_path = "requires_source_specific_truth_or_new_clean_holdout"
            marine_note += (
                " Candidate-rule audit also fails: "
                f"{marine_candidate.get('diagnostic_unique_epithet_rescues', 'NA')}; "
                f"decision={marine_candidate.get('decision', 'NA')}."
            )
        if (
            marine_source_readmap.get("promotion_decision")
            == "do_not_promote_marine_source_readmap_from_local_ids"
        ):
            action_type = "source_mapping_needed"
            release_upgrade_path = "requires_contig_or_otu_to_assembly_mapping_or_new_clean_holdout"
            marine_note += (
                " Source-readmap feasibility audit sampled "
                f"{marine_source_readmap.get('sampled_readmap_rows', 'NA')} rows and "
                f"{marine_source_readmap.get('unresolved_transfer_rows_sampled', 'NA')} unresolved rows, "
                "but found no direct source-sequence-to-GTDB assembly matches; "
                f"decision={marine_source_readmap.get('decision', 'NA')}."
            )
        if (
            marine_source_inventory.get("promotion_decision")
            == "do_not_promote_marine_from_local_inventory"
        ):
            action_type = "source_mapping_needed"
            release_upgrade_path = "requires_external_source_mapping_or_new_clean_holdout"
            marine_note += (
                " Local source-mapping inventory also fails: "
                f"{marine_source_inventory.get('local_missing', 'NA')}; "
                f"{marine_source_inventory.get('archive_status', 'NA')}; "
                f"decision={marine_source_inventory.get('decision', 'NA')}."
            )
        if marine_setup_metadata.get("promotion_decision"):
            if (
                marine_setup_metadata.get("promotion_decision")
                in {
                    "strict_setup_source_rule_clears_threshold_profiles_needed",
                    "strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed",
                }
            ):
                action_type = "profile_rescore_needed"
                release_upgrade_path = "requires_same_namespace_selected_default_profiles"
                marine_note += (
                    " Setup metadata plus local assembly-summary strict source rule clears the truth threshold, "
                    "but same-namespace selected-default profiles are still needed: "
                    f"{marine_setup_metadata.get('strict_rule', 'NA')}; "
                    f"decision={marine_setup_metadata.get('decision', 'NA')}."
                )
            else:
                action_type = "stronger_source_mapping_needed"
                release_upgrade_path = "requires_stronger_source_mapping_or_new_clean_holdout"
                marine_note += (
                    " Extracted setup metadata is present but still insufficient under "
                    "the release-grade strict source-name rule: "
                    f"{marine_setup_metadata.get('setup_inputs', 'NA')}; "
                    f"{marine_setup_metadata.get('match_summary', 'NA')}; "
                    f"{marine_setup_metadata.get('strict_rule', 'NA')}; "
                    f"{marine_setup_metadata.get('diagnostic_rule', 'NA')}; "
                    f"decision={marine_setup_metadata.get('decision', 'NA')}."
                )
        if (
            marine_profile_rescore.get("promotion_decision")
            == "selected_default_profile_scored_review_metrics"
        ):
            if marine_profile_rescore.get("comparison_decision") == "minco_higher_F1":
                action_type = "profile_rescore_review"
                release_upgrade_path = "requires_review_before_default_claim"
                marine_note += (
                    " Selected-default same-namespace profiles are now scored and MinCO has higher F1: "
                    f"{marine_profile_rescore.get('comparison', 'NA')}; "
                    f"decision={marine_profile_rescore.get('decision', 'NA')}."
                )
            else:
                priority = 65
                action_type = "completed_negative_holdout_review"
                release_upgrade_path = "completed_marine_profile_rescore_does_not_support_default_promotion"
                marine_note += (
                    " Selected-default same-namespace profiles are now scored and do not support promotion: "
                    f"{marine_profile_rescore.get('comparison', 'NA')}; "
                    f"decision={marine_profile_rescore.get('decision', 'NA')}."
                )
        if action_type == "profile_rescore_needed":
            next_action = (
                "Marine truth mapping now clears the >=95% mapped GTDB species threshold "
                "under the strict setup source-name or local assembly-summary rule. "
                "The next useful action is selected-default same-namespace profile scoring "
                f"and comparison, not another local threshold sweep.{marine_note}"
            )
        else:
            next_action = (
                "Marine is useful for domain breadth but current truth transfer is still partial. "
                "Do not promote until mapped GTDB species truth reaches >=95% and selected-default "
                "profiles are available in the same namespace. The next useful action is source-specific "
                f"truth mapping or a different clean holdout, not another local threshold sweep.{marine_note}"
            )
        if action_type == "completed_negative_holdout_review":
            next_action = (
                "Marine truth mapping and selected-default profile scoring are complete for "
                "the available same-namespace samples, and the result is negative for default "
                "promotion. Do not repeat this marine route unless the strategy or input "
                f"profiles change.{marine_note}"
            )
    elif "HMP_airskin" in family and grade == "release":
        priority = 80
        next_action = (
            "Release panel is refreshed to the 24-sample same-release r232 source-abundance run. "
            "Use it as a strong counterexample/accuracy stress panel, not as a Sylph-beating claim."
        )
    elif "HMP_gastrooral" in family and grade == "release":
        priority = 85
        next_action = (
            "Keep as a small but high-quality same-release source-abundance counterexample; "
            "use it to test abundance improvements before promotion."
        )
    elif "toy_mouse" in family and grade == "release":
        priority = 75
        next_action = (
            "Keep as clean GTDB species release evidence; useful for confirming F1 regressions "
            "when candidate call/rescue rules change."
        )

    return {
        "priority_rank": priority,
        "panel_id": panel_id,
        "dataset_family": family,
        "samples": samples,
        "sample_count": sample_count(samples),
        "audit_grade": grade,
        "truth_namespace": manifest.get("truth_namespace", ""),
        "truth_status": truth_status,
        "blocking_reason": blocker,
        "action_type": action_type,
        "release_upgrade_path": release_upgrade_path,
        "next_action": next_action,
        "metric_source": manifest.get("metric_source", ""),
    }


def build_plan() -> list[dict[str, object]]:
    manifest_by_panel = table_by(MANIFEST, "panel_id")
    audit_by_panel = table_by(AUDIT, "panel_id")
    cami3_scope = table_by(CAMI3_SCOPE_AUDIT, "metric") if CAMI3_SCOPE_AUDIT.exists() else {}
    cami3_resolver = build_cami3_resolver_summary()
    cami3_extension = build_cami3_extension_summary()
    marine_binomial = build_marine_binomial_summary()
    marine_candidate = build_marine_candidate_summary()
    marine_source_readmap = build_marine_source_readmap_summary()
    marine_source_inventory = build_marine_source_mapping_inventory_summary()
    marine_setup_metadata = build_marine_setup_metadata_summary()
    marine_profile_rescore = build_marine_profile_rescore_summary()
    rows = [
        action_for_panel(
            manifest_by_panel[panel_id],
            audit_by_panel.get(panel_id, {}),
            cami3_scope,
            cami3_resolver,
            cami3_extension,
            marine_binomial,
            marine_candidate,
            marine_source_readmap,
            marine_source_inventory,
            marine_setup_metadata,
            marine_profile_rescore,
        )
        for panel_id in manifest_by_panel
    ]
    rows.sort(key=lambda row: (int(row["priority_rank"]), row["panel_id"]))
    return rows


def build_cami3_resolver_summary() -> dict[str, str]:
    out: dict[str, str] = {}
    if not CAMI3_RESOLVER_SUMMARY.exists():
        return out
    rows = read_tsv(CAMI3_RESOLVER_SUMMARY)
    projected_ready = [
        row["sample"]
        for row in rows
        if str(row.get("projected_profile_scope_release_ready", "")).lower() == "true"
    ]
    recommended_rows = sum(int(float(row.get("recommended_resolvable_rows", "0") or 0)) for row in rows)
    projected_pct = ";".join(
        "sample{}={:.3f}%".format(
            row["sample"],
            float(row.get("projected_profile_scope_mapped_pct_after_recommended", "0") or 0.0),
        )
        for row in rows
    )
    out.update(
        {
            "projected_ready_samples": ",".join(projected_ready),
            "recommended_rows": str(recommended_rows),
            "projected_pct": projected_pct,
        }
    )
    if CAMI3_BINOMIAL_RESCORING_AUDIT.exists():
        audit = table_by(CAMI3_BINOMIAL_RESCORING_AUDIT, "metric")
        out.update(
            {
                "rescore_ready_samples": audit.get(
                    "profile_scope_release_ready_after_fallback", {}
                ).get("value", ""),
                "rescore_f1": audit.get("minco_vs_sylph_mean_F1_after_rescore", {}).get("value", ""),
                "rescore_l1": audit.get("minco_vs_sylph_mean_L1_after_rescore", {}).get("value", ""),
            }
        )
    if CAMI3_TRUTH_POLICY.exists():
        policy = table_by(CAMI3_TRUTH_POLICY, "metric")
        out["truth_policy_accepted"] = str(
            policy.get("truth_policy_decision", {}).get("value", "")
            == "accept_exact_binomial_fallback_for_release_truth"
        ).lower()
    return out


def build_cami3_extension_summary() -> dict[str, str]:
    out: dict[str, str] = {}
    if CAMI3_EXTENSION_CACHE_AUDIT.exists():
        audit = table_by(CAMI3_EXTENSION_CACHE_AUDIT, "metric")
        out.update(
            {
                "accepted_ready": audit.get("accepted_samples_release_ready", {}).get("value", ""),
                "extension_ready": audit.get("extension_samples_release_ready", {}).get("value", ""),
                "extension_blockers": audit.get("extension_blockers", {}).get("value", ""),
            }
        )
    if CAMI3_EXTENSION_AFTER_RECOVERY_AUDIT.exists():
        after = table_by(CAMI3_EXTENSION_AFTER_RECOVERY_AUDIT, "metric")
        if after.get("promotion_decision", {}).get("value") == "post_recovery_extension_scored":
            out.update(
                {
                    "after_recovery_scored": "true",
                    "candidate_vs_sylph": after.get("candidate_minco_vs_sylph", {}).get(
                        "value", ""
                    ),
                    "allocator_effect": after.get("refined_allocator_effect", {}).get(
                        "value", ""
                    ),
                }
            )
    return out


def build_marine_binomial_summary() -> dict[str, str]:
    if not MARINE_BINOMIAL_AUDIT.exists():
        return {}
    audit = table_by(MARINE_BINOMIAL_AUDIT, "metric")
    return {
        "min_all": audit.get("all_samples_min_mapped_pct_bacteria_archaea", {}).get("value", ""),
        "min_scored": audit.get("scored_samples_min_mapped_pct_bacteria_archaea", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
    }


def build_marine_candidate_summary() -> dict[str, str]:
    if not MARINE_CANDIDATE_AUDIT.exists():
        return {}
    audit = table_by(MARINE_CANDIDATE_AUDIT, "metric")
    return {
        "promotion_decision": audit.get("promotion_decision", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
        "diagnostic_unique_epithet_rescues": audit.get(
            "diagnostic_unique_epithet_rescues", {}
        ).get("value", ""),
    }


def build_marine_source_readmap_summary() -> dict[str, str]:
    if not MARINE_SOURCE_READMAP_AUDIT.exists():
        return {}
    audit = table_by(MARINE_SOURCE_READMAP_AUDIT, "metric")
    return {
        "promotion_decision": audit.get("promotion_decision", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
        "sampled_readmap_rows": audit.get("sampled_readmap_rows", {}).get("value", ""),
        "unresolved_transfer_rows_sampled": audit.get(
            "unresolved_transfer_rows_sampled", {}
        ).get("value", ""),
        "direct_sequence_to_gtdb_matches": audit.get(
            "direct_sequence_to_gtdb_matches", {}
        ).get("value", ""),
    }


def build_marine_source_mapping_inventory_summary() -> dict[str, str]:
    if not MARINE_SOURCE_MAPPING_INVENTORY_AUDIT.exists():
        return {}
    audit = table_by(MARINE_SOURCE_MAPPING_INVENTORY_AUDIT, "metric")
    return {
        "promotion_decision": audit.get("promotion_decision", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
        "local_missing": audit.get("local_marine_crosswalk_files_missing", {}).get(
            "value", ""
        ),
        "archive_status": audit.get("archive_relevant_member_status", {}).get("value", ""),
    }


def build_marine_setup_metadata_summary() -> dict[str, str]:
    if not MARINE_SETUP_METADATA_AUDIT.exists():
        return {}
    audit = table_by(MARINE_SETUP_METADATA_AUDIT, "metric")
    return {
        "promotion_decision": audit.get("promotion_decision", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
        "setup_inputs": audit.get("setup_metadata_inputs", {}).get("value", ""),
        "match_summary": audit.get("exact_source_name_match_summary", {}).get("value", ""),
        "assembly_summary": audit.get("assembly_summary_match_summary", {}).get("value", ""),
        "strict_rule": audit.get(
            "strict_source_name_or_assembly_sources_unique_summary", {}
        ).get("value", ""),
        "diagnostic_rule": audit.get(
            "diagnostic_partial_setup_unique_source_name_summary", {}
        ).get("value", ""),
    }


def build_marine_profile_rescore_summary() -> dict[str, str]:
    if not MARINE_PROFILE_RESCORE_AUDIT.exists():
        return {}
    audit = table_by(MARINE_PROFILE_RESCORE_AUDIT, "metric")
    return {
        "promotion_decision": audit.get("promotion_decision", {}).get("value", ""),
        "decision": audit.get("promotion_decision", {}).get("decision", ""),
        "comparison": audit.get("minco_vs_sylph_summary", {}).get("value", ""),
        "comparison_decision": audit.get("minco_vs_sylph_summary", {}).get("decision", ""),
    }


def build_audit(plan: list[dict[str, object]]) -> list[dict[str, object]]:
    summary = table_by(SUMMARY, "metric")
    gate = table_by(RELEASE_GATE, "gate_id")
    snapshot = table_by(DEFAULT_SNAPSHOT, "metric")
    release = [row for row in plan if row["audit_grade"] == "release"]
    diagnostic = [row for row in plan if row["audit_grade"] == "diagnostic"]
    nonrelease = [row for row in plan if row["audit_grade"] == "nonrelease"]
    top = [row["panel_id"] for row in plan[:3]]
    return [
        {
            "metric": "holdout_bundle_status",
            "value": (
                f"release={summary.get('release_grade_panels', {}).get('value', 'NA')};"
                f"diagnostic={summary.get('diagnostic_grade_panels', {}).get('value', 'NA')};"
                f"all_release={summary.get('release_grade_all_panels', {}).get('value', 'NA')}"
            ),
            "evidence": "holdout_bundle_summary.tsv",
            "decision": "expected_release_gap",
        },
        {
            "metric": "panel_counts_by_grade",
            "value": f"release={len(release)};diagnostic={len(diagnostic)};nonrelease={len(nonrelease)}",
            "evidence": "holdout_gap_action_plan.tsv",
            "decision": "action_plan_scope",
        },
        {
            "metric": "top_next_actions",
            "value": ",".join(top),
            "evidence": "holdout_gap_action_plan.tsv",
            "decision": "prioritized_by_release_upgrade_value",
        },
        {
            "metric": "default_claim_boundary",
            "value": snapshot.get("release_decision", {}).get("value", "NA"),
            "evidence": "default_release_snapshot_audit.tsv",
            "decision": snapshot.get("release_decision", {}).get("decision", "NA"),
        },
        {
            "metric": "release_gate_status",
            "value": gate.get("summary", {}).get("observed", "NA"),
            "evidence": "universal_strategy_release_gate.tsv",
            "decision": gate.get("summary", {}).get("decision", "NA"),
        },
        {
            "metric": "abundance_gate_status",
            "value": gate.get("abundance_release_claim", {}).get("observed", "NA"),
            "evidence": "universal_strategy_release_gate.tsv;ABUNDANCE_RELEASE_BLOCKER.md",
            "decision": gate.get("abundance_release_claim", {}).get("decision", "NA"),
        },
    ]


def write_markdown(plan: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Holdout Gap Action Plan",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note turns the holdout component of the release-gate",
        "expected gaps into concrete next actions. It uses cached manifests and",
        "summaries only.",
        "",
        "## Summary",
        "",
        f"- Holdout status: `{audit_by_metric['holdout_bundle_status']['value']}`",
        f"- Panel counts: `{audit_by_metric['panel_counts_by_grade']['value']}`",
        f"- Release gate: `{audit_by_metric['release_gate_status']['value']}`",
        f"- Abundance gate: `{audit_by_metric['abundance_gate_status']['value']}` "
        f"({audit_by_metric['abundance_gate_status']['decision']})",
        "",
        "## Priority Actions",
        "",
        "| Rank | Panel | Grade | Samples | Action |",
        "|---:|---|---|---:|---|",
    ]
    for row in plan:
        lines.append(
            "| {rank} | {panel} | {grade} | {n} | {action} |".format(
                rank=int(row["priority_rank"]),
                panel=row["panel_id"],
                grade=row["audit_grade"],
                n=int(row["sample_count"]),
                action=str(row["next_action"]).replace("|", "/"),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- The current default remains the candidate preset.",
            "- The abundance expected gap is tracked separately in",
            "  `ABUNDANCE_RELEASE_BLOCKER.md`; this plan only ranks",
            "  holdout/truth/profile coverage actions.",
            "- The next evidence target is not another local threshold sweep; it is",
            "  cleaner release-grade GTDB holdout coverage in panel families that",
            "  still lack high-coverage GTDB truth or selected-default profiles.",
            "- The 24-sample HMP airskin release panel is now reflected in the",
            "  manifest, but it is one dataset family and remains a counterexample",
            "  where the external baseline is stronger for F1 and abundance.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_PLAN.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    plan = build_plan()
    audit = build_audit(plan)
    write_tsv(
        OUT_PLAN,
        plan,
        [
            "priority_rank",
            "panel_id",
            "dataset_family",
            "samples",
            "sample_count",
            "audit_grade",
            "truth_namespace",
            "truth_status",
            "blocking_reason",
            "action_type",
            "release_upgrade_path",
            "next_action",
            "metric_source",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(plan, audit)
    print("\n".join(f"{row['priority_rank']}\t{row['panel_id']}\t{row['audit_grade']}" for row in plan))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
