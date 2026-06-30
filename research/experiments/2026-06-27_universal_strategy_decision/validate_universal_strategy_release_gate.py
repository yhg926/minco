#!/usr/bin/env python3
"""Validate the current universal-strategy decision gates.

This script does not rerun profiling. It checks the selected default, cached
benchmark summaries, documentation boundaries, and known release blockers. The
expected outcome is a supported current MinCO default plus a conservative claim
boundary, not broad Sylph-beating release readiness.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"

sys.path.insert(0, str(ROOT))
from scripts import minco_profile_default as default_launcher  # noqa: E402
from scripts import minco_profile_calibrated as calibrated  # noqa: E402


DEFAULT_SNAPSHOT = RESULTS / "default_release_snapshot_audit.tsv"
DEFAULT_MANIFEST = RESULTS / "current_default_strategy_manifest.tsv"
HOLDOUT_SUMMARY = RESULTS / "holdout_bundle_summary.tsv"
ALLOCATOR_EXTERNAL = RESULTS / "feature_allocator_external_exactsplit_audit.tsv"
ABUNDANCE_BLOCKER = RESULTS / "abundance_release_blocker_audit.tsv"
DEFAULT_CONTRACT_TEST = ROOT / "tests/test_minco_profile_calibrated_auto_exact.py"
CAMI3_AFTER_RECOVERY_AUDIT = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
CAMI3_AFTER_RECOVERY_RUNNER = EXP / "score_cami3_source_readmap_extension_after_recovery.py"
CAMI3_AFTER_RECOVERY_TEST = ROOT / "tests/test_cami3_extension_after_recovery.py"
NEXT_ROUTES_AUDIT = RESULTS / "universal_strategy_next_evidence_routes_audit.tsv"
NEXT_ROUTES_TEST = ROOT / "tests/test_universal_strategy_next_evidence_routes.py"

OUT_TSV = RESULTS / "universal_strategy_release_gate.tsv"
OUT_MD = EXP / "UNIVERSAL_STRATEGY_RELEASE_GATE.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def metric_table(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def status(ok: bool, expected_gap: bool = False) -> str:
    if ok:
        return "pass"
    return "expected_gap" if expected_gap else "fail"


def gate(
    gate_id: str,
    requirement: str,
    observed: str,
    ok: bool,
    evidence: str,
    decision: str,
    expected_gap: bool = False,
) -> dict[str, object]:
    return {
        "gate_id": gate_id,
        "requirement": requirement,
        "observed": observed,
        "status": status(ok, expected_gap),
        "expected_gap": bool(expected_gap),
        "evidence": evidence,
        "decision": decision,
    }


def contains_all(text: str, terms: Iterable[str]) -> bool:
    normalized = " ".join(text.split())
    return all(" ".join(term.split()) in normalized for term in terms)


def hmp_route_blocker_matches(blocker: str, scratch_decision: str) -> bool:
    return (
        "minco_profile=" in blocker
        and "sylph_profile=" in blocker
        and "read_input=" in blocker
        and blocker.endswith(f"scratch_headroom={scratch_decision}")
    )


def split_csv_set(value: str) -> set[str]:
    if value in {"", "none", "NA"}:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def build_rows() -> list[dict[str, object]]:
    snapshot = metric_table(DEFAULT_SNAPSHOT, "metric")
    manifest = metric_table(DEFAULT_MANIFEST, "item")
    holdout = metric_table(HOLDOUT_SUMMARY, "metric")
    allocator = metric_table(ALLOCATOR_EXTERNAL, "metric")
    abundance = metric_table(ABUNDANCE_BLOCKER, "metric")
    cami3_after_recovery = metric_table(CAMI3_AFTER_RECOVERY_AUDIT, "metric")
    next_routes = metric_table(NEXT_ROUTES_AUDIT, "metric")

    readme = read_text(ROOT / "README.md")
    manual = read_text(ROOT / "docs/USER_MANUAL.md")
    changelog = read_text(ROOT / "CHANGELOG.md")
    default_contract_test = read_text(DEFAULT_CONTRACT_TEST)
    cami3_after_recovery_test = read_text(CAMI3_AFTER_RECOVERY_TEST)
    next_routes_test = read_text(NEXT_ROUTES_TEST)

    rows: list[dict[str, object]] = []

    rows.append(
        gate(
            "default_launcher_entrypoint",
            "recommended no-manual-strategy entrypoint exists",
            str(ROOT / "scripts/minco_profile"),
            (ROOT / "scripts/minco_profile").is_file(),
            "scripts/minco_profile;scripts/minco_profile_default.py",
            "entrypoint_present",
        )
    )
    rows.append(
        gate(
            "default_launcher_strategy",
            "default launcher delegates to universal-auto-exact",
            default_launcher.DEFAULT_STRATEGY,
            default_launcher.DEFAULT_STRATEGY == "universal-auto-exact",
            "scripts/minco_profile_default.py",
            "f1_priority_strategy_selected",
        )
    )
    rows.append(
        gate(
            "default_profile_preset",
            "selected user default preset is candidate",
            default_launcher.DEFAULT_PRESET,
            default_launcher.DEFAULT_PRESET == default_launcher.CANDIDATE_PRESET == "candidate",
            "scripts/minco_profile_default.py",
            "candidate_preset_selected",
        )
    )
    rows.append(
        gate(
            "default_manifest",
            "manifest agrees that candidate is the selected default",
            manifest.get("selected_profile_preset", {}).get("value", "NA"),
            manifest.get("selected_profile_preset", {}).get("value") == "candidate",
            "results/current_default_strategy_manifest.tsv",
            "manifest_default_matches_code",
        )
    )
    contract_terms = [
        "test_default_launcher_constants_match_current_strategy_manifest",
        "current_default_strategy_manifest.tsv",
        "selected_profile_preset",
        "default_wrapper.DEFAULT_PRESET",
        "default_wrapper.CANDIDATE_PRESET",
        "default_wrapper.DEFAULT_STRATEGY",
        "wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70",
        "wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70",
        "wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2",
        "wrapper.DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA == 0.0",
        'wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF == "off"',
        'wrapper.ADAPTIVE_CALL_FILTER_SWITCH_OFF == "off"',
    ]
    rows.append(
        gate(
            "default_contract_test",
            "test suite locks selected default against current manifest",
            "test_default_launcher_constants_match_current_strategy_manifest",
            contains_all(default_contract_test, contract_terms),
            "tests/test_minco_profile_calibrated_auto_exact.py",
            "default_contract_regression_covered",
        )
    )
    rows.append(
        gate(
            "candidate_vs_previous_minco",
            "candidate default has no matched-panel regressions versus previous MinCO",
            snapshot.get("candidate_vs_previous_minco", {}).get("value", "NA"),
            snapshot.get("candidate_vs_previous_minco", {}).get("decision")
            == "candidate_dominates_previous_cached_default",
            "results/default_release_snapshot_audit.tsv",
            "selected_default_supported_vs_previous_minco",
        )
    )
    rows.append(
        gate(
            "sylph_claim_boundary",
            "current evidence must not claim broad Sylph beating",
            snapshot.get("candidate_vs_sylph", {}).get("value", "NA"),
            snapshot.get("candidate_vs_sylph", {}).get("decision") == "not_broad_sylph_beating",
            "results/default_release_snapshot_audit.tsv",
            "broad_sylph_beating_claim_rejected",
        )
    )
    rows.append(
        gate(
            "experimental_allocator_default",
            "guarded feature allocator remains experimental/off by default",
            allocator.get("promotion_decision", {}).get("value", "NA"),
            allocator.get("promotion_decision", {}).get("value")
            == "keep_feature_allocator_experimental_off_by_default",
            "results/feature_allocator_external_exactsplit_audit.tsv",
            "allocator_not_promoted",
        )
    )
    rows.append(
        gate(
            "calibrated_allocator_default",
            "calibrated wrapper abundance feature allocator default is off",
            calibrated.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF,
            calibrated.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF == "off",
            "scripts/minco_profile_calibrated.py",
            "experimental_allocator_opt_in_only",
        )
    )
    rows.append(
        gate(
            "ani_reporting",
            "diagnostic ANI reporting uses reported_ani/ZIP-AAF boundary",
            manifest.get("ani_reporting", {}).get("status", "NA"),
            manifest.get("ani_reporting", {}).get("status") == "reporting_ready_not_claim_win"
            and hasattr(calibrated, "reported_ani_values"),
            "results/current_default_strategy_manifest.tsv;scripts/minco_profile_calibrated.py",
            "ani_reporting_documented_not_overclaimed",
        )
    )
    docs_ok = (
        contains_all(
            readme,
            [
                "scripts/minco_profile",
                "scripts/minco_profile_default.py",
                "--profile-preset current",
                "not be read as a general Sylph-beating",
                "reported_ani",
            ],
        )
        and contains_all(
            manual,
            [
                "scripts/minco_profile",
                "scripts/minco_profile_default.py",
                "--profile-preset current",
                "not a general Sylph-beating",
                "reported_ani",
            ],
        )
        and "universal-auto-exact" in changelog
    )
    rows.append(
        gate(
            "docs_claim_boundary",
            "README/manual document default use, fallback preset, ANI, and claim boundary",
            "README/manual/changelog checked",
            docs_ok,
            "README.md;docs/USER_MANUAL.md;CHANGELOG.md",
            "user_facing_default_boundary_documented",
        )
    )
    rows.append(
        gate(
            "release_grade_bundle",
            "all holdout panels are release-grade GTDB species evidence",
            "release_grade_all_panels="
            + holdout.get("release_grade_all_panels", {}).get("value", "NA")
            + ";release_grade_panels="
            + holdout.get("release_grade_panels", {}).get("value", "NA"),
            holdout.get("release_grade_all_panels", {}).get("value") == "true",
            "results/holdout_bundle_summary.tsv",
            "not_release_complete_until_more_clean_gtdb_holdouts",
            expected_gap=True,
        )
    )
    cami3_test_ok = contains_all(
        cami3_after_recovery_test,
        [
            "test_missing_readmap_audit_is_clean_blocker",
            "test_exact_binomial_fallback_accepts_unmapped_source",
            "test_exact_binomial_fallback_rejects_priority_conflict",
            "truth_policy_used",
            "exact_binomial_fallback",
            "direct_source_readmap",
        ],
    )
    cami3_runner_blocked_cleanly = (
        CAMI3_AFTER_RECOVERY_RUNNER.is_file()
        and CAMI3_AFTER_RECOVERY_TEST.is_file()
        and cami3_after_recovery.get("post_recovery_input_status", {}).get("value")
        == "readmaps=0/3;selected_default=3/3;raw_tables=3/3;sylph=3/3"
        and cami3_after_recovery.get("post_recovery_input_status", {}).get("decision")
        == "blocked_until_readmaps_restored"
        and cami3_after_recovery.get("blocking_reason", {}).get("value") == "missing_reads_mapping"
        and cami3_test_ok
    )
    cami3_runner_scored_cleanly = (
        CAMI3_AFTER_RECOVERY_RUNNER.is_file()
        and CAMI3_AFTER_RECOVERY_TEST.is_file()
        and cami3_after_recovery.get("post_recovery_input_status", {}).get("value")
        == "readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3"
        and cami3_after_recovery.get("post_recovery_input_status", {}).get("decision")
        == "ready_to_score"
        and cami3_after_recovery.get("promotion_decision", {}).get("value")
        == "post_recovery_extension_scored"
        and cami3_after_recovery.get("candidate_minco_vs_sylph", {}).get("value")
        == "minco_F1_wins=0/3;minco_L1_wins=0/3"
        and cami3_after_recovery.get("refined_allocator_effect", {}).get("decision")
        == "do_not_promote_from_extension_result"
        and cami3_test_ok
    )
    cami3_runner_ok = cami3_runner_blocked_cleanly or cami3_runner_scored_cleanly
    rows.append(
        gate(
            "cami3_extension_runner_contract",
            "CAMI3 post-recovery extension scorer is ready or scored, and regression-tested",
            (
                cami3_after_recovery.get("post_recovery_input_status", {}).get("value", "NA")
                + ";blocker="
                + cami3_after_recovery.get("blocking_reason", {}).get("value", "NA")
                + ";promotion="
                + cami3_after_recovery.get("promotion_decision", {}).get("value", "NA")
            ),
            cami3_runner_ok,
            (
                "results/cami3_source_readmap_extension_after_recovery_audit.tsv;"
                "score_cami3_source_readmap_extension_after_recovery.py;"
                "tests/test_cami3_extension_after_recovery.py"
            ),
            "post_recovery_runner_scored_or_cleanly_blocked",
        )
    )
    next_route_blocker = next_routes.get("top_route_blocker_detail", {}).get("value", "")
    next_route_decision = next_routes.get("top_next_evidence_route", {}).get("decision")
    completed_route_ids = split_csv_set(
        next_routes.get("completed_evidence_routes", {}).get("value", "")
    )
    ready_route_ids = split_csv_set(
        next_routes.get("routes_ready_without_external_restore", {}).get("value", "")
    )
    external_route_ids = split_csv_set(
        next_routes.get("routes_requiring_external_restore", {}).get("value", "")
    )
    hmp_route_completed = "hmp_airskin_omitted_samples_2_8_12_26_27" in completed_route_ids
    marine_route_completed = "marine_truth_upgrade" in completed_route_ids
    next_route_ready_or_headroom_blocked = (
        (
            next_route_decision == "not_local_without_reads_or_headroom"
            and hmp_route_blocker_matches(
                next_route_blocker,
                "free_tmp_space_before_next_sample",
            )
        )
        or (
            next_route_decision == "not_local_without_reads_or_profiles"
            and hmp_route_blocker_matches(next_route_blocker, "run_next_sample")
        )
    )
    marine_setup_metadata_incomplete_ok = (
        next_route_decision == "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout"
        and "release-grade strict setup-source rule remains below 95%" in next_route_blocker
        and "partial source-name rule is diagnostic only" in next_route_blocker
        and "selected-default same-namespace profile set incomplete" in next_route_blocker
    )
    marine_profiles_needed_ok = (
        next_route_decision == "truth_upgrade_profiles_needed"
        and "strict setup source-name or local assembly-summary rule reaches >=95%" in next_route_blocker
        and "selected-default same-namespace profiles still need scoring" in next_route_blocker
    )
    next_route_progress_state_ok = (
        (
            not hmp_route_completed
            and next_routes.get("top_next_evidence_route", {}).get("value")
            == "hmp_airskin_omitted_samples_2_8_12_26_27"
            and next_route_ready_or_headroom_blocked
            and external_route_ids
            == {
                "hmp_airskin_omitted_samples_2_8_12_26_27",
                "marine_truth_upgrade",
                "refined_allocator_candidate",
            }
        )
        or (
            hmp_route_completed
            and not marine_route_completed
            and next_routes.get("top_next_evidence_route", {}).get("value")
            == "marine_truth_upgrade"
            and next_route_decision
            in {
                "truth_upgrade_needed_before_release_use",
                "truth_upgrade_needs_external_truth_source",
                "truth_upgrade_needs_source_mapping",
                "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout",
                "truth_upgrade_profiles_needed",
            }
            and (
                next_route_decision != "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout"
                or marine_setup_metadata_incomplete_ok
            )
            and (
                next_route_decision != "truth_upgrade_profiles_needed"
                or marine_profiles_needed_ok
            )
            and external_route_ids == {"marine_truth_upgrade", "refined_allocator_candidate"}
        )
        or (
            hmp_route_completed
            and marine_route_completed
            and next_routes.get("top_next_evidence_route", {}).get("value")
            == "refined_allocator_candidate"
            and next_route_decision == "keep_opt_in_until_independent_release_holdout_passes"
            and external_route_ids == {"refined_allocator_candidate"}
        )
    )
    next_routes_ok = (
        "cami3_source_readmap_samples3_5_completed" in completed_route_ids
        and next_routes.get("completed_evidence_routes", {}).get("decision")
        == "completed_routes_scored_negative"
        and next_route_progress_state_ok
        and {"cami3_source_readmap_samples3_5_completed", "current_default_candidate"}.issubset(
            ready_route_ids
        )
        and next_routes.get("routes_requiring_external_restore", {}).get("decision")
        == "expected_gap_requires_input_recovery_not_threshold_tuning"
        and next_routes.get("default_decision", {}).get("decision") == "keep_goal_active"
        and NEXT_ROUTES_TEST.is_file()
        and contains_all(
            next_routes_test,
            [
                "test_next_evidence_routes_rank_cami3_recovery_first",
                "cami3_source_readmap_samples3_5_completed",
                "hmp_airskin_omitted_samples_2_8_12_26_27",
                "test_next_evidence_routes_completed_hmp_omitted_route_moves_to_marine",
                "test_next_evidence_routes_marine_source_readmaps_need_mapping",
                "completed_routes_scored_negative",
                "marine_truth_upgrade",
                "truth_upgrade_needed_before_release_use",
                "truth_upgrade_needs_external_truth_source",
                "truth_upgrade_needs_source_mapping",
                "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout",
                "truth_upgrade_profiles_needed",
                "keep_opt_in_until_independent_release_holdout_passes",
                "marine_setup_truth_profile_rescore",
                "truth_transfer_source_mapping_absent_locally",
                "truth_transfer_setup_metadata_incomplete",
                "truth_transfer_setup_metadata_ready_profiles_needed",
                "test_next_evidence_routes_marine_setup_metadata_still_incomplete",
                "test_next_evidence_routes_marine_setup_metadata_threshold_passes_profiles_needed",
                "test_next_evidence_routes_marine_selected_default_scored_negative",
                "release-grade strict setup-source rule remains below 95",
                "selected-default same-namespace profiles still need scoring",
                "partial source-name rule is diagnostic only",
                "readmap_only_archives=3",
                "not_local_without_reads_or_headroom",
                "scratch_headroom=free_tmp_space_before_next_sample",
                "not_local_without_reads_or_profiles",
                "scratch_headroom=run_next_sample",
                "minco_profile=",
                "sylph_profile=",
                "read_input=",
                "expected_gap_requires_input_recovery_not_threshold_tuning",
                "keep_goal_active",
            ],
        )
    )
    rows.append(
        gate(
            "next_evidence_route_contract",
            "next evidence route is explicit, non-circular, and regression-tested",
            (
                "top="
                + next_routes.get("top_next_evidence_route", {}).get("value", "NA")
                + ";completed="
                + next_routes.get("completed_evidence_routes", {}).get("value", "NA")
                + ";ready_without_external_restore="
                + next_routes.get("routes_ready_without_external_restore", {}).get("value", "NA")
            ),
            next_routes_ok,
            (
                "results/universal_strategy_next_evidence_routes_audit.tsv;"
                "audit_universal_strategy_next_evidence_routes.py;"
                "tests/test_universal_strategy_next_evidence_routes.py"
            ),
            "next_route_ready_no_threshold_sweep_release_claim",
        )
    )
    rows.append(
        gate(
            "abundance_release_claim",
            "selected default should not claim broad abundance superiority",
            abundance.get("abundance_release_decision", {}).get("value", "NA"),
            abundance.get("abundance_release_decision", {}).get("decision")
            != "expected_gap_keep_default_no_abundance_claim",
            "results/abundance_release_blocker_audit.tsv",
            "abundance_expected_gap_not_release_claim",
            expected_gap=True,
        )
    )
    rows.append(
        gate(
            "runtime_claim_boundary",
            "runtime snapshot records speed/memory tradeoff instead of faster-than-Sylph claim",
            snapshot.get("runtime_tradeoff", {}).get("value", "NA"),
            snapshot.get("runtime_tradeoff", {}).get("decision")
            == "minco_slower_lower_memory_in_cached_timed_runs",
            "results/default_release_snapshot_audit.tsv",
            "speed_claim_rejected_memory_advantage_recorded",
        )
    )
    rows.append(
        gate(
            "current_release_decision",
            "release decision is current stable default, not final universal claim",
            snapshot.get("release_decision", {}).get("value", "NA"),
            snapshot.get("release_decision", {}).get("decision")
            == "stable_default_candidate_not_final_universal_claim",
            "results/default_release_snapshot_audit.tsv",
            "goal_progress_not_goal_complete",
        )
    )
    return rows


def audit_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    failures = [row for row in rows if row["status"] == "fail"]
    expected_gaps = [row for row in rows if row["status"] == "expected_gap"]
    passes = [row for row in rows if row["status"] == "pass"]
    return {
        "gate_id": "summary",
        "requirement": "universal strategy release gate summary",
        "observed": (
            f"pass={len(passes)};fail={len(failures)};"
            f"expected_gap={len(expected_gaps)}"
        ),
        "status": "pass" if not failures else "fail",
        "expected_gap": bool(expected_gaps),
        "evidence": str(OUT_TSV.relative_to(EXP)),
        "decision": (
            "current_default_supported_with_expected_release_gaps"
            if not failures
            else "review_failed_release_gates"
        ),
    }


def write_markdown(rows: list[dict[str, object]]) -> None:
    summary = rows[-1]
    lines = [
        "# Universal Strategy Release Gate",
        "",
        "Date: 2026-06-29",
        "",
        "This generated gate validates the current default-strategy decision from",
        "repo-local code, documentation, and cached benchmark summaries. It does",
        "not rerun raw profiling jobs.",
        "",
        "## Summary",
        "",
        f"- Status: `{summary['status']}`",
        f"- Observed: `{summary['observed']}`",
        f"- Decision: `{summary['decision']}`",
        "",
        "The expected gaps are intentional: the current default is supported as",
        "the best cached MinCO default, but the release-grade holdout bundle and",
        "broad abundance claim are not complete enough for a broad universal or",
        "Sylph-beating claim.",
        "",
        "## Gates",
        "",
        "| Gate | Status | Expected gap | Observed | Decision |",
        "|---|---|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            "| {gate_id} | {status} | {gap} | {observed} | {decision} |".format(
                gate_id=row["gate_id"],
                status=row["status"],
                gap=str(row["expected_gap"]),
                observed=str(row["observed"]).replace("|", "/"),
                decision=row["decision"],
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_TSV.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows = build_rows()
    rows.append(audit_summary(rows))
    write_tsv(
        OUT_TSV,
        rows,
        ["gate_id", "requirement", "observed", "status", "expected_gap", "evidence", "decision"],
    )
    write_markdown(rows)
    print("\n".join(f"{row['gate_id']}\t{row['status']}\t{row['decision']}" for row in rows))
    failures = [row for row in rows if row["status"] == "fail"]
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
