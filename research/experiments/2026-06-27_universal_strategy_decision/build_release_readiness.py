#!/usr/bin/env python3
"""Build a conservative release-readiness gate for the universal strategy.

This is intentionally stricter than the stage-decision summary. A strategy can
be the best current MinCO default while still failing release-grade evidence
requirements for a broad universal claim.
"""

from __future__ import annotations

import csv
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
RESOURCE_AUDIT_SUMMARY = RESULTS / "holdout_resource_audit_summary.tsv"
SYLPH_PREFLIGHT = RESULTS / "gtdb232_sylph_db_preflight.tsv"
RUNTIME_READINESS = RESULTS / "runtime_readiness.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows}


def status_pass(status: str) -> bool:
    return status == "pass" or status.startswith("implemented") or status == "supported"


def main() -> int:
    checks = by_key(read_tsv(RESULTS / "decision_checks.tsv"), "check")
    objective = by_key(read_tsv(RESULTS / "objective_audit.tsv"), "requirement")
    holdout = by_key(read_tsv(RESULTS / "holdout_bundle_summary.tsv"), "metric")
    resources = by_key(read_tsv(RESOURCE_AUDIT_SUMMARY), "metric") if RESOURCE_AUDIT_SUMMARY.exists() else {}
    preflight = by_key(read_tsv(SYLPH_PREFLIGHT), "metric") if SYLPH_PREFLIGHT.exists() else {}
    runtime = by_key(read_tsv(RUNTIME_READINESS), "metric") if RUNTIME_READINESS.exists() else {}

    release_grade_panels = int(float(holdout["release_grade_panels"]["value"]))
    diagnostic_grade_panels = int(float(holdout.get("diagnostic_grade_panels", {"value": "0"})["value"]))
    panels_total = int(float(holdout["panels_total"]["value"]))
    all_release_grade = holdout["release_grade_all_panels"]["value"].lower() == "true"
    same_release_hmp_ready = resources.get("hmp_airskin_same_release_sylph_ready", {}).get("value", "")
    resource_blocker = resources.get("main_blocker", {}).get("value", "")
    r232_preflight_status = preflight.get("preflight_status", {}).get("value", "")
    r232_ready_to_build = preflight.get("ready_to_build", {}).get("value", "")
    r232_chunked_ready = preflight.get("chunked_ready_to_build", {}).get("value", "")
    r232_chunked_complete = preflight.get("chunked_output_complete", {}).get("value", "")
    r232_chunk_count = preflight.get("chunk_expected_count", {}).get("value", "")
    default_speed_supported = (
        runtime.get("default_faster_than_sylph_claim_supported", {}).get("value", "").lower()
        == "true"
    )
    default_speed_evidence = runtime.get(
        "timed_default_faster_than_sylph_cases", {}
    ).get("value", "")
    default_memory_evidence = runtime.get(
        "timed_default_lower_memory_than_sylph_cases", {}
    ).get("value", "")
    sidecar_policy = runtime.get("sidecar_default_promotion_decision", {}).get("value", "")
    block_exact_reuse = runtime.get("block_exact_reuse_decision", {}).get("value", "")
    candidate_exact = runtime.get("candidate_restricted_exact_decision", {}).get("value", "")

    rows: list[dict[str, object]] = [
        {
            "gate": "single_default_entrypoint",
            "status": "pass"
            if status_pass(checks["default_launcher_forces_universal_auto_exact"]["status"])
            and status_pass(checks["default_launcher_help_explains_species_boundary"]["status"])
            and status_pass(checks["docs_recommend_default_launcher"]["status"])
            and status_pass(checks["default_launcher_discovers_packaged_sidecars"]["status"])
            else "fail",
            "evidence": "scripts/minco_profile_default.py help; README.md; docs/USER_MANUAL.md",
            "decision": "required_for_default",
        },
        {
            "gate": "c_profile_boundary_documented",
            "status": checks["c_profile_boundary_documented_as_conservative_direct"]["status"],
            "evidence": checks["c_profile_boundary_documented_as_conservative_direct"]["evidence"],
            "decision": "prevents treating conservative C profile as calibrated default",
        },
        {
            "gate": "f1_priority_beats_previous_minco_gate",
            "status": checks["mixed_panel_F1_beats_legacy_probability"]["status"],
            "evidence": checks["mixed_panel_F1_beats_legacy_probability"]["evidence"],
            "decision": "supports_current_minco_default",
        },
        {
            "gate": "unsafe_rescues_rejected",
            "status": "pass"
            if status_pass(checks["near_split_rejected_by_F1"]["status"])
            and status_pass(checks["loose_cami3_source_readmap_split_rescue_rejected_by_crossval_F1"]["status"])
            and status_pass(checks["hmp_raw_candidate_rescue_sweep_not_default"]["status"])
            and status_pass(checks["adaptive_call_filter_external_exactsplit_rejected_as_default"]["status"])
            else "fail",
            "evidence": "near-split, loose CAMI3 split-rescue, HMP raw-candidate rescue sweep, and adaptive call-filter external exact-split stress checks",
            "decision": "required_for_stability",
        },
        {
            "gate": "cross_panel_candidate_rescue_default_ready",
            "status": "known_gap"
            if checks["cross_panel_candidate_rescue_requires_wrapper_validation"]["status"]
            == "candidate_requires_wrapper_validation"
            else "pass",
            "evidence": checks["cross_panel_candidate_rescue_requires_wrapper_validation"]["value"],
            "decision": "wrapper_validation_subsumed_by_selected_candidate_preset"
            if checks["cross_panel_candidate_rescue_requires_wrapper_validation"]["status"] == "pass"
            else "promising_F1_candidate_needs_wrapper_raw_side_channel_and_abundance_policy",
        },
        {
            "gate": "candidate_preset_ergonomics",
            "status": "pass"
            if status_pass(checks["candidate_preset_launcher_available"]["status"])
            and status_pass(checks["candidate_preset_replay_matches_validated_calls"]["status"])
            else "fail",
            "evidence": (
                f"{checks['candidate_preset_launcher_available']['evidence']}; "
                f"{checks['candidate_preset_replay_matches_validated_calls']['value']}"
            ),
            "decision": "candidate_preset_is_selected_default",
        },
        {
            "gate": "abundance_default_not_replaced_by_local_candidate",
            "status": "pass"
            if status_pass(checks["max_su_zip1_abundance_improves_broad_fixed_call_L1"]["status"])
            and status_pass(checks["zip025_abundance_rejected_by_broad_fixed_call_L1"]["status"])
            and status_pass(checks["abundance_variant_safety_no_strict_fixed_call_replacement"]["status"])
            and status_pass(checks["abundance_oracle_bounds_fixed_calls_not_sufficient_all_panels"]["status"])
            and status_pass(checks["adaptive_abundance_switch_rejected_by_lopo"]["status"])
            and status_pass(checks["supervised_abundance_calibrator_not_default"]["status"])
            and status_pass(checks["supervised_abundance_sample28_holdout_not_default"]["status"])
            and status_pass(checks["supervised_abundance_external_exactsplit_rejected_as_default"]["status"])
            and status_pass(checks["edge_em_default_not_promoted"]["status"])
            else "fail",
            "evidence": "broad fixed-call abundance sweep, four-panel GTDB fixed-call variant sweep, fixed-call sample/panel safety audit, truth-aware fixed-call oracle bounds, output-feature adaptive switch LOPO audit, supervised abundance calibrator audit, sample28 same-panel holdout, external exact-split stress test, candidate normalized-depth abundance cross-panel wrapper validation, and cross-domain edge-EM policy check found no additional abundance replacement beyond the selected candidate preset; experimental genus-XnY blend, supervised calibrator, and edge-EM remain off by default",
            "decision": "keeps_selected_candidate_preset_abundance_default",
        },
        {
            "gate": "diagnostic_ani_reporting_available",
            "status": checks["reported_ani_zip_aaf_output_documented"]["status"],
            "evidence": checks["reported_ani_zip_aaf_output_documented"]["evidence"],
            "decision": "reporting_ready_not_claim_win",
        },
        {
            "gate": "default_faster_than_sylph_timed_claim",
            "status": "pass" if default_speed_supported else "known_gap",
            "evidence": (
                f"default_faster_cases={default_speed_evidence or 'NA'}; "
                f"default_lower_memory_cases={default_memory_evidence or 'NA'}"
                + "; sample22_rawread_timing=MinCO_3:28.88_vs_Sylph_sketch_profile_3:59.34_profile_only_3:10.78"
                + "; sample5_rawread_timing=MinCO_3:32.57_vs_Sylph_sketch_profile_2:52.33"
                + "; sample0_rawread_timing=MinCO_3:33.54_vs_Sylph_sketch_profile_3:18.89"
                + "; sample21_rawread_timing=MinCO_3:32.88_vs_Sylph_sketch_profile_3:30.07_profile_only_2:41.01"
                + "; sample18_rawread_timing=MinCO_3:31.30_vs_Sylph_sketch_profile_2:15.22_profile_only_1:26.17"
                + "; sample13_rawread_timing=MinCO_1:40.59_vs_Sylph_sketch_profile_3:20.23_profile_only_2:29.62"
                + "; sample25_rawread_timing=MinCO_3:34.91_vs_Sylph_sketch_profile_3:10.13_profile_only_2:19.41"
                + "; sample3_rawread_timing=MinCO_3:30.24_vs_Sylph_sketch_profile_2:23.51_profile_only_1:33.43"
                + "; sample1_rawread_timing=MinCO_3:33.32_vs_Sylph_sketch_profile_3:10.56_profile_only_2:20.30"
                + "; sample17_rawread_timing=MinCO_3:31.93_vs_Sylph_sketch_profile_2:59.16_profile_only_2:08.31"
                + "; sample24_rawread_timing=MinCO_3:30.53_vs_Sylph_sketch_profile_2:33.08_profile_only_1:42.73"
                + "; sample16_rawread_timing=MinCO_3:29.27_vs_Sylph_sketch_profile_2:55.53_profile_only_2:05.55"
                + "; sample15_rawread_timing=MinCO_3:30.25_vs_Sylph_sketch_profile_2:54.23_profile_only_2:04.17"
                + "; sample4_rawread_timing=MinCO_3:36.19_vs_Sylph_sketch_profile_2:25.11_profile_only_1:33.79"
                + "; sample23_rawread_timing=MinCO_3:33.17_vs_Sylph_sketch_profile_3:03.23_profile_only_2:12.73"
                + "; sample20_rawread_timing=MinCO_3:35.20_vs_Sylph_sketch_profile_3:02.95_profile_only_2:12.07"
                + "; sample7_rawread_timing=MinCO_3:34.80_vs_Sylph_sketch_profile_2:23.57_profile_only_1:32.90"
                + "; sample9_rawread_timing=MinCO_3:33.90_vs_Sylph_sketch_profile_3:03.06_profile_only_2:12.62"
                + "; sample10_rawread_timing=MinCO_3:36.69_vs_Sylph_sketch_profile_2:21.70_profile_only_1:30.85"
                + "; sample19_rawread_timing=MinCO_3:37.65_vs_Sylph_sketch_profile_3:13.38_profile_only_2:21.57"
                + "; sample14_rawread_timing=MinCO_3:31.61_vs_Sylph_sketch_profile_3:07.66_profile_only_2:17.34"
                + (f"; sidecar_policy={sidecar_policy}" if sidecar_policy else "")
                + (f"; block_exact_reuse={block_exact_reuse}" if block_exact_reuse else "")
                + (f"; candidate_restricted_exact={candidate_exact}" if candidate_exact else "")
            ),
            "decision": (
                "supports_default_speed_claim"
                if default_speed_supported
                else "do_not_claim_default_faster_than_sylph_yet"
            ),
        },
        {
            "gate": "extended_nonrelease_panels_audited",
            "status": "pass"
            if checks["marine_domain_recipe_F1_beats_sylph"]["status"] == "pass"
            and checks["marine_gtdb_transfer_F1_beats_sylph"]["status"] == "pass"
            and checks["marine_gtdb_transfer_abundance_L1_beats_sylph"]["status"] == "pass"
            and checks["marine_gtdb_transfer_truth_quality_release_grade"]["status"] == "known_gap"
            and checks["plant_calibrated_gate_F1_beats_sylph_slightly"]["status"] == "pass"
            and checks["plant_current_exactsplit_F1_beats_sylph_slightly"]["status"] == "pass"
            and checks["plant_current_exactsplit_is_not_release_grade"]["status"] == "known_gap"
            and checks["plant_gtdb_transfer_truth_quality_release_grade"]["status"] == "known_gap"
            and checks["strain_current_exactsplit_F1_beats_sylph"]["status"] == "pass"
            and checks["strain_current_exactsplit_L1_beats_sylph"]["status"] == "known_gap"
            and checks["strain_current_exactsplit_is_not_release_grade"]["status"] == "known_gap"
            and checks["strain_gtdb_transfer_truth_quality_release_grade"]["status"] == "known_gap"
            and checks["cached_exactsplit_diagnostic_audited"]["status"] == "pass"
            and checks["hmp_pilot_sylph_F1_beats_calibrated"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_gastrooral_current_universal_improves_old_calibrated_F1"]["status"]
            == "pass"
            and checks["hmp_gastrooral_current_universal_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_truth_quality_diagnostic"]["status"] == "pass"
            and checks["hmp_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_3sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_4sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin5_source_abundance_F1_beats_sylph"]["status"]
            == "pass"
            and checks["hmp_airskin0_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin1_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_5sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_6sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin13_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin16_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin17_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin18_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin21_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin24_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin25_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin3_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin4_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin23_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin20_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin7_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin9_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin10_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_airskin14_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_7sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_8sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_9sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_10sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_11sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_12sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_13sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_14sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_15sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_16sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_17sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_18sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_19sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_20sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_21sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_22sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_23sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_source_abundance_24sample_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["hmp_gastrooral_source_abundance_truth_quality_release_grade"]["status"]
            == "pass"
            and checks["hmp_gastrooral_source_abundance_F1_beats_sylph"]["status"]
            == "contradicts_broad_sylph_beating_claim"
            and checks["cami3_c_profile_direct_F1_beats_sylph"]["status"]
            == "documented_direct_profile_gap"
            else "fail",
            "evidence": "marine local plus GTDB-transfer partial evidence, plant legacy/current-exactsplit plus plant GTDB-transfer feasibility, strain current-exactsplit plus strain GTDB-transfer feasibility, cached exact-split diagnostic panel, current HMP gastrooral, HMP gastrooral source-abundance, HMP airskin source-abundance through samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28, and CAMI3 C-direct panels are represented in decision_checks.tsv",
            "decision": "prevents cherry-picking domain-specific wins",
        },
        {
            "gate": "c_direct_profile_gap_documented",
            "status": "pass"
            if checks["plant_c_profile_direct_recall_gap"]["status"] == "documented_direct_profile_gap"
            and checks["cami3_c_profile_direct_below_prior_offline_gate"]["status"]
            == "documented_direct_profile_gap"
            and status_pass(checks["default_launcher_help_explains_species_boundary"]["status"])
            else "fail",
            "evidence": "plant/CAMI3 C-direct-profile gaps plus default launcher help boundary",
            "decision": "C direct profile is documented as separate from the calibrated species default",
        },
        {
            "gate": "clean_release_grade_holdout_bundle",
            "status": "pass" if all_release_grade else "fail",
            "evidence": (
                f"{release_grade_panels}/{panels_total} panels audited release-grade; "
                f"{diagnostic_grade_panels} high-quality diagnostic panel(s)"
                + (f"; HMP same-release-ready={same_release_hmp_ready}" if same_release_hmp_ready else "")
                + (
                    f"; GTDB r232 Sylph preflight={r232_preflight_status}"
                    if r232_preflight_status
                    else ""
                )
                + (
                    f"; chunked_output_complete={r232_chunked_complete}"
                    if r232_chunked_complete
                    else ""
                )
                + (
                    f"; ready_to_build={r232_ready_to_build}"
                    if r232_ready_to_build and r232_chunked_complete != "true"
                    else ""
                )
                + (
                    f"; chunked_ready_to_build={r232_chunked_ready}"
                    f" ({r232_chunk_count} chunks)"
                    if r232_chunked_ready
                    else ""
                )
                + (f"; {resource_blocker}" if resource_blocker else "")
            ),
            "decision": "blocks_broad_universal_claim",
        },
        {
            "gate": "abundance_beats_sylph_on_key_panels",
            "status": "pass"
            if status_pass(checks["cami3_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["mouse57_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_3sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_4sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin5_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin0_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin1_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_5sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_6sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin13_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin15_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin16_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin17_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin18_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin21_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin24_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin25_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin3_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin4_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin23_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin20_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin7_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin9_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin10_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_airskin14_source_abundance_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_7sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_8sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_9sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_10sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_11sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_12sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_13sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_14sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_15sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_16sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_17sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_18sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_19sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_20sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_21sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_22sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_23sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_source_abundance_24sample_L1_beats_sylph"]["status"])
            and status_pass(checks["hmp_gastrooral_source_abundance_L1_beats_sylph"]["status"])
            else "known_gap",
            "evidence": "CAMI3, CAMI2 Toy Mouse, HMP gastrooral source-abundance, and HMP airskin source-abundance checks including samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28",
            "decision": "blocks_abundance_win_claim",
        },
        {
            "gate": "broad_sylph_beating_claim_supported",
            "status": "fail"
            if objective["prove broad Sylph-beating strategy"]["status"] == "not_supported"
            else "pass",
            "evidence": objective["prove broad Sylph-beating strategy"]["evidence"],
            "decision": "claim_boundary_must_stay_conservative",
        },
    ]

    blocking = [row["gate"] for row in rows if row["status"] in {"fail", "known_gap"}]
    final_status = "pre_release_candidate" if blocking else "release_ready"
    rows.append(
        {
            "gate": "final_decision",
            "status": final_status,
            "evidence": ";".join(blocking) if blocking else "all gates pass",
            "decision": (
                "Use universal-auto-exact plus the candidate preset as the current MinCO default, "
                "but do not claim broad universal or Sylph-beating release readiness."
                if blocking
                else "Release-grade universal default evidence is complete."
            ),
        }
    )

    write_tsv(RESULTS / "release_readiness.tsv", rows, ["gate", "status", "evidence", "decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
