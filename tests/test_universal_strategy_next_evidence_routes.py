#!/usr/bin/env python3
"""Regression checks for the universal-strategy next-evidence route audit."""

from __future__ import annotations

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
sys.path.insert(0, str(EXP))

import audit_universal_strategy_next_evidence_routes as routes  # noqa: E402


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def patch_paths(monkeypatch, tmp_path: Path) -> Path:
    results = tmp_path / "results"
    monkeypatch.setattr(routes, "EXP", tmp_path)
    monkeypatch.setattr(routes, "RESULTS", results)
    for name in [
        "READMAP_INPUTS",
        "READMAP_AFTER_RECOVERY",
        "HOLDOUT_INVENTORY",
        "HOLDOUT_INVENTORY_DETAIL",
        "HOLDOUT_PLAN",
        "MARINE_TRANSFER",
        "MARINE_UPGRADE_CANDIDATES",
        "MARINE_SOURCE_READMAP",
        "MARINE_SOURCE_MAPPING_INVENTORY",
        "MARINE_SETUP_METADATA",
        "MARINE_PROFILE_RESCORE",
        "ABUNDANCE",
        "RELEASE_GATE",
        "TMP_HEADROOM",
        "NEXT_ACTIONS",
        "OUT_TSV",
        "OUT_AUDIT",
    ]:
        monkeypatch.setattr(routes, name, results / getattr(routes, name).name)
    monkeypatch.setattr(routes, "OUT_MD", tmp_path / routes.OUT_MD.name)
    return results


def seed_route_inputs(
    results: Path,
    after_scored: bool = False,
    headroom_pass: bool = False,
    omitted_ready_pairs: int = 0,
) -> None:
    write_tsv(
        results / "cami3_source_readmap_recovery_inputs_audit.tsv",
        [
            {
                "metric": "rescoring_artifacts",
                "value": "candidate_profiles=3/3;raw_tables=3/3;sylph_profiles=3/3",
                "evidence": "fixture",
                "decision": "rescoring_artifacts_ready",
            },
            {
                "metric": "recovery_blockers",
                "value": "missing_reads_mapping",
                "evidence": "fixture",
                "decision": "requires_readmap_restore_or_stream_extract",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )
    write_tsv(
        results / "cami3_source_readmap_extension_after_recovery_audit.tsv",
        (
            [
                {
                    "metric": "post_recovery_input_status",
                    "value": "readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3",
                    "evidence": "fixture",
                    "decision": "ready_to_score",
                },
                {
                    "metric": "candidate_minco_vs_sylph",
                    "value": "minco_F1_wins=0/3;minco_L1_wins=0/3",
                    "evidence": "fixture",
                    "decision": "independent_holdout_comparison",
                },
                {
                    "metric": "refined_allocator_effect",
                    "value": "improved=0;worsened=0;mean_L1_delta_pp=0.000000000",
                    "evidence": "fixture",
                    "decision": "do_not_promote_from_extension_result",
                },
                {
                    "metric": "promotion_decision",
                    "value": "post_recovery_extension_scored",
                    "evidence": "fixture",
                    "decision": "review_holdout_result_before_default_change",
                },
            ]
            if after_scored
            else [
            {
                "metric": "post_recovery_input_status",
                "value": "readmaps=0/3;selected_default=3/3;raw_tables=3/3;sylph=3/3",
                "evidence": "fixture",
                "decision": "blocked_until_readmaps_restored",
            }
            ]
        ),
        ["metric", "value", "evidence", "decision"],
    )
    write_tsv(
        results / "feature_allocator_refined_independent_holdout_inventory_audit.tsv",
        [
            {
                "metric": "hmp_omitted_release_candidate_inputs",
                "value": (
                    "samples=2,8,12,26,27;truth=5/5;"
                    f"minco_profile={omitted_ready_pairs}/5;"
                    f"sylph_profile={omitted_ready_pairs}/5;"
                    f"ready_pairs={omitted_ready_pairs}/5"
                ),
                "evidence": "fixture",
                "decision": "missing_cached_profile_pairs",
            },
            {
                "metric": "promotion_decision",
                "value": "release_grade_independent_holdout_missing",
                "evidence": "fixture",
                "decision": "keep_refined_allocator_opt_in",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )
    inventory_rows = []
    omitted_samples = [2, 8, 12, 26, 27]
    ready_samples = set(omitted_samples[:omitted_ready_pairs])
    for sample in omitted_samples:
        ready = sample in ready_samples
        inventory_rows.append(
            {
                "candidate_set": f"hmp_airskin_omitted_sample{sample}",
                "audit_grade": "candidate_release_if_profiles_exist",
                "truth_namespace": "GTDB_species_source_abundance",
                "truth_status": "truth_file_present",
                "samples": str(sample),
                "sample_count": "1",
                "feature_panel": "hmp_airskin_gtdb_source_abundance",
                "selected_cache_overlap_samples": "",
                "selected_cache_overlap_n": "0",
                "independent_manifest_samples": str(sample),
                "independent_manifest_n": "1",
                "independent_ready_n": "1" if ready else "0",
                "decision": "independent_release_ready_profile_pair_available"
                if ready
                else "missing_minco_or_sylph_profile_for_release_grade_replay",
                "evidence": (
                    "truth=/tmp/truth.tsv;"
                    f"minco={'/tmp/minco.tsv' if ready else 'missing'};"
                    f"sylph={'/tmp/sylph.tsv' if ready else 'missing'};"
                    f"read_input={'/tmp/reads.fq.gz' if ready else 'missing'}"
                ),
            }
        )
    write_tsv(
        results / "feature_allocator_refined_independent_holdout_inventory.tsv",
        inventory_rows,
        [
            "candidate_set",
            "audit_grade",
            "truth_namespace",
            "truth_status",
            "samples",
            "sample_count",
            "feature_panel",
            "selected_cache_overlap_samples",
            "selected_cache_overlap_n",
            "independent_manifest_samples",
            "independent_manifest_n",
            "independent_ready_n",
            "decision",
            "evidence",
        ],
    )
    if omitted_ready_pairs == 5:
        write_tsv(
            results
            / (
                "hmp_current_refresh_r232_source_abundance_sample12_"
                "0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_summary.tsv"
            ),
            [
                {
                    "method": "sylph_gtdb_source_abundance",
                    "samples": "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28",
                    "mean_F1": "0.9580666992769648",
                    "pooled_TP": "1331",
                    "pooled_FP": "61",
                    "pooled_FN": "45",
                    "pooled_F1": "0.9617052023121387",
                    "mean_L1_union_pp": "6.5063842487957375",
                    "mean_L1_truth_only_pp": "5.554687588033204",
                    "mean_Pearson_union": "0.9970001808564978",
                    "mean_Pearson_truth_only": "0.9976435020126806",
                },
                {
                    "method": "minco_current_default_gtdb_source_abundance_sample12_0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28",
                    "samples": "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28",
                    "mean_F1": "0.8640773343135554",
                    "pooled_TP": "1293",
                    "pooled_FP": "407",
                    "pooled_FN": "83",
                    "pooled_F1": "0.840702210663199",
                    "mean_L1_union_pp": "23.811634317988677",
                    "mean_L1_truth_only_pp": "17.90396277873434",
                    "mean_Pearson_union": "0.9788938103023659",
                    "mean_Pearson_truth_only": "0.9855452376219913",
                },
            ],
            [
                "method",
                "samples",
                "mean_F1",
                "pooled_TP",
                "pooled_FP",
                "pooled_FN",
                "pooled_F1",
                "mean_L1_union_pp",
                "mean_L1_truth_only_pp",
                "mean_Pearson_union",
                "mean_Pearson_truth_only",
            ],
        )
    write_tsv(
        results / "holdout_gap_action_plan_audit.tsv",
        [
            {
                "metric": "release_gate_status",
                "value": "pass=15;fail=0;expected_gap=2",
                "evidence": "fixture",
                "decision": "current_default_supported_with_expected_release_gaps",
            }
        ],
        ["metric", "value", "evidence", "decision"],
    )
    write_tsv(
        results / "marine_binomial_transfer_audit.tsv",
        [
            {
                "metric": "all_samples_min_mapped_pct_bacteria_archaea",
                "value": "91.005806",
                "evidence": "fixture",
                "decision": "below_release_threshold",
            }
        ],
        ["metric", "value", "evidence", "decision"],
    )
    write_tsv(
        results / "abundance_release_blocker_audit.tsv",
        [
            {
                "metric": "feature_allocator_refined_score_replay",
                "value": "samples=38;switched=33;improved=33;worsened=0",
                "evidence": "fixture",
                "decision": "refined_guard_score_replay_pass_candidate_needs_independent_holdout",
            },
            {
                "metric": "feature_allocator_refined_marine_diagnostic",
                "value": "samples=2;switched=2;improved=2;worsened=0",
                "evidence": "fixture",
                "decision": "diagnostic_supports_refined_allocator_holdout_candidate",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )
    write_tsv(
        results / "universal_strategy_release_gate.tsv",
        [
            {
                "gate_id": "summary",
                "requirement": "summary",
                "observed": "pass=15;fail=0;expected_gap=2",
                "status": "pass",
                "expected_gap": "True",
                "evidence": "fixture",
                "decision": "current_default_supported_with_expected_release_gaps",
            }
        ],
        ["gate_id", "requirement", "observed", "status", "expected_gap", "evidence", "decision"],
    )
    write_tsv(
        results / "tmp_headroom_audit.tsv",
        [
            {"metric": "next_sample", "value": "26", "unit": "sample_id"},
            {
                "metric": "tmp_free_GiB",
                "value": "16.862" if headroom_pass else "2.403",
                "unit": "GiB",
            },
            {"metric": "estimated_required_free_GiB", "value": "6.503", "unit": "GiB"},
            {
                "metric": "headroom_pass",
                "value": "True" if headroom_pass else "False",
                "unit": "boolean",
            },
            {"metric": "alternate_free_GiB", "value": "10.217", "unit": "GiB"},
            {"metric": "alternate_headroom_pass", "value": "False", "unit": "boolean"},
            {
                "metric": "recommended_action",
                "value": "run_next_sample" if headroom_pass else "free_tmp_space_before_next_sample",
                "unit": "decision",
            },
        ],
        ["metric", "value", "unit"],
    )
    write_tsv(
        results / "next_release_grade_actions.tsv",
        [
            {
                "rank": "1",
                "action_id": "restore_profile_hmp_airskin_sample26",
                "action_type": "release_grade_holdout_expansion",
                "status": "ready_requires_network"
                if headroom_pass
                else "needs_tmp_cleanup_before_network_run",
                "evidence": "fixture",
                "reason": "fixture",
                "artifact": "results/hmp_airskin_next_release_sample_commands.sh",
            }
        ],
        ["rank", "action_id", "action_type", "status", "evidence", "reason", "artifact"],
    )


def seed_marine_candidate_audit(results: Path) -> None:
    write_tsv(
        results / "marine_truth_upgrade_candidate_audit.tsv",
        [
            {
                "metric": "diagnostic_unique_epithet_rescues",
                "value": (
                    "rescued_rows=51;rescued_mass_pct_all=3.146300;"
                    "min_mapped_pct_bacteria_archaea=91.292431;ready_samples=2/10"
                ),
                "evidence": "fixture",
                "decision": "candidate_rules_do_not_clear_release_threshold",
            },
            {
                "metric": "promotion_decision",
                "value": "do_not_promote_marine_truth_upgrade",
                "evidence": "fixture",
                "decision": "truth_ambiguity_not_resolved_by_current_candidate_rules",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def seed_marine_source_readmap_audit(results: Path) -> None:
    write_tsv(
        results / "marine_source_readmap_feasibility_audit.tsv",
        [
            {
                "metric": "sampled_readmap_rows",
                "value": "6000000",
                "evidence": "fixture",
                "decision": "source_readmaps_available_in_local_archives",
            },
            {
                "metric": "direct_sequence_to_gtdb_matches",
                "value": "sequence_accession_rows=0;core_rows=0",
                "evidence": "fixture",
                "decision": "no_direct_contig_to_gtdb_assembly_mapping",
            },
            {
                "metric": "unresolved_transfer_rows_sampled",
                "value": "314110",
                "evidence": "fixture",
                "decision": "readmaps_cover_unresolved_taxids_but_need_source_mapping",
            },
            {
                "metric": "promotion_decision",
                "value": "do_not_promote_marine_source_readmap_from_local_ids",
                "evidence": "fixture",
                "decision": "contig_or_otu_to_assembly_mapping_missing",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def seed_marine_source_mapping_inventory(results: Path) -> None:
    write_tsv(
        results / "marine_source_mapping_local_inventory_audit.tsv",
        [
            {
                "metric": "local_marine_crosswalk_files_missing",
                "value": (
                    "local_marine_genome_to_id.tsv,"
                    "local_marine_gsa_pooled_mapping.tsv.gz,local_marine_metadata.tsv"
                ),
                "evidence": "fixture",
                "decision": "source_crosswalk_not_found",
            },
            {
                "metric": "archive_relevant_member_status",
                "value": "readmap_only_archives=3",
                "evidence": "fixture",
                "decision": "archives_do_not_embed_source_crosswalk",
            },
            {
                "metric": "promotion_decision",
                "value": "do_not_promote_marine_from_local_inventory",
                "evidence": "fixture",
                "decision": "local_source_mapping_absent",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def seed_marine_setup_metadata_audit(results: Path) -> None:
    write_tsv(
        results / "marine_setup_metadata_truth_upgrade_audit.tsv",
        [
            {
                "metric": "setup_metadata_inputs",
                "value": (
                    "genome_to_id_exists=True;metadata_exists=True;setup_genomes=977"
                ),
                "evidence": "fixture",
                "decision": "setup_metadata_available",
            },
            {
                "metric": "exact_source_name_match_summary",
                "value": (
                    "setup_genomes=977;exact_unique=328;exact_ambiguous=8;no_match=641"
                ),
                "evidence": "fixture",
                "decision": "source_names_partially_resolve_gtdb_species",
            },
            {
                "metric": "strict_all_setup_sources_unique_summary",
                "value": (
                    "min_mapped_pct_bacteria_archaea=93.043000;ready_samples=3/10"
                ),
                "evidence": "fixture",
                "decision": "release_rule_fails_threshold",
            },
            {
                "metric": "diagnostic_partial_setup_unique_source_name_summary",
                "value": (
                    "min_mapped_pct_bacteria_archaea=94.810000;ready_samples=9/10"
                ),
                "evidence": "fixture",
                "decision": "diagnostic_not_release_rule",
            },
            {
                "metric": "promotion_decision",
                "value": "do_not_promote_marine_setup_metadata_truth_upgrade",
                "evidence": "fixture",
                "decision": "strict_setup_source_rule_below_release_threshold",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def seed_marine_setup_metadata_profiles_needed_audit(results: Path) -> None:
    write_tsv(
        results / "marine_setup_metadata_truth_upgrade_audit.tsv",
        [
            {
                "metric": "setup_metadata_inputs",
                "value": (
                    "genome_to_id_exists=True;metadata_exists=True;setup_genomes=977"
                ),
                "evidence": "fixture",
                "decision": "setup_metadata_available",
            },
            {
                "metric": "exact_source_name_match_summary",
                "value": (
                    "setup_genomes=977;exact_unique=300;exact_ambiguous=8;no_match=669"
                ),
                "evidence": "fixture",
                "decision": "source_names_partially_resolve_gtdb_species",
            },
            {
                "metric": "assembly_summary_match_summary",
                "value": (
                    "setup_genomes=977;assembly_unique=472;assembly_ambiguous=6;assembly_no_match=499"
                ),
                "evidence": "fixture",
                "decision": "assembly_summary_resolves_additional_sources",
            },
            {
                "metric": "strict_source_name_or_assembly_sources_unique_summary",
                "value": (
                    "min_mapped_pct_bacteria_archaea=95.029253;ready_samples=10/10"
                ),
                "evidence": "fixture",
                "decision": "release_rule_passes_threshold",
            },
            {
                "metric": "diagnostic_partial_setup_unique_source_name_summary",
                "value": (
                    "min_mapped_pct_bacteria_archaea=94.809534;ready_samples=9/10"
                ),
                "evidence": "fixture",
                "decision": "diagnostic_not_release_rule",
            },
            {
                "metric": "promotion_decision",
                "value": "strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed",
                "evidence": "fixture",
                "decision": "truth_mapping_threshold_cleared_but_profiles_need_same_namespace_rescore",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def seed_marine_profile_rescore_negative_audit(results: Path) -> None:
    write_tsv(
        results / "marine_setup_truth_profile_rescore_audit.tsv",
        [
            {
                "metric": "scoring_scope",
                "value": (
                    "samples=3,4,5;rule_id=strict_source_name_or_assembly_sources_unique;"
                    "profile_set_label=selected_default_same_namespace"
                ),
                "evidence": "fixture",
                "decision": "same_namespace_truth_rule_applied",
            },
            {
                "metric": "truth_rule_quality",
                "value": "ready_samples=3/3;min_mapped_pct_bacteria_archaea=96.463685",
                "evidence": "fixture",
                "decision": "truth_rule_passes_threshold",
            },
            {
                "metric": "minco_vs_sylph_summary",
                "value": (
                    "minco_mean_F1=0.781358;sylph_mean_F1=0.838716;"
                    "minco_mean_L1_union_pp=41.299694;sylph_mean_L1_union_pp=40.997305"
                ),
                "evidence": "fixture",
                "decision": "sylph_higher_F1",
            },
            {
                "metric": "promotion_decision",
                "value": "selected_default_profile_scored_review_metrics",
                "evidence": "fixture",
                "decision": "route_profile_evidence_ready_for_review",
            },
        ],
        ["metric", "value", "evidence", "decision"],
    )


def test_next_evidence_routes_rank_cami3_recovery_first(monkeypatch, tmp_path: Path) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert [row["route_id"] for row in built] == [
        "cami3_source_readmap_samples3_5",
        "hmp_airskin_omitted_samples_2_8_12_26_27",
        "marine_truth_upgrade",
        "refined_allocator_candidate",
        "current_default_candidate",
    ]
    assert built[0]["decision"] == "top_route_requires_readmap_restore_only"
    assert built[0]["ready_without_external_restore"] == "false"
    assert audit["top_next_evidence_route"]["value"] == "cami3_source_readmap_samples3_5"
    assert audit["routes_ready_without_external_restore"]["value"] == "current_default_candidate"
    assert audit["routes_ready_without_external_restore"]["decision"] == (
        "only_current_default_audit_ready_locally"
    )
    assert audit["routes_requiring_external_restore"]["decision"] == (
        "expected_gap_requires_input_recovery_not_threshold_tuning"
    )


def test_next_evidence_routes_after_headroom_cleanup_ready_for_network(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=True)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["route_id"] == "hmp_airskin_omitted_samples_2_8_12_26_27"
    assert built[0]["local_status"] == "truth_present_profiles_missing"
    assert built[0]["decision"] == "not_local_without_reads_or_profiles"
    assert built[0]["missing_evidence"] == (
        "minco_profile=0/5;sylph_profile=0/5;read_input=0/5;"
        "scratch_headroom=run_next_sample"
    )
    assert "headroom_pass=True" in str(built[0]["current_evidence"])
    assert audit["top_route_blocker_detail"]["value"] == (
        "minco_profile=0/5;sylph_profile=0/5;read_input=0/5;"
        "scratch_headroom=run_next_sample"
    )
    assert audit["top_route_blocker_detail"]["decision"] == (
        "not_local_without_reads_or_profiles"
    )
    assert audit["default_decision"]["decision"] == "keep_goal_active"


def test_next_evidence_routes_counts_restored_omitted_sample(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=True, omitted_ready_pairs=1)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["missing_evidence"] == (
        "minco_profile=1/5;sylph_profile=1/5;read_input=1/5;"
        "scratch_headroom=run_next_sample"
    )
    assert "minco_profile=1/5;sylph_profile=1/5;ready_pairs=1/5" in str(
        built[0]["current_evidence"]
    )
    assert audit["top_route_blocker_detail"]["value"] == (
        "minco_profile=1/5;sylph_profile=1/5;read_input=1/5;"
        "scratch_headroom=run_next_sample"
    )


def test_next_evidence_routes_after_cami3_completion_moves_to_hmp(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert [row["route_id"] for row in built] == [
        "hmp_airskin_omitted_samples_2_8_12_26_27",
        "marine_truth_upgrade",
        "refined_allocator_candidate",
        "cami3_source_readmap_samples3_5_completed",
        "current_default_candidate",
    ]
    assert audit["completed_evidence_routes"]["value"] == (
        "cami3_source_readmap_samples3_5_completed"
    )
    assert audit["completed_evidence_routes"]["decision"] == "completed_routes_scored_negative"
    assert audit["top_next_evidence_route"]["value"] == (
        "hmp_airskin_omitted_samples_2_8_12_26_27"
    )
    assert audit["top_next_evidence_route"]["decision"] == "not_local_without_reads_or_headroom"
    assert audit["top_route_blocker_detail"]["value"] == (
        "minco_profile=0/5;sylph_profile=0/5;read_input=0/5;"
        "scratch_headroom=free_tmp_space_before_next_sample"
    )
    assert audit["top_route_blocker_detail"]["decision"] == (
        "not_local_without_reads_or_headroom"
    )
    assert audit["routes_ready_without_external_restore"]["value"] == (
        "cami3_source_readmap_samples3_5_completed,current_default_candidate"
    )
    assert audit["routes_requiring_external_restore"]["value"] == (
        "hmp_airskin_omitted_samples_2_8_12_26_27,marine_truth_upgrade,refined_allocator_candidate"
    )
    assert audit["routes_requiring_external_restore"]["decision"] == (
        "expected_gap_requires_input_recovery_not_threshold_tuning"
    )


def test_next_evidence_routes_completed_hmp_omitted_route_moves_to_marine(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert [row["route_id"] for row in built] == [
        "marine_truth_upgrade",
        "refined_allocator_candidate",
        "hmp_airskin_omitted_samples_2_8_12_26_27",
        "cami3_source_readmap_samples3_5_completed",
        "current_default_candidate",
    ]
    completed = audit["completed_evidence_routes"]["value"]
    assert completed == (
        "hmp_airskin_omitted_samples_2_8_12_26_27,"
        "cami3_source_readmap_samples3_5_completed"
    )
    assert audit["completed_evidence_routes"]["decision"] == "completed_routes_scored_negative"
    assert audit["top_next_evidence_route"]["value"] == "marine_truth_upgrade"
    assert audit["top_next_evidence_route"]["decision"] == "truth_upgrade_needed_before_release_use"
    assert audit["routes_ready_without_external_restore"]["value"] == (
        "hmp_airskin_omitted_samples_2_8_12_26_27,"
        "cami3_source_readmap_samples3_5_completed,current_default_candidate"
    )
    assert audit["routes_requiring_external_restore"]["value"] == (
        "marine_truth_upgrade,refined_allocator_candidate"
    )


def test_next_evidence_routes_marine_candidate_rules_need_external_truth(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)
    seed_marine_candidate_audit(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["route_id"] == "marine_truth_upgrade"
    assert built[0]["local_status"] == "truth_transfer_candidate_rules_negative"
    assert built[0]["decision"] == "truth_upgrade_needs_external_truth_source"
    assert "min_mapped_pct_bacteria_archaea=91.292431" in str(
        built[0]["current_evidence"]
    )
    assert "requires accession-level/source-specific truth mapping" in str(
        built[0]["missing_evidence"]
    )
    assert audit["top_next_evidence_route"]["decision"] == (
        "truth_upgrade_needs_external_truth_source"
    )


def test_next_evidence_routes_marine_source_readmaps_need_mapping(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)
    seed_marine_candidate_audit(results)
    seed_marine_source_readmap_audit(results)
    seed_marine_source_mapping_inventory(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["route_id"] == "marine_truth_upgrade"
    assert built[0]["local_status"] == (
        "truth_transfer_source_mapping_absent_locally"
    )
    assert built[0]["decision"] == "truth_upgrade_needs_source_mapping"
    assert "source_readmap_rows=6000000" in str(built[0]["current_evidence"])
    assert "unresolved_rows_sampled=314110" in str(built[0]["current_evidence"])
    assert "readmap_only_archives=3" in str(built[0]["current_evidence"])
    assert "source IDs do not directly map to GTDB assemblies" in str(
        built[0]["missing_evidence"]
    )
    assert "contig/OTU-to-assembly metadata" in str(built[0]["missing_evidence"])
    assert "read archives are readmap-only" in str(built[0]["missing_evidence"])
    assert audit["top_next_evidence_route"]["decision"] == (
        "truth_upgrade_needs_source_mapping"
    )


def test_next_evidence_routes_marine_setup_metadata_still_incomplete(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)
    seed_marine_candidate_audit(results)
    seed_marine_source_readmap_audit(results)
    seed_marine_source_mapping_inventory(results)
    seed_marine_setup_metadata_audit(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["route_id"] == "marine_truth_upgrade"
    assert built[0]["local_status"] == "truth_transfer_setup_metadata_incomplete"
    assert built[0]["decision"] == (
        "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout"
    )
    assert "setup_genomes=977;exact_unique=328" in str(built[0]["current_evidence"])
    assert "release-grade strict setup-source rule remains below 95%" in str(
        built[0]["missing_evidence"]
    )
    assert "partial source-name rule is diagnostic only" in str(
        built[0]["missing_evidence"]
    )
    assert audit["top_next_evidence_route"]["decision"] == (
        "truth_upgrade_needs_stronger_source_mapping_or_clean_holdout"
    )


def test_next_evidence_routes_marine_setup_metadata_threshold_passes_profiles_needed(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)
    seed_marine_candidate_audit(results)
    seed_marine_source_readmap_audit(results)
    seed_marine_source_mapping_inventory(results)
    seed_marine_setup_metadata_profiles_needed_audit(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    assert built[0]["route_id"] == "marine_truth_upgrade"
    assert built[0]["local_status"] == "truth_transfer_setup_metadata_ready_profiles_needed"
    assert built[0]["decision"] == "truth_upgrade_profiles_needed"
    assert "assembly_unique=472" in str(built[0]["current_evidence"])
    assert "min_mapped_pct_bacteria_archaea=95.029253" in str(
        built[0]["current_evidence"]
    )
    assert "selected-default same-namespace profiles still need scoring" in str(
        built[0]["missing_evidence"]
    )
    assert audit["top_next_evidence_route"]["decision"] == "truth_upgrade_profiles_needed"


def test_next_evidence_routes_marine_selected_default_scored_negative(
    monkeypatch, tmp_path: Path
) -> None:
    results = patch_paths(monkeypatch, tmp_path)
    seed_route_inputs(results, after_scored=True, headroom_pass=False, omitted_ready_pairs=5)
    seed_marine_candidate_audit(results)
    seed_marine_source_readmap_audit(results)
    seed_marine_source_mapping_inventory(results)
    seed_marine_setup_metadata_profiles_needed_audit(results)
    seed_marine_profile_rescore_negative_audit(results)

    built = routes.build_routes()
    audit = {row["metric"]: row for row in routes.build_audit(built)}

    marine = next(row for row in built if row["route_id"] == "marine_truth_upgrade")
    assert marine["local_status"] == "scored_negative"
    assert marine["decision"] == "completed_negative_holdout_result"
    assert marine["ready_without_external_restore"] == "true"
    assert "minco_mean_F1=0.781358" in str(marine["current_evidence"])
    assert marine["missing_evidence"] == "none"
    assert audit["top_next_evidence_route"]["value"] == "refined_allocator_candidate"
    assert audit["top_next_evidence_route"]["decision"] == (
        "keep_opt_in_until_independent_release_holdout_passes"
    )
    assert audit["routes_requiring_external_restore"]["value"] == "refined_allocator_candidate"
