# Universal Strategy Release Gate

Date: 2026-06-29

This generated gate validates the current default-strategy decision from
repo-local code, documentation, and cached benchmark summaries. It does
not rerun raw profiling jobs.

## Summary

- Status: `pass`
- Observed: `pass=15;fail=0;expected_gap=2`
- Decision: `current_default_supported_with_expected_release_gaps`

The expected gaps are intentional: the current default is supported as
the best cached MinCO default, but the release-grade holdout bundle and
broad abundance claim are not complete enough for a broad universal or
Sylph-beating claim.

## Gates

| Gate | Status | Expected gap | Observed | Decision |
|---|---|---:|---|---|
| default_launcher_entrypoint | pass | False | /home/ubuntu/yihuiguang/tools/KSSD3mini/scripts/minco_profile | entrypoint_present |
| default_launcher_strategy | pass | False | universal-auto-exact | f1_priority_strategy_selected |
| default_profile_preset | pass | False | candidate | candidate_preset_selected |
| default_manifest | pass | False | candidate | manifest_default_matches_code |
| default_contract_test | pass | False | test_default_launcher_constants_match_current_strategy_manifest | default_contract_regression_covered |
| candidate_vs_previous_minco | pass | False | F1_regressions=0;L1_regressions=0;Pearson_regressions=0 | selected_default_supported_vs_previous_minco |
| sylph_claim_boundary | pass | False | F1_wins=1/4;L1_wins=0/4;Pearson_wins=0/4 | broad_sylph_beating_claim_rejected |
| experimental_allocator_default | pass | False | keep_feature_allocator_experimental_off_by_default | allocator_not_promoted |
| calibrated_allocator_default | pass | False | off | experimental_allocator_opt_in_only |
| ani_reporting | pass | False | reporting_ready_not_claim_win | ani_reporting_documented_not_overclaimed |
| docs_claim_boundary | pass | False | README/manual/changelog checked | user_facing_default_boundary_documented |
| release_grade_bundle | expected_gap | True | release_grade_all_panels=false;release_grade_panels=4 | not_release_complete_until_more_clean_gtdb_holdouts |
| cami3_extension_runner_contract | pass | False | readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3;blocker=NA;promotion=post_recovery_extension_scored | post_recovery_runner_scored_or_cleanly_blocked |
| next_evidence_route_contract | pass | False | top=refined_allocator_candidate;completed=marine_truth_upgrade,hmp_airskin_omitted_samples_2_8_12_26_27,cami3_source_readmap_samples3_5_completed;ready_without_external_restore=marine_truth_upgrade,hmp_airskin_omitted_samples_2_8_12_26_27,cami3_source_readmap_samples3_5_completed,current_default_candidate | next_route_ready_no_threshold_sweep_release_claim |
| abundance_release_claim | expected_gap | True | selected_default_L1_losses=4/4;posthoc_blend=do_not_promote_abundance_blend;oracle=prioritize_allocator_plus_call_recovery;feature_allocator=do_not_promote_until_external_sample_safe;refined_score_replay=refined_guard_score_replay_pass_candidate_needs_independent_holdout;marine_diagnostic=diagnostic_supports_refined_allocator_holdout_candidate;holdout_inventory=release_grade_independent_holdout_missing;recovery_plan=do_not_promote_refined_allocator_without_independent_release_holdout;source_profile_extension=keep_refined_allocator_opt_in;source_readmap_inputs=post_recovery_scoring_ready | abundance_expected_gap_not_release_claim |
| runtime_claim_boundary | pass | False | samples=7;seconds_ratio=2.515;rss_ratio=0.191 | speed_claim_rejected_memory_advantage_recorded |
| current_release_decision | pass | False | selected_candidate_default;experimental_allocators_off;do_not_claim_broad_sylph_beating | goal_progress_not_goal_complete |
| summary | pass | True | pass=15;fail=0;expected_gap=2 | current_default_supported_with_expected_release_gaps |

## Outputs

- `results/universal_strategy_release_gate.tsv`
