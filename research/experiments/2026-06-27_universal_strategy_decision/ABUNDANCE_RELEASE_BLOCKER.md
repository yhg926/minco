# Abundance Release Blocker

Date: 2026-06-29

This generated audit validates the abundance claim boundary for the
selected default. It uses cached score and allocator-audit TSVs only.

## Decision

- Candidate versus previous MinCO: `F1=0/4;L1=0/4;Pearson=0/4`.
- Candidate versus external baseline abundance: `L1_losses=4/4;Pearson_losses=4/4`.
- Posthoc allocator status: `do_not_promote_abundance_blend;improved=26;worsened=6;max_worse=0.9934199183847312`.
- Oracle next target: `allocation_only=31/32;call_recovery=1/32;residual=0/32;headroom_pp=23.237189576924983;detected_truth_pct=98.42068215947246`.
- Feature allocator candidate: `samples=38;switched=36;improved=34;worsened=2;mean_L1_delta_pp=-0.269171;max_worse_L1_delta_pp=0.000249;mean_F1_delta=0` (do_not_promote_until_external_sample_safe).
- Refined feature allocator guard: `s_xny_median>=230.3;mean_L1_delta_pp=-0.262955;selected_mean_L1_delta_pp=-0.311474;external_mean_L1_delta_pp=-0.004186;cached_gain_frac=0.976818` (refined_guard_score_replay_pass_candidate_needs_independent_holdout).
- Refined feature allocator score replay: `samples=38;switched=33;improved=33;worsened=0;mean_L1_delta_pp=-0.262955;selected_mean_L1_delta_pp=-0.311474;external_mean_L1_delta_pp=-0.004186;max_worse_L1_delta_pp=0.000000` (pass;refined_guard_score_replay_pass_candidate_needs_independent_holdout).
- Refined feature allocator marine diagnostic: `samples=2;switched=2;improved=2;worsened=0;mean_L1_delta=-0.000521634;max_worse_L1=-0.000503001;mean_F1_delta=0.000000000` (diagnostic_supports_refined_allocator_holdout_candidate).
- Refined feature allocator holdout inventory: `release_panels=4;release_samples=32;selected_cache_overlap=32;independent_release_ready_samples=0;samples=2,8,12,26,27;truth=5/5;minco_profile=3/5;sylph_profile=3/5;ready_pairs=3/5` (release_grade_independent_holdout_missing).
- Refined feature allocator recovery plan: `top_route=cami3_source_readmap_extension_3_5;completed=post_recovery_scored=True;candidate_vs_sylph=minco_F1_wins=0/3;minco_L1_wins=0/3;refined_effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000;remaining=extension_result_does_not_support_default_promotion;post_recovery_runner=readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3;blocker=none;promotion=post_recovery_extension_scored;decision=review_holdout_result_before_default_change` (do_not_promote_refined_allocator_without_independent_release_holdout).
- CAMI3 source-profile extension diagnostic: `samples=3,4,5;min_mapped_profile_scope_pct=62.728832;ready_samples=;improved=0;worsened=0;mean_L1_delta_pp=-0.000000000` (keep_refined_allocator_opt_in).
- CAMI3 source-readmap recovery inputs: `;candidate_profiles=3/3;raw_tables=3/3;sylph_profiles=3/3` (post_recovery_scoring_ready).
- Release decision: `expected_gap_keep_default_no_abundance_claim`.

## Panel Deltas

| Panel | Candidate - Previous L1 pp | Candidate - External L1 pp | Candidate - External F1 | Status |
|---|---:|---:|---:|---|
| cami2_toy_mouse_gut | -0.237646 | 6.798669 | -0.064928 | abundance_gap |
| cami3_toy_human_gut_gtdb_source_readmap | -4.984795 | 25.114394 | 0.082115 | abundance_gap |
| hmp_airskin_gtdb_source_abundance | -0.125430 | 17.056351 | -0.112163 | abundance_gap |
| hmp_gastrooral_gtdb_source_abundance | -2.629049 | 41.294813 | -0.028484 | abundance_gap |

## Outputs

- `results/abundance_release_blocker_panel.tsv`
- `results/abundance_release_blocker_audit.tsv`
