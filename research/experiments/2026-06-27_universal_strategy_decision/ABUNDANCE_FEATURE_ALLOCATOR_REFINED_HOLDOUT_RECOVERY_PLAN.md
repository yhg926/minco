# Refined Allocator Release Holdout Recovery Plan

Date: 2026-06-29

This generated note ranks local recovery routes for the refined guarded
abundance allocator. It does not promote the allocator or change the
selected default.

## Decision

- Top route: `cami3_source_readmap_extension_3_5` (completed_negative_holdout_result).
- Completed local step: `post_recovery_scored=True;candidate_vs_sylph=minco_F1_wins=0/3;minco_L1_wins=0/3;refined_effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000`.
- Remaining blocker: `extension_result_does_not_support_default_promotion`.
- Source-readmap recovery inputs: `post_recovery_scored=True;candidate_vs_sylph=minco_F1_wins=0/3;minco_L1_wins=0/3;refined_effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000`.
- Post-recovery scorer: `readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3;blocker=none;promotion=post_recovery_extension_scored;decision=review_holdout_result_before_default_change`.
- Tested local substitute: `samples=3,4,5;min_mapped_profile_scope_pct=62.728832;ready_samples=;improved=0;worsened=0;mean_L1_delta_pp=-0.000000000` (keep_refined_allocator_opt_in).
- Replay runtime: `samples=3;exit_statuses=0,0,0;max_rss_kb=1601816`.
- Default decision: `do_not_promote_refined_allocator_without_independent_release_holdout`.

## Ranked Routes

| Rank | Route | Samples | Status | Next Action |
|---:|---|---|---|---|
| 10 | `cami3_source_readmap_extension_3_5` | 3,4,5 | post_recovery_scored=True;candidate_vs_sylph=minco_F1_wins=0/3;minco_L1_wins=0/3;refined_effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000 | Do not promote the refined allocator from this route; the independent extension result has no refined-allocator gain and favors the external baseline on F1 and abundance. Move to the next independent holdout route. |
| 20 | `hmp_airskin_omitted_profile_pair_recovery` | 2,8,12,26,27 | truth=5/5;minco=0/5;sylph=0/5;read_input=0/5 | Recover read inputs, then run same-release selected-default MinCO and baseline profiles for the omitted sample IDs. |
| 25 | `cami3_source_profile_binomial_extension` | 3,4,5 | samples=3,4,5;min_mapped_profile_scope_pct=62.728832;ready_samples=;improved=0;worsened=0;mean_L1_delta_pp=-0.000000000 | Do not use the taxonomic-profile source rows as a release substitute; continue with per-read source-readmap recovery. |
| 30 | `marine_diagnostic_truth_upgrade` | marine0,marine2 | samples=marine0,marine2;min_truth_mass_mapped_pct_bacteria_archaea=74.120670;samples=2;switched=2;improved=2;worsened=0;mean_L1_delta=-0.000521634;max_worse_L1=-0.000503001;mean_F1_delta=0.000000000 | Improve the GTDB species truth transfer before using marine as a release-grade abundance holdout. |
| 40 | `new_clean_gtdb_holdout` |  | release_grade_independent_holdout_missing | Add a new clean GTDB truth panel only if CAMI3 or HMP recovery cannot produce independent release-grade profile pairs. |

## Outputs

- `results/feature_allocator_refined_holdout_recovery_plan.tsv`
- `results/feature_allocator_refined_holdout_recovery_runtime.tsv`
- `results/feature_allocator_refined_holdout_recovery_audit.tsv`
