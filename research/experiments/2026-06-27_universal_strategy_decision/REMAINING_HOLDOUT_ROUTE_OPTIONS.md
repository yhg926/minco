# Remaining Holdout Route Options

Date: 2026-06-29

This audit ranks remaining local evidence routes after the HMP omitted
route and CAMI3 samples3-5 route were scored negative. It incorporates
the current setup-metadata truth mapping state for the marine route.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `manifest_release_status` | release=4;nonrelease=4;total=8 | release_gap_remains |
| `local_viable_uncompleted_routes` | none | no_local_clean_holdout_ready |
| `locally_exhausted_routes` | hmp_omitted_samples_completed,marine_selected_default_profile_scored,cami3_samples3_5_completed,plant_local_transfer,strain_local_transfer,mixed_readiness_panel | do_not_repeat_without_new_truth_or_strategy_change |
| `top_remaining_route` | refined_allocator_independent_release_holdout | not_promotable_without_independent_release_holdout |
| `promotion_decision` | keep_current_default_with_expected_gaps | new_strategy_or_new_holdout_needed |

## Routes

| Rank | Route | Local Status | Best Truth % | Decision |
|---:|---|---|---:|---|
| 20 | `refined_allocator_independent_release_holdout` | opt_in_validated_cached_but_not_independent_release | NA | not_promotable_without_independent_release_holdout |
| 30 | `hmp_omitted_samples_completed` | scored_negative | high_mapped_source_abundance_same_release | completed_negative_do_not_repeat_without_strategy_change |
| 35 | `marine_selected_default_profile_scored` | scored_negative | min_mapped_pct_bacteria_archaea=95.029253;ready_samples=10/10 | completed_negative_do_not_repeat_without_strategy_change |
| 40 | `cami3_samples3_5_completed` | scored_negative | source_readmap_binomial_fallback_completed | completed_negative_do_not_repeat_without_strategy_change |
| 50 | `plant_local_transfer` | truth_transfer_far_below_release_threshold | min=19.083472;samples=3;release_ready=0/3;mean=20.565480 | not_viable_release_route_locally |
| 60 | `strain_local_transfer` | truth_transfer_far_below_release_threshold | min=0.477515;samples=3;release_ready=0/3;mean=1.575529 | not_viable_release_route_locally |
| 70 | `mixed_readiness_panel` | mixed_truth_namespaces | NA | not_viable_release_route_locally |
| 80 | `cami3_older_taxid_transfer` | completed_or_superseded_by_source_readmap | samples0-2 promoted; samples3-5 completed negative | superseded_by_source_readmap_evidence |

## Decision

No threshold-only sweep can close the release-holdout gap. If the
top route is `marine_selected_default_profile_rescore`, local truth
mapping has cleared the threshold and the next evidence step is
same-namespace selected-default profile scoring. Otherwise the next
step remains stronger source mapping or a different clean holdout.
The current default therefore stays selected, with the refined
abundance allocator opt-in.
