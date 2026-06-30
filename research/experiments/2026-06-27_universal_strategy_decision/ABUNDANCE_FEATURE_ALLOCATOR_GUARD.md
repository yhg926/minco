# Selected Call-Set Feature Allocator Guard Audit

Date: 2026-06-29

This tracked note is generated from cached selected-candidate profiles
and the selected-call feature allocator sweep. It tests whether
output-derived guards can apply a feature allocator only on samples
where it is stable.
Cached score columns (`selected_F1`, `selected_L1_pp`, and
`selected_Pearson`) are excluded from guard rule features.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `tested_guard_rules` | 5687 | truth_used_only_for_scoring |
| `sample_safe_guards` | 1990 | strict_sample_gate |
| `panel_safe_guards` | 5601 | panel_gate |
| `best_mean_guard` | genus_raw_mean_breadth_a0.02;single;max_species_per_genus>=5;mean_delta=-0.323744;worsened_samples=2;max_worse=0.065742;switched=31 | tradeoff |
| `best_sample_safe_guard` | genus_raw_hit_breadth_a0.02;single;base_mass_multi_genus_frac>=0.375015;mean_delta=-0.312324;switched=28 | candidate |
| `lopo_guard_result` | mean_delta=-0.146115;holdouts_with_mean_regression=0;holdouts_with_sample_regression=0;max_sample_worse=0.000000 | leave_one_panel_out_validation |
| `promotion_decision` | candidate_guard_needs_wrapper_and_independent_holdout | candidate_guard_needs_wrapper_and_independent_holdout |

## Top Rules

| Method | Rule | Mean sample delta pp | Worsened samples | Max worse pp | Switched |
|---|---|---:|---:|---:|---:|
| genus_raw_mean_breadth_a0.02 | n_candidate_calls>=0 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | candidate_call_frac>=0 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | candidate_raw_frac>=0 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | u_zip_median<=1 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | probability_min>=0 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | reported_ani_median<=1 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | n_candidate_calls>=0 AND candidate_call_frac>=0 | -0.328375 | 2 | 0.065742 | 32 |
| genus_raw_mean_breadth_a0.02 | n_candidate_calls>=0 AND candidate_raw_frac>=0 | -0.328375 | 2 | 0.065742 | 32 |

## LOPO

| Holdout | Method | Selection | Mean delta pp | Worsened samples | Switched |
|---|---|---|---:|---:|---:|
| cami2_toy_mouse_gut | genus_raw_hit_breadth_a0.02 | strict_train_sample_safe | -0.084767 | 0 | 1 |
| cami3_toy_human_gut_gtdb_source_readmap | genus_raw_hit_breadth_a0.02 | strict_train_sample_safe | 0.000000 | 0 | 0 |
| hmp_airskin_gtdb_source_abundance | genus_raw_hit_breadth_a0.02 | strict_train_sample_safe | -0.373669 | 0 | 22 |
| hmp_gastrooral_gtdb_source_abundance | genus_raw_hit_breadth_a0.02 | strict_train_sample_safe | -0.126023 | 0 | 2 |

## Decision

- This is a diagnostic guard search, not a promoted allocator.
- A guard is promotable only if it remains sample-safe in
  leave-one-panel-out validation of top-ranked rules and then survives
  independent holdouts.

## Outputs

- `results/selected_call_feature_allocator_guard_features.tsv`
- `results/selected_call_feature_allocator_guard_rules.tsv` top 1000 scored rules
- `results/selected_call_feature_allocator_guard_lopo.tsv`
- `results/selected_call_feature_allocator_guard_audit.tsv`
