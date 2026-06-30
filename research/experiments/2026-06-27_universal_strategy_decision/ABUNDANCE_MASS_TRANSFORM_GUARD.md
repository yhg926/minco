# Selected Call-Set Mass Transform Guard Audit

Date: 2026-06-29

This tracked note is generated from cached selected-candidate profiles and
the selected-call mass-transform sweep. It tests whether output-derived
guards can apply a simple transform only on samples where it is stable.
Cached score columns (`selected_F1`, `selected_L1_pp`, and
`selected_Pearson`) are excluded from guard rule features.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `tested_guard_rules` | 5659 | truth_used_only_for_scoring |
| `sample_safe_guards` | 822 | strict_sample_gate |
| `panel_safe_guards` | 3600 | panel_gate |
| `best_mean_guard` | base_global_xny_a0.10;single;max_species_per_genus>=5;mean_delta=-0.738736;worsened_samples=4;max_worse=0.825999;switched=31 | tradeoff |
| `best_sample_safe_guard` | base_global_xny_a0.10;single;n_calls>=71.4;mean_delta=-0.491859;switched=13 | candidate |
| `lopo_guard_result` | mean_delta=-0.129226;holdouts_with_mean_regression=1;holdouts_with_sample_regression=2;max_sample_worse=0.825999 | leave_one_panel_out_validation |
| `promotion_decision` | sample_safe_in_panel_but_fails_lopo | sample_safe_in_panel_but_fails_lopo |

## Top Rules

| Method | Rule | Mean sample delta pp | Worsened samples | Max worse pp | Switched |
|---|---|---:|---:|---:|---:|
| base_global_xny_a0.10 | max_species_per_genus>=5 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND n_candidate_calls>=0 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND candidate_call_frac>=0 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND candidate_raw_frac>=0 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND u_zip_median<=1 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND probability_min>=0 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | max_species_per_genus>=5 AND reported_ani_median<=1 | -0.738736 | 4 | 0.825999 | 31 |
| base_global_xny_a0.10 | n_candidate_calls>=0 | -0.728598 | 5 | 0.825999 | 32 |

## LOPO

| Holdout | Method | Selection | Mean delta pp | Worsened samples | Switched |
|---|---|---|---:|---:|---:|
| cami2_toy_mouse_gut | base_global_xny_a0.10 | strict_train_sample_safe | 0.207545 | 1 | 3 |
| cami3_toy_human_gut_gtdb_source_readmap | base_global_xny_a0.10 | strict_train_sample_safe | 0.000000 | 0 | 0 |
| hmp_airskin_gtdb_source_abundance | base_global_xny_a0.10 | strict_train_sample_safe | -0.724449 | 2 | 19 |
| hmp_gastrooral_gtdb_source_abundance | base_global_xny_a0.10 | strict_train_sample_safe | 0.000000 | 0 | 0 |

## Decision

- This is a diagnostic guard search, not a promoted allocator.
- A guard is promotable only if it remains sample-safe in
  leave-one-panel-out validation of top-ranked rules and then survives
  independent holdouts.

## Outputs

- `results/selected_call_mass_transform_guard_features.tsv`
- `results/selected_call_mass_transform_guard_rules.tsv` top 1000 scored rules
- `results/selected_call_mass_transform_guard_lopo.tsv`
- `results/selected_call_mass_transform_guard_audit.tsv`
