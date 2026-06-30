# Selected Call-Set Mass Transform Sweep

Date: 2026-06-29

This tracked note is generated from cached selected-candidate profile
TSVs. It keeps calls fixed, preserves candidate-added row mass, and
sweeps small base-called-row abundance transforms that are simple enough
to implement in the wrapper.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `variants_tested` | 14 | cached_selected_candidate_callset |
| `baseline_validation_max_abs_delta` | 1.14930287509e-12 | pass |
| `sample_safe_variants` | 0 | none |
| `panel_safe_variants` | 5 | panel_mean_only |
| `best_ranked_variant` | base_genus_xny_a0.10;mean_sample_delta=-0.382481;worsened_samples=4;max_worse=0.161705 | informational |
| `promotion_decision` | do_not_promote_selected_mass_transform | do_not_promote_selected_mass_transform |

## Top Variants

| Method | Mean sample L1 delta pp | Worsened samples | Max sample worse pp | Panel mean delta pp | Worsened panels |
|---|---:|---:|---:|---:|---:|
| base_genus_xny_a0.10 | -0.382481 | 4 | 0.161705 | -0.259783 | 0 |
| base_genus_xny_a0.05 | -0.205512 | 4 | 0.061795 | -0.141597 | 0 |
| base_genus_xny_a0.02 | -0.085030 | 4 | 0.018825 | -0.058311 | 0 |
| base_global_xny_a0.10 | -0.728598 | 5 | 0.825999 | -0.554313 | 1 |
| base_global_xny_a0.05 | -0.383641 | 5 | 0.274367 | -0.308378 | 0 |
| base_global_xny_a0.02 | -0.161870 | 5 | 0.063192 | -0.140380 | 0 |
| base_global_probability_a0.02 | -0.079584 | 5 | 0.088479 | -0.055133 | 1 |
| base_global_probability_a0.05 | -0.192356 | 6 | 0.230394 | -0.129249 | 1 |

## Decision

- Do not promote these simple mass transforms into the default.
- The best ranked transform, `base_genus_xny_a0.10`, improves mean
  L1 but regresses 4/32 samples.
- This supports the current next target: a more sample-aware
  matched-call allocator, not another unguarded panel-mean transform.

## Outputs

- `results/selected_call_mass_transform_scores.tsv`
- `results/selected_call_mass_transform_panel_delta.tsv`
- `results/selected_call_mass_transform_overall.tsv`
- `results/selected_call_mass_transform_audit.tsv`
