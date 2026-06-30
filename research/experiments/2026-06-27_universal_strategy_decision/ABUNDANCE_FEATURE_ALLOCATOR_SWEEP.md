# Selected Call-Set Feature Allocator Sweep

Date: 2026-06-29

This tracked note is generated from cached selected-candidate profile
TSVs. It keeps calls fixed, preserves candidate-added row mass, and
tests whether existing MinCO depth/quality output features can serve as
matched-call abundance targets.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `variants_tested` | 264 | cached_selected_candidate_callset |
| `baseline_validation_max_abs_delta` | 1.14930287509e-12 | pass |
| `sample_safe_variants` | 0 | none |
| `panel_safe_variants` | 46 | panel_mean_only |
| `best_ranked_variant` | genus_raw_hit_breadth_a0.02;mean_sample_delta=-0.318875;worsened_samples=1;max_worse=0.040714 | informational |
| `promotion_decision` | do_not_promote_feature_allocator | do_not_promote_feature_allocator |

## Top Variants

| Method | Mean sample L1 delta pp | Worsened samples | Max sample worse pp | Panel mean delta pp | Worsened panels |
|---|---:|---:|---:|---:|---:|
| genus_raw_hit_breadth_a0.02 | -0.318875 | 1 | 0.040714 | -0.193549 | 0 |
| genus_raw_hit_breadth_a0.01 | -0.161540 | 1 | 0.020357 | -0.100199 | 0 |
| genus_raw_hit_breadth_a0.005 | -0.081471 | 1 | 0.010179 | -0.051952 | 0 |
| genus_raw_mean_breadth_a0.02 | -0.328375 | 2 | 0.065742 | -0.187482 | 0 |
| genus_raw_norm_depth_a0.02 | -0.318771 | 2 | 0.062728 | -0.188789 | 0 |
| genus_raw_mean_breadth_a0.01 | -0.168165 | 2 | 0.032871 | -0.097579 | 0 |
| genus_raw_norm_depth_a0.01 | -0.161968 | 2 | 0.031364 | -0.097754 | 0 |
| genus_raw_mean_breadth_a0.005 | -0.084961 | 2 | 0.016436 | -0.050680 | 0 |
| genus_raw_norm_depth_a0.005 | -0.081780 | 2 | 0.015682 | -0.050762 | 0 |
| genus_raw_hit_breadth_a0.1 | -1.392514 | 3 | 0.492959 | -0.804090 | 0 |
| genus_raw_hit_breadth_a0.05 | -0.759858 | 3 | 0.112533 | -0.445332 | 0 |
| genus_raw_realaf_a0.01 | -0.043472 | 3 | 0.009413 | -0.029475 | 0 |

## Decision

- Do not promote a feature-derived allocator unless it is sample-safe and
  then survives external holdouts.
- If a sample-safe variant appears here, treat it as a candidate for
  independent replay rather than a default change.
- The per-sample score table saves the baseline plus the top
  32 ranked methods; the full sweep is
  summarized in the overall and audit TSVs.

## Outputs

- `results/selected_call_feature_allocator_scores.tsv`
- `results/selected_call_feature_allocator_panel_delta.tsv`
- `results/selected_call_feature_allocator_overall.tsv`
- `results/selected_call_feature_allocator_audit.tsv`
