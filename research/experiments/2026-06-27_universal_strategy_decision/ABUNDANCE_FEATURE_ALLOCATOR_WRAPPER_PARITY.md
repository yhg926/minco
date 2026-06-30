# Feature Allocator Wrapper Parity

Date: 2026-06-29

This note validates the implemented experimental switch
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002`
on cached selected-candidate profiles. It is not a raw-read rerun and
does not change the default.
The validator passes accession/taxmap labels to the allocator for genus
grouping when the cached profile records a candidate-surface taxmap sidecar.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `baseline_validation_max_abs_delta` | 1.14930287509e-12 | pass |
| `switched_samples` | 30 | wrapper_guard_applied |
| `mean_L1_delta_pp` | -0.318866 | abundance_delta |
| `worsened_samples` | 0 | sample_safety |
| `max_worse_L1_delta_pp` | 0.000000 | sample_safety |
| `promotion_decision` | candidate_needs_independent_holdout | candidate_needs_independent_holdout |

## Panel Summary

| Panel | Samples | Switched | Mean L1 delta pp | Max worse pp | Improved | Worsened |
|---|---:|---:|---:|---:|---:|---:|
| cami2_toy_mouse_gut | 3 | 3 | -0.187921 | -0.030594 | 3 | 0 |
| cami3_toy_human_gut_gtdb_source_readmap | 3 | 1 | -0.078847 | 0.000000 | 1 | 0 |
| hmp_airskin_gtdb_source_abundance | 24 | 24 | -0.381306 | -0.041632 | 24 | 0 |
| hmp_gastrooral_gtdb_source_abundance | 2 | 2 | -0.126023 | -0.012398 | 2 | 0 |
| all | 32 | 30 | -0.318866 | 0.000000 | 30 | 0 |

## Decision

- Keep the selected default unchanged.
- Treat this implemented switch as an independent-holdout candidate.
- The next validation step is a same-namespace holdout replay before
  considering any default promotion.

## Outputs

- `results/feature_allocator_wrapper_parity_scores.tsv`
- `results/feature_allocator_wrapper_parity_summary.tsv`
- `results/feature_allocator_wrapper_parity_audit.tsv`
