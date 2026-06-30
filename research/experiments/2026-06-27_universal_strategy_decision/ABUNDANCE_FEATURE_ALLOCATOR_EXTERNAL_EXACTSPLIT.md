# Feature Allocator External Exact-Split Stress Test

Date: 2026-06-29

This cached-profile replay applies the implemented experimental
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002`
to plant-associated samples 3-5 and strainmadness samples 0-2.
The call set is fixed; only abundance mass is adjusted.

These panels are diagnostic, not release-grade GTDB holdouts, because
the scorers use local bacteria-scope CAMI truth namespaces.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `baseline_validation` | counts=0;F1=5.55e-17;L1=2.22e-16;Pearson=2.22e-16 | matches_cached_exactsplit_scores |
| `candidate_effect` | mean_delta_l1=-0.004134;improved_datasets=2;worsened_datasets=0;improved_samples=4;worsened_samples=2;max_sample_worse=0.000249;switched_samples=6;adjusted_rows=130 | external_stress_rejects_default |
| `promotion_decision` | keep_feature_allocator_experimental_off_by_default | current_default_unchanged |

## Summary

| Dataset | Method | Samples | Mean F1 | Mean L1 | Delta L1 | Worsened L1 samples | Switched | Adjusted rows |
|---|---|---|---:|---:|---:|---:|---:|---:|
| plant_holdout | current_exactsplit_abundance | 3,4,5 | 0.612016 | 0.822986 | 0.000000 | 0 | 0 | 0 |
| plant_holdout | guarded_feature_allocator_exactsplit | 3,4,5 | 0.612016 | 0.822966 | -0.000021 | 2 | 3 | 103 |
| strainmadness | current_exactsplit_abundance | 0,1,2 | 0.628737 | 0.759779 | 0.000000 | 0 | 0 | 0 |
| strainmadness | guarded_feature_allocator_exactsplit | 0,1,2 | 0.628737 | 0.751531 | -0.008248 | 0 | 3 | 27 |

## Decision

- Keep the current default unchanged.
- Leave the guarded feature allocator experimental/off by default.
- It remains useful as a candidate for a future same-namespace
  release-grade holdout because this replay tests only diagnostic
  external exact-split panels.

## Outputs

- `results/feature_allocator_external_exactsplit_scores.tsv`
- `results/feature_allocator_external_exactsplit_summary.tsv`
- `results/feature_allocator_external_exactsplit_validation.tsv`
- `results/feature_allocator_external_exactsplit_audit.tsv`
