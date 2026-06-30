# Refined Feature Allocator Marine Diagnostic

Date: 2026-06-29

This cached-profile replay evaluates the opt-in
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002-xny230` on CAMI II marine
exact-split profiles that were not used to choose the refined guard.
The truth namespace is conservative GTDB taxid transfer, so this is
diagnostic evidence and not release-grade evidence.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `diagnostic_scope` | samples=marine0,marine2;min_truth_mass_mapped_pct_bacteria_archaea=74.120670 | diagnostic_nonrelease_truth_transfer |
| `baseline_internal_validation` | max_L1_delta_raw_vs_calibrated=6.63358257214e-15 | pass |
| `refined_allocator_marine_effect` | samples=2;switched=2;improved=2;worsened=0;mean_L1_delta=-0.000521634;max_worse_L1=-0.000503001;mean_F1_delta=0.000000000 | independent_diagnostic_effect |
| `promotion_decision` | diagnostic_supports_refined_allocator_holdout_candidate | current_default_unchanged |

## Summary

| Method | Samples | Mean F1 | Mean L1 | Mean Pearson | Mean L1 delta | Worsened L1 samples | Switched |
|---|---|---:|---:|---:|---:|---:|---:|
| marine_exactsplit_current_gtdb_transfer | 0,2 | 0.919315 | 0.140791 | 0.991722 | 0.000000000 | 0 | 0 |
| marine_exactsplit_refined_allocator_xny230 | 0,2 | 0.919315 | 0.140270 | 0.991759 | -0.000521634 | 0 | 2 |
| refined_vs_current | 0,2 | 0.919315 | 0.140270 | 0.991759 | -0.000521634 | 0 | 2 |

## Decision

- Keep the current default unchanged.
- Treat this as independent diagnostic stress support only.
- Release-grade promotion still needs clean GTDB holdout validation.

## Outputs

- `results/feature_allocator_refined_marine_diagnostic_scores.tsv`
- `results/feature_allocator_refined_marine_diagnostic_summary.tsv`
- `results/feature_allocator_refined_marine_diagnostic_validation.tsv`
- `results/feature_allocator_refined_marine_diagnostic_audit.tsv`
