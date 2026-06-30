# Refined Feature Allocator Score Replay

Date: 2026-06-29

This cached-profile replay validates the implemented opt-in switch
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002-xny230`. The call set is fixed;
only abundance mass is adjusted, and the default remains unchanged.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `baseline_validation` | counts=0;F1=1.11e-16;L1=1.15e-12;Pearson=5.55e-16 | pass |
| `score_replay_effect` | samples=38;switched=33;improved=33;worsened=0;mean_L1_delta_pp=-0.262955;selected_mean_L1_delta_pp=-0.311474;external_mean_L1_delta_pp=-0.004186;max_worse_L1_delta_pp=0.000000 | fixed_call_abundance_replay |
| `score_replay_matches_guard_estimate` | max_abs_mean_delta=3.10619365229e-07 | pass |
| `promotion_decision` | refined_guard_score_replay_pass_candidate_needs_independent_holdout | current_default_unchanged |

## Summary

| Evidence group | Profiles | Switched | Mean F1 delta | Mean L1 delta | Max worse L1 | Improved L1 | Worsened L1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| external_exactsplit_diagnostic | 6 | 4 | 0.000000 | -0.004186 | 0.000000 | 4 | 0 |
| selected_cached_wrapper | 32 | 29 | 0.000000 | -0.311474 | 0.000000 | 29 | 0 |
| all | 38 | 33 | 0.000000 | -0.262955 | 0.000000 | 33 | 0 |

## Decision

- Keep the current default unchanged.
- Treat the refined allocator as the strongest opt-in abundance candidate.
- Require independent holdout validation before any default promotion.

## Outputs

- `results/feature_allocator_refined_score_replay_scores.tsv`
- `results/feature_allocator_refined_score_replay_summary.tsv`
- `results/feature_allocator_refined_score_replay_validation.tsv`
- `results/feature_allocator_refined_score_replay_audit.tsv`
