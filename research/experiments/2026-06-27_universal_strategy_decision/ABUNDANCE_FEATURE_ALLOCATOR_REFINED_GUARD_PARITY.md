# Refined Feature Allocator Guard Parity

Date: 2026-06-29

This generated note validates the implemented refined allocator switch
against the cached combined guard-refinement rule. It does not rerun raw
profiling jobs.

## Decision

- Profiles: `38`.
- Rule parity: `mismatches=0;max_s_xny_median_delta=0` (pass).
- Applied profiles: `total=33;selected_cached=29;external=4`.
- Promotion decision: `guard_parity_only_not_default`.

## Outputs

- `results/feature_allocator_refined_guard_parity_apply.tsv`
- `results/feature_allocator_refined_guard_parity_audit.tsv`
