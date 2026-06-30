# Refined Allocator Independent Holdout Inventory

Date: 2026-06-29

This generated audit checks whether the current local cache contains an
unused release-grade GTDB holdout profile set for the refined allocator.
It does not rerun raw profiling.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `release_grade_manifest_overlap` | release_panels=4;release_samples=32;selected_cache_overlap=32;independent_release_ready_samples=0 | no_unused_release_grade_cache |
| `hmp_omitted_release_candidate_inputs` | samples=2,8,12,26,27;truth=5/5;minco_profile=5/5;sylph_profile=5/5;ready_pairs=5/5 | missing_cached_profile_pairs |
| `independent_diagnostic_support` | diagnostic_supports_refined_allocator_holdout_candidate | diagnostic_only_not_release_grade |
| `promotion_decision` | release_grade_independent_holdout_missing | keep_refined_allocator_opt_in |

## Decision

- The existing release-grade manifest samples are all already present in
  the refined guard-selection cache.
- HMP omitted sample IDs have truth files, but no cached MinCO/Sylph
  profile pairs were found locally.
- Marine exact-split replay is independent diagnostic support, not
  release-grade promotion evidence.

## Outputs

- `results/feature_allocator_refined_independent_holdout_inventory.tsv`
- `results/feature_allocator_refined_independent_holdout_inventory_audit.tsv`
