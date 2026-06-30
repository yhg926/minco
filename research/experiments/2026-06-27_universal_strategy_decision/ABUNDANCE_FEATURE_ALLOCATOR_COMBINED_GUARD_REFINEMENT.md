# Feature Allocator Combined Guard Refinement

Date: 2026-06-29

This generated note searches simple output-only guard refinements for
`guarded-genus-hit-breadth-a002` using cached selected-default profiles
plus external exact-split diagnostic profiles. It does not rerun raw
profiling jobs.

## Decision

- Current guard: `samples=38;improved=34;worsened=2;mean_L1_delta_pp=-0.269171;selected_mean_L1_delta_pp=-0.318866;external_worsened=2`.
- Tested rules: `1192`.
- Strict safe rules: `352`.
- Best strict rule: `s_xny_median>=230.3;mean_L1_delta_pp=-0.262955;selected_mean_L1_delta_pp=-0.311474;external_mean_L1_delta_pp=-0.004186;cached_gain_frac=0.976818`.
- Promotion decision: `refined_guard_candidate_needs_wrapper_validation`.

## Top Rules

| Rule | Mean L1 delta pp | Worsened | External worsened | Cached gain frac | Strict |
|---|---:|---:|---:|---:|---|
| s_xny_median>=230.3 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_mean>=345.779 | -0.262955 | 0 | 0 | 0.976818 | True |
| xny_max_median>=230.3 | -0.262955 | 0 | 0 | 0.976818 | True |
| xny_max_mean>=345.779 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_depth_median>=0.225311 | -0.262955 | 0 | 0 | 0.976818 | True |
| u_depth_median>=0.12905 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_median>=230.3 AND s_xny_mean>=345.779 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_median>=230.3 AND xny_max_median>=230.3 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_median>=230.3 AND xny_max_mean>=345.779 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_median>=230.3 AND s_depth_median>=0.225311 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_median>=230.3 AND u_depth_median>=0.12905 | -0.262955 | 0 | 0 | 0.976818 | True |
| s_xny_mean>=345.779 AND xny_max_median>=230.3 | -0.262955 | 0 | 0 | 0.976818 | True |

## Outputs

- `results/feature_allocator_combined_guard_features.tsv`
- `results/feature_allocator_combined_guard_rules.tsv`
- `results/feature_allocator_combined_guard_audit.tsv`
