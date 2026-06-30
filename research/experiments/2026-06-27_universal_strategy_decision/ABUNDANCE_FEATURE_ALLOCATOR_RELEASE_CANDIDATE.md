# Feature Allocator Release Candidate Audit

Date: 2026-06-29

This generated note combines cached wrapper-parity evidence with the
external exact-split stress test for the guarded feature allocator. It
does not rerun raw profiling jobs.

## Decision

- Cached wrapper evidence: `samples=32;switched=30;mean_L1_delta_pp=-0.3188657861319407;worsened=0;max_worse=0.0`.
- External stress evidence: `mean_delta_l1=-0.004134;improved_datasets=2;worsened_datasets=0;improved_samples=4;worsened_samples=2;max_sample_worse=0.000249;switched_samples=6;adjusted_rows=130`.
- Combined effect: `samples=38;switched=36;improved=34;worsened=2;mean_L1_delta_pp=-0.269171;max_worse_L1_delta_pp=0.000249;mean_F1_delta=0`.
- External regression detail: `plant_holdout3:delta_L1=0.000248910;base_multi=0.794009;rows=40;plant_holdout4:delta_L1=0.000060987;base_multi=0.762779;rows=33`.
- Promotion decision: `do_not_promote_until_external_sample_safe`.

## Panel Evidence

| Evidence group | Panel | Samples | Mean L1 delta pp | Max worse pp | Improved | Worsened | Decision |
|---|---|---:|---:|---:|---:|---:|---|
| selected_cached_wrapper | cami2_toy_mouse_gut | 3 | -0.187921 | -0.030594 | 3 | 0 | cached_wrapper_sample_safe |
| selected_cached_wrapper | cami3_toy_human_gut_gtdb_source_readmap | 3 | -0.078847 | 0.000000 | 1 | 0 | cached_wrapper_sample_safe |
| selected_cached_wrapper | hmp_airskin_gtdb_source_abundance | 24 | -0.381306 | -0.041632 | 24 | 0 | cached_wrapper_sample_safe |
| selected_cached_wrapper | hmp_gastrooral_gtdb_source_abundance | 2 | -0.126023 | -0.012398 | 2 | 0 | cached_wrapper_sample_safe |
| external_exactsplit_diagnostic | plant_holdout | 3 | -0.000021 | 0.000249 | 1 | 2 | external_sample_regression |
| external_exactsplit_diagnostic | strainmadness | 3 | -0.008248 | -0.003423 | 3 | 0 | external_sample_safe |

## Interpretation

- The allocator remains a useful experimental candidate because it
  improves the cached selected-default panels and both external panel
  means without changing calls.
- It is not promoted into the default because strict sample safety fails
  on the external exact-split diagnostic profiles.
- The next abundance step should change the output-only guard or target
  feature so independent-profile sample regressions are zero before
  default promotion is reconsidered.

## Outputs

- `results/feature_allocator_release_candidate_panel.tsv`
- `results/feature_allocator_release_candidate_audit.tsv`
