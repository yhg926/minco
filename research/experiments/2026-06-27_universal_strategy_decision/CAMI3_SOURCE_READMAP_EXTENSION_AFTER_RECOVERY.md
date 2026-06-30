# CAMI3 Source-Readmap Extension After Recovery

Date: 2026-06-29

This helper records whether CAMI3 samples3-5 are ready for the
source-readmap release-extension score. It does not download data. When
the required `reads_mapping.tsv.gz` files are present, it scores the
selected candidate MinCO default, the opt-in refined abundance allocator,
and Sylph in the same GTDB species namespace. The exact-binomial
fallback truth adjustment is used only if all policy checks pass.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `post_recovery_input_status` | readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3 | ready_to_score |
| `direct_readmap_truth_quality` | min_mapped_read_rows_pct=74.826330;ready_samples= | direct_readmap_truth_needs_fallback_review |
| `exact_binomial_fallback_policy` | accepted_samples=3/3;min_adjusted_profile_scope_pct=98.825642 | fallback_truth_policy_pass |
| `candidate_minco_vs_sylph` | minco_F1_wins=0/3;minco_L1_wins=0/3 | independent_holdout_comparison |
| `refined_allocator_effect` | improved=0;worsened=0;mean_L1_delta_pp=0.000000000 | do_not_promote_from_extension_result |
| `promotion_decision` | post_recovery_extension_scored | review_holdout_result_before_default_change |

## Inputs

| Sample | Readmap | Selected Default | Raw Tables | Sylph | Ready | Blocker |
|---:|---|---|---|---|---|---|
| 3 | True | True | True | True | True |  |
| 4 | True | True | True | True | True |  |
| 5 | True | True | True | True | True |  |

## Truth Quality

| Sample | Read rows | Mapped % | GTDB species |
|---:|---:|---:|---:|
| 3 | 33294108 | 77.968090 | 107 |
| 4 | 33285864 | 80.993019 | 86 |
| 5 | 33286454 | 74.826330 | 99 |

## Score Summary

| Method | Mean F1 | Pooled F1 | Mean L1 pp | Mean Pearson |
|---|---:|---:|---:|---:|
| sylph_source_readmap_extension | 0.784196 | 0.783582 | 21.448926 | 0.981150 |
| minco_candidate_source_readmap_extension | 0.706093 | 0.709193 | 43.530195 | 0.922077 |
| minco_refined_allocator_source_readmap_extension | 0.706093 | 0.709193 | 43.530195 | 0.922077 |

## Refined Allocator Delta

| Sample | Delta F1 | Delta L1 pp | Delta Pearson | Guard passed |
|---:|---:|---:|---:|---|
| 3 | 0.000000 | -0.000000 | 0.000000 | False |
| 4 | 0.000000 | 0.000000 | -0.000000 | False |
| 5 | 0.000000 | 0.000000 | 0.000000 | False |

## Outputs

- `results/cami3_source_readmap_extension_after_recovery_inputs.tsv`
- `results/cami3_source_readmap_extension_after_recovery_audit.tsv`
- `results/cami3_source_readmap_extension_after_recovery_truth.tsv`
- `results/cami3_source_readmap_extension_after_recovery_source_genomes.tsv`
- `results/cami3_source_readmap_extension_after_recovery_fallback_detail.tsv`
- `results/cami3_source_readmap_extension_after_recovery_quality.tsv`
- `results/cami3_source_readmap_extension_after_recovery_scores.tsv`
- `results/cami3_source_readmap_extension_after_recovery_summary.tsv`
- `results/cami3_source_readmap_extension_after_recovery_delta.tsv`
