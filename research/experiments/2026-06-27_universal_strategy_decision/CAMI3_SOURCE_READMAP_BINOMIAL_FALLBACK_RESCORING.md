# CAMI3 Source-Readmap Binomial Fallback Rescoring

Date: 2026-06-29

This cached-only audit adds the resolver-candidate exact GTDB-binomial
fallback rows to the existing source-readmap truth table and scores the
same cached MinCO/Sylph profile outputs against the adjusted truth. It
does not change MinCO calls and does not reread input data.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `fallback_added_rows` | 1732020 | cached_source_level_truth_expansion |
| `profile_scope_release_ready_after_fallback` | 0,1,2 | candidate_truth_transfer_reaches_coverage_threshold |
| `minco_vs_sylph_mean_F1_after_rescore` | 0.853058 vs 0.762024 | minco_higher_F1 |
| `minco_vs_sylph_mean_L1_after_rescore` | 61.884168 vs 31.103234 | sylph_lower_L1 |
| `mean_F1_delta_vs_previous_truth` | minco=0.037102;sylph=0.020502 | truth_expansion_changes_metric_denominator |
| `promotion_decision` | diagnostic_rescore_not_default_change | update_holdout_status_only_after_manual_truth_rule_review |

## Coverage

| Sample | Added rows | Base profile-scope mapped % | Adjusted profile-scope mapped % | Ready |
|---:|---:|---:|---:|---|
| 0 | 397912 | 93.827 | 95.941 | true |
| 1 | 521696 | 94.554 | 96.806 | true |
| 2 | 812412 | 96.870 | 99.908 | true |

## Summary Scores

| Method | Mean F1 | Pooled F1 | Mean L1 pp | Mean Pearson |
|---|---:|---:|---:|---:|
| minco_source_readmap_binomial_fallback_rescore | 0.853058 | 0.851675 | 61.884168 | 0.826285 |
| sylph_source_readmap_binomial_fallback_rescore | 0.762024 | 0.762617 | 31.103234 | 0.947203 |

## Delta Versus Previous Truth

| Sample | Method | Delta truth species | Delta F1 | Delta L1 pp | Delta Pearson |
|---:|---|---:|---:|---:|---:|
| 0 | minco_source_readmap_binomial_fallback_rescore | 9 | 0.043326 | -0.440127 | -0.000057 |
| 0 | sylph_source_readmap_binomial_fallback_rescore | 9 | 0.005005 | -0.615671 | -0.000063 |
| 1 | minco_source_readmap_binomial_fallback_rescore | 6 | 0.023523 | 0.572035 | 0.001342 |
| 1 | sylph_source_readmap_binomial_fallback_rescore | 6 | 0.018003 | -0.528251 | 0.000066 |
| 2 | minco_source_readmap_binomial_fallback_rescore | 8 | 0.044456 | -2.675491 | 0.005825 |
| 2 | sylph_source_readmap_binomial_fallback_rescore | 8 | 0.038499 | -3.444899 | 0.004481 |

## Decision

- The fallback reaches the profile-scope coverage threshold for samples
  0, 1, and 2 in this cached source-level rescore.
- This updates the CAMI3 evidence path, but it is still a truth-transfer
  audit. It does not justify changing the MinCO default strategy.
- Before promoting CAMI3 to release-grade evidence, review whether exact
  GTDB-binomial fallback is acceptable for the benchmark truth policy.

## Outputs

- `results/cami3_source_readmap_binomial_fallback_truth.tsv`
- `results/cami3_source_readmap_binomial_fallback_quality.tsv`
- `results/cami3_source_readmap_binomial_fallback_scores.tsv`
- `results/cami3_source_readmap_binomial_fallback_summary.tsv`
- `results/cami3_source_readmap_binomial_fallback_delta.tsv`
- `results/cami3_source_readmap_binomial_fallback_audit.tsv`
