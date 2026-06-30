# Default Release Snapshot

Date: 2026-06-29

This generated note consolidates cached evidence for the current MinCO
default. It separates the selected user default from broader release
claims. No raw profiling jobs are run by this snapshot.

## Decision

- Selected default: `candidate` preset via `scripts/minco_profile_default.py`.
- Keep experimental abundance allocators off by default.
- Current evidence supports the candidate preset as the best cached MinCO
  default versus the previous MinCO default, but not a broad Sylph-beating
  abundance claim.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `selected_default` | candidate | user_default_entrypoint |
| `matched_panel_scope` | panels=4;release=4;diagnostic=0 | cached_gtdb_default_evidence |
| `candidate_vs_previous_minco` | F1_regressions=0;L1_regressions=0;Pearson_regressions=0 | candidate_dominates_previous_cached_default |
| `candidate_vs_sylph` | F1_wins=1/4;L1_wins=0/4;Pearson_wins=0/4 | not_broad_sylph_beating |
| `external_exactsplit_diagnostic` | datasets=5;F1_win_datasets=3;L1_win_datasets=2 | diagnostic_only_nonrelease |
| `guarded_allocator_external_status` | mean_delta_l1=-0.004134;improved_datasets=2;worsened_datasets=0;improved_samples=4;worsened_samples=2;max_sample_worse=0.000249;switched_samples=6;adjusted_rows=130 | keep_feature_allocator_experimental_off_by_default |
| `runtime_tradeoff` | samples=7;seconds_ratio=2.515;rss_ratio=0.191 | minco_slower_lower_memory_in_cached_timed_runs |
| `release_decision` | selected_candidate_default;experimental_allocators_off;do_not_claim_broad_sylph_beating | stable_default_candidate_not_final_universal_claim |

## Matched Panels

| Panel | Grade | Samples | Delta F1 vs Sylph | Delta L1 pp vs Sylph | Delta Pearson vs Sylph |
|---|---|---|---:|---:|---:|
| CAMI2 Toy Mouse | release | 5,6,7 | -0.064928 | 6.798669 | -0.005824 |
| CAMI3 ToyGut source-readmap | release | 0,1,2 | 0.098143 | 24.816110 | -0.080610 |
| HMP airskin source-abundance | release | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | -0.112163 | 17.056351 | -0.016855 |
| HMP gastrooral source-abundance | release | 0,6 | -0.028484 | 41.294813 | -0.102844 |

## Diagnostic Exact-Split Panels

| Dataset | Samples | Mean delta F1 | Mean delta L1 | Mean delta Pearson |
|---|---|---:|---:|---:|
| hmp_gastrooral | hmp_gastrooral6 | -0.080378 | 0.063045 | -0.005277 |
| plant | plant0,plant1,plant2 | 0.099117 | -0.418542 | 0.219428 |
| plant_holdout | plant_holdout3,plant_holdout4,plant_holdout5 | 0.012023 | -0.247014 | 0.311106 |
| strain | strain0,strain1,strain2 | 0.101363 | 0.556779 | -0.170184 |
| toy_gut | toy0 | -0.201740 | 1.395508 | -1.163403 |

## Runtime

| Samples | Mean MinCO seconds | Mean Sylph seconds | Seconds ratio | Mean MinCO RSS GiB | Mean Sylph RSS GiB | RSS ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 265.27 | 105.47 | 2.515 | 3.528 | 18.490 | 0.191 |

## Outputs

- `results/default_release_snapshot_panel_metrics.tsv`
- `results/default_release_snapshot_diagnostic_metrics.tsv`
- `results/default_release_snapshot_runtime.tsv`
- `results/default_release_snapshot_audit.tsv`
