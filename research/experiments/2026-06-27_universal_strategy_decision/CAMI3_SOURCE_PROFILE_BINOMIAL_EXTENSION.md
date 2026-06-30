# CAMI3 Source-Profile Binomial Extension

Date: 2026-06-29

This diagnostic scores CAMI3 samples3-5 using local taxonomic-profile
strain/source rows with deterministic GTDB transfer. It does not replace
the missing per-read source-readmap truth files.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `source_profile_truth_scope` | samples=3,4,5;min_mapped_profile_scope_pct=62.728832;ready_samples= | source_profile_truth_mapping_incomplete |
| `baseline_minco_vs_sylph` | minco_F1_wins=0/3;minco_L1_wins=0/3 | diagnostic_comparison_only |
| `refined_allocator_effect` | improved=0;worsened=0;mean_L1_delta_pp=-0.000000000 | do_not_promote_from_source_profile_diagnostic |
| `promotion_decision` | source_profile_extension_diagnostic_not_source_readmap_release_truth | keep_refined_allocator_opt_in |

## Truth Quality

| Sample | In-scope abundance % | Mapped in-scope % | GTDB species | Ready |
|---:|---:|---:|---:|---|
| 3 | 38.428900 | 62.728832 | 91 | false |
| 4 | 43.307200 | 77.004747 | 75 | false |
| 5 | 40.439500 | 67.403900 | 87 | false |

## Score Summary

| Method | Samples | Mean F1 | Mean L1 union pp | Mean Pearson union |
|---|---|---:|---:|---:|
| sylph_source_profile | 3,4,5 | 0.743395 | 84.470851 | 0.602895 |
| minco_candidate_refined_allocator_source_profile | 3,4,5 | 0.661335 | 97.931764 | 0.556331 |
| minco_candidate_source_profile | 3,4,5 | 0.661335 | 97.931764 | 0.556331 |

## Decision

- Keep this as diagnostic evidence until the truth policy is reviewed.
- The main release route remains recovering per-read source mapping for
  CAMI3 samples3-5.
- The refined allocator remains opt-in by default.

## Outputs

- `results/cami3_source_profile_binomial_extension_truth.tsv`
- `results/cami3_source_profile_binomial_extension_source_rows.tsv`
- `results/cami3_source_profile_binomial_extension_quality.tsv`
- `results/cami3_source_profile_binomial_extension_scores.tsv`
- `results/cami3_source_profile_binomial_extension_summary.tsv`
- `results/cami3_source_profile_binomial_extension_delta.tsv`
- `results/cami3_source_profile_binomial_extension_audit.tsv`

Release decision remains `keep_refined_allocator_opt_in`.
