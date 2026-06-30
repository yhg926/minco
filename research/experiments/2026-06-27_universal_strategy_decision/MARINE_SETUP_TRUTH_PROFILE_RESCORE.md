# Marine Setup Truth Profile Rescore

Date: 2026-06-30

This scorer applies the strict setup source-name or local assembly-summary
GTDB truth rule to marine profile outputs. It does not run profiling.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `scoring_scope` | samples=3,4,5;rule_id=strict_source_name_or_assembly_sources_unique;profile_set_label=selected_default_same_namespace | same_namespace_truth_rule_applied |
| `truth_rule_quality` | ready_samples=3/3;min_mapped_pct_bacteria_archaea=96.463685 | truth_rule_passes_threshold |
| `minco_vs_sylph_summary` | minco_mean_F1=0.781358;sylph_mean_F1=0.838716;minco_mean_L1_union_pp=41.299694;sylph_mean_L1_union_pp=40.997305;minco_mean_Pearson_union=0.886706;sylph_mean_Pearson_union=0.923667 | sylph_higher_F1 |
| `promotion_decision` | selected_default_profile_scored_review_metrics | route_profile_evidence_ready_for_review |

## Outputs

- `results/marine_setup_truth_profile_rescore_truth.tsv`
- `results/marine_setup_truth_profile_rescore_quality.tsv`
- `results/marine_setup_truth_profile_rescore_scores.tsv`
- `results/marine_setup_truth_profile_rescore_summary.tsv`
- `results/marine_setup_truth_profile_rescore_audit.tsv`
