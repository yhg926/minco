# CAMI3 Source-Readmap Scope Audit

Date: 2026-06-29

This cached audit separates CAMI3 source-readmap rows that are in scope
for GTDB species profiling from host, viral, fungal/eukaryotic, plasmid,
and other out-of-scope rows. It does not rerun profiling.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `all_rows_mapping_status` | sample0=79.832%;sample1=79.059%;sample2=84.914% | all_rows_include_out_of_scope_sources |
| `profile_scope_mapping_status` | sample0=93.827%;sample1=94.554%;sample2=96.870% | profile_scope_nearly_release_grade |
| `samples_below_95pct_profile_scope` | 0,1 | targeted_mapping_needed |
| `additional_rows_needed_for_95pct` | 324026 | feasible_if_top_sources_can_be_mapped |
| `top_in_scope_sources_to_review` | ASV579.0,ASV783.5,ASV134.7,ASV579.0,ASV134.7,ASV476.8 | manual_or_metadata_mapping_review |
| `release_upgrade_decision` | do_not_promote_yet | needs_targeted_in_scope_source_mapping_before_release_grade |

## Summary

| Sample | All rows mapped % | Profile-scope mapped % | Rows needed for 95% | Scope ready |
|---:|---:|---:|---:|---|
| 0 | 79.832 | 93.827 | 220749 | false |
| 1 | 79.059 | 94.554 | 103277 | false |
| 2 | 84.914 | 96.870 | 0 | true |

## Top In-Scope Unmapped Sources

| Sample | Source | Taxid | Species label | Unmapped rows |
|---:|---|---:|---|---:|
| 0 | ASV579.0 | 626940 | Phascolarctobacterium succinatutens | 763090 |
| 0 | ASV783.5 | 562 | Escherichia coli | 137208 |
| 0 | ASV134.7 | 817 | Bacteroides fragilis | 104128 |
| 0 | ASV476.8 | 1496 | Clostridioides difficile | 84300 |
| 0 | ASV607.0 | 29466 | Veillonella parvula | 29434 |
| 0 | ASV444.5 | 33038 | Mediterraneibacter gnavus | 28894 |
| 0 | ASV71.10 | 1681 | Bifidobacterium bifidum | 6846 |
| 0 | ASV74.2 | 28901 | Salmonella enterica | 1516 |
| 0 | ASV670.1 | 225992 | Comamonas kerstersii | 1136 |
| 0 | ASV156.5 | 329854 | Bacteroides intestinalis | 1094 |
| 0 | ASV577.0 | 187327 | Acidaminococcus intestini | 708 |
| 0 | ASV364.1 | 1776384 | Emergencia timonensis | 654 |
| 0 | ASV627.1 | 76857 | Fusobacterium polymorphum | 546 |
| 0 | ASV676.10 | 28901 | Salmonella enterica | 512 |
| 0 | ASV790.3 | 287 | Pseudomonas aeruginosa | 502 |
| 0 | ASV368.3 | 1496 | Clostridioides difficile | 328 |
| 0 | ASV709.10 | 562 | Escherichia coli | 194 |
| 0 | ASV737.3 | 562 | Escherichia coli | 194 |
| 0 | ASV592.0 | 907 | Megasphaera elsdenii | 188 |
| 0 | ASV348.5 | 562 | Escherichia coli | 184 |

## Decision

- CAMI3 should remain diagnostic for now.
- After excluding clear out-of-scope rows, sample 2 crosses 95% mapped
  profile-scope coverage; samples 0 and 1 are close but still below 95%.
- The next CAMI3 release-upgrade task is targeted mapping of the top
  in-scope unmapped sources, not another MinCO call-threshold sweep.

## Outputs

- `results/cami3_source_readmap_scope_summary.tsv`
- `results/cami3_source_readmap_unmapped_by_category.tsv`
- `results/cami3_source_readmap_top_unmapped_in_scope.tsv`
- `results/cami3_source_readmap_scope_audit.tsv`
