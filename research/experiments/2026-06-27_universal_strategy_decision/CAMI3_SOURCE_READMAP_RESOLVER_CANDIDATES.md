# CAMI3 Source-Readmap Resolver Candidates

Date: 2026-06-29

This cached-only audit checks deterministic metadata fallbacks for the
top in-scope CAMI3 source-readmap rows that remain unmapped after the
current WGS-prefix and unique-taxid transfer. It does not use raw
input data and does not rerun MinCO or Sylph.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `reviewed_top_sources_per_sample` | 20 | cached_only_no_profile_rerun |
| `reviewed_unmapped_rows` | 3258804 | top_in_scope_source_rows |
| `recommended_resolvable_rows` | 1732020 | candidate_truth_transfer_fallback |
| `recommended_resolvable_pct_of_reviewed` | 53.149 | large_enough_to_test_release_upgrade |
| `projected_samples_release_ready_after_candidate` | 0,1,2 | rescored_truth_test_is_worth_running |
| `promotion_decision` | do_not_promote_without_rescored_truth_audit | exact_binomial_fallback_is_candidate_not_current_truth |

## Summary

| Sample | Reviewed rows | Unique taxid rows | Unique NCBI-name rows | Unique GTDB-binomial rows | Recommended rows |
|---:|---:|---:|---:|---:|---:|
| 0 | 1161656 | 0 | 0 | 397912 | 397912 |
| 1 | 1261518 | 0 | 0 | 521696 | 521696 |
| 2 | 835630 | 0 | 0 | 812412 | 812412 |

## Release-Coverage Projection

| Sample | Current mapped % | Rows needed for 95% | Recommended rows | Projected mapped % | Projected ready |
|---:|---:|---:|---:|---:|---|
| 0 | 93.827 | 220749 | 397912 | 95.941 | true |
| 1 | 94.554 | 103277 | 521696 | 96.806 | true |
| 2 | 96.870 | 0 | 812412 | 99.908 | true |

## Top Candidate Resolutions

| Sample | Source | Species label | Unmapped rows | Recommended GTDB species | Rule | Ambiguity note |
|---:|---|---|---:|---|---|---|
| 0 | ASV579.0 | Phascolarctobacterium succinatutens | 763090 |  |  | taxid_candidates=17;ncbi_name_candidates=17 |
| 0 | ASV783.5 | Escherichia coli | 137208 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |
| 0 | ASV134.7 | Bacteroides fragilis | 104128 | s__Bacteroides fragilis | unique_gtdb_binomial | taxid_candidates=22;ncbi_name_candidates=22 |
| 0 | ASV476.8 | Clostridioides difficile | 84300 | s__Clostridioides difficile | unique_gtdb_binomial | taxid_candidates=9;ncbi_name_candidates=9 |
| 0 | ASV607.0 | Veillonella parvula | 29434 | s__Veillonella parvula | unique_gtdb_binomial | taxid_candidates=28;ncbi_name_candidates=28 |
| 0 | ASV444.5 | Mediterraneibacter gnavus | 28894 | s__Mediterraneibacter gnavus | unique_gtdb_binomial | taxid_candidates=4;ncbi_name_candidates=4 |
| 0 | ASV71.10 | Bifidobacterium bifidum | 6846 | s__Bifidobacterium bifidum | unique_gtdb_binomial | taxid_candidates=3;ncbi_name_candidates=3 |
| 0 | ASV74.2 | Salmonella enterica | 1516 | s__Salmonella enterica | unique_gtdb_binomial | taxid_candidates=11;ncbi_name_candidates=11 |
| 0 | ASV670.1 | Comamonas kerstersii | 1136 | s__Comamonas kerstersii | unique_gtdb_binomial | taxid_candidates=2;ncbi_name_candidates=2 |
| 0 | ASV156.5 | Bacteroides intestinalis | 1094 | s__Bacteroides intestinalis | unique_gtdb_binomial | taxid_candidates=7;ncbi_name_candidates=7 |
| 0 | ASV577.0 | Acidaminococcus intestini | 708 | s__Acidaminococcus intestini | unique_gtdb_binomial | taxid_candidates=5;ncbi_name_candidates=5 |
| 0 | ASV364.1 | Emergencia timonensis | 654 |  |  | taxid_candidates=2;ncbi_name_candidates=2 |
| 0 | ASV627.1 | Fusobacterium polymorphum | 546 | s__Fusobacterium polymorphum | unique_gtdb_binomial | taxid_candidates=3;ncbi_name_candidates=3 |
| 0 | ASV676.10 | Salmonella enterica | 512 | s__Salmonella enterica | unique_gtdb_binomial | taxid_candidates=11;ncbi_name_candidates=11 |
| 0 | ASV790.3 | Pseudomonas aeruginosa | 502 | s__Pseudomonas aeruginosa | unique_gtdb_binomial | taxid_candidates=5;ncbi_name_candidates=5 |
| 0 | ASV368.3 | Clostridioides difficile | 328 | s__Clostridioides difficile | unique_gtdb_binomial | taxid_candidates=9;ncbi_name_candidates=9 |
| 0 | ASV709.10 | Escherichia coli | 194 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |
| 0 | ASV737.3 | Escherichia coli | 194 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |
| 0 | ASV592.0 | Megasphaera elsdenii | 188 | s__Megasphaera elsdenii | unique_gtdb_binomial | taxid_candidates=5;ncbi_name_candidates=5 |
| 0 | ASV348.5 | Escherichia coli | 184 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |
| 1 | ASV579.0 | Phascolarctobacterium succinatutens | 734058 |  |  | taxid_candidates=17;ncbi_name_candidates=17 |
| 1 | ASV134.7 | Bacteroides fragilis | 280814 | s__Bacteroides fragilis | unique_gtdb_binomial | taxid_candidates=22;ncbi_name_candidates=22 |
| 1 | ASV476.8 | Clostridioides difficile | 158700 | s__Clostridioides difficile | unique_gtdb_binomial | taxid_candidates=9;ncbi_name_candidates=9 |
| 1 | ASV444.5 | Mediterraneibacter gnavus | 24668 | s__Mediterraneibacter gnavus | unique_gtdb_binomial | taxid_candidates=4;ncbi_name_candidates=4 |
| 1 | ASV71.10 | Bifidobacterium bifidum | 22762 | s__Bifidobacterium bifidum | unique_gtdb_binomial | taxid_candidates=3;ncbi_name_candidates=3 |
| 1 | ASV815.4 | Clostridioides difficile | 8132 | s__Clostridioides difficile | unique_gtdb_binomial | taxid_candidates=9;ncbi_name_candidates=9 |
| 1 | ASV670.1 | Comamonas kerstersii | 7870 | s__Comamonas kerstersii | unique_gtdb_binomial | taxid_candidates=2;ncbi_name_candidates=2 |
| 1 | ASV364.1 | Emergencia timonensis | 5764 |  |  | taxid_candidates=2;ncbi_name_candidates=2 |
| 1 | ASV709.10 | Escherichia coli | 4324 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |
| 1 | ASV783.5 | Escherichia coli | 4308 | s__Escherichia coli | unique_gtdb_binomial | taxid_candidates=51;ncbi_name_candidates=51 |

## Decision

- Do not promote CAMI3 source-readmap to release-grade evidence yet.
- The exact GTDB-binomial fallback is promising enough to test in the
  scorer because it can resolve many top in-scope rows, but it needs
  a rescored truth audit before it becomes evidence.
- This is a metadata/truth-transfer task, not a MinCO threshold task.

## Outputs

- `results/cami3_source_readmap_resolver_candidates.tsv`
- `results/cami3_source_readmap_resolver_candidate_summary.tsv`
- `results/cami3_source_readmap_resolver_candidate_audit.tsv`
