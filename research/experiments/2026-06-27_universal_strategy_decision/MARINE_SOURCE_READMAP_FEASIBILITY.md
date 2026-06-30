# Marine Source Readmap Feasibility

Date: 2026-06-29

This audit samples source read maps embedded in local marine archives.
It checks whether local read IDs can directly resolve unresolved marine
truth rows to GTDB assemblies. It does not extract FASTQ files.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `sampled_readmap_rows` | 6000000 | source_readmaps_available_in_local_archives |
| `direct_sequence_to_gtdb_matches` | sequence_accession_rows=0;core_rows=0 | no_direct_contig_to_gtdb_assembly_mapping |
| `unresolved_transfer_rows_sampled` | 314110 | readmaps_cover_unresolved_taxids_but_need_source_mapping |
| `promotion_decision` | do_not_promote_marine_source_readmap_from_local_ids | contig_or_otu_to_assembly_mapping_missing |

## Samples

| Sample | Rows Read | Genome IDs | Source Seq IDs | Unresolved Rows | Direct GTDB Matches |
|---:|---:|---:|---:|---:|---:|
| 3 | 2000000 | 789 | 12016 | 117704 | 0 |
| 4 | 2000000 | 581 | 9397 | 87118 | 0 |
| 5 | 2000000 | 656 | 9570 | 109288 | 0 |

## Decision

The read maps are useful, but the local IDs are sequence/contig IDs or
CAMI genome IDs. They do not directly match GTDB assembly accessions in
the local metadata. Marine still needs contig/OTU-to-assembly metadata,
source FASTA provenance, or another clean holdout before release use.

## Outputs

- `results/marine_source_readmap_feasibility_summary.tsv`
- `results/marine_source_readmap_unresolved_top.tsv`
- `results/marine_source_readmap_feasibility_audit.tsv`
