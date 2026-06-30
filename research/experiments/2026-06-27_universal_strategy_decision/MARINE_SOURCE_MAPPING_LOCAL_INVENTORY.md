# Marine Source Mapping Local Inventory

Date: 2026-06-30

This audit checks whether the local marine workspace contains the source
crosswalk needed to use marine readmaps as release-grade GTDB truth.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `local_marine_setup_metadata_status` | local_marine_genome_to_id.tsv,local_marine_metadata.tsv | setup_metadata_present |
| `local_marine_crosswalk_files_missing` | local_marine_gsa_pooled_mapping.tsv.gz | pooled_mapping_missing_but_setup_metadata_present |
| `archive_relevant_member_status` | readmap_only_archives=3 | archives_do_not_embed_source_crosswalk |
| `promotion_decision` | do_not_promote_marine_from_inventory_alone | setup_metadata_present_requires_truth_upgrade_audit |

## Inventory

| Item | Status | Value | Decision |
|---|---|---|---|
| `local_marine_workdir` | true | /tmp/cami2_marine_samples3_5_20260625 | workdir_available |
| `local_marine_genome_to_id.tsv` | present | 1 | local_setup_metadata_present |
| `local_marine_gsa_pooled_mapping.tsv.gz` | missing | 0 | source_crosswalk_not_found |
| `local_marine_metadata.tsv` | present | 1 | local_setup_metadata_present |
| `sample3_archive_relevant_members` | present | simulation_short_read/2018.08.15_09.49.32_sample_3/reads/reads_mapping.tsv.gz | readmap_only_no_source_crosswalk |
| `sample4_archive_relevant_members` | present | simulation_short_read/2018.08.15_09.49.32_sample_4/reads/reads_mapping.tsv.gz | readmap_only_no_source_crosswalk |
| `sample5_archive_relevant_members` | present | simulation_short_read/2018.08.15_09.49.32_sample_5/reads/reads_mapping.tsv.gz | readmap_only_no_source_crosswalk |
| `sibling_panel_crosswalk_files` | present | /mnt/new3T/minco_cami2_plant_20260621/simulation_short_read/genome_to_id.tsv;/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read/metadata.tsv;/mnt/new3T/minco_cami2_strain_20260621/short_read/genome_to_id.tsv;/mnt/new3T/minco_cami2_strain_20260621/short_read/metadata.tsv | sibling_crosswalks_exist_but_do_not_resolve_marine |

## Decision

The local marine workspace now contains extracted setup metadata
(`genome_to_id.tsv` and `metadata.tsv`). That improves source-name
auditing but is not, by itself, proof that the truth transfer reaches
the release threshold. Use `audit_marine_setup_metadata_truth_upgrade.py`
for the threshold decision.
