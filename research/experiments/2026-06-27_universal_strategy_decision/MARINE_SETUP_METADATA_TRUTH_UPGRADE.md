# Marine Setup Metadata Truth Upgrade

Date: 2026-06-30

This audit tests whether extracted marine setup metadata is enough to
upgrade GTDB truth transfer. It optionally uses a local RefSeq assembly
summary as source metadata. It does not rerun profilers and does not
change the accepted marine scorer.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `setup_metadata_inputs` | genome_to_id_exists=True;metadata_exists=True;setup_genomes=977 | setup_metadata_available |
| `assembly_summary_inputs` | assembly_summary_exists=True;path=/mnt/new3T/All_pathogen/metadata/assembly_summary_refseq.txt | assembly_summary_available |
| `exact_source_name_match_summary` | setup_genomes=977;exact_unique=300;exact_ambiguous=8;no_match=669 | source_names_partially_resolve_gtdb_species |
| `assembly_summary_match_summary` | setup_genomes=977;assembly_unique=472;assembly_ambiguous=6;assembly_no_match=499 | assembly_summary_resolves_additional_sources |
| `accepted_exact_binomial_summary` | min_mapped_pct_bacteria_archaea=91.005806;ready_samples=0/10 | release_rule_fails_threshold |
| `strict_all_setup_sources_unique_summary` | min_mapped_pct_bacteria_archaea=93.043167;ready_samples=4/10 | release_rule_fails_threshold |
| `strict_assembly_summary_sources_unique_summary` | min_mapped_pct_bacteria_archaea=93.978664;ready_samples=9/10 | release_rule_fails_threshold |
| `strict_source_name_or_assembly_sources_unique_summary` | min_mapped_pct_bacteria_archaea=95.029253;ready_samples=10/10 | release_rule_passes_threshold |
| `diagnostic_partial_setup_unique_source_name_summary` | min_mapped_pct_bacteria_archaea=94.809534;ready_samples=9/10 | diagnostic_not_release_rule |
| `strict_all_setup_sources_unique_rescues` | rescued_rows=83;rescued_mass_pct_all=4.583400 | rescues_recorded_for_rule |
| `strict_assembly_summary_sources_unique_rescues` | rescued_rows=175;rescued_mass_pct_all=10.237800 | rescues_recorded_for_rule |
| `strict_source_name_or_assembly_sources_unique_rescues` | rescued_rows=186;rescued_mass_pct_all=11.537300 | rescues_recorded_for_rule |
| `diagnostic_partial_setup_unique_source_name_rescues` | rescued_rows=93;rescued_mass_pct_all=10.019100 | rescues_recorded_for_rule |
| `promotion_decision` | strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed | truth_mapping_threshold_cleared_but_profiles_need_same_namespace_rescore |

## Decision

The setup metadata is present and source names partially resolve GTDB
species. Adding local RefSeq assembly-summary strain/isolate matches
clears the 95% mapped in-scope truth threshold under the strict
source-specific rule, but this still requires selected-default
profiles to be scored in the same namespace before the marine route
can close the holdout gap.

## Outputs

- `results/marine_setup_metadata_source_matches.tsv`
- `results/marine_setup_metadata_truth_upgrade_rules.tsv`
- `results/marine_setup_metadata_truth_upgrade_rescues.tsv`
- `results/marine_setup_metadata_truth_upgrade_audit.tsv`

Promotion decision: `strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed`.
