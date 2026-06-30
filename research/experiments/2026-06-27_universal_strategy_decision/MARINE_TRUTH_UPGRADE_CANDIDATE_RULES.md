# Marine Truth Upgrade Candidate Rules

Date: 2026-06-29

This cached audit tests whether deterministic-looking fallback rules can
raise marine GTDB truth mapping to the 95% release threshold. It does not
rerun profilers and does not change the accepted release truth policy.

## Audit

| Metric | Value | Decision |
|---|---|---|
| `candidate_rules_tested` | accepted_exact_binomial,diagnostic_unique_ambiguous_epithet,unsafe_all_ambiguous_assigned | diagnostic_rules_only |
| `accepted_exact_binomial_summary` | min_mapped_pct_bacteria_archaea=91.005806;ready_samples=0/10 | release_rule_fails_threshold |
| `diagnostic_unique_ambiguous_epithet_summary` | min_mapped_pct_bacteria_archaea=91.292431;ready_samples=2/10 | diagnostic_not_release_rule |
| `unsafe_all_ambiguous_assigned_summary` | min_mapped_pct_bacteria_archaea=99.411153;ready_samples=10/10 | diagnostic_not_release_rule |
| `diagnostic_unique_epithet_rescues` | rescued_rows=51;rescued_mass_pct_all=3.146300;min_mapped_pct_bacteria_archaea=91.292431;ready_samples=2/10 | candidate_rules_do_not_clear_release_threshold |
| `promotion_decision` | do_not_promote_marine_truth_upgrade | truth_ambiguity_not_resolved_by_current_candidate_rules |

## Main Remaining Blockers

| Method | Taxid | Name | Total mass | Samples | Candidate count |
|---|---:|---|---:|---:|---:|
| ambiguous_taxid | 1263979 | Candidatus Endolissoclinum faulkneri | 5.4357 | 10 | 2 |
| ambiguous_taxid | 467094 | Ilumatobacter coccineus | 2.8633 | 10 | 2 |
| ambiguous_taxid | 1219 | Prochlorococcus marinus | 2.3293 | 10 | 55 |
| ambiguous_taxid | 2296 | Desulfobacterium autotrophicum | 1.5470 | 9 | 2 |
| ambiguous_taxid | 1491 | Clostridium botulinum | 1.2305 | 10 | 13 |
| ambiguous_taxid | 28229 | Colwellia psychrerythraea | 0.9454 | 10 | 4 |
| ambiguous_taxid | 73141 | Thiocystis violascens | 0.6399 | 10 | 2 |
| ambiguous_taxid | 316279 | Synechococcus sp. CC9902 | 0.5486 | 10 | 2 |
| ambiguous_taxid | 669041 | Tenacibaculum dicentrarchi | 0.5314 | 8 | 2 |
| ambiguous_taxid | 324767 | Bacillus infantis | 0.4745 | 10 | 5 |
| unmapped | 1250153 | Maribacter sp. MAR_2009_60 | 0.4391 | 8 | 0 |
| ambiguous_taxid | 316 | Pseudomonas stutzeri | 0.4253 | 10 | 51 |
| ambiguous_taxid | 154 | Spirochaeta thermophila | 0.3922 | 7 | 2 |
| ambiguous_taxid | 89184 | Ruegeria pomeroyi | 0.3705 | 9 | 2 |
| ambiguous_taxid | 713887 | Candidatus Atelocyanobacterium thalassa | 0.3451 | 10 | 3 |
| ambiguous_taxid | 2743 | Marinobacter hydrocarbonoclasticus | 0.3291 | 9 | 4 |
| unmapped | 400668 | Marinomonas sp. MWYL1 | 0.2912 | 7 | 0 |
| ambiguous_taxid | 1336795 | Formosa sp. Hel3_A1_48 | 0.2846 | 9 | 2 |
| ambiguous_taxid | 111501 | Muricauda ruestringensis | 0.2502 | 7 | 2 |
| ambiguous_taxid | 1155739 | Moorea producens | 0.2461 | 8 | 2 |

## Decision

The current marine route should not be promoted from cached truth alone.
The relaxed epithet fallback is diagnostic only and still leaves most
samples below 95% mapped Bacteria/Archaea truth mass.

## Outputs

- `results/marine_truth_upgrade_candidate_rules.tsv`
- `results/marine_truth_upgrade_candidate_rescues.tsv`
- `results/marine_truth_upgrade_remaining_blockers.tsv`
- `results/marine_truth_upgrade_candidate_audit.tsv`
