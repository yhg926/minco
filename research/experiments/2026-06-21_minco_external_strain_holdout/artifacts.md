# Artifacts

## Repo-tracked

- `manifest_strain0_2.tsv`: external strain holdout manifest.
- `scripts/train_current_test_external_calls.py`: train-on-current-panel, test-on-external-panel scorer.
- `summary.tsv`: per-sample call metrics.
- `mean_summary.tsv`: mean metrics by method and train/test split.
- `rf_feature_importance.tsv`: RF feature importances from the nine-sample training panel.
- `simple_filter_diagnostics.tsv`: raw versus generic-name-filtered simple rules and model calls.
- `u_or_s_direct_validtax_missed_gold.tsv`: strain gold candidate rows missed by the best simple rule.
- `u_or_s_direct_validtax_false_positives.tsv`: remaining strain FPs from the best simple rule.
- `u_or_s_direct_validtax_fp_rank_match.tsv`: FP lineage context using direct rank-name matching.
- `genus_capped_rule_diagnostics.tsv`: diagnostic genus-capped relaxed rescue rules.
- `commands.sh`: rerun commands.

## External storage

- CAMI II strain-madness working directory: `/mnt/new3T/minco_cami2_strain_20260621`
- Raw strain FASTQs:
  - `/mnt/new3T/minco_cami2_strain_20260621/sample_0/short_read/2018.09.07_11.43.52_sample_0/reads/anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_strain_20260621/sample_1/short_read/2018.09.07_11.43.52_sample_1/reads/anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_strain_20260621/sample_2/short_read/2018.09.07_11.43.52_sample_2/reads/anonymous_reads.fq.gz`
- Gold profiles: `/mnt/new3T/minco_cami2_strain_20260621/short_read/taxonomic_profile_{0,1,2}.txt`
- minco unfiltered outputs: `/mnt/new3T/minco_cami2_strain_20260621/minco_outputs/strain_sample*_s1000_gtdb_*_zip_unfiltered.tsv`
- Sylph outputs: `/mnt/new3T/minco_cami2_strain_20260621/sylph_sample*/profile.tsv`
- Full external scoring output: `/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2`

## Reference Data

- minco GTDB S1000 reference: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
- Sylph GTDB database: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`
- CAMI taxmap used for species scoring: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
