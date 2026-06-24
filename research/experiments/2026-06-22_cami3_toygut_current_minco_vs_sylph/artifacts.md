# Artifacts

## Inputs

- CAMI3 ToyGut sample0 reads: `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`
- CAMI3 ToyGut sample1 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz`
- CAMI3 ToyGut sample2 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz`
- CAMI3 taxonomic truth profiles: `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_{0,1,2}.txt`
- Current best MinCO S2000 global context markerdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Sylph GTDB r226 database: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`
- Existing Sylph sample0 profile: `/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv`

## Outputs

- Run directory: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622`
- MinCO sample0 output: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/minco_s0_current_product0.tsv`
- MinCO sample1 output: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/minco_s1_current_product0.tsv`
- MinCO sample2 output: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/minco_s2_current_product0.tsv`
- Sylph sample1 profile: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv`
- Sylph sample2 profile: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample2/profile.tsv`
- Score table: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/cami3_toygut_minco_current_vs_sylph_scores.tsv`
- Abundance table: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/cami3_toygut_minco_current_vs_sylph_abundance.tsv`
- FP/FN detail table: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/cami3_toygut_minco_current_vs_sylph_details.tsv`
- Collapsed selected species table: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/cami3_toygut_minco_current_vs_sylph_selected_species.tsv`
- Repo copy of score summary: `summary.tsv`
- Repo copy of abundance summary: `abundance.tsv`
- Scoring script: `score_cami3_toygut.py`
- Command log: `commands.sh`
