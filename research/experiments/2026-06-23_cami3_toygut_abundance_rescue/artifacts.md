# Artifacts

## Repo-Tracked Experiment Files

- `NOTE.md`: experiment question, method, results, and conclusion.
- `commands.sh`: exact rerun recipe for MinCO profiles and scoring.
- `score_cami3_abundance_rescue.py`: scorer for old MinCO abundance, robust effective depth, robust depth plus intra-genus rescue, and Sylph.
- `summary.tsv`: mean presence and renormalized abundance metrics joined by method.
- `mean_presence.tsv`: mean presence/absence metrics over samples 0-2.
- `mean_abundance_renorm.tsv`: mean abundance metrics with predictions renormalized over reported bacterial species.
- `sample_presence.tsv`: per-sample presence/absence metrics.
- `sample_abundance.tsv`: per-sample abundance metrics.
- `selected_species.tsv`: collapsed selected species table used for scoring.
- `fp_fn_details.tsv`: false-positive and false-negative details by sample and method.
- `rescued_species.tsv`: rescued rows from the strict intra-genus rescue rule; header only because CAMI3 had zero rescues.

## Large Temporary Outputs

- Run directory: `/tmp/cami3_toygut_abundance_rescue_20260623`.
- Ctx-marker MinCO profiles:
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s0_ctxmarker_current_binary.tsv`
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s1_ctxmarker_current_binary.tsv`
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s2_ctxmarker_current_binary.tsv`
- Ctx+obj MinCO profiles:
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s0_ctxobj_current_binary.tsv`
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s1_ctxobj_current_binary.tsv`
  - `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s2_ctxobj_current_binary.tsv`
- Time logs: `/tmp/cami3_toygut_abundance_rescue_20260623/*.time.log`.

## External Inputs

- CAMI3 ToyGut sample0 reads: `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`.
- CAMI3 ToyGut sample1 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz`.
- CAMI3 ToyGut sample2 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz`.
- CAMI3 truth profiles: `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_{0,1,2}.txt`.
- Sylph sample0 profile: `/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv`.
- Sylph sample1 profile: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv`.
- Sylph sample2 profile: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample2/profile.tsv`.
- Ctx-marker refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`.
- Ctx+obj marker refdb: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`.
