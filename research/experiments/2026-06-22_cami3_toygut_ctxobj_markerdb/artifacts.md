# Artifacts

## Repository files

- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/NOTE.md`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/commands.sh`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/score_cami3_ctxobj_markerdb.py`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/mean_summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/total_summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/abundance.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/abundance_renorm_mean.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/abundance_raw_mean.tsv`

## Inputs

- Ctx+obj S2000 markerdb: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`
- Ctx-only S2000 baseline markerdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Previous CAMI3 ctx-only and Sylph run directory: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622`
- CAMI3 ToyGut sample0 reads: `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`
- CAMI3 ToyGut sample1 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz`
- CAMI3 ToyGut sample2 reads: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz`
- CAMI3 ToyGut truth profiles: `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt`, `taxonomic_profile_1.txt`, `taxonomic_profile_2.txt`
- Sylph database from previous comparison: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`

## Generated outputs

- Run directory: `/tmp/cami3_toygut_ctxobj_markerdb_20260622`
- Ctx+obj sample0 profile: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s0_ctxobj_product0.tsv`
- Ctx+obj sample1 profile: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s1_ctxobj_product0.tsv`
- Ctx+obj sample2 profile: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s2_ctxobj_product0.tsv`
- Combined score table: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/cami3_toygut_ctxobj_markerdb_scores.tsv`
- Combined abundance table: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/cami3_toygut_ctxobj_markerdb_abundance.tsv`
- FP/FN detail table: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/cami3_toygut_ctxobj_markerdb_details.tsv`
- Selected species table: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/cami3_toygut_ctxobj_markerdb_selected_species.tsv`

## Logs

- Ctx+obj sample0 time log: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s0_ctxobj_product0.time.log`
- Ctx+obj sample1 time log: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s1_ctxobj_product0.time.log`
- Ctx+obj sample2 time log: `/tmp/cami3_toygut_ctxobj_markerdb_20260622/minco_s2_ctxobj_product0.time.log`

## Storage note

Large generated files are under `/tmp`. Recreate them with `commands.sh` if `/tmp` is cleaned.
