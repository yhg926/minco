# Artifacts

## Repository files

- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/NOTE.md`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/commands.sh`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/score_ctxobj_markerdb.py`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/summary.tsv`

## Inputs

- Deduplicated S2000 MinCO sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Current best ctx-only markerdb baseline: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Toy Mouse sample0 reads: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- GTDB truth species profile: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/mouse0_gtdb_species_profile.tsv`
- Current best ctx-only profile table: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv`

## Generated outputs

- Ctx+obj markerdb sketch: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`
- Ctx+obj markerdb psmp table: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/markerdb_ctxobj.psmp.tsv`
- Ctx+obj Toy Mouse readwise profile: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.tsv`
- GTDB score table: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/ctxobj_markerdb_gtdb_scores.tsv`
- FP/FN detail table: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/ctxobj_markerdb_gtdb_details.tsv`
- Abundance comparison table: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/ctxobj_markerdb_gtdb_abundance.tsv`
- Gate sweep table: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/ctxobj_markerdb_gate_sweep.tsv`

## Logs

- Markerdb build log: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/markerdb_ctxobj.log`
- Markerdb build time log: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/markerdb_ctxobj.time.log`
- Index log: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/index_ctxobjmarker.log`
- Index time log: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/index_ctxobjmarker.time.log`
- Toy Mouse run time log: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.time.log`

## Storage note

The large generated artifacts are under `/tmp`. Recreate them with `commands.sh` if `/tmp` is cleaned.
