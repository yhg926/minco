#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
WORK=/tmp/gtdb232_s2000_dedup_marker.qKJofv
SRC=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz

cd "$REPO"

/usr/bin/time -v -o "$WORK/downsample_s2000.time.log" \
  bin/minco set --downsample -S 2000 \
    -o "$WORK/sketch_T_S2000_anno" \
    "$SRC" \
  > "$WORK/downsample_s2000.log" 2>&1

/usr/bin/time -v -o "$WORK/index_s2000.time.log" \
  bin/minco sketch -i "$WORK/sketch_T_S2000_anno" \
  > "$WORK/index_s2000.log" 2>&1

/usr/bin/time -v -o "$WORK/dedup_aaf003.time.log" \
  bin/minco matrix --format dedup-plan -m aaf --cut 0.03 -p8 \
    --keep-out "$WORK/keep.aaf003.txt" \
    --remove-out "$WORK/remove.aaf003.txt" \
    -o "$WORK/dedup_plan.aaf003.tsv" \
    "$WORK/sketch_T_S2000_anno" \
  > "$WORK/dedup_aaf003.log" 2>&1

/usr/bin/time -v -o "$WORK/remove_dups_s2000.time.log" \
  bin/minco sketch --remove "$WORK/remove.aaf003.txt" \
    -o "$WORK/sketch_T_S2000_aaf003_dedup" \
    "$WORK/sketch_T_S2000_anno" \
  > "$WORK/remove_dups_s2000.log" 2>&1

/usr/bin/time -v -o "$WORK/markerdb_ctx_s2000_dedup.time.log" \
  bin/minco set --uniq_union --markerdb-ctx \
    -o "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    "$WORK/sketch_T_S2000_aaf003_dedup" \
  > "$WORK/markerdb_ctx_s2000_dedup.log" 2>&1

/usr/bin/time -v -o "$WORK/index_ctxmarker_s2000_dedup.time.log" \
  bin/minco sketch -i "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
  > "$WORK/index_ctxmarker_s2000_dedup.log" 2>&1

/usr/bin/time -v -o "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip.time.log" \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-unique \
    --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip_unfiltered.tsv" \
  > "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip.log" 2>&1

/usr/bin/time -v -o "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_split_zip.time.log" \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_split_zip_unfiltered.tsv" \
  > "$WORK/toymouse_sample0_s2000_dedup_ctxmarker_split_zip.log" 2>&1

python3 research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/score_sourceaware_direct.py

MINCO_READWISE_TRACE_REF=GCF_018987235 \
MINCO_READWISE_TRACE_OUT=research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/l_crispatus_s2000_marker_split_readwise_assignment_trace.tsv \
/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_perread_trace_lcri.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani zip-aaf \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_zip_unfiltered_perread_trace_lcri.tsv

python3 research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/join_l_crispatus_trace_to_truth.py

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_unfiltered.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_unfiltered.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani zip-aaf \
    --readwise-ctx-filter poisson-depth \
    --readwise-fake-threshold 1.30103 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter poisson-depth \
    --readwise-fake-threshold 1.30103 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter poisson-product \
    --readwise-fake-threshold 1.30103 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_nb_p005.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter product-nb \
    --readwise-fake-threshold 1.30103 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_nb_p005.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_nb_reliable_abundance_p005.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter product-nb \
    --readwise-fake-threshold 1.30103 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_nb_reliable_abundance_p005.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_topfrac025.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter product-topfrac \
    --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_topfrac025.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median \
    --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv

/usr/bin/time -v \
  -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product1_topfrac_median025.time.log \
  bin/minco ani -p16 \
    -r "$WORK/sketch_T_S2000_aaf003_dedup_ctxmarker" \
    --qraw "$READS" \
    --query-density ref \
    --abundance-est depth \
    --readwise-profile-only \
    --readwise-assign best-diff-split \
    --readwise-ani naive \
    --readwise-ctx-filter product1-topfrac-median \
    --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/toymouse_sample0_s2000_marker_split_naive_product1_topfrac_median025.tsv

# Authoritative benchmark truth and scoring. Source-aware scoring below is diagnostic.
python3 research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/build_and_score_gtdb_ground_truth.py

python3 research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/score_reliable_truncated_depth_gate.py

python3 research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/score_sourceaware_direct.py
