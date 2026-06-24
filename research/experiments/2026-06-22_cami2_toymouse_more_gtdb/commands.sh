#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
BASE=/mnt/new3T/minco_cami2_toymouse_20260621
RUN=/tmp/cami2_toymouse_more_gtdb_20260622
CTX=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
CTXOBJ=/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb

cd "$REPO"
mkdir -p "$RUN"

# Download the additional CAMI II Toy Mouse Gut short-read samples.
for i in 1 2; do
  curl -L --fail -C - --retry 8 --retry-all-errors --retry-delay 5 \
    -o "$BASE/2017.12.29_11.37.26_sample_${i}_reads.tar" \
    "https://frl.publisso.de/data/frl%3A6421672/dataset/2017.12.29_11.37.26_sample_${i}_reads.tar"
  tar -xf "$BASE/2017.12.29_11.37.26_sample_${i}_reads.tar" \
    -C "$BASE" --one-top-level="sample_${i}"
done

# Build GTDB truth files and show missing profile outputs before profiling.
python3 research/experiments/2026-06-22_cami2_toymouse_more_gtdb/score_more_toymouse_gtdb.py

# Run current best ctx-only markerdb and the ctx+obj markerdb for new samples.
for i in 1 2; do
  READS="$BASE/sample_${i}/2017.12.29_11.37.26_sample_${i}/reads/anonymous_reads.fq.gz"
  /usr/bin/time -v -o "$RUN/toymouse_sample${i}_ctxmarker_split_naive_product_topfrac_median025.time.log" \
    bin/minco ani -p16 -r "$CTX" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/toymouse_sample${i}_ctxmarker_split_naive_product_topfrac_median025.tsv"
  /usr/bin/time -v -o "$RUN/toymouse_sample${i}_ctxobjmarker_split_naive_product_topfrac_median025.time.log" \
    bin/minco ani -p16 -r "$CTXOBJ" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/toymouse_sample${i}_ctxobjmarker_split_naive_product_topfrac_median025.tsv"
done

# Run Sylph for the same two new samples.
for i in 1 2; do
  READS="$BASE/sample_${i}/2017.12.29_11.37.26_sample_${i}/reads/anonymous_reads.fq.gz"
  mkdir -p "$BASE/sylph_sample${i}"
  /usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 \
    -r "$READS" -d "$BASE/sylph_sample${i}"
  /usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph profile -t 16 \
    "$SYLPH_DB" "$BASE/sylph_sample${i}/anonymous_reads.fq.gz.sylsp" \
    -o "$BASE/sylph_sample${i}/profile.tsv"
done

# Final GTDB-species scoring.
python3 research/experiments/2026-06-22_cami2_toymouse_more_gtdb/score_more_toymouse_gtdb.py
