#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
RUN=/tmp/cami2_toymouse_density_block_20260622
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
BASE=/mnt/new3T/minco_cami2_toymouse_20260621

cd "$REPO"
mkdir -p "$RUN"

# Sample0: compare N=0 and N=1 directly.
for n in 0 1; do
  /usr/bin/time -v -o "$RUN/sample0_ctxmarker_block${n}.time.log" \
    bin/minco ani -p16 -r "$REF" \
      --qraw "$BASE/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      --density-block-ctx "$n" \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/sample0_ctxmarker_block${n}.tsv"
done

sha256sum "$RUN/sample0_ctxmarker_block0.tsv" "$RUN/sample0_ctxmarker_block1.tsv"
cmp -s "$RUN/sample0_ctxmarker_block0.tsv" "$RUN/sample0_ctxmarker_block1.tsv"

# Samples1-2: N=0 per-read benchmark. N=1 is the same code path as N=0.
for sample in 1 2; do
  /usr/bin/time -v -o "$RUN/sample${sample}_ctxmarker_block0.time.log" \
    bin/minco ani -p16 -r "$REF" \
      --qraw "$BASE/sample_${sample}/2017.12.29_11.37.26_sample_${sample}/reads/anonymous_reads.fq.gz" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      --density-block-ctx 0 \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/sample${sample}_ctxmarker_block0.tsv"
done

# Scoring used the helper functions from:
# research/experiments/2026-06-22_cami2_toymouse_more_gtdb/score_more_toymouse_gtdb.py
