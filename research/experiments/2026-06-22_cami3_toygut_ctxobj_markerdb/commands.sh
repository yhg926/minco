#!/usr/bin/env bash
set -euo pipefail

RUN=/tmp/cami3_toygut_ctxobj_markerdb_20260622
REF=/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker

mkdir -p "$RUN"

/usr/bin/time -v -o "$RUN/minco_s0_ctxobj_product0.time.log" \
  bin/minco ani -p16 -r "$REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s0_ctxobj_product0.tsv"

/usr/bin/time -v -o "$RUN/minco_s1_ctxobj_product0.time.log" \
  bin/minco ani -p16 -r "$REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s1_ctxobj_product0.tsv"

/usr/bin/time -v -o "$RUN/minco_s2_ctxobj_product0.time.log" \
  bin/minco ani -p16 -r "$REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s2_ctxobj_product0.tsv"

python3 research/experiments/2026-06-22_cami3_toygut_ctxobj_markerdb/score_cami3_ctxobj_markerdb.py
