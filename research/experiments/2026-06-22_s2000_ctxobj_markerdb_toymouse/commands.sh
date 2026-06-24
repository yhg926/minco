#!/usr/bin/env bash
set -euo pipefail

RUN=/tmp/gtdb232_s2000_ctxobj_marker_20260622
SRC=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
OUT=$RUN/sketch_T_S2000_aaf003_dedup_ctxobjmarker
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz

mkdir -p "$RUN"

/usr/bin/time -v -o "$RUN/markerdb_ctxobj.time.log" \
  bin/minco set --uniq_union --markerdb -p 8 \
    -o "$OUT" "$SRC" \
  > "$RUN/markerdb_ctxobj.log" 2>&1

bin/minco sketch --psmp "$OUT" > "$RUN/markerdb_ctxobj.psmp.tsv"

/usr/bin/time -v -o "$RUN/index_ctxobjmarker.time.log" \
  bin/minco sketch -i "$OUT" \
  > "$RUN/index_ctxobjmarker.log" 2>&1

/usr/bin/time -v -o "$RUN/toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.time.log" \
  bin/minco ani -p16 -r "$OUT" \
    --qraw "$READS" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.tsv"

python3 research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/score_ctxobj_markerdb.py
