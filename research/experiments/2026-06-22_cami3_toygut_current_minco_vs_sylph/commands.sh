#!/usr/bin/env bash
set -euo pipefail

RUN=/tmp/cami3_toygut_current_minco_vs_sylph_20260622
MINCO_REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb

mkdir -p "$RUN"/sylph_sample{0,1,2}

/usr/bin/time -v -o "$RUN/minco_s0_current_product0.time.log" \
  bin/minco ani -p16 -r "$MINCO_REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s0_current_product0.tsv"

/usr/bin/time -v -o "$RUN/minco_s1_current_product0.time.log" \
  bin/minco ani -p16 -r "$MINCO_REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s1_current_product0.tsv"

/usr/bin/time -v -o "$RUN/minco_s2_current_product0.time.log" \
  bin/minco ani -p16 -r "$MINCO_REF" \
    --qraw /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/minco_s2_current_product0.tsv"

/usr/bin/time -v -o "$RUN/sylph_sample1/sketch.time.log" \
  sylph sketch -t 16 -d "$RUN/sylph_sample1" --reads \
    /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz

/usr/bin/time -v -o "$RUN/sylph_sample1/profile.time.log" \
  sylph profile -t 16 -o "$RUN/sylph_sample1/profile.tsv" \
    "$SYLPH_DB" "$RUN/sylph_sample1/anonymous_reads.fq.gz.sylsp"

/usr/bin/time -v -o "$RUN/sylph_sample2/sketch.time.log" \
  sylph sketch -t 16 -d "$RUN/sylph_sample2" --reads \
    /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz

/usr/bin/time -v -o "$RUN/sylph_sample2/profile.time.log" \
  sylph profile -t 16 -o "$RUN/sylph_sample2/profile.tsv" \
    "$SYLPH_DB" "$RUN/sylph_sample2/anonymous_reads.fq.gz.sylsp"

python3 research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph/score_cami3_toygut.py
