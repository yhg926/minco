#!/usr/bin/env bash
set -euo pipefail

REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
FASTQ=/tmp/cami_marine_sample0_reads.fq.gz
OUT=/tmp/minco_readwise_assign_20260620

make -C minco_core -j8

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r "$REF" --qraw "$FASTQ" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 --readwise-assign best-diff \
  -o "$OUT/s10000_bestdiff_unfiltered.tsv"

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r "$REF" --qraw "$FASTQ" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 --readwise-assign best-diff-split \
  -o "$OUT/s10000_bestdiffsplit_unfiltered.tsv"

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r "$REF" --qraw "$FASTQ" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 --readwise-assign best-diff-unique \
  --readwise-ani naive \
  -o "$OUT/s10000_bestdiffunique_unfiltered.tsv"

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r "$REF" --qraw "$FASTQ" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 \
  -o "$OUT/s10000_unique_zipaaf_unfiltered_c.tsv"

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r "$REF" --qraw "$FASTQ" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  -m0 -f0.05 -n0.94 -t10 \
  -o "$OUT/s10000_unique_zipaaf_f0.05_n0.94_t10.tsv"

python3 research/experiments/2026-06-20_minco_readwise_assignment_fix/scripts/grid_profile_thresholds.py \
  --minco "$OUT/s10000_unique_zipaaf_unfiltered_c.tsv" \
  --label s10000_unique_zipaaf_c \
  --outdir "$OUT"
