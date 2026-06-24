#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
RUN=/tmp/cami3_toygut_abundance_rescue_20260623
CTX=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
CTXOBJ=/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker

mkdir -p "$RUN"
cd "$REPO"

for sample in 0 1 2; do
  if [ "$sample" = 0 ]; then
    READS=/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz
  else
    READS=/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_${sample}_reads/anonymous_reads.fq.gz
  fi

  /usr/bin/time -v -o "$RUN/minco_s${sample}_ctxmarker_current_binary.time.log" \
    minco_core/bin/minco ani -p16 -r "$CTX" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/minco_s${sample}_ctxmarker_current_binary.tsv"

  /usr/bin/time -v -o "$RUN/minco_s${sample}_ctxobj_current_binary.time.log" \
    minco_core/bin/minco ani -p16 -r "$CTXOBJ" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$RUN/minco_s${sample}_ctxobj_current_binary.tsv"
done

python3 research/experiments/2026-06-23_cami3_toygut_abundance_rescue/score_cami3_abundance_rescue.py
