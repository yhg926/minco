#!/usr/bin/env bash
set -euo pipefail

EXP=research/experiments/2026-06-23_cami3_toygut_coden15_formula_af
MINCO=/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/bin_coden15/minco
REF=/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker

mkdir -p "$EXP/logs"

for sample in 0 1 2; do
  case "$sample" in
    0) READS=/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz ;;
    1) READS=/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz ;;
    2) READS=/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz ;;
  esac

  OUT=$EXP/cami3_sample${sample}_coden15_ctxmarker_split_naive_product.tsv
  /usr/bin/time -v -o "$EXP/logs/cami3_sample${sample}_coden15.time.log" \
    "$MINCO" ani -p16 -r "$REF" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 -o "$OUT" \
    > "$EXP/logs/cami3_sample${sample}_coden15.stdout.log" \
    2> "$EXP/logs/cami3_sample${sample}_coden15.stderr.log"
done

python3 "$EXP/score_cami3_coden15_formula_af.py"
