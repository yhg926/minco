#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ubuntu/yihuiguang/tools/KSSD3mini
EXP=$ROOT/research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse
WORK=/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623
MINCO=$ROOT/minco_core/bin_coden15/minco
CFLAGS_CODEN15='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1'

make -C "$ROOT/minco_core" BINDIR=bin_coden15 CFLAGS="$CFLAGS_CODEN15" clean
make -C "$ROOT/minco_core" BINDIR=bin_coden15 CFLAGS="$CFLAGS_CODEN15"

awk 'substr($0,1,1)=="/"{print}' \
  "$EXP/gtdb232_200709_paths.list" \
  > "$EXP/gtdb232_199924_available_abs_paths.list"
perl -ne 'chomp; print "$_\n" unless -e $_' \
  "$EXP/gtdb232_199924_available_abs_paths.list" \
  > "$EXP/gtdb232_199924_available_abs_paths.missing"
head -32 "$EXP/gtdb232_199924_available_abs_paths.list" \
  > "$EXP/gtdb232_first32_available_abs_paths.list"

rm -rf "$WORK/test32_p1" "$WORK/test32_p4"
/usr/bin/time -v -o "$WORK/logs/test32_p1.time.log" \
  "$MINCO" sketch -p1 -S 2000 \
  -l "$EXP/gtdb232_first32_available_abs_paths.list" \
  -o "$WORK/test32_p1"
/usr/bin/time -v -o "$WORK/logs/test32_p4.time.log" \
  "$MINCO" sketch -p4 -S 2000 \
  -l "$EXP/gtdb232_first32_available_abs_paths.list" \
  -o "$WORK/test32_p4"
cmp "$WORK/test32_p1/minco.ctxobj96" "$WORK/test32_p4/minco.ctxobj96"
cmp "$WORK/test32_p1/minco.ctxobj96.offsets" "$WORK/test32_p4/minco.ctxobj96.offsets"
cmp "$WORK/test32_p1/minco.ctxmeta" "$WORK/test32_p4/minco.ctxmeta"

REF=$WORK/sketch_T15_S2000_199924_native
rm -rf "$REF"
/usr/bin/time -v -o "$WORK/logs/sketch_T15_S2000_199924_native.time.log" \
  "$MINCO" sketch -p8 -S 2000 \
  -l "$EXP/gtdb232_199924_available_abs_paths.list" \
  -o "$REF" \
  > "$WORK/logs/sketch_T15_S2000_199924_native.stdout.log" \
  2> "$WORK/logs/sketch_T15_S2000_199924_native.stderr.log"

/usr/bin/time -v -o "$WORK/logs/index_T15_S2000_199924_native.time.log" \
  "$MINCO" sketch -i "$REF" \
  > "$WORK/logs/index_T15_S2000_199924_native.stdout.log" \
  2> "$WORK/logs/index_T15_S2000_199924_native.stderr.log"

MARKER=$WORK/sketch_T15_S2000_199924_ctxmarker
rm -rf "$MARKER"
/usr/bin/time -v -o "$WORK/logs/markerdb_ctx_T15_S2000_199924.fixed_order.time.log" \
  "$MINCO" set --uniq_union --markerdb-ctx -p8 -o "$MARKER" "$REF" \
  > "$WORK/logs/markerdb_ctx_T15_S2000_199924.fixed_order.stdout.log" \
  2> "$WORK/logs/markerdb_ctx_T15_S2000_199924.fixed_order.stderr.log"

/usr/bin/time -v -o "$WORK/logs/index_T15_S2000_199924_ctxmarker.time.log" \
  "$MINCO" sketch -i "$MARKER" \
  > "$WORK/logs/index_T15_S2000_199924_ctxmarker.stdout.log" \
  2> "$WORK/logs/index_T15_S2000_199924_ctxmarker.stderr.log"

for sample in 0 1 2; do
  READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_${sample}/2017.12.29_11.37.26_sample_${sample}/reads/anonymous_reads.fq.gz
  OUT=$WORK/toymouse_sample${sample}_T15_s2000_ctxmarker_split_naive_product.tsv
  /usr/bin/time -v -o "$WORK/logs/readwise_sample${sample}_T15_s2000_ctxmarker.time.log" \
    "$MINCO" ani -p16 -r "$MARKER" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 -o "$OUT" \
    > "$WORK/logs/readwise_sample${sample}_T15_s2000_ctxmarker.stdout.log" \
    2> "$WORK/logs/readwise_sample${sample}_T15_s2000_ctxmarker.stderr.log"
done

python3 "$EXP/score_coden15_toymouse.py"
python3 "$EXP/score_coden15_toymouse_three_samples.py"
python3 "$EXP/score_coden15_abundance_rescue.py"
python3 "$EXP/analyze_coden15_mechanism.py"
python3 "$EXP/optimize_coden15_gate.py"
python3 "$EXP/optimize_coden15_params_fast.py"
python3 "$EXP/score_coden15_formula_af_gate.py"
