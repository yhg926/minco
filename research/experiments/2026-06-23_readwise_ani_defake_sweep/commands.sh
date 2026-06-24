#!/usr/bin/env bash
set -euo pipefail

# Build adjusted-ANI coden11 binary.
make -C minco_core \
  BINDIR=bin_adjani \
  OBJDIR=obj_adjani \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DMINCO_REPORT_FILTERED_READWISE_ANI=1' \
  -j4

# Build adjusted-ANI coden15 binary.
make -C minco_core \
  BINDIR=bin_coden15_adjani \
  OBJDIR=obj_coden15_adjani \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1 -DMINCO_REPORT_FILTERED_READWISE_ANI=1' \
  -j4

# Run product0 top-fraction median replacement for p75, p80, p85, p90.
# The actual run was started interactively so progress could be monitored.

EXP=research/experiments/2026-06-23_readwise_ani_defake_sweep
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
REF11=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
REF15=/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker

for frac in 0.25 0.20 0.15 0.10; do
  label=$(python3 -c "print(f'{int(round((1-$frac)*100)):03d}')")
  /usr/bin/time -v -o "$EXP/logs/c11_p${label}.time.log" \
    minco_core/bin_adjani/minco ani -p8 \
      -r "$REF11" \
      --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold "$frac" \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/results/toymouse_sample0_c11_adjusted_ani_p${label}.tsv" \
    > "$EXP/logs/c11_p${label}.stdout.log" \
    2> "$EXP/logs/c11_p${label}.stderr.log"
done

# c15 product0 is much slower. Only p75 was kept as a baseline; the planned
# remaining product0 c15 loop was interrupted after p75 finished.
/usr/bin/time -v -o "$EXP/logs/c15_p075.time.log" \
  minco_core/bin_coden15_adjani/minco ani -p16 \
    -r "$REF15" \
    --qraw "$READS" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$EXP/results/toymouse_sample0_c15_adjusted_ani_p075.tsv" \
  > "$EXP/logs/c15_p075.stdout.log" \
  2> "$EXP/logs/c15_p075.stderr.log"

python3 "$EXP/score_defake_sweep.py"

# Targeted product1 follow-up. This ranks contexts by (best_diff + 1) * depth,
# so high-depth exact-hit contexts can also be median-replaced. This is the
# plausible path when product0 leaves readwise naive ANI saturated.
for frac in 0.25 0.20 0.15 0.10; do
  label=$(python3 -c "print(f'{int(round((1-$frac)*100)):03d}')")
  /usr/bin/time -v -o "$EXP/logs/c11_product1_p${label}.time.log" \
    minco_core/bin_adjani/minco ani -p8 \
      -r "$REF11" \
      --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product1-topfrac-median --readwise-fake-threshold "$frac" \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/results/toymouse_sample0_c11_adjusted_ani_product1_p${label}.tsv" \
    > "$EXP/logs/c11_product1_p${label}.stdout.log" \
    2> "$EXP/logs/c11_product1_p${label}.stderr.log"
done

for item in 0.20:080 0.15:085 0.10:090; do
  frac=${item%%:*}
  label=${item##*:}
  /usr/bin/time -v -o "$EXP/logs/c15_product1_p${label}.time.log" \
    minco_core/bin_coden15_adjani/minco ani -p16 \
      -r "$REF15" \
      --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product1-topfrac-median --readwise-fake-threshold "$frac" \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/results/toymouse_sample0_c15_adjusted_ani_product1_p${label}.tsv" \
    > "$EXP/logs/c15_product1_p${label}.stdout.log" \
    2> "$EXP/logs/c15_product1_p${label}.stderr.log"
done

python3 "$EXP/score_defake_sweep.py"
