#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
EXP=$REPO/research/experiments/2026-06-22_minco_sylph_abundance_model
BASE=/mnt/new3T/minco_cami2_toymouse_20260621
CTX=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
CTXOBJ=/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker
S1000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
S2000_FULL=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup

cd "$REPO"

# Inspect Sylph source behavior: first-pass stats, winner-map reassignment,
# second-pass stats, and abundance as Eff_cov / sum(Eff_cov).
sed -n '260,330p' /home/ubuntu/yihuiguang/tools/sylph/src/contain.rs
sed -n '349,455p' /home/ubuntu/yihuiguang/tools/sylph/src/contain.rs
sed -n '680,760p' /home/ubuntu/yihuiguang/tools/sylph/src/contain.rs

# Build MinCO after adding reliable median and robust effective abundance columns.
make -C minco_core -j4

# Full-sketch diagnostics. These did not beat markerdb because current MinCO
# gate allowed many weak false rows.
for i in 0 1 2; do
  READS="$BASE/sample_${i}/2017.12.29_11.37.26_sample_${i}/reads/anonymous_reads.fq.gz"
  /usr/bin/time -v -o "$EXP/toymouse_sample${i}_s1000_full_split_naive_product_topfrac_median025.time.log" \
    bin/minco ani -p8 -r "$S1000" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/toymouse_sample${i}_s1000_full_split_naive_product_topfrac_median025.tsv"

  /usr/bin/time -v -o "$EXP/toymouse_sample${i}_s2000_dedup_full_split_naive_product_topfrac_median025.time.log" \
    bin/minco ani -p8 -r "$S2000_FULL" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/toymouse_sample${i}_s2000_dedup_full_split_naive_product_topfrac_median025.tsv"
done

# Full/shared-context best-diff-unique diagnostics. This is closer to keeping
# shared contexts in the refdb, but it did not beat markerdb plus rescue here.
for i in 0 1 2; do
  READS="$BASE/sample_${i}/2017.12.29_11.37.26_sample_${i}/reads/anonymous_reads.fq.gz"
  /usr/bin/time -v -o "$EXP/toymouse_sample${i}_s2000_dedup_full_unique_naive_product_topfrac_median025.time.log" \
    minco_core/bin/minco ani -p8 -r "$S2000_FULL" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-unique --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/toymouse_sample${i}_s2000_dedup_full_unique_naive_product_topfrac_median025.tsv"
done

# Rebuilt markerdb profiles with Reliable_Ref_hit_median_depth.
for i in 0 1 2; do
  READS="$BASE/sample_${i}/2017.12.29_11.37.26_sample_${i}/reads/anonymous_reads.fq.gz"
  /usr/bin/time -v -o "$EXP/toymouse_sample${i}_ctxmarker_median_col.time.log" \
    minco_core/bin/minco ani -p8 -r "$CTX" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/toymouse_sample${i}_ctxmarker_median_col.tsv"

  /usr/bin/time -v -o "$EXP/toymouse_sample${i}_ctxobjmarker_median_col.time.log" \
    minco_core/bin/minco ani -p8 -r "$CTXOBJ" --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$EXP/toymouse_sample${i}_ctxobjmarker_median_col.tsv"
done

# Final smoke-test with output-level robust effective abundance columns.
/usr/bin/time -v -o "$EXP/toymouse_sample0_ctxmarker_effective_abundance.time.log" \
  minco_core/bin/minco ani -p8 -r "$CTX" \
    --qraw "$BASE/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o "$EXP/toymouse_sample0_ctxmarker_effective_abundance.tsv"

python3 "$EXP/evaluate_abundance_estimators.py"
python3 "$EXP/score_effective_abundance.py"
python3 "$EXP/score_intragenus_winner_rescue.py"
