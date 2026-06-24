#!/usr/bin/env bash
set -euo pipefail

# Experiment: Hybrid Full Sketch Unique Marker Fallback
# Date: 2026-06-23
# Project: KSSD3mini
cd /home/ubuntu/yihuiguang/tools/KSSD3mini

# 1. Syntax-check the hybrid scoring scripts.
python3 -m py_compile \
  research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_mouse.py \
  research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_mouse_baseline_plus.py \
  research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_cami3.py

# 2. Run the exploratory mouse size-cutoff sweep.
/usr/bin/time -f 'elapsed=%E maxrss=%MKB' \
  python3 research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_mouse.py \
  > research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_mouse.stdout.tsv \
  2> research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_mouse.stderr.txt

# 3. Generate CAMI3 full S2000 readwise outputs used by the hybrid fallback.
mkdir -p /tmp/cami3_toygut_hybrid_full_20260623
for sample in 0 1 2; do
  if [ "$sample" = 0 ]; then
    READS=/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz
  else
    READS=/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_${sample}_reads/anonymous_reads.fq.gz
  fi

  /usr/bin/time -v -o /tmp/cami3_toygut_hybrid_full_20260623/minco_s${sample}_full_product0.time.log \
    minco_core/bin/minco ani -p8 \
      -r /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup \
      --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o /tmp/cami3_toygut_hybrid_full_20260623/minco_s${sample}_full_product0.tsv
done

# 4. Run the final mouse baseline-plus hybrid scorer.
/usr/bin/time -f 'elapsed=%E maxrss=%MKB' \
  python3 research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_mouse_baseline_plus.py \
  > research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_mouse_baseline_plus.stdout.tsv \
  2> research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_mouse_baseline_plus.stderr.txt

# 5. Run the CAMI3 baseline-plus hybrid scorer.
/usr/bin/time -f 'elapsed=%E maxrss=%MKB' \
  python3 research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/score_hybrid_cami3.py \
  > research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_cami3.stdout.tsv \
  2> research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/results/score_hybrid_cami3.stderr.txt

# 6. Main result tables:
# - research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/hybrid_mouse_baseline_plus_mean_metrics.tsv
# - research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback/hybrid_cami3_summary.tsv
