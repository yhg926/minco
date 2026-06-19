#!/usr/bin/env bash
set -euo pipefail

# Working directory:
# /home/ubuntu/yihuiguang/tools/minco

command -v mash
mash --version

python3 -m py_compile \
  research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/scripts/compare_mash_minco_anim.py

/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/scripts/compare_mash_minco_anim.py \
  --sample-pairs research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv \
  --predictions research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv \
  --outdir research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100 \
  --sketch-size 10000 \
  --per-bin 100 \
  --jobs 16

column -t \
  research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/mash_minco_anim_metrics.tsv
