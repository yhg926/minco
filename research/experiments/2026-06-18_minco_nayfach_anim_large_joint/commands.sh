#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ubuntu/yihuiguang/tools/KSSD3mini
EXP=$ROOT/research/experiments/2026-06-18_minco_nayfach_anim_large_joint
RUN=$EXP/large_10k_perbin20000

cd "$ROOT"

python3 -m py_compile \
  research/experiments/2026-06-18_minco_nayfach_anim_large_joint/scripts/train_joint_moe_hgb.py \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/optimize_moe_model.py \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py

/usr/bin/time -v -o "$RUN/time_sample_features_train_initial.txt" python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  --outdir "$RUN" \
  --sketch-size 10000 \
  --per-bin 20000 \
  --jobs 16 \
  --minco-threads 1 \
  --seed 20260618

/usr/bin/time -v -o "$RUN/time_optimize_moe_long.txt" python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/optimize_moe_model.py \
  --features "$RUN/minco_features.tsv" \
  --outdir "$RUN/moe_opt_long" \
  --maxiter 35 \
  --popsize 8 \
  --polish-iter 160 \
  --seed 20260618

/usr/bin/time -v -o "$RUN/time_train_joint_hgb.txt" python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_large_joint/scripts/train_joint_moe_hgb.py \
  --features "$RUN/minco_features.tsv" \
  --moe-params "$RUN/moe_opt_long/optimized_moe_params.tsv" \
  --outdir "$RUN/joint_hgb" \
  --seed 20260618

/usr/bin/time -v -o "$RUN/time_validate_current_binary_s1.txt" python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py \
  --sample-pairs "$RUN/sample_pairs.tsv" \
  --predictions "$RUN/joint_hgb/joint_model_predictions.tsv" \
  --out "$RUN/current_binary_validation_s1_best_alltest.tsv" \
  --per-bin 6000 \
  --jobs 16 \
  --sketch-size 10000 \
  --slmetrics 1

cp minco_core/src/model_refaf_hgb.h "$EXP/model_refaf_hgb.before_large_hgb.h"
cp "$RUN/model_refaf_hgb.minco_generated.h" minco_core/src/model_refaf_hgb.h

# The installed header comment was edited to identify the large-run calibration.
/usr/bin/time -v -o "$EXP/time_make_test_after_large_hgb.txt" make test

/usr/bin/time -v -o "$RUN/time_validate_installed_large_hgb_s1.txt" python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py \
  --sample-pairs "$RUN/sample_pairs.tsv" \
  --predictions "$RUN/model_predictions.tsv" \
  --out "$RUN/installed_large_hgb_validation_s1_best_alltest.tsv" \
  --per-bin 6000 \
  --jobs 16 \
  --sketch-size 10000 \
  --slmetrics 1
