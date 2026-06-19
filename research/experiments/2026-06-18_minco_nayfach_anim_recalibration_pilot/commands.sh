#!/usr/bin/env bash
set -euo pipefail

# Working directory for all commands:
# /home/ubuntu/yihuiguang/tools/KSSD3mini

python3 -m py_compile \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py

# 30-pair smoke test: 5 pairs per ANI band.
/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin5 \
  --sketch-size 10000 \
  --per-bin 5 \
  --jobs 4 \
  --minco-threads 1

# 3,000-pair pilot: 500 pairs per ANI band.
/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin500 \
  --sketch-size 10000 \
  --per-bin 500 \
  --jobs 8 \
  --minco-threads 1

# Initial 12,000-pair pilot before MoE replacement: 2,000 pairs per ANI band.
/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000 \
  --sketch-size 10000 \
  --per-bin 2000 \
  --jobs 16 \
  --minco-threads 1

# Fixed-denominator MoE coefficient refit from the initial feature table.
python3 -m py_compile \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/train_fixedp_moe.py

python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/train_fixedp_moe.py \
  --features research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/minco_features.tsv \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000

# Optimized MoE denominator bases and linear coefficients.
python3 -m py_compile \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/optimize_moe_model.py

/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/optimize_moe_model.py \
  --features research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/minco_features.tsv \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/moe_opt_long \
  --maxiter 35 \
  --popsize 8 \
  --polish-iter 160

# Preserve old model_ani.h, then patch NUM_CODENS==11 constants from:
# research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/moe_opt_long/optimized_moe_model_ani_arrays.hfrag
cp minco_core/src/model_ani.h \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_ani.before_moe_recalibration.h

make test

# Regenerate the same 12,000-pair feature table after MoE replacement and
# retrain/export the HGB layer on top of the new raw MoE score.
mkdir -p research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated
cp research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/sample_pairs.tsv \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv
cp research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/sample_bin_counts.tsv \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_bin_counts.tsv

/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py \
  --outdir research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated \
  --steps features,train \
  --sketch-size 10000 \
  --jobs 16 \
  --minco-threads 1

# Preserve HGB header before replacing it with the model retrained on new MoE.
cp minco_core/src/model_refaf_hgb.h \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_refaf_hgb.before_moe_recalibrated_hgb.h

cp research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_refaf_hgb.minco_generated.h \
  minco_core/src/model_refaf_hgb.h

make test

python3 -m py_compile \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py

# Check compiled final default Best output on all 3,600 held-out pairs.
/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py \
  --sample-pairs research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv \
  --predictions research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv \
  --out research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s1_best_alltest.tsv \
  --per-bin 600 \
  --jobs 16 \
  --sketch-size 10000 \
  --slmetrics 1

# Check compiled final raw CtxMoE output on all 3,600 held-out pairs.
/usr/bin/time -v python3 \
  research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py \
  --sample-pairs research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv \
  --predictions research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv \
  --out research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s3_raw_alltest.tsv \
  --per-bin 600 \
  --jobs 16 \
  --sketch-size 10000 \
  --slmetrics 3

./bin/minco ani --help | grep -A4 'CtxMoE is'
