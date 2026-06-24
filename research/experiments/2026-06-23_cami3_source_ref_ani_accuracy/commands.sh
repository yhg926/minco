#!/usr/bin/env bash
set -euo pipefail

python3 -m py_compile \
  research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/score_source_ref_ani.py \
  research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/score_alternative_read_ani.py

python3 research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/score_source_ref_ani.py
python3 research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/score_alternative_read_ani.py
