#!/usr/bin/env bash
set -euo pipefail

python3 -m py_compile \
  research/experiments/2026-06-23_cami3_toygut_gtdb_ani_accuracy/score_gtdb_ani_accuracy.py

python3 \
  research/experiments/2026-06-23_cami3_toygut_gtdb_ani_accuracy/score_gtdb_ani_accuracy.py
