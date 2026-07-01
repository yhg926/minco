#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

PYTHONDONTWRITEBYTECODE=1 python3 -B \
  research/experiments/2026-07-01_profile_default_rescue_productization_validation/run_default_entrypoint_validation.py

PYTHONDONTWRITEBYTECODE=1 python3 -B \
  research/experiments/2026-07-01_profile_default_rescue_productization_validation/run_raw_default_smoke.py

make minco

bash tests/smoke.sh

bash tests/full_cli.sh
