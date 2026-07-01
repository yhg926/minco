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

scripts/minco_profile --check-ref \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --reads /tmp/minco_default_raw_smoke_20260701/reads.fq \
  > research/experiments/2026-07-01_profile_default_rescue_productization_validation/results/packaged_ref_preflight.tsv
