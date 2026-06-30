#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini
test -s /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample6/profile.tsv
test -s /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample11/profile.tsv
python3 -B /home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py --sylph-source r232 --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232
