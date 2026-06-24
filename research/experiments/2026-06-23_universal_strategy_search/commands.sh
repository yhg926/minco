#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

python3 research/experiments/2026-06-23_universal_strategy_search/score_universal_rules.py \
  > research/experiments/2026-06-23_universal_strategy_search/results/score_universal_rules.stdout.tsv

MINCO_MODEL_NAMES=logreg,hgb \
python3 research/experiments/2026-06-23_universal_strategy_search/score_universal_models.py \
  > research/experiments/2026-06-23_universal_strategy_search/results/score_universal_models_light.stdout.tsv

python3 research/experiments/2026-06-23_universal_strategy_search/score_previous_rf_hgb_threshold.py \
  > research/experiments/2026-06-23_universal_strategy_search/results/score_previous_rf_hgb_threshold.stdout.tsv

curl -fsSL https://cami-challenge.org/static/examples/gs_marine_short.profile \
  -o /tmp/gs_marine_short.profile

python3 research/experiments/2026-06-23_universal_strategy_search/score_universal_abundance.py \
  > research/experiments/2026-06-23_universal_strategy_search/results/score_universal_abundance.stdout.tsv
