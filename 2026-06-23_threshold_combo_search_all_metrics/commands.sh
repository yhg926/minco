#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

python3 -m py_compile 2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py

/usr/bin/time -v \
  python3 2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py \
    --max-combos 1000 \
    --ani-estimator aaf \
  > /tmp/minco_threshold_combo_1000.stdout \
  2> /tmp/minco_threshold_combo_1000.stderr

/usr/bin/time -v \
  python3 2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py \
    --preset focused \
    --max-combos 5000 \
    --ani-estimator aaf \
  > /tmp/minco_threshold_combo_focused5000.stdout \
  2> /tmp/minco_threshold_combo_focused5000.stderr

/usr/bin/time -v \
  python3 2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py \
    --ref-mode ctxobj \
    --preset focused \
    --max-combos 5000 \
    --ani-estimator aaf \
  > /tmp/minco_threshold_combo_ctxobj_focused5000.stdout \
  2> /tmp/minco_threshold_combo_ctxobj_focused5000.stderr

# The following two exploratory ensemble grids were run as inline Python
# snippets in the working directory. They import search_threshold_combos.py,
# precompute ctx-marker/full and ctx+obj values, then score scaled add-back.
# Outputs:
#   2026-06-23_threshold_combo_search_all_metrics/results/ensemble_scaled_grid.tsv
#   2026-06-23_threshold_combo_search_all_metrics/results/ensemble_scaled_highcami_grid.tsv
