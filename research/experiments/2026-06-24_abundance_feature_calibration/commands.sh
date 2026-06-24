#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

# Code repository: /home/ubuntu/yihuiguang/tools/KSSD3mini
# Git metadata directory: /home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git
# Commit: b55b4d3f4c986de29098d2f1a092ee28251008aa
# Ref/describe: main, b55b4d3-dirty
# Working tree: dirty

python3 -m py_compile research/experiments/2026-06-24_abundance_feature_calibration/search_feature_calibration.py

/usr/bin/time -v python3 research/experiments/2026-06-24_abundance_feature_calibration/search_feature_calibration.py --max-models 1 --save-feature-table > /tmp/minco_feature_calibration_features.stdout 2> /tmp/minco_feature_calibration_features.stderr

/usr/bin/time -v python3 research/experiments/2026-06-24_abundance_feature_calibration/search_feature_calibration.py --max-models 120 --feature-table research/experiments/2026-06-24_abundance_feature_calibration/results/selected_call_features.tsv --truth-table research/experiments/2026-06-24_abundance_feature_calibration/results/truth_profiles.tsv > /tmp/minco_feature_calibration_120.stdout 2> /tmp/minco_feature_calibration_120.stderr
