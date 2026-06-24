#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

python3 -m py_compile \
  research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/assemble_full_benchmark.py

/usr/bin/time -v \
  python3 research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/assemble_full_benchmark.py

column -t -s $'\t' \
  research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/results/benchmark_decision.tsv

column -t -s $'\t' \
  research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/results/f1_abundance_comparison.tsv

column -t -s $'\t' \
  research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/results/ani_accuracy_comparison.tsv

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git status --short \
  > research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/provenance/code_status.txt

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff --stat \
  > research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/provenance/code_diffstat.txt
