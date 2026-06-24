#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

python3 -m py_compile research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/benchmark_virtual_markers.py

/usr/bin/time -v python3 research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/benchmark_virtual_markers.py \
  --chunk-records 10000000 \
  --outdir research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results \
  > /tmp/virtual_s2000_marker_scan_fixed.stdout \
  2> /tmp/virtual_s2000_marker_scan_fixed.stderr

du -sh /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup \
  /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker

du -b /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup \
  /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
