#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

make -C minco_core -j4
make -C minco_core \
  BINDIR=bin_coden15 \
  OBJDIR=obj_coden15 \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1' \
  -j4

/usr/bin/time -v \
  -o research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_ctxmarker_rawani_report.time.log \
  minco_core/bin/minco ani -p8 \
    -r /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker \
    --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_ctxmarker_rawani_report.tsv \
  > research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_ctxmarker_rawani_report.stdout.log \
  2> research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_ctxmarker_rawani_report.stderr.log

/usr/bin/time -v \
  -o research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_T15_ctxmarker_rawani_report.time.log \
  minco_core/bin_coden15/minco ani -p16 \
    -r /mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker \
    --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    -m0 -f0 -n0 -t0 \
    -o research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_T15_ctxmarker_rawani_report.tsv \
  > research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_T15_ctxmarker_rawani_report.stdout.log \
  2> research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_T15_ctxmarker_rawani_report.stderr.log

python3 -m py_compile \
  research/experiments/2026-06-23_toymouse_source_rep_ani/build_source_rep_ani_table.py

python3 \
  research/experiments/2026-06-23_toymouse_source_rep_ani/build_source_rep_ani_table.py

cat research/experiments/2026-06-23_toymouse_source_rep_ani/summary.tsv

MINCO_PROFILE=research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_T15_ctxmarker_rawani_report.tsv \
MINCO_PROFILE_LABEL='coden15 ctx-marker raw-ANI-report' \
MINCO_CODEN_LEN=15 \
MINCO_STORAGE=ctxobj96 \
MINCO_REFERENCE_KIND='coden15 s2000 ctx-marker' \
OUT_TABLE=toymouse_sample0_source_rep_ani_coden15.tsv \
OUT_SUMMARY=summary_coden15.tsv \
python3 research/experiments/2026-06-23_toymouse_source_rep_ani/build_source_rep_ani_table.py

python3 - <<'PY'
import pandas as pd
from pathlib import Path

base = Path("research/experiments/2026-06-23_toymouse_source_rep_ani")
rows = []
for method, fname in [("coden11", "summary.tsv"), ("coden15", "summary_coden15.tsv")]:
    df = pd.read_csv(base / fname, sep="\t")
    df.insert(0, "method", method)
    rows.append(df)
pd.concat(rows, ignore_index=True).to_csv(
    base / "summary_coden11_vs_coden15.tsv", sep="\t", index=False
)
PY

python3 -m py_compile \
  research/experiments/2026-06-23_toymouse_source_rep_ani/run_source_ref_assembly_minco_naive.py

python3 \
  research/experiments/2026-06-23_toymouse_source_rep_ani/run_source_ref_assembly_minco_naive.py
