#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

EXP=research/experiments/2026-06-24_context_level_ambiguity_em_prototype
RUN=/tmp/minco_context_em_20260624
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
READS1=/mnt/new3T/minco_cami2_toymouse_20260621/sample_1/2017.12.29_11.37.26_sample_1/reads/anonymous_reads.fq.gz
READS2=/mnt/new3T/minco_cami2_toymouse_20260621/sample_2/2017.12.29_11.37.26_sample_2/reads/anonymous_reads.fq.gz
MINCO=./minco_core/bin/minco
mkdir -p "$RUN" "$EXP/results" "$EXP/provenance"

make -C minco_core

zcat "$READS" | head -4000 | gzip -c > /tmp/minco_edge_smoke_1k.fq.gz

MINCO_READWISE_EDGE_OUT=/tmp/minco_edge_smoke_edges.tsv \
MINCO_READWISE_EDGE_MAX=100000 \
MINCO_READWISE_EDGE_MAX_CANDIDATES=64 \
MINCO_READWISE_EDGE_AMBIGUOUS_ONLY=1 \
MINCO_READWISE_EDGE_SELECTED_ONLY=0 \
/usr/bin/time -v -o /tmp/minco_edge_smoke.time.log \
  "$MINCO" ani -p2 \
    -r "$REF" \
    --qraw /tmp/minco_edge_smoke_1k.fq.gz \
    --query-density ref --abundance-est depth --readwise-profile-only --readwise-dual-evidence \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 \
    -o /tmp/minco_edge_smoke_profile.tsv \
    > /tmp/minco_edge_smoke.stdout.log \
    2> /tmp/minco_edge_smoke.stderr.log

MINCO_READWISE_EDGE_OUT="$RUN/mouse_s0_edges.tsv" \
MINCO_READWISE_EDGE_MAX=5000000 \
MINCO_READWISE_EDGE_MAX_CANDIDATES=64 \
MINCO_READWISE_EDGE_AMBIGUOUS_ONLY=1 \
MINCO_READWISE_EDGE_SELECTED_ONLY=0 \
/usr/bin/time -v -o "$RUN/mouse_s0_edge.time.log" \
  "$MINCO" ani -p4 \
    -r "$REF" \
    --qraw "$READS" \
    --query-density ref --abundance-est depth --readwise-profile-only --readwise-dual-evidence \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/mouse_s0_profile.tsv" \
    > "$RUN/mouse_s0_edge.stdout.log" \
    2> "$RUN/mouse_s0_edge.stderr.log"

python3 -m py_compile "$EXP/score_mouse0_edge_em.py"

/usr/bin/time -v -o "$EXP/results/score_mouse0_edge_em.time.log" \
  python3 "$EXP/score_mouse0_edge_em.py" \
  > "$EXP/results/score_mouse0_edge_em.stdout.tsv" \
  2> "$EXP/results/score_mouse0_edge_em.stderr.txt"

python3 -m py_compile "$EXP/search_mouse0_filtered_edge_em.py"
python3 "$EXP/search_mouse0_filtered_edge_em.py" \
  > "$EXP/results/search_mouse0_filtered_edge_em.stdout.tsv" \
  2> "$EXP/results/search_mouse0_filtered_edge_em.stderr.txt"

python3 -m py_compile "$EXP/refine_mouse0_filtered_edge_em_beta.py"
python3 "$EXP/refine_mouse0_filtered_edge_em_beta.py" \
  > "$EXP/results/refine_mouse0_filtered_edge_em_beta.stdout.tsv" \
  2> "$EXP/results/refine_mouse0_filtered_edge_em_beta.stderr.txt"

MINCO_READWISE_EDGE_OUT="$RUN/mouse_s1_edges.tsv" \
MINCO_READWISE_EDGE_MAX=20000000 \
MINCO_READWISE_EDGE_MAX_CANDIDATES=64 \
MINCO_READWISE_EDGE_AMBIGUOUS_ONLY=1 \
MINCO_READWISE_EDGE_SELECTED_ONLY=0 \
/usr/bin/time -v -o "$RUN/mouse_s1_edge.time.log" \
  "$MINCO" ani -p4 \
    -r "$REF" \
    --qraw "$READS1" \
    --query-density ref --abundance-est depth --readwise-profile-only --readwise-dual-evidence \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/mouse_s1_profile.tsv" \
    > "$RUN/mouse_s1_edge.stdout.log" \
    2> "$RUN/mouse_s1_edge.stderr.log"

MINCO_READWISE_EDGE_OUT="$RUN/mouse_s2_edges.tsv" \
MINCO_READWISE_EDGE_MAX=20000000 \
MINCO_READWISE_EDGE_MAX_CANDIDATES=64 \
MINCO_READWISE_EDGE_AMBIGUOUS_ONLY=1 \
MINCO_READWISE_EDGE_SELECTED_ONLY=0 \
/usr/bin/time -v -o "$RUN/mouse_s2_edge.time.log" \
  "$MINCO" ani -p4 \
    -r "$REF" \
    --qraw "$READS2" \
    --query-density ref --abundance-est depth --readwise-profile-only --readwise-dual-evidence \
    --readwise-assign best-diff-split --readwise-ani naive \
    --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 \
    -o "$RUN/mouse_s2_profile.tsv" \
    > "$RUN/mouse_s2_edge.stdout.log" \
    2> "$RUN/mouse_s2_edge.stderr.log"

python3 -m py_compile "$EXP/score_mouse012_filtered_edge_em.py"
/usr/bin/time -v -o "$EXP/results/score_mouse012_filtered_edge_em.time.log" \
  python3 "$EXP/score_mouse012_filtered_edge_em.py" \
  > "$EXP/results/score_mouse012_filtered_edge_em.stdout.tsv" \
  2> "$EXP/results/score_mouse012_filtered_edge_em.stderr.txt"

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git status --short > "$EXP/provenance/code_status.txt" || true

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff -- minco_core/src/command_ani.c "$EXP" > "$EXP/provenance/code_diff.patch" || true

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff --stat > "$EXP/provenance/code_diffstat.txt" || true
