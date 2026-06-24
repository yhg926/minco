#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

EXP=research/experiments/2026-06-24_integrated_dual_evidence_benchmark
RUN=/tmp/minco_dual_s2000_20260624
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
MINCO=minco_core/bin/minco
mkdir -p "$RUN" "$EXP/results" "$EXP/provenance"

python3 -m py_compile "$EXP/score_integrated_dual_evidence.py"

run_dual() {
  local label="$1"
  local reads="$2"
  local out="$RUN/${label}_dual.tsv"
  if [ -s "$out" ]; then
    return 0
  fi
  /usr/bin/time -v -o "$RUN/${label}_dual.time.log" \
    "$MINCO" ani -p4 \
      -r "$REF" \
      --qraw "$reads" \
      --query-density ref --abundance-est depth --readwise-profile-only --readwise-dual-evidence \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -m0 -f0 -n0 -t0 \
      -o "$out" \
      > "$RUN/${label}_dual.stdout.log" \
      2> "$RUN/${label}_dual.stderr.log"
}

run_dual mouse_s0 /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
run_dual mouse_s1 /mnt/new3T/minco_cami2_toymouse_20260621/sample_1/2017.12.29_11.37.26_sample_1/reads/anonymous_reads.fq.gz
run_dual mouse_s2 /mnt/new3T/minco_cami2_toymouse_20260621/sample_2/2017.12.29_11.37.26_sample_2/reads/anonymous_reads.fq.gz

run_dual cami3_s0 /mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz
run_dual cami3_s1 /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz
run_dual cami3_s2 /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz

/usr/bin/time -v -o "$EXP/results/score_integrated_dual_evidence.time.log" \
  python3 "$EXP/score_integrated_dual_evidence.py" \
  > "$EXP/results/score_integrated_dual_evidence.stdout.tsv" \
  2> "$EXP/results/score_integrated_dual_evidence.stderr.txt"

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git status --short > "$EXP/provenance/code_status.txt" || true

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff --stat > "$EXP/provenance/code_diffstat.txt" || true
