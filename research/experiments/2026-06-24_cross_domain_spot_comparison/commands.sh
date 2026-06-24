#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

EXP=research/experiments/2026-06-24_cross_domain_spot_comparison
RUN=/tmp/minco_context_em_external_20260624
FAIR=/tmp/minco_fair_l1_edge_adaptive_20260624
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
REF_S1000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
MINCO=./minco_core/bin/minco

mkdir -p "$RUN" "$FAIR/marine_current_best" "$FAIR/marine_sylph" "$EXP/results" "$EXP/provenance"

# Cached comparison table assembled from already-existing benchmark TSV files.
python3 - <<'PY' > research/experiments/2026-06-24_cross_domain_spot_comparison/results/summary_stdout.tsv
import pandas as pd

rows = pd.read_csv(
    "research/experiments/2026-06-24_cross_domain_spot_comparison/summary.tsv",
    sep="\t",
)
print(rows.to_csv(sep="\t", index=False), end="")
PY

# Marine sample0 read cache used by the external edge-EM run.
curl --retry 8 --retry-all-errors --retry-delay 5 -fL \
  -o "$RUN/marmgCAMI2_sample_0_reads.tar.gz" \
  "https://frl.publisso.de/data/frl%3A6425521/marine/short_read/marmgCAMI2_sample_0_reads.tar.gz"
tar -xOzf "$RUN/marmgCAMI2_sample_0_reads.tar.gz" \
  simulation_short_read/2018.08.15_09.49.32_sample_0/reads/anonymous_reads.fq.gz \
  > "$RUN/cami_marine_sample0_reads.fq.gz"

# Shared edge-dump settings. The cap was not hit by any of the three samples.
export MINCO_READWISE_EDGE_MAX=20000000
export MINCO_READWISE_EDGE_MAX_CANDIDATES=16
export MINCO_READWISE_EDGE_AMBIGUOUS_ONLY=0
export MINCO_READWISE_EDGE_SELECTED_ONLY=0

export MINCO_READWISE_EDGE_OUT="$RUN/cami3_s0_edges.tsv"
/usr/bin/time -v -o "$RUN/cami3_s0_edge.time.log" \
  "$MINCO" ani -p4 -r "$REF" \
  --qraw /mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz \
  --query-density ref --abundance-est depth \
  --readwise-profile-only --readwise-dual-evidence \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  --density-block-ctx 0 -m0 -f0 -n0 -t0 \
  -o "$RUN/cami3_s0_profile.tsv"

export MINCO_READWISE_EDGE_OUT="$RUN/marine_s0_edges.tsv"
/usr/bin/time -v -o "$RUN/marine_s0_edge.time.log" \
  "$MINCO" ani -p4 -r "$REF" \
  --qraw "$RUN/cami_marine_sample0_reads.fq.gz" \
  --query-density ref --abundance-est depth \
  --readwise-profile-only --readwise-dual-evidence \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  --density-block-ctx 0 -m0 -f0 -n0 -t0 \
  -o "$RUN/marine_s0_profile.tsv"

export MINCO_READWISE_EDGE_OUT="$RUN/strain_s0_edges.tsv"
/usr/bin/time -v -o "$RUN/strain_s0_edge.time.log" \
  "$MINCO" ani -p4 -r "$REF" \
  --qraw /mnt/new3T/minco_cami2_strain_20260621/sample_0/short_read/2018.09.07_11.43.52_sample_0/reads/anonymous_reads.fq.gz \
  --query-density ref --abundance-est depth \
  --readwise-profile-only --readwise-dual-evidence \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  --density-block-ctx 0 -m0 -f0 -n0 -t0 \
  -o "$RUN/strain_s0_profile.tsv"

/usr/bin/time -v -o "$EXP/results/score_external_edge_em.time.log" \
  python3 "$EXP/score_external_edge_em.py" \
  > "$EXP/results/score_external_edge_em.stdout.tsv" \
  2> "$EXP/results/score_external_edge_em.stderr.txt"

# Regenerate missing marine current-best and Sylph abundance profiles for the
# fair full-profile L1 comparison.
/usr/bin/time -v -o "$FAIR/marine_current_best/minco_s1000_unique_zipaaf.time.log" \
  "$MINCO" ani -p16 -r "$REF_S1000" \
  --qraw "$RUN/cami_marine_sample0_reads.fq.gz" \
  --query-density ref --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o "$FAIR/marine_current_best/minco_s1000_unique_zipaaf_sample0.tsv"

/usr/bin/time -v -o "$FAIR/marine_sylph/sketch.time.log" \
  /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 \
  -r "$RUN/cami_marine_sample0_reads.fq.gz" \
  -d "$FAIR/marine_sylph"

/usr/bin/time -v -o "$FAIR/marine_sylph/profile.time.log" \
  /home/ubuntu/yihuiguang/bin/sylph profile -t 16 "$SYLPH_DB" \
  "$FAIR/marine_sylph/cami_marine_sample0_reads.fq.gz.sylsp" \
  -o "$FAIR/marine_sylph/profile.tsv"

/usr/bin/time -v -o "$EXP/results/score_fair_full_l1_and_adaptive_em.time.log" \
  python3 "$EXP/score_fair_full_l1_and_adaptive_em.py" \
  > "$EXP/results/score_fair_full_l1_and_adaptive_em.stdout.tsv" \
  2> "$EXP/results/score_fair_full_l1_and_adaptive_em.stderr.txt"

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git status --short > "$EXP/provenance/code_status.txt" || true

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff -- "$EXP" > "$EXP/provenance/code_diff.patch" || true

GIT_DIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git \
GIT_WORK_TREE=/home/ubuntu/yihuiguang/tools/KSSD3mini \
  git diff --stat > "$EXP/provenance/code_diffstat.txt" || true
