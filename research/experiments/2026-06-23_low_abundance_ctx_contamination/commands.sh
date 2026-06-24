#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

EXP=research/experiments/2026-06-23_low_abundance_ctx_contamination
READ_TAR=/mnt/new3T/minco_cami2_toymouse_20260621/2017.12.29_11.37.26_sample_0_reads.tar
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
MINCO=/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/bin/minco

mkdir -p "$EXP/traces_perread" "$EXP/logs_perread"

tar -xOf "$READ_TAR" 2017.12.29_11.37.26_sample_0/reads/reads_mapping.tsv.gz > /tmp/sample0_reads_mapping.tsv.gz

run_trace() {
  local label="$1"
  local ref_match="$2"
  /usr/bin/time -v -o "$EXP/logs_perread/${label}.time.log" \
    env MINCO_READWISE_TRACE_REF="$ref_match" \
        MINCO_READWISE_TRACE_OUT="$PWD/$EXP/traces_perread/${label}.trace.tsv" \
    "$MINCO" ani \
      -r "$REF" \
      --qraw "$READS" \
      --query-density ref --abundance-est depth --readwise-profile-only \
      --density-block-ctx 0 \
      --readwise-assign best-diff-split --readwise-ani naive \
      --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
      -p 8 \
      -o "$EXP/traces_perread/${label}.profile.tsv" \
    > "$EXP/logs_perread/${label}.stdout.log" \
    2> "$EXP/logs_perread/${label}.stderr.log"
}

run_trace kinnaridis_low_bad GCA_946487795.1
run_trace johnsonii_low_bad GCF_000159355.1
run_trace crispatus_low_bad GCF_018987235.1
run_trace rhamnosus_low_multisource GCF_900636965.1
run_trace paracasei_high_control GCF_000829035.1

python3 "$EXP/analyze_ctx_trace_sources.py"
python3 "$EXP/check_origin_marker_membership.py"
python3 "$EXP/check_crispatus_helveticus_wholegenome.py"  # ANIm table only; S10000000 sketch membership is not comparable to S2000.
python3 "$EXP/check_crispatus_helveticus_s2000_context_scan.py"
