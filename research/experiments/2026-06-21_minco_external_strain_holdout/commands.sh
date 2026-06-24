#!/usr/bin/env bash
set -euo pipefail

# Data source checked from CAMI datasets page:
# https://frl.publisso.de/data/frl:6425521/strain/short_read/

BASE=/mnt/new3T/minco_cami2_strain_20260621
OUT=$BASE/minco_outputs
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv

mkdir -p "$OUT"
for i in 0 1 2; do
  Q=$BASE/sample_${i}/short_read/2018.09.07_11.43.52_sample_${i}/reads/anonymous_reads.fq.gz
  for assign in best-diff-unique best-diff-split; do
    label=unique
    if [ "$assign" = "best-diff-split" ]; then label=split; fi
    /usr/bin/time -v ./minco_core/bin/minco ani -p16 \
      -r "$REF" \
      --qraw "$Q" \
      --query-density ref \
      --abundance-est depth \
      --readwise-profile-only \
      --readwise-assign "$assign" \
      --readwise-ani zip-aaf \
      -m0 -f0 -n0 -t0 \
      -o "$OUT/strain_sample${i}_s1000_gtdb_${label}_zip_unfiltered.tsv" \
      > "$OUT/strain_sample${i}_${label}.time.log" 2>&1
  done
done

for i in 0 1 2; do
  Q=$BASE/sample_${i}/short_read/2018.09.07_11.43.52_sample_${i}/reads/anonymous_reads.fq.gz
  D=$BASE/sylph_sample${i}
  mkdir -p "$D"
  /usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 -r "$Q" -d "$D" > "$D/sketch.log" 2>&1
  SP=$(find "$D" -name "*.sylsp" -type f | sort | head -1)
  /usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph profile -t 16 "$SYLPH_DB" "$SP" -o "$D/profile.tsv" > "$D/profile.log" 2>&1
done

python3 research/experiments/2026-06-21_minco_external_strain_holdout/scripts/train_current_test_external_calls.py \
  --train-manifest research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2_plant0_2.tsv \
  --test-manifest research/experiments/2026-06-21_minco_external_strain_holdout/manifest_strain0_2.tsv \
  --taxmap "$TAXMAP" \
  --outdir "$BASE/external_test_train9_strain0_2"

cp "$BASE/external_test_train9_strain0_2/mean_summary.tsv" research/experiments/2026-06-21_minco_external_strain_holdout/mean_summary.tsv
cp "$BASE/external_test_train9_strain0_2/summary.tsv" research/experiments/2026-06-21_minco_external_strain_holdout/summary.tsv
cp "$BASE/external_test_train9_strain0_2/rf_feature_importance.tsv" research/experiments/2026-06-21_minco_external_strain_holdout/rf_feature_importance.tsv
