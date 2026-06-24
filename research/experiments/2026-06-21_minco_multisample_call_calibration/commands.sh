#!/usr/bin/env bash
set -euo pipefail

MINCO=./minco_core/bin/minco
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
OUT=/tmp/minco_multisample_calibration_20260621
TOY_EXTRA=/mnt/new3T/minco_cami3_toygut_extra_20260621
SYLPH=/home/ubuntu/yihuiguang/bin/sylph
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb

mkdir -p "$OUT"

# CAMI III Toy Human Gut sample1/sample2 downloads.
for i in 1 2; do
  curl -L --fail -C - \
    -o "$TOY_EXTRA/sample_${i}_reads.tar.gz" \
    "https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_${i}_reads.tar.gz"
  tar -xzf "$TOY_EXTRA/sample_${i}_reads.tar.gz" -C "$TOY_EXTRA"
done

# Extract all Toy Human Gut profiles from the existing tarball.
tar -xzf /mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles.tar.gz \
  -C /mnt/new3T/minco_cami3_toygut_20260620

# Generate missing unfiltered marine candidate tables.
for sample in 1 2; do
  /usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
    --qraw "/tmp/cami_marine_sample${sample}_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-unique --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$OUT/marine_sample${sample}_s1000_gtdb_unique_zip_unfiltered.tsv"

  /usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
    --qraw "/tmp/cami_marine_sample${sample}_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$OUT/marine_sample${sample}_s1000_gtdb_split_zip_unfiltered.tsv"
done

# Generate Toy Human Gut sample1/sample2 minco candidate tables.
for sample in 1 2; do
  /usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
    --qraw "$TOY_EXTRA/sample_${sample}_reads/anonymous_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-unique --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$OUT/toy_sample${sample}_s1000_gtdb_unique_zip_unfiltered.tsv"

  /usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
    --qraw "$TOY_EXTRA/sample_${sample}_reads/anonymous_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 \
    -o "$OUT/toy_sample${sample}_s1000_gtdb_split_zip_unfiltered.tsv"
done

# Generate Toy Human Gut sample1/sample2 Sylph baselines.
for sample in 1 2; do
  mkdir -p "$OUT/sylph_toy_sample${sample}"
  /usr/bin/time -v "$SYLPH" sketch -t 16 \
    -r "$TOY_EXTRA/sample_${sample}_reads/anonymous_reads.fq.gz" \
    -d "$OUT/sylph_toy_sample${sample}"
  /usr/bin/time -v "$SYLPH" profile -t 16 "$SYLPH_DB" \
    "$OUT/sylph_toy_sample${sample}/anonymous_reads.fq.gz.sylsp" \
    -o "$OUT/sylph_toy_sample${sample}/profile.tsv"
done

# Multi-sample model and threshold calibration.
python3 research/experiments/2026-06-21_minco_multisample_call_calibration/scripts/calibrate_multisample_calls.py \
  --manifest research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2.tsv \
  --taxmap "$TAXMAP" \
  --outdir "$OUT/model_marine_toy0_2"

python3 research/experiments/2026-06-21_minco_multisample_call_calibration/scripts/grid_loso_thresholds.py \
  --manifest research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2.tsv \
  --taxmap "$TAXMAP" \
  --outdir "$OUT/threshold_marine_toy0_2"
