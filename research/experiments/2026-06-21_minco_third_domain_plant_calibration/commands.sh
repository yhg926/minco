#!/usr/bin/env bash
set -euo pipefail

# Download CAMI II plant-associated setup and read samples 0-2.
BASE=/mnt/new3T/minco_cami2_plant_20260621
mkdir -p "$BASE"
curl -L --fail -C - \
  -o "$BASE/rhimgCAMI2_setup.tar.gz" \
  "https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_setup.tar.gz"
tar -xzf "$BASE/rhimgCAMI2_setup.tar.gz" -C "$BASE"

for i in 0 1 2; do
  curl -L --fail -C - \
    -o "$BASE/rhimgCAMI2_sample_${i}_reads.tar.gz" \
    "https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_sample_${i}_reads.tar.gz"
  mkdir -p "$BASE/sample_${i}"
  tar -xzf "$BASE/rhimgCAMI2_sample_${i}_reads.tar.gz" -C "$BASE/sample_${i}"
done

# Generate minco unfiltered readwise candidate tables.
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
OUT="$BASE/minco_outputs"
mkdir -p "$OUT"

for i in 0 1 2; do
  Q="$BASE/sample_${i}/simulation_short_read/2019.09.27_13.59.10_sample_${i}/reads/anonymous_reads.fq.gz"
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
      -o "$OUT/plant_sample${i}_s1000_gtdb_${label}_zip_unfiltered.tsv" \
      > "$OUT/plant_sample${i}_${label}.time.log" 2>&1
  done
done

# Generate Sylph baseline.
SYLPH=/home/ubuntu/yihuiguang/bin/sylph
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
for i in 0 1 2; do
  Q="$BASE/sample_${i}/simulation_short_read/2019.09.27_13.59.10_sample_${i}/reads/anonymous_reads.fq.gz"
  D="$BASE/sylph_sample${i}"
  mkdir -p "$D"
  /usr/bin/time -v "$SYLPH" sketch -t 16 -r "$Q" -d "$D" > "$D/sketch.log" 2>&1
  SP=$(find "$D" -name "*.sylsp" -type f | sort | head -1)
  /usr/bin/time -v "$SYLPH" profile -t 16 "$SYLPH_DB" "$SP" -o "$D/profile.tsv" > "$D/profile.log" 2>&1
done

# Run no-leak, rule-feature, RF/HGB ensemble LOSO calibration.
MAN=research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2_plant0_2.tsv
TAX=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
CAL="$BASE/calibration_marine_toy_plant_20260621"
mkdir -p "$CAL/model_noleak_rules_ensemble" "$CAL/threshold"

/usr/bin/time -v python3 \
  research/experiments/2026-06-21_minco_multisample_call_calibration/scripts/calibrate_multisample_calls.py \
  --manifest "$MAN" \
  --taxmap "$TAX" \
  --outdir "$CAL/model_noleak_rules_ensemble" \
  > "$CAL/model_noleak_rules_ensemble.log" 2>&1

python3 \
  research/experiments/2026-06-21_minco_multisample_call_calibration/scripts/grid_loso_thresholds.py \
  --manifest "$MAN" \
  --taxmap "$TAX" \
  --outdir "$CAL/threshold"
