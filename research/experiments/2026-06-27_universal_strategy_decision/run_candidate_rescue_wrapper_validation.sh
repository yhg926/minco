#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

OUT=/tmp/minco_candidate_rescue_validation_20260629
TRAIN_FEATURES=/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2
MODEL_CACHE=${MINCO_PROFILE_MODEL_CACHE:-/tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib}
TOY_TAXMAP=/tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv
GTDB_TAXMAP=/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv
RESCUE_SWITCH=emitted-ani90-xny100-br01-af70

mkdir -p "$OUT"

MODEL_ARGS=()
if [[ -s "$MODEL_CACHE" ]]; then
  MODEL_ARGS=(--model-cache "$MODEL_CACHE")
fi

for sample in 5 6 7; do
  split=/tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_split.unfiltered.tsv
  if [[ "$sample" == "5" ]]; then
    split=/tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_split.exact.unfiltered.tsv
  fi
  python3 -B scripts/minco_profile_calibrated.py \
    --unique-table /tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_unique.unfiltered.tsv \
    --split-table "$split" \
    --taxmap "$TOY_TAXMAP" \
    --train-features "$TRAIN_FEATURES" \
    "${MODEL_ARGS[@]}" \
    --train-pool train12 \
    --scope bacteria \
    --strategy universal-auto-exact \
    --candidate-rescue-switch "$RESCUE_SWITCH" \
    --report-all \
    -o "$OUT/toymouse_sample${sample}.tsv"
done

python3 -B scripts/minco_profile_calibrated.py \
  --unique-table /tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample6_unique_zip_unfiltered.tsv \
  --split-table /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_split_exact_unfiltered.tsv \
  --taxmap "$GTDB_TAXMAP" \
  --train-features "$TRAIN_FEATURES" \
  "${MODEL_ARGS[@]}" \
  --train-pool train12 \
  --scope bacteria \
  --strategy universal-auto-exact \
  --candidate-rescue-switch "$RESCUE_SWITCH" \
  --report-all \
  -o "$OUT/hmp_airskin_sample6.tsv"

python3 -B scripts/minco_profile_calibrated.py \
  --unique-table /tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample11_unique_zip_unfiltered.tsv \
  --split-table /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_split_exact_unfiltered.tsv \
  --taxmap "$GTDB_TAXMAP" \
  --train-features "$TRAIN_FEATURES" \
  "${MODEL_ARGS[@]}" \
  --train-pool train12 \
  --scope bacteria \
  --strategy universal-auto-exact \
  --candidate-rescue-switch "$RESCUE_SWITCH" \
  --report-all \
  -o "$OUT/hmp_airskin_sample11.tsv"

for sample in 0 6; do
  unique=/tmp/minco_current_code_hmp_gastrooral_20260627/sample${sample}_work/minco.best_diff_unique.unfiltered.tsv
  split=/tmp/minco_current_code_hmp_gastrooral_20260627/sample${sample}_work/minco.best_diff_split.unfiltered.tsv
  if [[ "$sample" == "0" ]]; then
    split=/tmp/minco_current_code_hmp_gastrooral_20260627/sample${sample}_work/minco.best_diff_split.exact.unfiltered.tsv
  fi
  python3 -B scripts/minco_profile_calibrated.py \
    --unique-table "$unique" \
    --split-table "$split" \
    --taxmap "$GTDB_TAXMAP" \
    --train-features "$TRAIN_FEATURES" \
    "${MODEL_ARGS[@]}" \
    --train-pool train12 \
    --scope bacteria \
    --strategy universal-auto-exact \
    --candidate-rescue-switch "$RESCUE_SWITCH" \
    --report-all \
    -o "$OUT/hmp_gastrooral_sample${sample}.tsv"
done

for sample in 0 1 2; do
  python3 -B scripts/minco_profile_calibrated.py \
    --unique-table /tmp/cami3_toy_human_gut_20260626/run/sample${sample}_autoexact/minco.best_diff_unique.unfiltered.tsv \
    --split-table /tmp/cami3_toy_human_gut_20260626/run/sample${sample}_autoexact/minco.best_diff_split.unfiltered.tsv \
    --taxmap "$GTDB_TAXMAP" \
    --train-features "$TRAIN_FEATURES" \
    "${MODEL_ARGS[@]}" \
    --train-pool train12 \
    --scope bacteria \
    --strategy universal-auto-exact \
    --candidate-rescue-switch "$RESCUE_SWITCH" \
    --report-all \
    -o "$OUT/cami3_sample${sample}.tsv"
done
