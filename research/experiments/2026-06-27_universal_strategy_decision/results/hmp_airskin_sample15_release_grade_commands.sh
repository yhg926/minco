#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

SAMPLE=15
WORK_ROOT="${MINCO_HMP_WORK_ROOT:-/tmp}"
WORK="$WORK_ROOT/minco_current_code_hmp_airskin15_20260627"
HMP_WORK="$WORK_ROOT/cami2_hmp_unseen_transfer_20260626"
TRUTH=/tmp/cami2_hmp_airskin_20260625/truth
FASTQ="$HMP_WORK/reads/airskinurogenital_sample${SAMPLE}.nonzero.fastq.gz"
BAM_LIST="$HMP_WORK/nonzero_bam_paths_sample${SAMPLE}.txt"
URL="https://frl.publisso.de/data/frl:6425518/airskinurogenital/sample_${SAMPLE}.tar.gz"
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
TAXMAP=/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv
MODEL_CACHE=/tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib
SYLPH=/home/ubuntu/yihuiguang/bin/sylph
SYLPH_CHUNKS=/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list
SYLPH_RUN="$HMP_WORK/run_sylph_r232/sylph_sample${SAMPLE}"
SYLPH_INPUTS="$HMP_WORK/run_sylph_r232/sylph_sample${SAMPLE}.profile_inputs.list"

mkdir -p "$HMP_WORK/reads" "$HMP_WORK/logs" "$SYLPH_RUN" "$WORK/work"

awk -F'\t' -v prefix="2017.12.04_18.56.22_sample_${SAMPLE}" '$2 > 0 {print prefix "/bam/" $1 ".bam"}' \
  "$TRUTH/abundance${SAMPLE}.tsv" > "$BAM_LIST"

/usr/bin/time -v -o "$HMP_WORK/logs/sample${SAMPLE}_stream_bam_to_fastq_20260627.time.log" \
  bash -c 'set -euo pipefail; curl -fL "$1" | tar -xzf - -T "$2" --to-command="samtools fastq -" | gzip -1 > "$3"' \
  _ "$URL" "$BAM_LIST" "$FASTQ"

gzip -t "$FASTQ"

/usr/bin/time -v -o "$WORK/minco_sample${SAMPLE}_current_default.time.log" \
  env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
    -r "$REF" \
    --reads "$FASTQ" \
    --taxmap "$TAXMAP" \
    --model-cache "$MODEL_CACHE" \
    --scope bacteria \
    --workdir "$WORK/work" \
    -p16 \
    -o "$WORK/minco_sample${SAMPLE}_current_default.tsv"

/usr/bin/time -v -o "$SYLPH_RUN/sketch.time.log" \
  "$SYLPH" sketch -t 16 -r "$FASTQ" -d "$SYLPH_RUN"

sed '$a '"$SYLPH_RUN/airskinurogenital_sample${SAMPLE}.nonzero.fastq.gz.sylsp" "$SYLPH_CHUNKS" > "$SYLPH_INPUTS"

/usr/bin/time -v -o "$SYLPH_RUN/profile.time.log" \
  "$SYLPH" profile -t 16 -l "$SYLPH_INPUTS" -o "$SYLPH_RUN/profile.tsv"

python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples "${SAMPLE}" \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample "${SAMPLE}=$WORK/minco_sample${SAMPLE}_current_default.tsv" \
  --minco-unique-sample "${SAMPLE}=$WORK/work/minco.best_diff_unique.unfiltered.tsv" \
  --minco-split-sample "${SAMPLE}=$WORK/work/minco.best_diff_split.unfiltered.tsv" \
  --sylph-sample "${SAMPLE}=$SYLPH_RUN/profile.tsv" \
  --minco-method "minco_current_default_gtdb_source_abundance_sample${SAMPLE}" \
  --output-prefix "hmp_airskin${SAMPLE}_r232_source_abundance"

python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples "${SAMPLE},0,1,3,5,6,11,13,16,17,18,21,22,24,25,28" \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample "${SAMPLE}=$WORK/minco_sample${SAMPLE}_current_default.tsv" \
  --minco-unique-sample "${SAMPLE}=$WORK/work/minco.best_diff_unique.unfiltered.tsv" \
  --minco-split-sample "${SAMPLE}=$WORK/work/minco.best_diff_split.unfiltered.tsv" \
  --sylph-sample "${SAMPLE}=$SYLPH_RUN/profile.tsv" \
  --minco-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/minco_sample18_current_default.tsv \
  --minco-unique-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 18=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample18/profile.tsv \
  --minco-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv \
  --minco-unique-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 21=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.tsv \
  --minco-sample 13=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/minco_sample13_current_default.tsv \
  --minco-unique-sample 13=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 13=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 13=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample13/profile.tsv \
  --minco-sample 17=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/minco_sample17_current_default.tsv \
  --minco-unique-sample 17=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 17=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 17=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample17/profile.tsv \
  --minco-sample 16=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/minco_sample16_current_default.tsv \
  --minco-unique-sample 16=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 16=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 16=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample16/profile.tsv \
  --minco-sample 24=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/minco_sample24_current_default.tsv \
  --minco-unique-sample 24=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 24=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 24=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample24/profile.tsv \
  --minco-sample 25=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/minco_sample25_current_default.tsv \
  --minco-unique-sample 25=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 25=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 25=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample25/profile.tsv \
  --minco-sample 3=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/minco_sample3_current_default.tsv \
  --minco-unique-sample 3=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 3=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 3=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample3/profile.tsv \
  --minco-sample 1=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/minco_sample1_current_default.tsv \
  --minco-unique-sample 1=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 1=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 1=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample1/profile.tsv \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method "minco_current_default_gtdb_source_abundance_sample${SAMPLE}_0_1_3_5_6_11_13_16_17_18_21_22_24_25_28" \
  --output-prefix "hmp_current_refresh_r232_source_abundance_sample${SAMPLE}_0_1_3_5_6_11_13_16_17_18_21_22_24_25_28"
