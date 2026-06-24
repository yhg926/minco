#!/usr/bin/env bash
set -euo pipefail

BASE=/mnt/new3T/minco_cami2_toymouse_20260621
OUT=$BASE/minco_outputs
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
Q=$BASE/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
EXP=research/experiments/2026-06-21_minco_toy_mouse_gut_sample0

mkdir -p "$OUT"

python3 "$EXP/scripts/build_species_profile_from_camisim.py" \
  --distribution "$BASE/setup/distributions/distribution_0.txt" \
  --metadata "$BASE/setup/internal/meta_data.tsv" \
  --nodes /mnt/new3T/gtdbr220/gtdbr226/taxdump/nodes.dmp \
  --names /mnt/new3T/gtdbr220/gtdbr226/taxdump/names.dmp \
  --sample-id mouse0 \
  --out-profile "$BASE/taxonomic_profile_mouse0_species.txt" \
  --out-summary "$BASE/taxonomic_profile_mouse0_species.summary.tsv"

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
    -o "$OUT/toymouse_sample0_s1000_gtdb_${label}_zip_unfiltered.tsv" \
    > "$OUT/toymouse_sample0_${label}.time.log" 2>&1
done

mkdir -p "$BASE/sylph_sample0"
/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 -r "$Q" -d "$BASE/sylph_sample0" \
  > "$BASE/sylph_sample0/sketch.log" 2>&1
/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph profile -t 16 \
  "$SYLPH_DB" "$BASE/sylph_sample0/anonymous_reads.fq.gz.sylsp" \
  -o "$BASE/sylph_sample0/profile.tsv" \
  > "$BASE/sylph_sample0/profile.log" 2>&1

python3 research/experiments/2026-06-21_minco_external_strain_holdout/scripts/train_current_test_external_calls.py \
  --train-manifest research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2_plant0_2.tsv \
  --test-manifest "$EXP/manifest_mouse0.tsv" \
  --taxmap "$TAXMAP" \
  --outdir "$EXP/train9_test_mouse0"

python3 research/experiments/2026-06-21_minco_external_strain_holdout/scripts/diagnose_simple_filters.py \
  --score-dir "$EXP/train9_test_mouse0" \
  --taxmap "$TAXMAP" \
  --outdir "$EXP/train9_test_mouse0"

WORK=$BASE/source_aware_crosswalk
mkdir -p "$WORK/source_positive"
python3 "$EXP/scripts/prepare_source_aware_inputs.py" \
  --setup-dir "$BASE/setup" \
  --unique-minco "$OUT/toymouse_sample0_s1000_gtdb_unique_zip_unfiltered.tsv" \
  --split-minco "$OUT/toymouse_sample0_s1000_gtdb_split_zip_unfiltered.tsv" \
  --sylph-profile "$BASE/sylph_sample0/profile.tsv" \
  --gtdb-genome-dir /mnt/new3T/gtdbr220/GTDBr226_genomes \
  --outdir "$WORK"
tar -xzf "$BASE/CAMISIM_setup.tar.gz" -C "$WORK/source_positive" -T "$WORK/positive_source_tar_paths.list"
find "$WORK/source_positive/source_genomes" -type f | sort > "$WORK/positive_source_paths.list"
skani dist --ql "$WORK/selected_called_ref_paths.list" \
  --rl "$WORK/positive_source_paths.list" \
  -t 16 --min-af 5 \
  -o "$WORK/skani_called_vs_sources.tsv"
python3 "$EXP/scripts/build_source_aware_taxmap.py" \
  --source-genomes "$WORK/positive_source_genomes.tsv" \
  --skani "$WORK/skani_called_vs_sources.tsv" \
  --base-taxmap "$TAXMAP" \
  --nodes /mnt/new3T/gtdbr220/gtdbr226/taxdump/nodes.dmp \
  --names /mnt/new3T/gtdbr220/gtdbr226/taxdump/names.dmp \
  --min-ani 95 \
  --min-af 50 \
  --af-mode min \
  --out-overrides "$WORK/mouse0_sourceaware_overrides_ani95_minaf50.tsv" \
  --out-taxmap "$WORK/mouse0_sourceaware_taxmap_ani95_minaf50.tsv" \
  --out-summary "$WORK/mouse0_sourceaware_taxmap_ani95_minaf50.summary.tsv"
python3 research/experiments/2026-06-21_minco_external_strain_holdout/scripts/train_current_test_external_calls.py \
  --train-manifest research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2_plant0_2.tsv \
  --test-manifest "$EXP/manifest_mouse0.tsv" \
  --taxmap "$WORK/mouse0_sourceaware_taxmap_ani95_minaf50.tsv" \
  --outdir "$EXP/train9_test_mouse0_sourceaware_ani95_minaf50"
python3 research/experiments/2026-06-21_minco_external_strain_holdout/scripts/diagnose_simple_filters.py \
  --score-dir "$EXP/train9_test_mouse0_sourceaware_ani95_minaf50" \
  --taxmap "$WORK/mouse0_sourceaware_taxmap_ani95_minaf50.tsv" \
  --outdir "$EXP/train9_test_mouse0_sourceaware_ani95_minaf50"

python3 -m py_compile \
  research/experiments/2026-06-20_minco_readwise_correction_research/scripts/analyze_readwise_corrections.py \
  "$EXP/scripts/build_species_profile_from_camisim.py" \
  "$EXP/scripts/prepare_source_aware_inputs.py" \
  "$EXP/scripts/build_source_aware_taxmap.py"

make -C minco_core test

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani naive \
  --density-block-ctx 0 \
  --readwise-ctx-filter fake-prob \
  --readwise-fake-threshold 3.0 \
  -m0 -f0 -n0 -t0 \
  -o research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_fakeprob_t30.tsv

/usr/bin/time -v ./minco_core/bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani naive \
  --density-block-ctx 0 \
  --readwise-ctx-filter fake-prob \
  --readwise-fake-threshold 3.5 \
  -m0 -f0 -n0 -t0 \
  -o research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_fakeprob_t35.tsv
