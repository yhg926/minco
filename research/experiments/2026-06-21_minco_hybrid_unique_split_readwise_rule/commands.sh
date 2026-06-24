#!/usr/bin/env bash
set -euo pipefail

MINCO=./minco_core/bin/minco
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
OUT=/tmp/minco_hybrid_20260621

mkdir -p "$OUT"

# Existing toy unique table:
# /mnt/new3T/minco_cami3_toygut_20260620/minco_s1000_gtdbonly_sample0_unfiltered.tsv

# Existing toy split table generated before this note:
# /tmp/minco_toy_s1000_gtdbonly_split_zip_unfiltered.tsv

# Existing marine sample0 reads and gold:
# /tmp/cami_marine_sample0_reads.fq.gz
# /tmp/gs_marine_short.profile

# Generate missing marine split table.
/usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0 \
  -o "$OUT/marine_s1000_gtdb_split_zip_sample0_unfiltered.tsv"

# Generate clean current-binary marine unique table.
/usr/bin/time -v "$MINCO" ani -p16 -r "$REF" \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0 \
  -o "$OUT/marine_s1000_gtdb_unique_zip_sample0_unfiltered.tsv"

# Join unique/split features and grid-search the hybrid rule family.
python3 research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/scripts/hybrid_unique_split_score.py \
  --taxmap "$TAXMAP" \
  --outdir "$OUT" \
  --toy-gold /mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt \
  --toy-unique /mnt/new3T/minco_cami3_toygut_20260620/minco_s1000_gtdbonly_sample0_unfiltered.tsv \
  --toy-split /tmp/minco_toy_s1000_gtdbonly_split_zip_unfiltered.tsv \
  --toy-sylph /mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv \
  --marine-gold /tmp/gs_marine_short.profile \
  --marine-unique "$OUT/marine_s1000_gtdb_unique_zip_sample0_unfiltered.tsv" \
  --marine-split "$OUT/marine_s1000_gtdb_split_zip_sample0_unfiltered.tsv" \
  --marine-sylph /tmp/sylph_marine_sample0/profile.tsv

cp "$OUT/summary.tsv" research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/summary.tsv
cp "$OUT/hybrid_best_params.tsv" research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/hybrid_best_params.tsv
cp "$OUT/hybrid_grid_best50.tsv" research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/hybrid_grid_best50.tsv
cp "$OUT/toy_bacteria_gold_taxids.tsv" research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/toy_bacteria_gold_taxids.tsv
cp "$OUT/marine_gold_taxids.tsv" research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/marine_gold_taxids.tsv

# Test whether a learned unique/split call model generalizes when holding out
# an entire sample.
python3 research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/scripts/train_hybrid_classifier.py \
  --toy-features "$OUT/toy_gtdb_bacteria.joined_features.tsv" \
  --toy-gold "$OUT/toy_bacteria_gold_taxids.tsv" \
  --marine-features "$OUT/marine_gtdb_species.joined_features.tsv" \
  --marine-gold "$OUT/marine_gold_taxids.tsv" \
  --out "$OUT/classifier_leave_one_dataset_out.tsv"

cp "$OUT/classifier_leave_one_dataset_out.tsv" \
  research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/classifier_leave_one_dataset_out.tsv
