#!/usr/bin/env bash
set -euo pipefail

# Main paths used in this experiment. Python one-off scripts expanded the
# taxmap, scored the threshold grid, and created post-hoc CAMI profiles from
# the unfiltered candidate table.

REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620
FASTQ=/tmp/cami_marine_sample0_reads.fq.gz
SPECIES_TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.tsv
FULL_TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
OUT=/tmp/minco_cami_lineage_20260620
SAMPLE=marmgCAMI2_short_read_sample_0

mkdir -p "$OUT"

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap "$FULL_TAXMAP" \
  --cami-profile "$OUT/s1000_plusvirus_default.full_lineage.profile" \
  --cami-sample-id "$SAMPLE" -m0 \
  -o "$OUT/s1000_plusvirus_default.full_lineage.tsv" \
  > "$OUT/s1000_plusvirus_default.full_lineage.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 \
  -o "$OUT/s1000_plusvirus_unfiltered.tsv" \
  > "$OUT/s1000_plusvirus_unfiltered.log" 2>&1

# Best sample-0 post-hoc threshold from the grid:
#   -f0.1 -n0.94 -t10
