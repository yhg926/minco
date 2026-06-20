#!/usr/bin/env bash
set -euo pipefail

REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
FASTQ=/tmp/cami_marine_sample0_reads.fq.gz
FULL_TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
OUT=/tmp/minco_cami_lineage_s10000_20260620
SAMPLE=marmgCAMI2_short_read_sample_0

mkdir -p "$OUT"

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap "$FULL_TAXMAP" \
  --cami-profile "$OUT/s10000_gtdb_default.full_lineage.profile" \
  --cami-sample-id "$SAMPLE" -m0 \
  -o "$OUT/s10000_gtdb_default.full_lineage.tsv" \
  > "$OUT/s10000_gtdb_default.full_lineage.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 \
  -o "$OUT/s10000_gtdb_unfiltered.tsv" \
  > "$OUT/s10000_gtdb_unfiltered.log" 2>&1

# Best sample-0 thresholds from the grid:
#   species-only:          -f0.05 -n0.94 -t10
#   species+genus combined: -f0.05 -n0.92 -t10
