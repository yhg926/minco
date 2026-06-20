#!/usr/bin/env bash
set -euo pipefail

OUT=/tmp/minco_cami_newani_20260620
mkdir -p "$OUT"

FASTQ=/tmp/cami_marine_sample0_reads.fq.gz
TAX=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.tsv
REF_S1000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
REF_S10000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
REF_S1000_PLUSVIRUS=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S1000" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 \
  -o "$OUT/s1000_gtdb_default.tsv" > "$OUT/s1000_gtdb_default.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S1000" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 \
  -o "$OUT/s1000_gtdb_unfiltered.tsv" > "$OUT/s1000_gtdb_unfiltered.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S10000" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 \
  -o "$OUT/s10000_gtdb_default.tsv" > "$OUT/s10000_gtdb_default.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S10000" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 \
  -o "$OUT/s10000_gtdb_unfiltered.tsv" > "$OUT/s10000_gtdb_unfiltered.log" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S1000_PLUSVIRUS" \
  --qraw "$FASTQ" --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap "$TAX" \
  --cami-profile "$OUT/s1000_plusvirus_default.profile" \
  --cami-sample-id marmgCAMI2_short_read_sample_0 -m0 \
  -o "$OUT/s1000_plusvirus_default.tsv" > "$OUT/s1000_plusvirus_default.log" 2>&1

objdump -d bin/minco_stage3_native | rg -n "popcnt|vpopcnt"
