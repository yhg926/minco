#!/usr/bin/env bash
set -euo pipefail

DATA_TAR=/mnt/new3T/gtdbr220/test/patmgCAMI2.tar.gz
WORK=/tmp/patmg_CAMI2_reads_20260620
APPENDED_REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620
VIRAL_REF=/mnt/new3T/IMG_VR_data/refseq_viro/refseq_virus_minco_S1000_anno_20260619
MINCO=/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/minco

mkdir -p "$WORK"

/usr/bin/time -v -o "$WORK/extract.time.txt" \
  tar -xzf "$DATA_TAR" -C "$WORK" \
  patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz \
  patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz

zcat "$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz" | wc -l
zcat "$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz" | wc -l

/usr/bin/time -v -o "$WORK/minco_appended215k_patmgCAMI2_n0.9_f0.01.time.txt" \
  bash -lc "set -o pipefail; gzip -cd '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz' '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz' | '$MINCO' ani -p16 -r '$APPENDED_REF' --qraw - --query-density ref --abundance-est depth -f0.01 -n0.9 -m0 -o '$WORK/minco_appended215k_patmgCAMI2_n0.9_f0.01.tsv'"

/usr/bin/time -v -o "$WORK/minco_viral_patmgCAMI2_n0.9_f0.01.time.txt" \
  bash -lc "set -o pipefail; gzip -cd '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz' '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz' | '$MINCO' ani -p16 -r '$VIRAL_REF' --qraw - --query-density ref --abundance-est depth -f0.01 -n0.9 -m0 -o '$WORK/minco_viral_patmgCAMI2_n0.9_f0.01.tsv'"

/usr/bin/time -v -o "$WORK/minco_viral_patmgCAMI2_permissive.time.txt" \
  bash -lc "set -o pipefail; gzip -cd '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz' '$WORK/patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz' | '$MINCO' ani -p16 -r '$VIRAL_REF' --qraw - --query-density ref --abundance-est depth -f0 -n0 -t1 -m0 -o '$WORK/minco_viral_patmgCAMI2_permissive.tsv'"

