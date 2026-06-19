#!/usr/bin/env bash
set -euo pipefail

KSSD=/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/kssd3mini_stage3_native
ASM_2009=/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz
ASM_GTDB=/mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_001272115.1_Salmonella_enterica_CVM_N46840_v1.0_genomic.fna.gz
SK_2009=/tmp/kssdmini_asm_2009K1545_0618
SK_READS_2009=/tmp/kssdmini_readsqc_kssdqcsample_counts_all_0618
SK_GTDB=/tmp/kssdmini_asm_GCF001272115_0618

fastANI -q "$ASM_2009" -r "$ASM_GTDB" -o /tmp/2009K1545_vs_GCF001272115_fastani.tsv

"$KSSD" sketch -p 4 --ctxmeta both -o "$SK_GTDB" "$ASM_GTDB"

"$KSSD" ani -r "$SK_GTDB" --qraw "$SK_2009" -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/kssdmini_2009asm_vs_GCF001272115asm_0618.tsv

"$KSSD" ani -r "$SK_GTDB" --qraw "$SK_READS_2009" -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/kssdmini_2009reads_vs_GCF001272115asm_0618.tsv

"$KSSD" ani -r "$SK_READS_2009" --qraw "$SK_GTDB" -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/kssdmini_GCF001272115asm_vs_2009reads_0618.tsv
