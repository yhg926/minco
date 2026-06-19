#!/usr/bin/env bash
set -euo pipefail

KSSD=/home/ubuntu/yihuiguang/tools/minco/bin/minco_stage3_native
ASM=/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz
READ_DIR=/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield

"$KSSD" sketch -p 4 --ctxmeta both \
  -o /tmp/minco_asm_2009K1545_0618 \
  "$ASM"

"$KSSD" sketch -p 4 --conflict --ctxmeta both \
  -o /tmp/minco_ont_yield_2009K1545_noqc_0618 \
  "$READ_DIR"/5m/2009K-1545__SRR19787631__5m.fastq.gz \
  "$READ_DIR"/10m/2009K-1545__SRR19787631__10m.fastq.gz \
  "$READ_DIR"/20m/2009K-1545__SRR19787631__20m.fastq.gz \
  "$READ_DIR"/50m/2009K-1545__SRR19787631__50m.fastq.gz \
  "$READ_DIR"/75m/2009K-1545__SRR19787631__75m.fastq.gz \
  "$READ_DIR"/100m/2009K-1545__SRR19787631__100m.fastq.gz \
  "$READ_DIR"/250m/2009K-1545__SRR19787631__250m.fastq.gz

"$KSSD" sketch -p 4 --conflict --readsQC --ctxmeta both \
  -o /tmp/minco_ont_yield_2009K1545_readsqc_0618 \
  "$READ_DIR"/5m/2009K-1545__SRR19787631__5m.fastq.gz \
  "$READ_DIR"/10m/2009K-1545__SRR19787631__10m.fastq.gz \
  "$READ_DIR"/20m/2009K-1545__SRR19787631__20m.fastq.gz \
  "$READ_DIR"/50m/2009K-1545__SRR19787631__50m.fastq.gz \
  "$READ_DIR"/75m/2009K-1545__SRR19787631__75m.fastq.gz \
  "$READ_DIR"/100m/2009K-1545__SRR19787631__100m.fastq.gz \
  "$READ_DIR"/250m/2009K-1545__SRR19787631__250m.fastq.gz

"$KSSD" ani -r /tmp/minco_asm_2009K1545_0618 \
  --qraw /tmp/minco_ont_yield_2009K1545_noqc_0618 \
  -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/minco_asm_vs_ont_yield_noqc_0618.tsv

"$KSSD" ani -r /tmp/minco_asm_2009K1545_0618 \
  --qraw /tmp/minco_ont_yield_2009K1545_readsqc_0618 \
  -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/minco_asm_vs_ont_yield_readsqc_0618.tsv

"$KSSD" sketch -p 4 --ctxmeta both \
  -o /tmp/minco_ont_yield_2009K1545_noconflict_noqc_0618 \
  "$READ_DIR"/5m/2009K-1545__SRR19787631__5m.fastq.gz \
  "$READ_DIR"/10m/2009K-1545__SRR19787631__10m.fastq.gz \
  "$READ_DIR"/20m/2009K-1545__SRR19787631__20m.fastq.gz \
  "$READ_DIR"/50m/2009K-1545__SRR19787631__50m.fastq.gz \
  "$READ_DIR"/75m/2009K-1545__SRR19787631__75m.fastq.gz \
  "$READ_DIR"/100m/2009K-1545__SRR19787631__100m.fastq.gz \
  "$READ_DIR"/250m/2009K-1545__SRR19787631__250m.fastq.gz

"$KSSD" sketch -p 4 --readsQC --ctxmeta both \
  -o /tmp/minco_ont_yield_2009K1545_noconflict_readsqc_0618 \
  "$READ_DIR"/5m/2009K-1545__SRR19787631__5m.fastq.gz \
  "$READ_DIR"/10m/2009K-1545__SRR19787631__10m.fastq.gz \
  "$READ_DIR"/20m/2009K-1545__SRR19787631__20m.fastq.gz \
  "$READ_DIR"/50m/2009K-1545__SRR19787631__50m.fastq.gz \
  "$READ_DIR"/75m/2009K-1545__SRR19787631__75m.fastq.gz \
  "$READ_DIR"/100m/2009K-1545__SRR19787631__100m.fastq.gz \
  "$READ_DIR"/250m/2009K-1545__SRR19787631__250m.fastq.gz

"$KSSD" ani -r /tmp/minco_asm_2009K1545_0618 \
  --qraw /tmp/minco_ont_yield_2009K1545_noconflict_noqc_0618 \
  -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/minco_asm_vs_ont_yield_noconflict_noqc_0618.tsv

"$KSSD" ani -r /tmp/minco_asm_2009K1545_0618 \
  --qraw /tmp/minco_ont_yield_2009K1545_noconflict_readsqc_0618 \
  -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/minco_asm_vs_ont_yield_noconflict_readsqc_0618.tsv

KSSD3_READSQC_DEBUG=1 "$KSSD" sketch -p 4 --conflict --readsQC --ctxmeta both \
  -o /tmp/minco_readsqc_kssdqcsample_counts_all_0618 \
  "$READ_DIR"/5m/2009K-1545__SRR19787631__5m.fastq.gz \
  "$READ_DIR"/10m/2009K-1545__SRR19787631__10m.fastq.gz \
  "$READ_DIR"/20m/2009K-1545__SRR19787631__20m.fastq.gz \
  "$READ_DIR"/50m/2009K-1545__SRR19787631__50m.fastq.gz \
  "$READ_DIR"/75m/2009K-1545__SRR19787631__75m.fastq.gz \
  "$READ_DIR"/100m/2009K-1545__SRR19787631__100m.fastq.gz \
  "$READ_DIR"/250m/2009K-1545__SRR19787631__250m.fastq.gz

"$KSSD" ani -r /tmp/minco_asm_2009K1545_0618 \
  --qraw /tmp/minco_readsqc_kssdqcsample_counts_all_0618 \
  -n0 -f0 -t0 -s4 -m0 \
  -o /tmp/minco_asm_vs_ont_readsqc_kssdqcsample_counts_0618.tsv
