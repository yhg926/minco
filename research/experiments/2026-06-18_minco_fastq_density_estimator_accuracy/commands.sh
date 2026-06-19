#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/minco

head -40000  /mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq > /tmp/minco_fastq_R1_10k.fq
head -200000 /mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq > /tmp/minco_fastq_R1_50k.fq
head -400000 /mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq > /tmp/minco_fastq_R1_100k.fq
head -800000 /mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq > /tmp/minco_fastq_R1_200k.fq

for n in 10k 50k 100k 200k; do
  /usr/bin/time -v bin/minco_stage3_native sketch \
    -p 4 \
    --conflict \
    --ctxmeta both \
    -o /tmp/minco_fastq_density_${n}_std_0618 \
    /tmp/minco_fastq_R1_${n}.fq

  /usr/bin/time -v bin/minco_exact20m sketch \
    -p 2 \
    --conflict \
    --ctxmeta both \
    -o /tmp/minco_fastq_density_${n}_exact_0618 \
    /tmp/minco_fastq_R1_${n}.fq
done

/usr/bin/time -v bin/minco_stage3_native sketch \
  -p 8 \
  --conflict \
  --ctxmeta both \
  -o /tmp/minco_fastq_full_R1_conflict_0618 \
  /mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq
