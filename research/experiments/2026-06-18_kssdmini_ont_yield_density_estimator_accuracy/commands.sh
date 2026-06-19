#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

BASE=/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield
SAMPLE=2009K-1545__SRR19787631

# Normal 10k KSSDmini density sketches.
for y in 5m 10m 20m 50m 75m 100m 250m; do
  f=$BASE/$y/${SAMPLE}__${y}.fastq.gz
  out=/tmp/kssdmini_ont_density_${y}_std_0618
  log=/tmp/kssdmini_ont_density_${y}_std_0618.time
  /usr/bin/time -v bin/kssd3mini_stage3_native sketch \
    -p 4 \
    --conflict \
    --ctxmeta both \
    -o "$out" \
    "$f" 2> "$log"
done

# Exact validation. 5m-20m fit under the existing 20M cap.
for y in 5m 10m 20m; do
  f=$BASE/$y/${SAMPLE}__${y}.fastq.gz
  out=/tmp/kssdmini_ont_density_${y}_exact20m_0618
  log=/tmp/kssdmini_ont_density_${y}_exact20m_0618.time
  /usr/bin/time -v bin/kssd3mini_exact20m sketch \
    -p 2 \
    --conflict \
    --ctxmeta both \
    -o "$out" \
    "$f" 2> "$log"
done

# 50m-250m required a larger cap; exact100m did not cap for this isolate.
make -C kssd3a_stage0 \
  PRONAME=kssd3mini_exact100m \
  BINDIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/bin \
  OBJDIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/kssd3a_stage0/obj_exact100m \
  CFLAGS="-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -march=native -mavx2 -mbmi2 -DKSSD3MINI_HASH_BOTTOMK=1 -DKSSD3MINI_KEEP_KSSD_FILTER=0 -DKSSD3MINI_STREAM_BOTTOMK=1 -DKSSD3MINI_SKETCH_SIZE=100000000"

for y in 50m 75m 100m 250m; do
  f=$BASE/$y/${SAMPLE}__${y}.fastq.gz
  out=/tmp/kssdmini_ont_density_${y}_exact100m_0618
  log=/tmp/kssdmini_ont_density_${y}_exact100m_0618.time
  /usr/bin/time -v bin/kssd3mini_exact100m sketch \
    -p 2 \
    --conflict \
    --ctxmeta both \
    -o "$out" \
    "$f" 2> "$log"
done
