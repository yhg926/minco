#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

head -50 /tmp/kssd3a_based_kssdmini_stages/genomes_50.list > /tmp/kssdmini_density_test50.list

make -C kssd3a_stage0 \
  PRONAME=kssd3mini_exact20m \
  BINDIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/bin \
  OBJDIR=/home/ubuntu/yihuiguang/tools/KSSD3mini/kssd3a_stage0/obj_exact20m \
  CFLAGS="-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -march=native -mavx2 -mbmi2 -DKSSD3MINI_HASH_BOTTOMK=1 -DKSSD3MINI_KEEP_KSSD_FILTER=0 -DKSSD3MINI_STREAM_BOTTOMK=1 -DKSSD3MINI_SKETCH_SIZE=20000000"

/usr/bin/time -v bin/kssd3mini_stage3_native sketch \
  -p 8 \
  --ctxmeta both \
  -l /tmp/kssdmini_density_test50.list \
  -o /tmp/kssdmini_density50_std_0618

/usr/bin/time -v bin/kssd3mini_exact20m sketch \
  -p 2 \
  --conflict \
  --ctxmeta both \
  -l /tmp/kssdmini_density_test50.list \
  -o /tmp/kssdmini_density50_exact_0618

research/experiments/2026-06-18_kssdmini_density_estimator_predicts_unique_contexts/analyze_density.py \
  --std-ctxmeta /tmp/kssdmini_density50_std_0618/kssdmini.ctxmeta.tsv \
  --exact-ctxmeta /tmp/kssdmini_density50_exact_0618/kssdmini.ctxmeta.tsv \
  --infilemeta /tmp/kssdmini_density50_std_0618/lcofiles.infilemeta \
  --summary research/experiments/2026-06-18_kssdmini_density_estimator_predicts_unique_contexts/summary.tsv \
  --per-sample research/experiments/2026-06-18_kssdmini_density_estimator_predicts_unique_contexts/per_sample.tsv
