#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/minco

head -50 /tmp/kssd3a_based_minco_stages/genomes_50.list > /tmp/minco_density_test50.list

make -C kssd3a_stage0 \
  PRONAME=minco_exact20m \
  BINDIR=/home/ubuntu/yihuiguang/tools/minco/bin \
  OBJDIR=/home/ubuntu/yihuiguang/tools/minco/kssd3a_stage0/obj_exact20m \
  CFLAGS="-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -march=native -mavx2 -mbmi2 -DMINCO_HASH_BOTTOMK=1 -DMINCO_KEEP_KSSD_FILTER=0 -DMINCO_STREAM_BOTTOMK=1 -DMINCO_SKETCH_SIZE=20000000"

/usr/bin/time -v bin/minco_stage3_native sketch \
  -p 8 \
  --ctxmeta both \
  -l /tmp/minco_density_test50.list \
  -o /tmp/minco_density50_std_0618

/usr/bin/time -v bin/minco_exact20m sketch \
  -p 2 \
  --conflict \
  --ctxmeta both \
  -l /tmp/minco_density_test50.list \
  -o /tmp/minco_density50_exact_0618

research/experiments/2026-06-18_minco_density_estimator_predicts_unique_contexts/analyze_density.py \
  --std-ctxmeta /tmp/minco_density50_std_0618/minco.ctxmeta.tsv \
  --exact-ctxmeta /tmp/minco_density50_exact_0618/minco.ctxmeta.tsv \
  --infilemeta /tmp/minco_density50_std_0618/lcofiles.infilemeta \
  --summary research/experiments/2026-06-18_minco_density_estimator_predicts_unique_contexts/summary.tsv \
  --per-sample research/experiments/2026-06-18_minco_density_estimator_predicts_unique_contexts/per_sample.tsv
