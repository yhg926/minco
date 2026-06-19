#!/usr/bin/env bash
set -euo pipefail

# Build and smoke-test minco.
make
make test

# minco benchmarks.
scripts/run_minco_benchmark.sh \
  /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_neisseria_n20/selected_genomes.tsv \
  /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_neisseria_n20/per_pair_errors_with_kssd3a_naive.tsv \
  /tmp/minco_neisseria_n20_final \
  10000 \
  8

scripts/run_minco_benchmark.sh \
  /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/selected_genomes.tsv \
  /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/per_pair_errors_with_naive.tsv \
  /tmp/minco_all_n20_final \
  10000 \
  8

# Compact uint64 ctxhash_obj comparison on the same five-species set.
mkdir -p /tmp/minco_all_n20_compact
awk -F'\t' 'NR > 1 { print $3 }' \
  /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/selected_genomes.tsv \
  > /tmp/minco_all_n20_compact/genomes.list

/usr/bin/time -v -o /tmp/minco_all_n20_compact/sketch.time.log \
  bin/minco_compact sketch -p 8 -s 10000 \
  -l /tmp/minco_all_n20_compact/genomes.list \
  -o /tmp/minco_all_n20_compact/minco_compact_s10000.kmini

/usr/bin/time -v -o /tmp/minco_all_n20_compact/ani.time.log \
  bin/minco_compact ani -p 8 \
  -r /tmp/minco_all_n20_compact/minco_compact_s10000.kmini \
  -q /tmp/minco_all_n20_compact/minco_compact_s10000.kmini \
  -o /tmp/minco_all_n20_compact/minco_compact_s10000_ani.tsv

python3 scripts/summarize_minco_benchmark.py \
  --minco /tmp/minco_all_n20_compact/minco_compact_s10000_ani.tsv \
  --ground /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/per_pair_errors.tsv \
  --per-pair-out /tmp/minco_all_n20_compact/per_pair_errors.tsv \
  --summary-out /tmp/minco_all_n20_compact/summary.tsv

# ntHash-style compact rolling-context hash comparison.
mkdir -p /tmp/minco_all_n20_compact_roll
cp /tmp/minco_all_n20_compact/genomes.list \
  /tmp/minco_all_n20_compact_roll/genomes.list

/usr/bin/time -v -o /tmp/minco_all_n20_compact_roll/sketch.time.log \
  bin/minco_compact_roll sketch -p 8 -s 10000 \
  -l /tmp/minco_all_n20_compact_roll/genomes.list \
  -o /tmp/minco_all_n20_compact_roll/minco_compact_roll_s10000.kmini

/usr/bin/time -v -o /tmp/minco_all_n20_compact_roll/ani.time.log \
  bin/minco_compact_roll ani -p 8 \
  -r /tmp/minco_all_n20_compact_roll/minco_compact_roll_s10000.kmini \
  -q /tmp/minco_all_n20_compact_roll/minco_compact_roll_s10000.kmini \
  -o /tmp/minco_all_n20_compact_roll/minco_compact_roll_s10000_ani.tsv

python3 scripts/summarize_minco_benchmark.py \
  --minco /tmp/minco_all_n20_compact_roll/minco_compact_roll_s10000_ani.tsv \
  --ground /mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/per_pair_errors.tsv \
  --per-pair-out /tmp/minco_all_n20_compact_roll/per_pair_errors.tsv \
  --summary-out /tmp/minco_all_n20_compact_roll/summary.tsv

# Neisseria n=20 timing baselines.
mkdir -p /tmp/minco_baselines_neisseria_n20
sed 's/\r$//' /tmp/minco_neisseria_n20_final/genomes.list \
  > /tmp/minco_baselines_neisseria_n20/genomes.clean.list

/usr/bin/time -v -o /tmp/minco_baselines_neisseria_n20/kssd3a_sketch.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a sketch \
  -f8 -p8 \
  -l /tmp/minco_neisseria_n20_final/genomes.list \
  -o /tmp/minco_baselines_neisseria_n20/kssd3a_sketch

/usr/bin/time -v -o /tmp/minco_baselines_neisseria_n20/kssd3a_ani.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a ani \
  -r /tmp/minco_baselines_neisseria_n20/kssd3a_sketch \
  -q /tmp/minco_baselines_neisseria_n20/kssd3a_sketch \
  -m0 -p8 -n0 -f0 \
  -o /tmp/minco_baselines_neisseria_n20/kssd3a_ani.tsv

/usr/bin/time -v -o /tmp/minco_baselines_neisseria_n20/skani_triangle.time.log \
  skani triangle -t 8 \
  -l /tmp/minco_baselines_neisseria_n20/genomes.clean.list \
  -o /tmp/minco_baselines_neisseria_n20/skani_triangle.tsv

/usr/bin/time -v -o /tmp/minco_baselines_neisseria_n20/mash_triangle.time.log \
mash triangle -p 8 -s 10000 -E \
  -l /tmp/minco_baselines_neisseria_n20/genomes.clean.list \
  > /tmp/minco_baselines_neisseria_n20/mash_triangle.tsv

# Five-species n=20 each timing baselines.
mkdir -p /tmp/minco_baselines_all_n20
sed 's/\r$//' /tmp/minco_all_n20_final/genomes.list \
  > /tmp/minco_baselines_all_n20/genomes.clean.list

/usr/bin/time -v -o /tmp/minco_baselines_all_n20/kssd3a_sketch.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a sketch \
  -f8 -p8 \
  -l /tmp/minco_all_n20_final/genomes.list \
  -o /tmp/minco_baselines_all_n20/kssd3a_sketch

/usr/bin/time -v -o /tmp/minco_baselines_all_n20/kssd3a_ani.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a ani \
  -r /tmp/minco_baselines_all_n20/kssd3a_sketch \
  -q /tmp/minco_baselines_all_n20/kssd3a_sketch \
  -m0 -p8 -n0 -f0 \
  -o /tmp/minco_baselines_all_n20/kssd3a_ani.tsv

/usr/bin/time -v -o /tmp/minco_baselines_all_n20/skani_triangle.time.log \
  skani triangle -t 8 \
  -l /tmp/minco_baselines_all_n20/genomes.clean.list \
  -o /tmp/minco_baselines_all_n20/skani_triangle.tsv

/usr/bin/time -v -o /tmp/minco_baselines_all_n20/mash_triangle.time.log \
  mash triangle -p 8 -s 10000 -E \
  -l /tmp/minco_baselines_all_n20/genomes.clean.list \
  > /tmp/minco_baselines_all_n20/mash_triangle.tsv

# GTDBr226 first-1000 speed benchmark.
mkdir -p /tmp/minco_1000_speed/baselines
head -1000 /mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list \
  > /tmp/minco_1000_speed/genomes.list

/usr/bin/time -v -o /tmp/minco_1000_speed/full.sketch.time.log \
  bin/minco sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/full.kmini

/usr/bin/time -v -o /tmp/minco_1000_speed/full.ani.time.log \
  bin/minco ani -p 8 \
  -r /tmp/minco_1000_speed/full.kmini \
  -q /tmp/minco_1000_speed/full.kmini \
  -o /tmp/minco_1000_speed/full.ani.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/full.index_m2.time.log \
  bin/minco ani -q /tmp/minco_1000_speed/full.kmini \
  -m2 -d -p8 \
  -o /tmp/minco_1000_speed/full.index_m2.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/compact.sketch.time.log \
  bin/minco_compact sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/compact.kmini

/usr/bin/time -v \
  bin/minco_compact_bmi2 sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/compact_bmi2_lut.kmini \
  > /tmp/minco_1000_speed/compact_bmi2_lut.sketch.stdout \
  2> /tmp/minco_1000_speed/compact_bmi2_lut.sketch.time.log

cmp -s /tmp/minco_1000_speed/compact.kmini \
  /tmp/minco_1000_speed/compact_bmi2_lut.kmini

/usr/bin/time -v \
  bin/minco_compact sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/compact_vec_radix.kmini \
  > /tmp/minco_1000_speed/compact_vec_radix.sketch.stdout \
  2> /tmp/minco_1000_speed/compact_vec_radix.sketch.time.log

cmp -s /tmp/minco_1000_speed/compact.kmini \
  /tmp/minco_1000_speed/compact_vec_radix.kmini

/usr/bin/time -v \
  bin/minco_compact_bmi2 sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/compact_bmi2_vec_radix.kmini \
  > /tmp/minco_1000_speed/compact_bmi2_vec_radix.sketch.stdout \
  2> /tmp/minco_1000_speed/compact_bmi2_vec_radix.sketch.time.log

cmp -s /tmp/minco_1000_speed/compact.kmini \
  /tmp/minco_1000_speed/compact_bmi2_vec_radix.kmini

/usr/bin/time -v -o /tmp/minco_1000_speed/compact.ani.time.log \
  bin/minco_compact ani -p 8 \
  -r /tmp/minco_1000_speed/compact.kmini \
  -q /tmp/minco_1000_speed/compact.kmini \
  -o /tmp/minco_1000_speed/compact.ani.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/compact.index_detail.time.log \
  bin/minco_compact ani --index -p8 \
  -r /tmp/minco_1000_speed/compact.kmini \
  -o /tmp/minco_1000_speed/compact.index_detail.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/compact.index_m2.time.log \
  bin/minco_compact ani -q /tmp/minco_1000_speed/compact.kmini \
  -m2 -d -p8 \
  -o /tmp/minco_1000_speed/compact.index_m2.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/compact_roll.sketch.time.log \
  bin/minco_compact_roll sketch -p 8 -s 10000 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/compact_roll.kmini

/usr/bin/time -v -o /tmp/minco_1000_speed/compact_roll.ani.time.log \
  bin/minco_compact_roll ani -p 8 \
  -r /tmp/minco_1000_speed/compact_roll.kmini \
  -q /tmp/minco_1000_speed/compact_roll.kmini \
  -o /tmp/minco_1000_speed/compact_roll.ani.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/compact_roll.index_m2.time.log \
  bin/minco_compact_roll ani -q /tmp/minco_1000_speed/compact_roll.kmini \
  -m2 -d -p8 \
  -o /tmp/minco_1000_speed/compact_roll.index_m2.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/kssd3a_sketch.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a sketch \
  -f8 -p8 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/baselines/kssd3a_sketch

/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/kssd3a_ani.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a ani \
  -r /tmp/minco_1000_speed/baselines/kssd3a_sketch \
  -q /tmp/minco_1000_speed/baselines/kssd3a_sketch \
  -m0 -p8 -n0 -f0 \
  -o /tmp/minco_1000_speed/baselines/kssd3a_ani.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/kssd3a_ani_m2.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a ani \
  -q /tmp/minco_1000_speed/baselines/kssd3a_sketch \
  -m2 -s -1 -d -p8 \
  -o /tmp/minco_1000_speed/baselines/kssd3a_ani_m2.tsv

KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 \
/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/kssd3a_ani_m2_index.time.log \
  /home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a ani \
  -q /tmp/minco_1000_speed/baselines/kssd3a_sketch \
  -m2 -s -1 -d -p8 \
  -o /tmp/minco_1000_speed/baselines/kssd3a_ani_m2_index.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/skani_triangle.time.log \
  skani triangle -t 8 \
  -l /tmp/minco_1000_speed/genomes.list \
  -o /tmp/minco_1000_speed/baselines/skani_triangle.tsv

/usr/bin/time -v -o /tmp/minco_1000_speed/baselines/mash_triangle.time.log \
  mash triangle -p 8 -s 10000 -E \
  -l /tmp/minco_1000_speed/genomes.list \
  > /tmp/minco_1000_speed/baselines/mash_triangle.tsv

# Cross-species guard check from all_n20 all-vs-all output.
python3 - <<'PY'
import csv
selected='/mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20/selected_genomes.tsv'
ani='/tmp/minco_all_n20_final/minco_s10000_ani.tsv'
sp={}
with open(selected, newline='') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        sp[r['genome_id'].strip()] = r['species'].strip()
counts={'same_pass':0,'same_total':0,'cross_pass':0,'cross_total':0}
with open(ani, newline='') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        q=r['Qry'].strip(); b=r['Ref'].strip()
        if q>=b or q not in sp or b not in sp: continue
        key = 'same' if sp[q] == sp[b] else 'cross'
        counts[key + '_total'] += 1
        if r['Pass_default'] == '1':
            counts[key + '_pass'] += 1
print(counts)
PY
