#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/minco

OUT=/tmp/minco_validation_20260618
mkdir -p "$OUT/genomes1000" "$OUT/fastq250m"

# Functional test suite.
/usr/bin/time -v make test > "$OUT/make_test.stdout" 2> "$OUT/make_test.time.log"

# 1000-genome GTDBr226 stress test.
head -1000 /mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list \
  > "$OUT/genomes1000/genomes1000.list"
awk 'BEGIN{bad=0} { if (system("test -r \"" $0 "\"") != 0) { print $0; bad++ } } END{ print "missing=" bad; exit bad ? 1 : 0 }' \
  "$OUT/genomes1000/genomes1000.list"
/usr/bin/time -v bin/minco sketch -p 8 --ctxmeta both \
  -l "$OUT/genomes1000/genomes1000.list" \
  -o "$OUT/genomes1000/sketch" \
  > "$OUT/genomes1000/sketch.stdout" \
  2> "$OUT/genomes1000/sketch.time.log"
/usr/bin/time -v bin/minco sketch -i "$OUT/genomes1000/sketch" \
  > "$OUT/genomes1000/index.stdout" \
  2> "$OUT/genomes1000/index.time.log"
/usr/bin/time -v bin/minco ani -q "$OUT/genomes1000/sketch" \
  -m2 -s -1 -d -p 8 \
  -o "$OUT/genomes1000/ani_triangle.tsv" \
  > "$OUT/genomes1000/ani.stdout" \
  2> "$OUT/genomes1000/ani.time.log"
/usr/bin/time -v bin/minco matrix --format triangle \
  -q "$OUT/genomes1000/sketch" -d -p 8 \
  -o "$OUT/genomes1000/matrix_triangle.tsv" \
  > "$OUT/genomes1000/matrix.stdout" \
  2> "$OUT/genomes1000/matrix.time.log"

# Four 250m ONT FASTQ stress test.
find /home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/250m \
  -maxdepth 1 -type f -name '*.fastq.gz' | sort | head -4 \
  > "$OUT/fastq250m/reads250m_4.list"
/usr/bin/time -v bin/minco sketch --conflict --readsQC --ctxmeta both -p 8 \
  -l "$OUT/fastq250m/reads250m_4.list" \
  -o "$OUT/fastq250m/reads_sketch" \
  > "$OUT/fastq250m/reads_sketch.stdout" \
  2> "$OUT/fastq250m/reads_sketch.time.log"
printf '%s\n' \
  /home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz \
  /home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1553.fa.gz \
  /home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1559.fa.gz \
  /home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1562.fa.gz \
  > "$OUT/fastq250m/assemblies4.list"
/usr/bin/time -v bin/minco sketch --ctxmeta both -p 4 \
  -l "$OUT/fastq250m/assemblies4.list" \
  -o "$OUT/fastq250m/asm_sketch" \
  > "$OUT/fastq250m/asm_sketch.stdout" \
  2> "$OUT/fastq250m/asm_sketch.time.log"
/usr/bin/time -v bin/minco sketch -i "$OUT/fastq250m/asm_sketch" \
  > "$OUT/fastq250m/asm_index.stdout" \
  2> "$OUT/fastq250m/asm_index.time.log"
/usr/bin/time -v bin/minco ani \
  -r "$OUT/fastq250m/asm_sketch" \
  --qraw "$OUT/fastq250m/reads_sketch" \
  -m 0 -f 0 -n 0 -t 0 -p 8 \
  -o "$OUT/fastq250m/qraw_ani.tsv" \
  > "$OUT/fastq250m/qraw_ani.stdout" \
  2> "$OUT/fastq250m/qraw_ani.time.log"
