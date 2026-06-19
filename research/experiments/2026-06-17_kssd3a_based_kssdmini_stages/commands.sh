#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/ubuntu/yihuiguang/tools/KSSD3mini"
KSSD3A="/home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a"
WORK="/tmp/kssd3a_based_kssdmini_stages"
LIST="/tmp/kssdmini_1000_speed/genomes.list"

mkdir -p "$WORK"

if [[ ! -s "$LIST" ]]; then
  head -1000 /mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list > "$LIST"
fi

cd "$ROOT"
make
make stage0 stage1 stage2 stage3_native

# Baseline A: true KSSD3A fixed-rate sketch.
/usr/bin/time -v "$KSSD3A" sketch -f8 -p8 \
  -l "$LIST" \
  -o "$WORK/kssd3a_f8" \
  > "$WORK/kssd3a_f8.sketch.stdout" \
  2> "$WORK/kssd3a_f8.sketch.time.log"

du -sb "$WORK/kssd3a_f8" > "$WORK/kssd3a_f8.du.tsv"
find "$WORK/kssd3a_f8" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$WORK/kssd3a_f8.files.tsv"

# Baseline B: current KSSDmini best compact sketch.
/usr/bin/time -v "$ROOT/bin/kssdmini_compact_bmi2" sketch -p 8 -s 10000 \
  -l "$LIST" \
  -o "$WORK/kssdmini_compact_bmi2.kmini" \
  > "$WORK/kssdmini_compact_bmi2.sketch.stdout" \
  2> "$WORK/kssdmini_compact_bmi2.sketch.time.log"

wc -c "$WORK/kssdmini_compact_bmi2.kmini" > "$WORK/kssdmini_compact_bmi2.size.tsv"

# Stage 0: copied KSSD3A source path, unchanged fixed-rate selection.
/usr/bin/time -v "$ROOT/bin/kssd3mini_stage0" sketch -f8 -p8 \
  -l "$LIST" \
  -o "$WORK/stage0_kssd3a_f8" \
  > "$WORK/stage0_kssd3a_f8.sketch.stdout" \
  2> "$WORK/stage0_kssd3a_f8.sketch.time.log"

du -sb "$WORK/stage0_kssd3a_f8" > "$WORK/stage0_kssd3a_f8.du.tsv"
find "$WORK/stage0_kssd3a_f8" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$WORK/stage0_kssd3a_f8.files.tsv"

# Stage 1: KSSD3A -f8 prefilter retained, then KSSDmini hash-bottom-k cap.
/usr/bin/time -v "$ROOT/bin/kssd3mini_stage1" sketch -f8 -p8 \
  -l "$LIST" \
  -o "$WORK/stage1_f8_hash_bottomk" \
  > "$WORK/stage1_f8_hash_bottomk.sketch.stdout" \
  2> "$WORK/stage1_f8_hash_bottomk.sketch.time.log"

du -sb "$WORK/stage1_f8_hash_bottomk" > "$WORK/stage1_f8_hash_bottomk.du.tsv"
find "$WORK/stage1_f8_hash_bottomk" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$WORK/stage1_f8_hash_bottomk.files.tsv"

# Stage 2: exact fixed-size MinHash after collecting all contexts.
# This is intentionally expensive; on the benchmark it used ~28 GB RSS.
/usr/bin/time -v "$ROOT/bin/kssd3mini_stage2" sketch -f8 -p8 \
  -l "$LIST" \
  -o "$WORK/stage2_no_filter_hash_bottomk_1000" \
  > "$WORK/stage2_no_filter_hash_bottomk_1000.sketch.stdout" \
  2> "$WORK/stage2_no_filter_hash_bottomk_1000.sketch.time.log"

du -sb "$WORK/stage2_no_filter_hash_bottomk_1000" > "$WORK/stage2_no_filter_hash_bottomk_1000.du.tsv"
find "$WORK/stage2_no_filter_hash_bottomk_1000" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$WORK/stage2_no_filter_hash_bottomk_1000.files.tsv"

# Stage 3: selected streaming exact fixed-size MinHash implementation.
/usr/bin/time -v "$ROOT/bin/kssd3mini_stage3_native" sketch -f8 -p8 \
  -l "$LIST" \
  -o "$WORK/stage3_native_default_m2_1000" \
  > "$WORK/stage3_native_default_m2_1000.sketch.stdout" \
  2> "$WORK/stage3_native_default_m2_1000.sketch.time.log"

du -sb "$WORK/stage3_native_default_m2_1000" > "$WORK/stage3_native_default_m2_1000.du.tsv"
find "$WORK/stage3_native_default_m2_1000" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$WORK/stage3_native_default_m2_1000.files.tsv"

cmp -s "$WORK/stage3_native_default_m2_1000/comblco" "$WORK/stage2_no_filter_hash_bottomk_1000/comblco"
