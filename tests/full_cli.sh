#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/minco_full_cli.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

BIN="$ROOT/bin/minco"

make_fasta() {
  local out="$1"
  local name="$2"
  local seed="$3"
  local mutate="$4"

  awk -v name="$name" -v seed="$seed" -v mutate="$mutate" '
    BEGIN {
      bases = "ACGT";
      x = seed;
      print ">" name;
      for (i = 1; i <= 36000; i++) {
        x = (1103515245 * x + 12345) % 2147483648;
        b = (x % 4) + 1;
        if (mutate && i % 997 == 0) b = (b % 4) + 1;
        printf "%s", substr(bases, b, 1);
        if (i % 80 == 0) printf "\n";
      }
      if ((i - 1) % 80 != 0) printf "\n";
    }
  ' > "$out"
}

make_fastq() {
  local out="$1"
  local name="$2"

  awk -v name="$name" '
    BEGIN {
      bases = "ACGT";
      for (r = 1; r <= 1200; r++) {
        printf "@%s_%04d\n", name, r;
        for (i = 1; i <= 180; i++) {
          b = ((r * 37 + i * 17 + int(i / 11)) % 4) + 1;
          printf "%s", substr(bases, b, 1);
        }
        printf "\n+\n";
        for (i = 1; i <= 180; i++) printf "I";
        printf "\n";
      }
    }
  ' > "$out"
}

assert_nonempty() {
  test -s "$1"
}

make_fasta "$WORK/a.fna" a 17 0
make_fasta "$WORK/b.fna" b 17 1
cp "$WORK/a.fna" "$WORK/c.fna"
make_fastq "$WORK/reads.fq" reads
printf '%s\n%s\n%s\n' "$WORK/a.fna" "$WORK/b.fna" "$WORK/c.fna" > "$WORK/genomes.list"

"$BIN" examples > "$WORK/examples.txt"
"$BIN" doctor > "$WORK/doctor.txt"

"$BIN" sketch -p 2 --anno --position -l "$WORK/genomes.list" -o "$WORK/base" > "$WORK/base.sketch.log" 2>&1
assert_nonempty "$WORK/base/comblco"
assert_nonempty "$WORK/base/lcofiles.stat"
assert_nonempty "$WORK/base/lcofiles.anno"
assert_nonempty "$WORK/base/lcofiles.infilemeta"
assert_nonempty "$WORK/base/comblco.position"

"$BIN" sketch -i "$WORK/base" > "$WORK/base.index.log" 2>&1
assert_nonempty "$WORK/base/comblco.index"
assert_nonempty "$WORK/base/sortedcomb_ctxgid64obj32"

"$BIN" sketch --psmp "$WORK/base" > "$WORK/base.psmp.tsv"
"$BIN" sketch --psketch "$WORK/base" > "$WORK/base.psketch.tsv"
"$BIN" sketch --pindex "$WORK/base" > "$WORK/base.pindex.tsv"
"$BIN" sketch --ppos "$WORK/base" > "$WORK/base.ppos.tsv"
test "$(wc -l < "$WORK/base.psmp.tsv")" -eq 3
assert_nonempty "$WORK/base.psketch.tsv"
assert_nonempty "$WORK/base.pindex.tsv"
assert_nonempty "$WORK/base.ppos.tsv"

awk 'NR == 1 { print $2 }' "$WORK/base.psmp.tsv" > "$WORK/keep.txt"
awk 'NR == 2 { print $2 }' "$WORK/base.psmp.tsv" > "$WORK/remove.txt"
"$BIN" sketch --keep "$WORK/keep.txt" -o "$WORK/kept" "$WORK/base" > "$WORK/keep.log"
"$BIN" sketch --remove "$WORK/remove.txt" -o "$WORK/removed" "$WORK/base" > "$WORK/remove.log"
"$BIN" sketch --keep "$WORK/keep.txt" --drop-position -o "$WORK/kept_drop" "$WORK/base" > "$WORK/keep_drop.log"
assert_nonempty "$WORK/kept/comblco"
assert_nonempty "$WORK/removed/comblco"
assert_nonempty "$WORK/kept_drop/comblco"
test ! -e "$WORK/kept_drop/comblco.position"

"$BIN" sketch --append -o "$WORK/appended" "$WORK/kept" "$WORK/removed" > "$WORK/append.log"
assert_nonempty "$WORK/appended/comblco"

"$BIN" sketch --dedup 0.001 --metric ctx-naive -o "$WORK/dedup" "$WORK/base" > "$WORK/dedup.log"
"$BIN" sketch --dedup 0.001 --metric ctx-naive --dedup-index -o "$WORK/dedup_index" "$WORK/base" > "$WORK/dedup_index.log"
assert_nonempty "$WORK/dedup/comblco"
assert_nonempty "$WORK/dedup_index/comblco"

for mode in none preconflict postconflict both; do
  "$BIN" sketch -p 2 --ctxmeta "$mode" -o "$WORK/ctx_$mode" "$WORK/a.fna" "$WORK/b.fna" > "$WORK/ctx_$mode.log" 2>&1
  assert_nonempty "$WORK/ctx_$mode/comblco"
  if [[ "$mode" == "none" ]]; then
    test ! -e "$WORK/ctx_$mode/minco.ctxmeta.tsv"
  else
    assert_nonempty "$WORK/ctx_$mode/minco.ctxmeta.tsv"
  fi
done

{
  sed 's/^>a/>multi_a/' "$WORK/a.fna"
  sed 's/^>b/>multi_b/' "$WORK/b.fna"
} > "$WORK/multi.fna"
"$BIN" sketch --splitmfa -p 2 -o "$WORK/splitmfa" "$WORK/multi.fna" > "$WORK/splitmfa.log" 2>&1
"$BIN" sketch --psmp "$WORK/splitmfa" > "$WORK/splitmfa.psmp.tsv"
test "$(wc -l < "$WORK/splitmfa.psmp.tsv")" -eq 2

"$BIN" sketch --asone -p 2 -o "$WORK/asone" "$WORK/a.fna" "$WORK/b.fna" > "$WORK/asone.log" 2>&1
"$BIN" sketch --psmp "$WORK/asone" > "$WORK/asone.psmp.tsv"
test "$(wc -l < "$WORK/asone.psmp.tsv")" -eq 1

cat "$WORK/a.fna" | "$BIN" sketch -p 1 -o "$WORK/stdin" - > "$WORK/stdin.log" 2>&1
"$BIN" sketch --pipecmd 'cat {}' -p 1 -o "$WORK/pipecmd" "$WORK/a.fna" > "$WORK/pipecmd.log" 2>&1
assert_nonempty "$WORK/stdin/comblco"
assert_nonempty "$WORK/pipecmd/comblco"

"$BIN" sketch --conflict --readsQC --qc-hash-target 2000 -p 2 --ctxmeta both -o "$WORK/reads" "$WORK/reads.fq" > "$WORK/reads.log" 2>&1
assert_nonempty "$WORK/reads/comblco"
assert_nonempty "$WORK/reads/lcofiles.qc"
assert_nonempty "$WORK/reads/minco.ctxmeta.tsv"

"$BIN" sketch -A --conflict --readsQC --qc-hash-target 2000 -p 2 -o "$WORK/reads_abund" "$WORK/reads.fq" > "$WORK/reads_abund.log" 2>&1
"$BIN" sketch --sketchQC -o "$WORK/reads_abund_qc" "$WORK/reads_abund" > "$WORK/reads_abund_qc.log" 2>&1
assert_nonempty "$WORK/reads_abund/comblco.a"
assert_nonempty "$WORK/reads_abund_qc/comblco"

"$BIN" ani -r "$WORK/base" -q "$WORK/base" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_detail.tsv"
"$BIN" ani -q "$WORK/base" -m 1 -s -1 -d -p 2 -o "$WORK/ani_matrix.tsv"
"$BIN" ani -q "$WORK/base" -m 2 -s -1 -d -p 2 -o "$WORK/ani_triangle.tsv"
"$BIN" ani --pair -f 0 -n 0 -t 0 -o "$WORK/ani_pair.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -f 0 -n 0 -t 0 -o "$WORK/ani_positional.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -S 5000 -f 0 -n 0 -t 0 -o "$WORK/ani_positional_5k.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -r "$WORK/base" -q "$WORK/b.fna" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_seq_query.tsv"
"$BIN" ani -r "$WORK/base" --qraw "$WORK/reads" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_qraw.tsv"
assert_nonempty "$WORK/ani_detail.tsv"
assert_nonempty "$WORK/ani_matrix.tsv"
assert_nonempty "$WORK/ani_triangle.tsv"
assert_nonempty "$WORK/ani_pair.tsv"
assert_nonempty "$WORK/ani_positional.tsv"
assert_nonempty "$WORK/ani_positional_5k.tsv"
assert_nonempty "$WORK/ani_seq_query.tsv"
assert_nonempty "$WORK/ani_qraw.tsv"

printf '%s\n' "$WORK/a.fna" "$WORK/b.fna" > "$WORK/ref.list"
printf '%s\n' "$WORK/b.fna" "$WORK/c.fna" > "$WORK/qry.list"
"$BIN" ani --reflist "$WORK/ref.list" --qrylist "$WORK/qry.list" -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_lists.tsv"
assert_nonempty "$WORK/ani_lists.tsv"

"$BIN" matrix --format full -q "$WORK/base" -d -o "$WORK/matrix_full.tsv"
"$BIN" matrix --format triangle -q "$WORK/base" -d -o "$WORK/matrix_triangle.tsv"
"$BIN" matrix -r "$WORK/kept" -q "$WORK/removed" --format full -o "$WORK/matrix_rect.tsv"
"$BIN" matrix --format edges --cut 0.2 -q "$WORK/base" -o "$WORK/matrix_edges.tsv"
"$BIN" matrix --format clusters --cut 0.2 -q "$WORK/base" -o "$WORK/matrix_clusters.tsv"
"$BIN" matrix --format dedup-plan --cut 0.001 --keep-out "$WORK/plan_keep.txt" \
  --remove-out "$WORK/plan_remove.txt" --edge-out "$WORK/plan_edges.tsv" \
  --keep-matrix-out "$WORK/plan_keep_matrix.tsv" -q "$WORK/base" -o "$WORK/matrix_dedup_plan.tsv"
"$BIN" matrix --format full --matrix-format phylip --matrix-idmap "$WORK/matrix.idmap.tsv" \
  -q "$WORK/base" -d -o "$WORK/matrix.phy"
assert_nonempty "$WORK/matrix_full.tsv"
assert_nonempty "$WORK/matrix_triangle.tsv"
assert_nonempty "$WORK/matrix_rect.tsv"
assert_nonempty "$WORK/matrix_edges.tsv"
assert_nonempty "$WORK/matrix_clusters.tsv"
assert_nonempty "$WORK/matrix_dedup_plan.tsv"
assert_nonempty "$WORK/matrix.phy"
assert_nonempty "$WORK/matrix.idmap.tsv"
