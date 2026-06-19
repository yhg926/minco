#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/minco_smoke.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

BIN="$ROOT/bin/minco"

make_fasta() {
  local out="$1"
  local name="$2"
  local mutate="$3"

  awk -v name="$name" -v mutate="$mutate" '
    BEGIN {
      bases = "ACGT";
      x = 17;
      print ">" name;
      for (i = 1; i <= 24000; i++) {
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

make_fasta "$WORK/a.fna" a 0
make_fasta "$WORK/b.fna" b 1

"$BIN" --help > "$WORK/help.txt"
"$BIN" examples > "$WORK/examples.txt"
"$BIN" doctor > "$WORK/doctor.txt"

for hidden in shuffle dist reverse; do
  if "$BIN" "$hidden" --help > "$WORK/${hidden}.txt" 2>&1; then
    echo "hidden command unexpectedly succeeded: $hidden" >&2
    exit 1
  fi
done
"$BIN" set --help > "$WORK/set_help.txt"
grep -q -- '--downsample' "$WORK/set_help.txt"

"$BIN" sketch -p 2 --ctxmeta both -o "$WORK/pair.minco" "$WORK/a.fna" "$WORK/b.fna" > "$WORK/sketch.log" 2>&1
test -s "$WORK/pair.minco/minco.ctxobj64"
test -s "$WORK/pair.minco/minco.ctxobj64.offsets"
test -s "$WORK/pair.minco/minco.stat"
test -s "$WORK/pair.minco/minco.ctxmeta"
test ! -e "$WORK/pair.minco/minco.ctxmeta.tsv"
test ! -e "$WORK/pair.minco/minco.ctxsetmeta.tsv"
"$BIN" sketch --pctxmeta "$WORK/pair.minco" > "$WORK/pair.ctxmeta.tsv"
"$BIN" sketch --pctxsetmeta "$WORK/pair.minco" > "$WORK/pair.ctxsetmeta.tsv"
test "$(wc -l < "$WORK/pair.ctxmeta.tsv")" -eq 3
grep -q $'^universal_density_policy\tlargest_sample_density$' "$WORK/pair.ctxsetmeta.tsv"

"$BIN" set --downsample -S 5 -o "$WORK/pair.S5.minco" "$WORK/pair.minco" > "$WORK/downsample.log" 2>&1
"$BIN" sketch --psmp "$WORK/pair.S5.minco" > "$WORK/pair.S5.psmp.tsv"
test "$(wc -l < "$WORK/pair.S5.psmp.tsv")" -eq 2
awk '$1 > 5 { exit 1 }' "$WORK/pair.S5.psmp.tsv"
"$BIN" sketch -i "$WORK/pair.S5.minco" > "$WORK/pair.S5.index.log" 2>&1
test -s "$WORK/pair.S5.minco/minco.refindex.ctxgid64obj32"

"$BIN" sketch -p 1 --sketch-size 5 -o "$WORK/size5.minco" "$WORK/a.fna" > "$WORK/size5.log" 2>&1
"$BIN" sketch --psmp "$WORK/size5.minco" > "$WORK/size5.psmp.tsv"
test "$(awk 'NR == 1 { print $1 }' "$WORK/size5.psmp.tsv")" -eq 5

mkdir "$WORK/legacy_names.minco"
cp "$WORK/size5.minco/minco.ctxobj64" "$WORK/legacy_names.minco/comblco"
cp "$WORK/size5.minco/minco.ctxobj64.offsets" "$WORK/legacy_names.minco/comblco.index"
cp "$WORK/size5.minco/minco.stat" "$WORK/legacy_names.minco/lcofiles.stat"
"$BIN" sketch --psmp "$WORK/legacy_names.minco" > "$WORK/legacy_names.psmp.tsv" 2> "$WORK/legacy_names.psmp.err"
grep -q "using legacy sketch filename" "$WORK/legacy_names.psmp.err"
test "$(awk 'NR == 1 { print $1 }' "$WORK/legacy_names.psmp.tsv")" -eq 5
"$BIN" sketch -i "$WORK/legacy_names.minco" > "$WORK/legacy_names.index.log" 2>&1
test -s "$WORK/legacy_names.minco/minco.refindex.ctxgid64obj32"
"$BIN" ani -r "$WORK/legacy_names.minco" -q "$WORK/size5.minco" -m0 -f0 -n0 -t0 -o "$WORK/legacy_names.ani.tsv" > "$WORK/legacy_names.ani.log" 2>&1
test -s "$WORK/legacy_names.ani.tsv"

"$BIN" sketch -i "$WORK/pair.minco" > "$WORK/index.log" 2>&1
test -s "$WORK/pair.minco/minco.refindex.ctxgid64obj32"

"$BIN" ani -q "$WORK/pair.minco" -m2 -s -1 -d -p 2 -o "$WORK/ani.tsv" > "$WORK/ani.log" 2>&1
test -s "$WORK/ani.tsv"
test "$(wc -l < "$WORK/ani.tsv")" -ge 2

"$BIN" matrix --format triangle -q "$WORK/pair.minco" --diagonal -o "$WORK/matrix.tsv" > "$WORK/matrix.log" 2>&1
test -s "$WORK/matrix.tsv"
