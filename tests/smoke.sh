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

for hidden in set shuffle dist reverse; do
  if "$BIN" "$hidden" --help > "$WORK/${hidden}.txt" 2>&1; then
    echo "hidden command unexpectedly succeeded: $hidden" >&2
    exit 1
  fi
done

"$BIN" sketch -p 2 --ctxmeta both -o "$WORK/pair.minco" "$WORK/a.fna" "$WORK/b.fna" > "$WORK/sketch.log" 2>&1
test -s "$WORK/pair.minco/comblco"
test -s "$WORK/pair.minco/lcofiles.stat"
test -s "$WORK/pair.minco/minco.ctxmeta.tsv"
test "$(wc -l < "$WORK/pair.minco/minco.ctxmeta.tsv")" -eq 3

"$BIN" sketch -p 1 --sketch-size 5 -o "$WORK/size5.minco" "$WORK/a.fna" > "$WORK/size5.log" 2>&1
"$BIN" sketch --psmp "$WORK/size5.minco" > "$WORK/size5.psmp.tsv"
test "$(awk 'NR == 1 { print $1 }' "$WORK/size5.psmp.tsv")" -eq 5

"$BIN" sketch -i "$WORK/pair.minco" > "$WORK/index.log" 2>&1
test -s "$WORK/pair.minco/comblco.index"

"$BIN" ani -q "$WORK/pair.minco" -m2 -s -1 -d -p 2 -o "$WORK/ani.tsv" > "$WORK/ani.log" 2>&1
test -s "$WORK/ani.tsv"
test "$(wc -l < "$WORK/ani.tsv")" -ge 2

"$BIN" matrix --format triangle -q "$WORK/pair.minco" --diagonal -o "$WORK/matrix.tsv" > "$WORK/matrix.log" 2>&1
test -s "$WORK/matrix.tsv"
