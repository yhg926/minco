#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 SELECTED_GENOMES.tsv OUTDIR [THREADS]" >&2
  exit 2
fi

SELECTED="$1"
OUTDIR="$2"
THREADS="${3:-8}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$OUTDIR"

awk -F'\t' 'NR > 1 { print $3 }' "$SELECTED" > "$OUTDIR/genomes.list"

/usr/bin/time -v "$ROOT/bin/minco" sketch \
  -p "$THREADS" \
  -l "$OUTDIR/genomes.list" \
  -o "$OUTDIR/minco" \
  > "$OUTDIR/sketch.stdout" \
  2> "$OUTDIR/sketch.time.log"

"$ROOT/bin/minco" sketch -i "$OUTDIR/minco" > "$OUTDIR/index.stdout" 2> "$OUTDIR/index.stderr"

/usr/bin/time -v "$ROOT/bin/minco" ani \
  -q "$OUTDIR/minco" \
  -m2 -s -1 -d \
  -p "$THREADS" \
  -o "$OUTDIR/minco_ani.tsv" \
  > "$OUTDIR/ani.stdout" \
  2> "$OUTDIR/ani.time.log"
