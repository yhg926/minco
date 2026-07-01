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
cp "$WORK/a.fna" "$WORK/a_dup.fna"

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
"$BIN" profile --help > "$WORK/profile_help.txt"
grep -q -- '--reads' "$WORK/profile_help.txt"
grep -q -- '--track' "$WORK/profile_help.txt"
grep -q 'conservative direct' "$WORK/profile_help.txt"
env PYTHONDONTWRITEBYTECODE=1 python3 -B "$ROOT/scripts/minco_profile" --help \
  > "$WORK/profile_default_help.txt" 2>&1
grep -q 'recommended no-manual-strategy entry point' "$WORK/profile_default_help.txt"
grep -q 'scripts/minco_profile -r ref.minco' "$WORK/profile_default_help.txt"
grep -q 'scripts/minco_profile --check-ref -r ref.minco' "$WORK/profile_default_help.txt"
grep -q 'sidecar beside --ref' "$WORK/profile_default_help.txt"
grep -q 'AMR/gene/virus/mixed-domain profiling' "$WORK/profile_default_help.txt"
env PYTHONDONTWRITEBYTECODE=1 python3 -B "$ROOT/tests/test_minco_profile_calibrated_auto_exact.py"

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

"$BIN" sketch -p 2 --sketch-size 31 -o "$WORK/dup.minco" "$WORK/a.fna" "$WORK/a_dup.fna" > "$WORK/dup.log" 2>&1
"$BIN" set --uniq_union --markerdb -o "$WORK/dup.marker.minco" "$WORK/dup.minco" > "$WORK/dup.marker.log" 2>&1
grep -q "WARNING: 2/2 refs have markerdb sketch size below 500" "$WORK/dup.marker.log"
test "$(wc -c < "$WORK/dup.marker.minco/minco.ctxobj64")" -eq 0
"$BIN" set --uniq_union --markerdb-ctx -o "$WORK/dup.ctxmarker.minco" "$WORK/dup.minco" > "$WORK/dup.ctxmarker.log" 2>&1
grep -q "WARNING: 2/2 refs have markerdb sketch size below 500" "$WORK/dup.ctxmarker.log"
grep -q $'low-marker-ref\t0\t0\t' "$WORK/dup.ctxmarker.log"
grep -q $'low-marker-ref\t1\t0\t' "$WORK/dup.ctxmarker.log"
test "$(wc -c < "$WORK/dup.ctxmarker.minco/minco.ctxobj64")" -eq 0
test "$(wc -c < "$WORK/dup.ctxmarker.minco/minco.ctxobj64.offsets")" -eq 24
"$BIN" sketch --psmp "$WORK/dup.ctxmarker.minco" > "$WORK/dup.ctxmarker.psmp.tsv"
awk '$1 != 0 { exit 1 }' "$WORK/dup.ctxmarker.psmp.tsv"
"$BIN" set --uniq_union --markerdb-ctx --markerdb-ctx-pairwise \
  --markerdb-ctx-min-xny 100 --markerdb-ctx-min-af 0.01 \
  -o "$WORK/dup.ctxmarker.pairkeep.minco" "$WORK/dup.minco" \
  > "$WORK/dup.ctxmarker.pairkeep.log" 2>&1
"$BIN" sketch --psmp "$WORK/dup.ctxmarker.pairkeep.minco" \
  > "$WORK/dup.ctxmarker.pairkeep.psmp.tsv"
test "$(awk '{s += $1} END { print s }' "$WORK/dup.ctxmarker.pairkeep.psmp.tsv")" -eq 62
grep -q "accepted_pairs=0" "$WORK/dup.ctxmarker.pairkeep.log"
"$BIN" set --uniq_union --markerdb-ctx --markerdb-ctx-pairwise \
  --markerdb-ctx-min-xny 10 --markerdb-ctx-min-af 0.01 \
  -o "$WORK/dup.ctxmarker.pairdrop.minco" "$WORK/dup.minco" \
  > "$WORK/dup.ctxmarker.pairdrop.log" 2>&1
"$BIN" sketch --psmp "$WORK/dup.ctxmarker.pairdrop.minco" \
  > "$WORK/dup.ctxmarker.pairdrop.psmp.tsv"
test "$(awk '{s += $1} END { print s }' "$WORK/dup.ctxmarker.pairdrop.psmp.tsv")" -eq 0
grep -q "accepted_pairs=1" "$WORK/dup.ctxmarker.pairdrop.log"

"$BIN" sketch -p 1 --sketch-size 5 -o "$WORK/size5.minco" "$WORK/a.fna" > "$WORK/size5.log" 2>&1
"$BIN" sketch --psmp "$WORK/size5.minco" > "$WORK/size5.psmp.tsv"
test "$(awk 'NR == 1 { print $1 }' "$WORK/size5.psmp.tsv")" -eq 5
"$BIN" set --uniq_union --markerdb -o "$WORK/size5.marker.minco" "$WORK/size5.minco" > "$WORK/size5.marker.log" 2>&1
! grep -q "only 1 sketch" "$WORK/size5.marker.log"
"$BIN" sketch --psmp "$WORK/size5.marker.minco" > "$WORK/size5.marker.psmp.tsv"
test "$(awk 'NR == 1 { print $1 }' "$WORK/size5.marker.psmp.tsv")" -eq 5

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

cat > "$WORK/gtdb.taxmap.tsv" <<EOF
ref_key	TAXID	RANK	TAXPATH	TAXPATHSN
a.fna	1	species	d__Bacteria|p__Test|c__Test|o__Test|f__Test|g__Minco|s__a	d__Bacteria|p__Test|c__Test|o__Test|f__Test|g__Minco|s__a
b.fna	2	species	d__Bacteria|p__Test|c__Test|o__Test|f__Test|g__Minco|s__b	d__Bacteria|p__Test|c__Test|o__Test|f__Test|g__Minco|s__b
EOF
"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --readwise-track "$WORK/read.track.tsv" \
  --readwise-track-summary "$WORK/read.track.summary.tsv" \
  --readwise-taxonomy gtdb --gtdb-taxmap "$WORK/gtdb.taxmap.tsv" \
  -m0 -f0 -n0 -t0 -o "$WORK/read.track.ani.tsv" \
  > "$WORK/read.track.log" 2>&1
test -s "$WORK/read.track.tsv"
test -s "$WORK/read.track.summary.tsv"
grep -q "forces --density-block-ctx 0" "$WORK/read.track.log"
grep -q $'^read_id\tread_ord\tread_len\tpossible_ctx\tdensity_ctx\tmatched_ctx' "$WORK/read.track.tsv"
awk 'NR > 1 && $8 >= 1 && $11 != "NA" && $15 != "NA" { found = 1 } END { exit found ? 0 : 1 }' "$WORK/read.track.tsv"
grep -q $'^total_reads\t1$' "$WORK/read.track.summary.tsv"
grep -q $'^reads_with_ref_hit\t1$' "$WORK/read.track.summary.tsv"
grep -q $'^density_positive_no_ref_hit_read_pct\t' "$WORK/read.track.summary.tsv"
grep -q $'^sketch_corrected_available\t1$' "$WORK/read.track.summary.tsv"
grep -q $'^estimated_unknown_reads_pct\t' "$WORK/read.track.summary.tsv"
grep -q $'^estimated_unknown_reads_basis\tsketch_corrected_ref_absent_contexts' "$WORK/read.track.summary.tsv"
! grep -q $'^estimated_ref_absent_ctx_pct\t' "$WORK/read.track.summary.tsv"

"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-edge-out "$WORK/read.edges.tsv" \
  --readwise-edge-max 10000 \
  -m0 -f0 -n0 -t0 -o "$WORK/read.edge.ani.tsv" \
  > "$WORK/read.edge.log" 2>&1
test -s "$WORK/read.edges.tsv"
grep -q "forces --density-block-ctx 0" "$WORK/read.edge.log"
grep -q $'^edge_id\tread_id\tunit_id\tqctx\tedge_rank\tref_begin\tgid\t' "$WORK/read.edges.tsv"

"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-unique-out "$WORK/read.unique.sidecar.tsv" \
  -m0 -f0 -n0 -t0 -o "$WORK/read.split.with_unique.tsv" \
  > "$WORK/read.split.with_unique.log" 2>&1
test -s "$WORK/read.unique.sidecar.tsv"
test -s "$WORK/read.split.with_unique.tsv"
grep -q "wrote best-diff-unique sidecar" "$WORK/read.split.with_unique.log"
"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique \
  -m0 -f0 -n0 -t0 -o "$WORK/read.unique.separate.tsv" \
  > "$WORK/read.unique.separate.log" 2>&1
python3 - "$WORK/read.unique.sidecar.tsv" "$WORK/read.unique.separate.tsv" <<'PY'
import csv
import sys

sidecar, separate = sys.argv[1:3]
cols = [
    "Qry", "Ref", "ANI", "XnY_ctx", "N_diff_obj", "N_diff_obj_section",
    "N_mut2_ctx", "Unique_ref_ctx_hit", "Raw_XnY_ctx",
    "Ref_breadth", "Ref_mean_depth", "Ref_hit_mean_depth", "Ref_zip_af",
    "Relative_abundance_depth", "Normalized_abundance_depth",
]

def read_rows(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in cols if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        return sorted(
            ({c: row[c] for c in cols} for row in reader),
            key=lambda row: (row["Qry"], row["Ref"]),
        )

a = read_rows(sidecar)
b = read_rows(separate)
if a != b:
    raise SystemExit(f"unique sidecar differs from separate unique pass: {sidecar} {separate}")
PY

"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --density-block-ctx 100 \
  --readwise-exact-split-out "$WORK/read.exact_split.sidecar.tsv" \
  -m0 -f0 -n0 -t0 -o "$WORK/read.split.with_exact_sidecar.tsv" \
  > "$WORK/read.split.with_exact_sidecar.log" 2>&1
test -s "$WORK/read.exact_split.sidecar.tsv"
grep -q "wrote exact best-diff-split sidecar" "$WORK/read.split.with_exact_sidecar.log"
"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --density-block-ctx 0 \
  -m0 -f0 -n0 -t0 -o "$WORK/read.exact_split.separate.tsv" \
  > "$WORK/read.exact_split.separate.log" 2>&1
python3 - "$WORK/read.exact_split.sidecar.tsv" "$WORK/read.exact_split.separate.tsv" <<'PY'
import csv
import sys

sidecar, separate = sys.argv[1:3]
cols = [
    "Qry", "Ref", "ANI", "XnY_ctx", "N_diff_obj", "N_diff_obj_section",
    "N_mut2_ctx", "Unique_ref_ctx_hit", "Raw_XnY_ctx",
    "Ref_breadth", "Ref_mean_depth", "Ref_hit_mean_depth", "Ref_zip_af",
    "Relative_abundance_depth", "Normalized_abundance_depth",
]

def read_rows(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in cols if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        return sorted(
            ({c: row[c] for c in cols} for row in reader),
            key=lambda row: (row["Qry"], row["Ref"]),
        )

a = read_rows(sidecar)
b = read_rows(separate)
if a != b:
    raise SystemExit(f"exact split sidecar differs from separate exact pass: {sidecar} {separate}")
PY

"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --density-block-ctx 100 \
  --readwise-density-cache-out "$WORK/read.density.cache.bin" \
  -m0 -f0 -n0 -t0 -o "$WORK/read.split.with_density_cache.tsv" \
  > "$WORK/read.split.with_density_cache.log" 2>&1
test -s "$WORK/read.density.cache.bin"
grep -q "wrote density cache" "$WORK/read.split.with_density_cache.log"
"$BIN" ani -r "$WORK/pair.minco" --qraw "$WORK/a.fna" --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --density-block-ctx 0 \
  --readwise-density-cache-in "$WORK/read.density.cache.bin" \
  -m0 -f0 -n0 -t0 -o "$WORK/read.exact_split.cache_replay.tsv" \
  > "$WORK/read.exact_split.cache_replay.log" 2>&1
grep -q "density cache replay active" "$WORK/read.exact_split.cache_replay.log"
python3 - "$WORK/read.exact_split.cache_replay.tsv" "$WORK/read.exact_split.separate.tsv" <<'PY'
import csv
import sys

cache_replay, separate = sys.argv[1:3]
cols = [
    "Qry", "Ref", "ANI", "XnY_ctx", "N_diff_obj", "N_diff_obj_section",
    "N_mut2_ctx", "Unique_ref_ctx_hit", "Raw_XnY_ctx",
    "Ref_breadth", "Ref_mean_depth", "Ref_hit_mean_depth", "Ref_zip_af",
    "Relative_abundance_depth", "Normalized_abundance_depth",
]

def read_rows(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in cols if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        return sorted(
            ({c: row[c] for c in cols} for row in reader),
            key=lambda row: (row["Qry"], row["Ref"]),
        )

a = read_rows(cache_replay)
b = read_rows(separate)
if a != b:
    raise SystemExit(f"density cache replay differs from separate exact pass: {cache_replay} {separate}")
PY

"$BIN" profile -r "$WORK/pair.minco" "$WORK/a.fna" -p 2 \
  -o "$WORK/profile.tsv" > "$WORK/profile.log" 2>&1
test -s "$WORK/profile.tsv"
grep -q 'Relative_abundance_depth' "$WORK/profile.tsv"
grep -q 'Default_call' "$WORK/profile.tsv"
grep -q 'major' "$WORK/profile.tsv"

"$BIN" profile -r "$WORK/pair.minco" "$WORK/a.fna" \
  --track "$WORK/profile.track.tsv" \
  --track-summary "$WORK/profile.track.summary.tsv" \
  --taxonomy gtdb --gtdb-taxmap "$WORK/gtdb.taxmap.tsv" \
  --report-all -o "$WORK/profile.track.profile.tsv" \
  > "$WORK/profile.track.log" 2>&1
test -s "$WORK/profile.track.tsv"
test -s "$WORK/profile.track.summary.tsv"
grep -q "forces --density-block-ctx 0" "$WORK/profile.track.log"
grep -q $'^read_id\tread_ord\tread_len\tpossible_ctx\tdensity_ctx\tmatched_ctx' "$WORK/profile.track.tsv"
grep -q $'^estimated_unknown_reads_basis\tsketch_corrected_ref_absent_contexts' "$WORK/profile.track.summary.tsv"

"$BIN" profile -r "$WORK/pair.minco" "$WORK/a.fna" \
  --edge-out "$WORK/profile.edges.tsv" \
  --edge-max 10000 \
  --report-all -o "$WORK/profile.edge.profile.tsv" \
  > "$WORK/profile.edge.log" 2>&1
test -s "$WORK/profile.edges.tsv"
grep -q "forces --density-block-ctx 0" "$WORK/profile.edge.log"
grep -q $'^edge_id\tread_id\tunit_id\tqctx\tedge_rank\tref_begin\tgid\t' "$WORK/profile.edges.tsv"

"$BIN" ani -q "$WORK/pair.minco" -m2 -s -1 -d -p 2 -o "$WORK/ani.tsv" > "$WORK/ani.log" 2>&1
test -s "$WORK/ani.tsv"
test "$(wc -l < "$WORK/ani.tsv")" -ge 2

"$BIN" matrix --format triangle -q "$WORK/pair.minco" --diagonal -o "$WORK/matrix.tsv" > "$WORK/matrix.log" 2>&1
test -s "$WORK/matrix.tsv"
