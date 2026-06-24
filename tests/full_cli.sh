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

assert_stat_target() {
  local stat_path="$1"
  local expected="$2"
  python3 - "$stat_path" "$expected" <<'PY'
import struct
import sys

path = sys.argv[1]
expected = int(sys.argv[2])
data = open(path, "rb").read()
prefix_size = 32
pathlen = 256
if len(data) < prefix_size:
    raise SystemExit(f"{path}: stat is too small")
infile_num = struct.unpack_from("<i", data, 28)[0]
compat_filter_shift = struct.unpack_from("<i", data, 24)[0]
base = prefix_size + infile_num * pathlen
if len(data) < base + 64:
    raise SystemExit(f"{path}: missing minco stat extension")
magic, version, ext_size, target = struct.unpack_from("<IHHI", data, base)
if magic != 0x4D434F53:
    raise SystemExit(f"{path}: bad minco stat magic {magic:#x}")
if version != 1 or ext_size < 64:
    raise SystemExit(f"{path}: unsupported minco stat extension v{version} size {ext_size}")
if compat_filter_shift != 0:
    raise SystemExit(f"{path}: stale compat_filter_shift stored in minco stat prefix: {compat_filter_shift}")
if target != expected:
    raise SystemExit(f"{path}: target sketch size {target}, expected {expected}")
PY
}

assert_stat_density_summary() {
  local stat_path="$1"
  local expected_valid="$2"
  python3 - "$stat_path" "$expected_valid" <<'PY'
import struct
import sys

path = sys.argv[1]
expected_valid = int(sys.argv[2])
data = open(path, "rb").read()
prefix_size = 32
pathlen = 256
infile_num = struct.unpack_from("<i", data, 28)[0]
base = prefix_size + infile_num * pathlen
if len(data) < base + 120:
    raise SystemExit(f"{path}: missing v1 density summary in minco stat")
magic, version, ext_size = struct.unpack_from("<IHH", data, base)
if magic != 0x4D434F53 or version != 1 or ext_size < 120:
    raise SystemExit(f"{path}: unsupported minco stat extension v{version} size {ext_size}")
min_thr, max_thr, universal_thr = struct.unpack_from("<QQQ", data, base + 64)
valid_count = struct.unpack_from("<I", data, base + 100)[0]
hash_bits = struct.unpack_from("<I", data, base + 104)[0]
policy = struct.unpack_from("<I", data, base + 108)[0]
if valid_count != expected_valid:
    raise SystemExit(f"{path}: density valid count {valid_count}, expected {expected_valid}")
if hash_bits == 0 or hash_bits > 64:
    raise SystemExit(f"{path}: invalid density hash_bits {hash_bits}")
if min_thr > max_thr or universal_thr != max_thr:
    raise SystemExit(f"{path}: invalid density thresholds min={min_thr} max={max_thr} universal={universal_thr}")
expected_policy = 1 if expected_valid == 1 else 2
if policy != expected_policy:
    raise SystemExit(f"{path}: density policy {policy}, expected {expected_policy}")
PY
}

make_fasta "$WORK/a.fna" a 17 0
make_fasta "$WORK/b.fna" b 17 1
cp "$WORK/a.fna" "$WORK/c.fna"
make_fastq "$WORK/reads.fq" reads
awk '
  /^>/ { next }
  { seq = seq $0 }
  END {
    n = length(seq);
    for (r = 1; r <= 80; r++) {
      start = 1 + ((r - 1) * 113) % (n - 220);
      read = substr(seq, start, 220);
      printf "@a_read_%04d\n%s\n+\n", r, read;
      for (i = 1; i <= length(read); i++) printf "I";
      printf "\n";
    }
  }
' "$WORK/a.fna" > "$WORK/a_reads.fq"
printf '%s\n%s\n%s\n' "$WORK/a.fna" "$WORK/b.fna" "$WORK/c.fna" > "$WORK/genomes.list"

"$BIN" examples > "$WORK/examples.txt"
"$BIN" doctor > "$WORK/doctor.txt"

"$BIN" sketch -p 2 --anno --position -l "$WORK/genomes.list" -o "$WORK/base" > "$WORK/base.sketch.log" 2>&1
assert_nonempty "$WORK/base/minco.ctxobj64"
assert_nonempty "$WORK/base/minco.stat"
assert_stat_target "$WORK/base/minco.stat" 10000
assert_nonempty "$WORK/base/minco.anno"
assert_nonempty "$WORK/base/minco.infilemeta"
assert_nonempty "$WORK/base/minco.ctxobj64.position"

"$BIN" sketch -i "$WORK/base" > "$WORK/base.index.log" 2>&1
assert_nonempty "$WORK/base/minco.ctxobj64.offsets"
assert_nonempty "$WORK/base/minco.refindex.ctxgid64obj32"

"$BIN" set --downsample -S 1000 -o "$WORK/base_1k" "$WORK/base" > "$WORK/base_1k.downsample.log" 2>&1
assert_nonempty "$WORK/base_1k/minco.ctxobj64"
assert_nonempty "$WORK/base_1k/minco.ctxobj64.offsets"
assert_nonempty "$WORK/base_1k/minco.stat"
assert_nonempty "$WORK/base_1k/minco.anno"
assert_nonempty "$WORK/base_1k/minco.infilemeta"
assert_nonempty "$WORK/base_1k/minco.ctxobj64.position"
assert_stat_target "$WORK/base_1k/minco.stat" 1000
assert_stat_density_summary "$WORK/base_1k/minco.stat" 3
"$BIN" sketch --psmp "$WORK/base_1k" > "$WORK/base_1k.psmp.tsv"
test "$(wc -l < "$WORK/base_1k.psmp.tsv")" -eq 3
awk '$1 > 1000 { exit 1 }' "$WORK/base_1k.psmp.tsv"
test ! -e "$WORK/base_1k/minco.refindex.ctxgid64obj32"
"$BIN" sketch -i "$WORK/base_1k" > "$WORK/base_1k.index.log" 2>&1
assert_nonempty "$WORK/base_1k/minco.refindex.ctxgid64obj32"

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
assert_nonempty "$WORK/kept/minco.ctxobj64"
assert_nonempty "$WORK/removed/minco.ctxobj64"
assert_nonempty "$WORK/kept_drop/minco.ctxobj64"
test ! -e "$WORK/kept_drop/minco.ctxobj64.position"

"$BIN" sketch --append -o "$WORK/appended" "$WORK/kept" "$WORK/removed" > "$WORK/append.log"
assert_nonempty "$WORK/appended/minco.ctxobj64"

"$BIN" sketch --dedup 0.001 --metric ctx-naive -o "$WORK/dedup" "$WORK/base" > "$WORK/dedup.log"
"$BIN" sketch --dedup 0.001 --metric ctx-naive --dedup-index -o "$WORK/dedup_index" "$WORK/base" > "$WORK/dedup_index.log"
assert_nonempty "$WORK/dedup/minco.ctxobj64"
assert_nonempty "$WORK/dedup_index/minco.ctxobj64"

for mode in none preconflict postconflict both; do
  "$BIN" sketch -p 2 --ctxmeta "$mode" -o "$WORK/ctx_$mode" "$WORK/a.fna" "$WORK/b.fna" > "$WORK/ctx_$mode.log" 2>&1
  assert_nonempty "$WORK/ctx_$mode/minco.ctxobj64"
  if [[ "$mode" == "none" ]]; then
    test ! -e "$WORK/ctx_$mode/minco.ctxmeta"
    test ! -e "$WORK/ctx_$mode/minco.ctxmeta.tsv"
    test ! -e "$WORK/ctx_$mode/minco.ctxsetmeta.tsv"
  else
    assert_nonempty "$WORK/ctx_$mode/minco.ctxmeta"
    test ! -e "$WORK/ctx_$mode/minco.ctxmeta.tsv"
    test ! -e "$WORK/ctx_$mode/minco.ctxsetmeta.tsv"
    "$BIN" sketch --pctxmeta "$WORK/ctx_$mode" > "$WORK/ctx_$mode.ctxmeta.tsv"
    "$BIN" sketch --pctxsetmeta "$WORK/ctx_$mode" > "$WORK/ctx_$mode.ctxsetmeta.tsv"
    test "$(wc -l < "$WORK/ctx_$mode.ctxmeta.tsv")" -eq 3
    grep -q $'^min_sample_density\t' "$WORK/ctx_$mode.ctxsetmeta.tsv"
    grep -q $'^largest_sample_density\t' "$WORK/ctx_$mode.ctxsetmeta.tsv"
    assert_stat_density_summary "$WORK/ctx_$mode/minco.stat" 2
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
assert_nonempty "$WORK/stdin/minco.ctxobj64"
assert_nonempty "$WORK/pipecmd/minco.ctxobj64"

"$BIN" sketch --conflict --readsQC --qc-hash-target 2000 -p 2 --ctxmeta both -o "$WORK/reads" "$WORK/reads.fq" > "$WORK/reads.log" 2>&1
assert_nonempty "$WORK/reads/minco.ctxobj64"
assert_nonempty "$WORK/reads/minco.qc"
assert_nonempty "$WORK/reads/minco.ctxmeta"
test ! -e "$WORK/reads/minco.ctxmeta.tsv"

"$BIN" sketch -A --conflict --readsQC --qc-hash-target 2000 -p 2 -o "$WORK/reads_abund" "$WORK/reads.fq" > "$WORK/reads_abund.log" 2>&1
"$BIN" sketch --sketchQC -o "$WORK/reads_abund_qc" "$WORK/reads_abund" > "$WORK/reads_abund_qc.log" 2>&1
assert_nonempty "$WORK/reads_abund/minco.ctxobj64.abund"
assert_nonempty "$WORK/reads_abund_qc/minco.ctxobj64"

"$BIN" ani -r "$WORK/base" -q "$WORK/base" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_detail.tsv"
"$BIN" ani -q "$WORK/base" -m 1 -s -1 -d -p 2 -o "$WORK/ani_matrix.tsv"
"$BIN" ani -q "$WORK/base" -m 2 -s -1 -d -p 2 -o "$WORK/ani_triangle.tsv"
"$BIN" ani --pair -f 0 -n 0 -t 0 -o "$WORK/ani_pair.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -f 0 -n 0 -t 0 -o "$WORK/ani_positional.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -S 5000 -f 0 -n 0 -t 0 -o "$WORK/ani_positional_5k.tsv" "$WORK/a.fna" "$WORK/b.fna"
"$BIN" ani -r "$WORK/base" -q "$WORK/b.fna" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_seq_query.tsv"
"$BIN" sketch -p 1 --ctxmeta both -o "$WORK/ctx_one" "$WORK/a.fna" > "$WORK/ctx_one.log" 2>&1
assert_stat_target "$WORK/ctx_one/minco.stat" 10000
assert_stat_density_summary "$WORK/ctx_one/minco.stat" 1
"$BIN" ani -r "$WORK/ctx_one" -q "$WORK/b.fna" --query-density ref --save-query-sketch "$WORK/saved_query_one" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_query_density_one.tsv"
assert_stat_target "$WORK/saved_query_one/minco.stat" 10000
assert_stat_density_summary "$WORK/saved_query_one/minco.stat" 1
"$BIN" sketch -p 1 --ctxmeta both --sketch-size 5000 -o "$WORK/ctx_one_5k" "$WORK/a.fna" > "$WORK/ctx_one_5k.log" 2>&1
assert_stat_target "$WORK/ctx_one_5k/minco.stat" 5000
assert_stat_density_summary "$WORK/ctx_one_5k/minco.stat" 1
"$BIN" ani -r "$WORK/ctx_one_5k" -q "$WORK/b.fna" --query-density ref --save-query-sketch "$WORK/saved_query_5k" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_query_density_5k.tsv"
assert_stat_target "$WORK/saved_query_5k/minco.stat" 5000
assert_stat_density_summary "$WORK/saved_query_5k/minco.stat" 1
"$BIN" sketch -p 1 --sketch-size 7000 -o "$WORK/query_7k" "$WORK/b.fna" > "$WORK/query_7k.log" 2>&1
assert_stat_target "$WORK/query_7k/minco.stat" 7000
if "$BIN" ani -r "$WORK/ctx_one_5k" -q "$WORK/query_7k" -m 0 -f 0 -n 0 -t 0 -o "$WORK/ani_mismatch.tsv" > "$WORK/ani_mismatch.log" 2>&1; then
  echo "mismatched sketch sizes unexpectedly succeeded" >&2
  exit 1
fi
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_readwise_density.tsv"
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref --abundance-est depth -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_readwise_abundance.tsv"
{
  printf '%s\t2\tsuperkingdom\t2\tBacteria\n' "$WORK/a.fna"
  printf '%s\t111\tspecies\t2|111\tBacteria|a_species\n' "$WORK/a.fna"
} > "$WORK/cami_taxmap.tsv"
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref --abundance-est depth \
  --cami-taxmap "$WORK/cami_taxmap.tsv" --cami-profile "$WORK/ani_readwise_abundance.profile" \
  --cami-sample-id full_cli_sample -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_readwise_abundance_cami.tsv"
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap "$WORK/cami_taxmap.tsv" \
  --cami-profile "$WORK/ani_readwise_abundance_profile_only.profile" \
  --cami-sample-id full_cli_profile_only -m 0 -f 0 -n 0 -t 0 -p 2 \
  -o "$WORK/ani_readwise_abundance_profile_only.tsv"
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref --abundance-est depth \
  --readwise-profile-only --readwise-ani naive \
  --readwise-ctx-filter poisson-depth --readwise-fake-threshold 1.30103 \
  -m 0 -f 0 -n 0 -t 0 -p 2 \
  -o "$WORK/ani_readwise_poisson_depth.tsv"
"$BIN" ani -r "$WORK/ctx_one" --qraw "$WORK/a_reads.fq" --query-density ref --abundance-est depth \
  --readwise-profile-only --readwise-ani naive \
  --readwise-ctx-filter poisson-product --readwise-fake-threshold 1.30103 \
  -m 0 -f 0 -n 0 -t 0 -p 2 \
  -o "$WORK/ani_readwise_poisson_product.tsv"
"$BIN" ani -r "$WORK/ctx_both" -q "$WORK/b.fna" --query-density ref -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_query_density_combined.tsv"
"$BIN" ani -r "$WORK/base" --qraw "$WORK/reads" -m 0 -f 0 -n 0 -t 0 -p 2 -o "$WORK/ani_qraw.tsv"
assert_nonempty "$WORK/ani_detail.tsv"
assert_nonempty "$WORK/ani_matrix.tsv"
assert_nonempty "$WORK/ani_triangle.tsv"
assert_nonempty "$WORK/ani_pair.tsv"
assert_nonempty "$WORK/ani_positional.tsv"
assert_nonempty "$WORK/ani_positional_5k.tsv"
assert_nonempty "$WORK/ani_seq_query.tsv"
assert_nonempty "$WORK/ani_query_density_one.tsv"
assert_nonempty "$WORK/ani_query_density_5k.tsv"
assert_nonempty "$WORK/ani_readwise_density.tsv"
assert_nonempty "$WORK/ani_readwise_abundance.tsv"
assert_nonempty "$WORK/ani_readwise_abundance_cami.tsv"
assert_nonempty "$WORK/ani_readwise_abundance.profile"
assert_nonempty "$WORK/ani_readwise_abundance_profile_only.tsv"
assert_nonempty "$WORK/ani_readwise_abundance_profile_only.profile"
assert_nonempty "$WORK/ani_readwise_poisson_depth.tsv"
assert_nonempty "$WORK/ani_readwise_poisson_product.tsv"
assert_nonempty "$WORK/ani_query_density_combined.tsv"
assert_nonempty "$WORK/ani_qraw.tsv"
grep -q 'Reads_with_ctx_match' "$WORK/ani_readwise_density.tsv"
grep -q 'readwise_coverage' "$WORK/ani_readwise_density.tsv"
grep -q 'Relative_abundance_depth' "$WORK/ani_readwise_abundance.tsv"
grep -q 'Normalized_abundance_depth' "$WORK/ani_readwise_abundance.tsv"
grep -q 'Ref_mean_depth' "$WORK/ani_readwise_abundance.tsv"
grep -q 'Rejected_ctx' "$WORK/ani_readwise_poisson_depth.tsv"
grep -q 'Fake_ctx_fraction' "$WORK/ani_readwise_poisson_depth.tsv"
grep -q 'Rejected_ctx' "$WORK/ani_readwise_poisson_product.tsv"
grep -q 'Fake_ctx_fraction' "$WORK/ani_readwise_poisson_product.tsv"
grep -q '^@SampleID:full_cli_sample$' "$WORK/ani_readwise_abundance.profile"
grep -q '^@@TAXID' "$WORK/ani_readwise_abundance.profile"
grep -q $'^2\tsuperkingdom\t2\tBacteria\t' "$WORK/ani_readwise_abundance.profile"
grep -q $'^111\tspecies\t2|111\tBacteria|a_species\t' "$WORK/ani_readwise_abundance.profile"
grep -q '^@SampleID:full_cli_profile_only$' "$WORK/ani_readwise_abundance_profile_only.profile"
grep -q '^@@TAXID' "$WORK/ani_readwise_abundance_profile_only.profile"
grep -q $'^2\tsuperkingdom\t2\tBacteria\t' "$WORK/ani_readwise_abundance_profile_only.profile"
grep -q $'^111\tspecies\t2|111\tBacteria|a_species\t' "$WORK/ani_readwise_abundance_profile_only.profile"
awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == "Normalized_abundance_depth") col = i
    if (!col) exit 2
    next
  }
  { sum += $col; n++ }
  END {
    if (n < 1 || sum < 0.999 || sum > 1.001) exit 1
  }
' "$WORK/ani_readwise_abundance.tsv"
awk -F '\t' '
  !/^@/ { sum[$2] += $5; n[$2]++ }
  END {
    if (!("superkingdom" in n) || !("species" in n)) exit 1
    for (rank in n) {
      if (sum[rank] < 99.999 || sum[rank] > 100.001) exit 1
    }
  }
' "$WORK/ani_readwise_abundance.profile"
assert_nonempty "$WORK/saved_query_one/minco.ctxobj64"
assert_nonempty "$WORK/saved_query_one/minco.ctxobj64.offsets"
assert_nonempty "$WORK/saved_query_one/minco.stat"
assert_nonempty "$WORK/saved_query_one/minco.ctxmeta"
assert_nonempty "$WORK/saved_query_5k/minco.stat"
assert_nonempty "$WORK/saved_query_5k/minco.ctxmeta"

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
  --keep-matrix-out "$WORK/plan_keep_matrix.tsv" -q "$WORK/base" \
  -o "$WORK/matrix_dedup_plan.tsv" > "$WORK/matrix_dedup_plan.log" 2>&1
grep -q "predicted context markerdb after dedup-plan" "$WORK/matrix_dedup_plan.log"
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
