#!/usr/bin/env bash
set -euo pipefail

REPO=/home/ubuntu/yihuiguang/tools/KSSD3mini
WORK=/tmp/minco_readwise_tracking_toymouse0_20260624
REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker
READS=/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
TAXMAP="$WORK/gtdb_r232_accession.taxmap.tsv"

cd "$REPO"
mkdir -p "$WORK"

python3 - <<'PY'
import csv, gzip, re
from pathlib import Path
work = Path("/tmp/minco_readwise_tracking_toymouse0_20260624")
out = work / "gtdb_r232_accession.taxmap.tsv"
meta = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]
acc_re = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
seen = set()
with out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["ref_key", "TAXID", "RANK", "TAXPATH", "TAXPATHSN"])
    for path in meta:
        with gzip.open(path, "rt") as gz:
            for row in csv.DictReader(gz, delimiter="\t"):
                tax = row.get("gtdb_taxonomy", "")
                if not tax or "s__" not in tax:
                    continue
                for field in ("accession", "ncbi_genbank_assembly_accession"):
                    match = acc_re.search(row.get(field, "") or "")
                    if not match:
                        continue
                    key = match.group(1)
                    if key in seen:
                        continue
                    seen.add(key)
                    taxpath = tax.replace(";", "|")
                    writer.writerow([key, key, "species", taxpath, taxpath])
PY

/usr/bin/time -v -o "$WORK/baseline_default.time.log" \
  bin/minco ani -p16 -r "$REF" --qraw "$READS" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    -m0 -f0 -n0 -t0 -o "$WORK/baseline_default.tsv" \
  > "$WORK/baseline_default.stdout.log" 2> "$WORK/baseline_default.stderr.log"

/usr/bin/time -v -o "$WORK/baseline_exact.time.log" \
  bin/minco ani -p16 -r "$REF" --qraw "$READS" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --density-block-ctx 0 \
    -m0 -f0 -n0 -t0 -o "$WORK/baseline_exact.tsv" \
  > "$WORK/baseline_exact.stdout.log" 2> "$WORK/baseline_exact.stderr.log"

/usr/bin/time -v -o "$WORK/track_none.time.log" \
  bin/minco ani -p16 -r "$REF" --qraw "$READS" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --readwise-track "$WORK/track_none.tsv" \
    --readwise-track-summary "$WORK/track_none.summary.tsv" \
    -m0 -f0 -n0 -t0 -o "$WORK/track_none.profile.tsv" \
  > "$WORK/track_none.stdout.log" 2> "$WORK/track_none.stderr.log"

/usr/bin/time -v -o "$WORK/track_gtdb.time.log" \
  bin/minco ani -p16 -r "$REF" --qraw "$READS" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --readwise-track "$WORK/track_gtdb.tsv" \
    --readwise-track-summary "$WORK/track_gtdb.summary.tsv" \
    --readwise-taxonomy gtdb --gtdb-taxmap "$TAXMAP" \
    -m0 -f0 -n0 -t0 -o "$WORK/track_gtdb.profile.tsv" \
  > "$WORK/track_gtdb.stdout.log" 2> "$WORK/track_gtdb.stderr.log"

/usr/bin/time -v -o "$WORK/sample_reads_0p003.time.log" \
  /home/ubuntu/miniconda3/bin/seqkit sample -p 0.003 -s 20260624 \
    "$READS" -o "$WORK/mouse0_reads_sample_p0003.fq.gz" \
  > "$WORK/sample_reads_0p003.stdout.log" \
  2> "$WORK/sample_reads_0p003.stderr.log"

awk -F '\t' 'NR > 1 {print $NF}' \
  research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/mouse0_gtdb_source_species_ground_truth.tsv \
  > "$WORK/whole_ctx_absence_sample_p0003.source_paths.txt"

bin/minco sketch --psmp "$REF" > "$WORK/full_s2000_psmp.tsv"

python3 - <<'PY'
import csv, gzip, re
from pathlib import Path
truth = Path("research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/mouse0_gtdb_source_species_ground_truth.tsv")
psmp = Path("/tmp/minco_readwise_tracking_toymouse0_20260624/full_s2000_psmp.tsv")
out = Path("/tmp/minco_readwise_tracking_toymouse0_20260624/gtdb_source_species_rep_paths.txt")
meta_paths = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]
acc_re = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
def acc(s):
    m = acc_re.search(s or "")
    return m.group(1) if m else ""
def core(a):
    a = acc(a)
    return a.split("_", 1)[1] if a else ""
by_acc = {}
by_core = {}
for mp in meta_paths:
    with gzip.open(mp, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rec = dict(row)
            for k in [acc(row.get("accession", "")), acc(row.get("ncbi_genbank_assembly_accession", ""))]:
                if k:
                    by_acc.setdefault(k, rec)
                    by_core.setdefault(core(k), []).append(rec)
def lookup(a):
    a = acc(a)
    if a in by_acc:
        return by_acc[a]
    vals = by_core.get(core(a), [])
    uniq = []
    seen = set()
    for r in vals:
        key = (r.get("accession", ""), r.get("gtdb_taxonomy", ""))
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    return uniq[0] if len(uniq) == 1 else None
path_by_acc = {}
with psmp.open() as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t", 1)
        if len(parts) == 2 and acc(parts[1]):
            path_by_acc.setdefault(acc(parts[1]), parts[1])
paths = []
seen = set()
with truth.open() as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rec = lookup(row.get("source_accession", ""))
        rep = acc(rec.get("gtdb_genome_representative", "")) if rec else ""
        path = path_by_acc.get(rep, "")
        if path and path not in seen:
            seen.add(path)
            paths.append(path)
with out.open("w") as fh:
    for path in paths:
        fh.write(path + "\n")
PY

SOURCE_EXACT="$WORK/source_exact_minco"
GTDB_REP_EXACT="$WORK/gtdb_rep_exact_minco"
mkdir -p "$SOURCE_EXACT" "$GTDB_REP_EXACT"

/usr/bin/time -v -o "$SOURCE_EXACT/sketch_source_exact.time.log" \
  bin/minco sketch -p16 --conflict -S 10000000 \
    -l "$WORK/whole_ctx_absence_sample_p0003.source_paths.txt" \
    -o "$SOURCE_EXACT/source_exact.minco" \
  > "$SOURCE_EXACT/sketch_source_exact.stdout.log" \
  2> "$SOURCE_EXACT/sketch_source_exact.stderr.log"

/usr/bin/time -v -o "$SOURCE_EXACT/index_source_exact.time.log" \
  bin/minco sketch -i "$SOURCE_EXACT/source_exact.minco" \
  > "$SOURCE_EXACT/index_source_exact.stdout.log" \
  2> "$SOURCE_EXACT/index_source_exact.stderr.log"

/usr/bin/time -v -o "$SOURCE_EXACT/sample_reads_vs_source_exact.time.log" \
  bin/minco ani -p16 -r "$SOURCE_EXACT/source_exact.minco" \
    --qraw "$WORK/mouse0_reads_sample_p0003.fq.gz" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --readwise-track /dev/null \
    --readwise-track-summary "$SOURCE_EXACT/sample_reads_vs_source_exact.summary.tsv" \
    -m0 -f0 -n0 -t0 -o "$SOURCE_EXACT/sample_reads_vs_source_exact.profile.tsv" \
  > "$SOURCE_EXACT/sample_reads_vs_source_exact.stdout.log" \
  2> "$SOURCE_EXACT/sample_reads_vs_source_exact.stderr.log"

/usr/bin/time -v -o "$GTDB_REP_EXACT/sketch_gtdb_rep_exact.time.log" \
  bin/minco sketch -p16 --conflict -S 10000000 \
    -l "$WORK/gtdb_source_species_rep_paths.txt" \
    -o "$GTDB_REP_EXACT/gtdb_rep_exact.minco" \
  > "$GTDB_REP_EXACT/sketch_gtdb_rep_exact.stdout.log" \
  2> "$GTDB_REP_EXACT/sketch_gtdb_rep_exact.stderr.log"

/usr/bin/time -v -o "$GTDB_REP_EXACT/index_gtdb_rep_exact.time.log" \
  bin/minco sketch -i "$GTDB_REP_EXACT/gtdb_rep_exact.minco" \
  > "$GTDB_REP_EXACT/index_gtdb_rep_exact.stdout.log" \
  2> "$GTDB_REP_EXACT/index_gtdb_rep_exact.stderr.log"

/usr/bin/time -v -o "$GTDB_REP_EXACT/sample_reads_vs_gtdb_rep_exact.time.log" \
  bin/minco ani -p16 -r "$GTDB_REP_EXACT/gtdb_rep_exact.minco" \
    --qraw "$WORK/mouse0_reads_sample_p0003.fq.gz" --query-density ref \
    --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --readwise-track /dev/null \
    --readwise-track-summary "$GTDB_REP_EXACT/sample_reads_vs_gtdb_rep_exact.summary.tsv" \
    -m0 -f0 -n0 -t0 -o "$GTDB_REP_EXACT/sample_reads_vs_gtdb_rep_exact.profile.tsv" \
  > "$GTDB_REP_EXACT/sample_reads_vs_gtdb_rep_exact.stdout.log" \
  2> "$GTDB_REP_EXACT/sample_reads_vs_gtdb_rep_exact.stderr.log"

FULL_S2000_REF=/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup

/usr/bin/time -v -o "$WORK/full_s2000_track_corrected.time.log" \
  bin/minco ani -p16 -r "$FULL_S2000_REF" --qraw "$READS" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-split --readwise-ani zip-aaf \
    --readwise-track /dev/null \
    --readwise-track-summary "$WORK/full_s2000_track_corrected.summary.tsv" \
    -m0 -f0 -n0 -t0 \
    -o "$WORK/full_s2000_track_corrected.profile.tsv"
