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
