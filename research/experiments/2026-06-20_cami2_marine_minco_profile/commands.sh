#!/usr/bin/env bash
set -euo pipefail

REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620
TAX=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.tsv
URL=https://frl.publisso.de/data/frl%3A6425521/marine/short_read/marmgCAMI2_sample_0_reads.tar.gz
MEMBER=simulation_short_read/2018.08.15_09.49.32_sample_0/reads/anonymous_reads.fq.gz
PROFILE=/tmp/minco_cami_marine_sample0_s1000.profileonly.profile
TSV=/tmp/minco_cami_marine_sample0_s1000.profileonly.tsv
LOG=/tmp/minco_cami_marine_sample0_s1000.profileonly.log

curl -fsSL https://cami-challenge.org/static/examples/gs_marine_short.profile \
  -o /tmp/gs_marine_short.profile

/usr/bin/time -v bash -c \
  'set -o pipefail; curl -fsSL "$0" | tar -xzO "$1" | ./bin/minco ani -p16 -r "$2" --qraw - --query-density ref --abundance-est depth --readwise-profile-only --cami-taxmap "$3" --cami-profile "$4" --cami-sample-id marmgCAMI2_short_read_sample_0 -m0 -o "$5"' \
  "$URL" "$MEMBER" "$REF" "$TAX" "$PROFILE" "$TSV" > "$LOG" 2>&1

# Full permissive candidate table used for S=1000 false-negative diagnosis and
# threshold optimization. This intentionally disables AF, ANI, and context
# reporting filters so post-hoc threshold grids can be tested.
UNFILTERED_TSV=/tmp/minco_cami_marine_sample0_s1000.profileonly.unfiltered.tsv
UNFILTERED_LOG=/tmp/minco_cami_marine_sample0_s1000.profileonly.unfiltered.log

/usr/bin/time -v bash -c \
  'set -o pipefail; curl -fsSL "$0" | tar -xzO "$1" | ./bin/minco ani -p16 -r "$2" --qraw - --query-density ref --abundance-est depth --readwise-profile-only -m0 -f0 -n0 -t0 -o "$3"' \
  "$URL" "$MEMBER" "$REF" "$UNFILTERED_TSV" > "$UNFILTERED_LOG" 2>&1

# S=10000 GTDB-only reference. The streaming form timed out because minco spent
# ~3 minutes loading the 23 GB refindex before consuming stdin, so cache only
# the inner FASTQ gzip first.
REF_S10000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
FASTQ_CACHE=/tmp/cami_marine_sample0_reads.fq.gz
PROFILE_S10000=/tmp/minco_cami_marine_sample0_s10000_gtdb.profileonly.profile
TSV_S10000=/tmp/minco_cami_marine_sample0_s10000_gtdb.profileonly.tsv
LOG_S10000=/tmp/minco_cami_marine_sample0_s10000_gtdb.profileonly.log
UNFILTERED_TSV_S10000=/tmp/minco_cami_marine_sample0_s10000_gtdb.profileonly.unfiltered.tsv
UNFILTERED_LOG_S10000=/tmp/minco_cami_marine_sample0_s10000_gtdb.profileonly.unfiltered.log

curl --retry 8 --retry-all-errors --retry-delay 5 -fsSL "$URL" \
  | tar -xzO "$MEMBER" > "$FASTQ_CACHE"

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S10000" \
  --qraw "$FASTQ_CACHE" --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap "$TAX" --cami-profile "$PROFILE_S10000" \
  --cami-sample-id marmgCAMI2_short_read_sample_0 -m0 -o "$TSV_S10000" \
  > "$LOG_S10000" 2>&1

/usr/bin/time -v ./bin/minco ani -p16 -r "$REF_S10000" \
  --qraw "$FASTQ_CACHE" --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 -o "$UNFILTERED_TSV_S10000" \
  > "$UNFILTERED_LOG_S10000" 2>&1

# CAMI upload/evaluation used browser-style curl with CSRF token/cookie.
# Uploaded submission: https://cami-challenge.org/submission/9089740f62614e468422/
# OPAL metric table: https://cami-challenge.org/opal_run/471/
