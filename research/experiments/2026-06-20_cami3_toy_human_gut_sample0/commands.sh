#!/usr/bin/env bash
set -euo pipefail

WORK=/mnt/new3T/minco_cami3_toygut_20260620
REF=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
GOLD=$WORK/taxonomic_profiles/taxonomic_profile_0.txt

mkdir -p "$WORK"

# Disk check before large files.
df -h /tmp /mnt/new3T

# Small truth archive.
curl --retry 5 --retry-all-errors --retry-delay 5 -fL \
  https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/taxonomic_profiles.tar.gz \
  -o "$WORK/taxonomic_profiles.tar.gz"
tar -xzf "$WORK/taxonomic_profiles.tar.gz" -C "$WORK" taxonomic_profiles/taxonomic_profile_0.txt

# Sample0 reads, archived. Object size observed with HEAD: 4,866,106,604 bytes.
/usr/bin/time -v curl --retry 5 --retry-all-errors --retry-delay 5 -fL \
  https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_0_reads.tar.gz \
  -o "$WORK/sample_0_reads.tar.gz"

tar -tzf "$WORK/sample_0_reads.tar.gz" | sed -n '1,40p'

# Current marine-tuned minco recipe.
/usr/bin/time -v bash -lc \
  'set -o pipefail; tar -xOzf "$0" sample_0_reads/anonymous_reads.fq.gz | ./minco_core/bin/minco ani -p16 -r "$1" --qraw - --query-density ref --abundance-est depth --readwise-profile-only --readwise-assign best-diff-unique --readwise-ani zip-aaf -m0 -f0.05 -n0.94 -t10 -o "$2"' \
  "$WORK/sample_0_reads.tar.gz" "$REF" "$WORK/minco_s1000_plusvirus_sample0_f0.05_n0.94_t10.tsv"

# Unfiltered minco candidates for threshold diagnosis.
/usr/bin/time -v bash -lc \
  'set -o pipefail; tar -xOzf "$0" sample_0_reads/anonymous_reads.fq.gz | ./minco_core/bin/minco ani -p16 -r "$1" --qraw - --query-density ref --abundance-est depth --readwise-profile-only --readwise-assign best-diff-unique --readwise-ani zip-aaf -m0 -f0 -n0 -t0 -o "$2"' \
  "$WORK/sample_0_reads.tar.gz" "$REF" "$WORK/minco_s1000_plusvirus_sample0_unfiltered.tsv"

# Sylph needs a physical FASTQ gzip.
/usr/bin/time -v tar -xOzf "$WORK/sample_0_reads.tar.gz" sample_0_reads/anonymous_reads.fq.gz \
  > "$WORK/sample_0_anonymous_reads.fq.gz"
mkdir -p "$WORK/sylph_sample0"
/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 \
  -r "$WORK/sample_0_anonymous_reads.fq.gz" -d "$WORK/sylph_sample0"
/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph profile -t 16 \
  /mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb \
  "$WORK/sylph_sample0/sample_0_anonymous_reads.fq.gz.sylsp" \
  -o "$WORK/sylph_sample0/profile.tsv"

# Scoring commands used local Python helpers and ad hoc scope/group scoring;
# exact score outputs are in /tmp/minco_cami3_toygut_20260620.
