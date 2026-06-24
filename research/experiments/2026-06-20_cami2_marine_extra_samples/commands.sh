#!/usr/bin/env bash
set -euo pipefail

MINCO=./minco_core/bin/minco
REF_S10000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
GOLD=/tmp/gs_marine_short.profile
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
SYLPH=/home/ubuntu/yihuiguang/bin/sylph
SYLPH_DB=/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
OUT=/tmp/minco_readwise_assign_20260620

# Cache sample1 and sample2 inner FASTQ gzip files from the CAMI II marine short-read archives.
/usr/bin/time -v bash -c 'set -o pipefail; curl --retry 8 --retry-all-errors --retry-delay 5 -fsSL https://frl.publisso.de/data/frl%3A6425521/marine/short_read/marmgCAMI2_sample_1_reads.tar.gz | tar -xzO simulation_short_read/2018.08.15_09.49.32_sample_1/reads/anonymous_reads.fq.gz > /tmp/cami_marine_sample1_reads.fq.gz'
/usr/bin/time -v bash -c 'set -o pipefail; curl --retry 8 --retry-all-errors --retry-delay 5 -fsSL https://frl.publisso.de/data/frl%3A6425521/marine/short_read/marmgCAMI2_sample_2_reads.tar.gz | tar -xzO simulation_short_read/2018.08.15_09.49.32_sample_2/reads/anonymous_reads.fq.gz > /tmp/cami_marine_sample2_reads.fq.gz'

# Fixed minco S10000 recipe tested on samples 1 and 2.
/usr/bin/time -v "$MINCO" ani -p16 -r "$REF_S10000" \
  --qraw /tmp/cami_marine_sample1_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o "$OUT/s10000_unique_zipaaf_sample1_f0.05_n0.94_t10.tsv"

/usr/bin/time -v "$MINCO" ani -p16 -r "$REF_S10000" \
  --qraw /tmp/cami_marine_sample2_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o "$OUT/s10000_unique_zipaaf_sample2_f0.05_n0.94_t10.tsv"

# Local Sylph baseline for samples 1 and 2.
mkdir -p /tmp/sylph_marine_sample1 /tmp/sylph_marine_sample2
/usr/bin/time -v "$SYLPH" sketch -t 16 -r /tmp/cami_marine_sample1_reads.fq.gz -d /tmp/sylph_marine_sample1
/usr/bin/time -v "$SYLPH" profile -t 16 "$SYLPH_DB" /tmp/sylph_marine_sample1/cami_marine_sample1_reads.fq.gz.sylsp -o /tmp/sylph_marine_sample1/profile.tsv
/usr/bin/time -v "$SYLPH" sketch -t 16 -r /tmp/cami_marine_sample2_reads.fq.gz -d /tmp/sylph_marine_sample2
/usr/bin/time -v "$SYLPH" profile -t 16 "$SYLPH_DB" /tmp/sylph_marine_sample2/cami_marine_sample2_reads.fq.gz.sylsp -o /tmp/sylph_marine_sample2/profile.tsv

# Score sample0, sample1, and sample2 using the same local species-taxid scorer.
python3 research/experiments/2026-06-20_cami2_marine_extra_samples/scripts/score_profiles.py \
  --sample-id marmgCAMI2_short_read_sample_0 \
  --minco "$OUT/s10000_unique_zipaaf_f0.05_n0.94_t10.tsv" \
  --sylph /tmp/sylph_marine_sample0/profile.tsv \
  --gold "$GOLD" --taxmap "$TAXMAP" \
  --out "$OUT/sample0_minco_sylph_score.tsv"

python3 research/experiments/2026-06-20_cami2_marine_extra_samples/scripts/score_profiles.py \
  --sample-id marmgCAMI2_short_read_sample_1 \
  --minco "$OUT/s10000_unique_zipaaf_sample1_f0.05_n0.94_t10.tsv" \
  --sylph /tmp/sylph_marine_sample1/profile.tsv \
  --gold "$GOLD" --taxmap "$TAXMAP" \
  --out "$OUT/sample1_minco_sylph_score.tsv"

python3 research/experiments/2026-06-20_cami2_marine_extra_samples/scripts/score_profiles.py \
  --sample-id marmgCAMI2_short_read_sample_2 \
  --minco "$OUT/s10000_unique_zipaaf_sample2_f0.05_n0.94_t10.tsv" \
  --sylph /tmp/sylph_marine_sample2/profile.tsv \
  --gold "$GOLD" --taxmap "$TAXMAP" \
  --out "$OUT/sample2_minco_sylph_score.tsv"
