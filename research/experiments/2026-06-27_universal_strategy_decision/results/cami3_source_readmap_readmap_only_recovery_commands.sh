#!/usr/bin/env bash
set -euo pipefail

# Stream-extract only reads_mapping.tsv.gz for CAMI3 samples3-5.
# This avoids storing the full sample read archives or extracting read FASTQ.

mkdir -p /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_3_reads
curl -fL https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_3_reads.tar.gz | tar -xz -C /mnt/new3T/minco_cami3_toygut_extra_20260621 sample_3_reads/reads_mapping.tsv.gz

mkdir -p /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_4_reads
curl -fL https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_4_reads.tar.gz | tar -xz -C /mnt/new3T/minco_cami3_toygut_extra_20260621 sample_4_reads/reads_mapping.tsv.gz

mkdir -p /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_5_reads
curl -fL https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_5_reads.tar.gz | tar -xz -C /mnt/new3T/minco_cami3_toygut_extra_20260621 sample_5_reads/reads_mapping.tsv.gz
