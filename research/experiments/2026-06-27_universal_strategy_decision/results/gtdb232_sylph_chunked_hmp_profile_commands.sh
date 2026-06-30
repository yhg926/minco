#!/usr/bin/env bash
set -euo pipefail

test -s /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list

mkdir -p /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample6 /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs
cat /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample6.profile_inputs.list
echo /tmp/cami2_hmp_unseen_transfer_20260626/run/sylph_sample6/airskinurogenital_sample6.nonzero.fastq.gz.sylsp >> /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample6.profile_inputs.list
/usr/bin/time -v -o /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample6/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample6.profile_inputs.list -o /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample6/profile.chunked.tsv

mkdir -p /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample11 /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs
cat /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample11.profile_inputs.list
echo /tmp/cami2_hmp_unseen_transfer_20260626/run/sylph_sample11/airskinurogenital_sample11.nonzero.fastq.gz.sylsp >> /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample11.profile_inputs.list
/usr/bin/time -v -o /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample11/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/profile_inputs/sample11.profile_inputs.list -o /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample11/profile.chunked.tsv
