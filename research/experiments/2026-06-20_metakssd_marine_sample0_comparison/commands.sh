#!/usr/bin/env bash
set -euo pipefail

# Fetch OPAL marine report.
curl -fsSL https://yhg926.github.io/KSSD2/OPAL/marine/ \
  -o /tmp/metakssd_opal_marine_index.html

# Fetch MetaKSSD README speed/memory screenshot.
curl -L -f https://github.com/user-attachments/assets/ea17814b-748f-4734-939c-bebf071d5487 \
  -o /tmp/metakssd_speed_memory.png

# Extracted tables were generated with ad hoc Python parsers from the embedded
# Bokeh JSON in /tmp/metakssd_opal_marine_index.html and CAMI profile files.
# See NOTE.md for metric definitions and caveats.

# Sylph local sample0 run.
/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 \
  -r /tmp/cami_marine_sample0_reads.fq.gz \
  -d /tmp/sylph_marine_sample0

/usr/bin/time -v /home/ubuntu/yihuiguang/bin/sylph profile -t 16 \
  /mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb \
  /tmp/sylph_marine_sample0/cami_marine_sample0_reads.fq.gz.sylsp \
  -o /tmp/sylph_marine_sample0/profile.tsv
