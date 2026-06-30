#!/usr/bin/env bash
set -euo pipefail

/usr/bin/time -v -o /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1.build.time.log /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 -c 200 -k 31 --gl /mnt/new3T/gtdbr220/gtdb232/manifests/r232.reps_fullpath.list -o /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1
