#!/usr/bin/env bash
set -euo pipefail

SYLPH=/home/ubuntu/yihuiguang/bin/sylph
REP_LIST=/mnt/new3T/gtdbr220/gtdb232/manifests/r232.reps_fullpath.list
CHUNK_ROOT=/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks
CHUNK_LIST_DIR=/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/lists
CHUNK_DB_DIR=/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/db
CHUNK_SIZE=10000
THREADS=16
C_VALUE=200
K_VALUE=31

mkdir -p "$CHUNK_LIST_DIR" "$CHUNK_DB_DIR"
if ! compgen -G "$CHUNK_LIST_DIR/chunk_*.list" > /dev/null; then
  split -l "$CHUNK_SIZE" -d -a 4 --additional-suffix=.list "$REP_LIST" "$CHUNK_LIST_DIR/chunk_"
fi

for list in "$CHUNK_LIST_DIR"/chunk_*.list; do
  name=$(basename "$list" .list)
  out="$CHUNK_DB_DIR/$name"
  if [[ -s "$out.syldb" ]]; then
    continue
  fi
  /usr/bin/time -v -o "$CHUNK_DB_DIR/$name.build.time.log" \
    "$SYLPH" sketch -t "$THREADS" -c "$C_VALUE" -k "$K_VALUE" --gl "$list" -o "$out"
done

find "$CHUNK_DB_DIR" -maxdepth 1 -name 'chunk_*.syldb' | sort > "$CHUNK_ROOT/gtdb-r232-c200-dbv1.chunk_syldb.list"
wc -l "$CHUNK_ROOT/gtdb-r232-c200-dbv1.chunk_syldb.list"
