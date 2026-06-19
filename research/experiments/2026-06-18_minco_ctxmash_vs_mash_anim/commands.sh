#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ubuntu/yihuiguang/tools/minco
OUT=$ROOT/research/experiments/2026-06-18_minco_ctxmash_vs_mash_anim
SRC=/mnt/new3T/gtdbr220/eval_runs/gtdb_random_genus_vibrio_20260609/vibrio_pair_predictions_with_anim.tsv

cd "$ROOT"
mkdir -p "$OUT"

# Pair selection was deterministic: sort available Vibrio pairs by ANIm and
# take 120 evenly spaced pairs from the source table, yielding 32 unique genomes.
python3 - <<'PY'
import csv, pathlib
src=pathlib.Path('/mnt/new3T/gtdbr220/eval_runs/gtdb_random_genus_vibrio_20260609/vibrio_pair_predictions_with_anim.tsv')
out=pathlib.Path('/home/ubuntu/yihuiguang/tools/minco/research/experiments/2026-06-18_minco_ctxmash_vs_mash_anim')
rows=[]
with src.open() as f:
    r=csv.DictReader(f, delimiter='\t')
    for row in r:
        if row.get('anim_identity_mean') and row.get('path_a') and row.get('path_b'):
            try:
                float(row['anim_identity_mean'])
            except ValueError:
                continue
            rows.append(row)
rows.sort(key=lambda x: float(x['anim_identity_mean']))
if len(rows)>120:
    idx=sorted(set(round(i*(len(rows)-1)/119) for i in range(120)))
    rows=[rows[i] for i in idx]
with (out/'pairs_truth.tsv').open('w', newline='') as g:
    w=csv.writer(g, delimiter='\t')
    w.writerow(['pair_id','sample_a','sample_b','same_species','path_a','path_b','anim'])
    for row in rows:
        w.writerow([row['pair_id'], row['sample_a'], row['sample_b'], row['same_species'], row['path_a'], row['path_b'], row['anim_identity_mean']])
seen=[]; seen_set=set()
for row in rows:
    for p in (row['path_a'], row['path_b']):
        if p not in seen_set:
            seen.append(p); seen_set.add(p)
with (out/'genomes.list').open('w') as g:
    for p in seen:
        g.write(p+'\n')
PY

/usr/bin/time -v ./bin/minco sketch -p8 -S 10000 -l "$OUT/genomes.list" -o "$OUT/genomes.minco"
/usr/bin/time -v ./bin/minco sketch -i "$OUT/genomes.minco"
/usr/bin/time -v mash sketch -p 8 -s 10000 -l -o "$OUT/genomes.mash_s10000" "$OUT/genomes.list"
/usr/bin/time -v mash sketch -p 8 -s 10000 -k 32 -l -o "$OUT/genomes.mash_k32_s10000" "$OUT/genomes.list"

/usr/bin/time -v ./bin/minco ani -p8 -f0 -n0 -t0 -s -1 -m0 -r "$OUT/genomes.minco" -q "$OUT/genomes.minco" -o "$OUT/minco_best_all.tsv"
/usr/bin/time -v ./bin/minco ani -p8 -f0 -n0 -t0 -s -5 -m0 -r "$OUT/genomes.minco" -q "$OUT/genomes.minco" -o "$OUT/minco_ctxmash_all.tsv"
./bin/minco ani -p8 -f0 -n0 -t0 -s -6 -m0 -r "$OUT/genomes.minco" -q "$OUT/genomes.minco" -o "$OUT/minco_aaf_all.tsv"

/usr/bin/time -v -o "$OUT/time_mash_dist.txt" mash dist -p 8 "$OUT/genomes.mash_s10000.msh" "$OUT/genomes.mash_s10000.msh" > "$OUT/mash_dist_all.tsv"
/usr/bin/time -v -o "$OUT/time_mash_k32_dist.txt" mash dist -p 8 "$OUT/genomes.mash_k32_s10000.msh" "$OUT/genomes.mash_k32_s10000.msh" > "$OUT/mash_k32_dist_all.tsv"

# The join/metrics step is recorded by the generated outputs:
# per_pair_predictions.with_mash_k32.tsv and metrics.with_mash_k32.tsv.
