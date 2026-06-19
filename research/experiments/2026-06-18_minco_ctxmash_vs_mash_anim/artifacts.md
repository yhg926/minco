# Artifacts

Experiment directory:

```text
/home/ubuntu/yihuiguang/tools/minco/research/experiments/2026-06-18_minco_ctxmash_vs_mash_anim
```

Inputs:

```text
/mnt/new3T/gtdbr220/eval_runs/gtdb_random_genus_vibrio_20260609/vibrio_pair_predictions_with_anim.tsv
```

Selected input tables:

```text
pairs_truth.tsv
genomes.list
```

Sketches and pairwise outputs:

```text
genomes.minco/
genomes.mash_s10000.msh
genomes.mash_k32_s10000.msh
minco_best_all.tsv
minco_ctxmash_all.tsv
minco_aaf_all.tsv
mash_dist_all.tsv
mash_k32_dist_all.tsv
```

Joined outputs:

```text
per_pair_predictions.tsv
per_pair_predictions.with_mash_k32.tsv
metrics.tsv
metrics.with_mash_k32.tsv
summary.tsv
```

Timing notes captured in this run:

```text
Minco sketch: 32 genomes, -S 10000, -p8, 0.26 s wall, 92.5 MB max RSS
Minco index: 0.01 s wall, 8.1 MB max RSS
Minco all-vs-all ani -s -1: 0.01 s wall, 10.4 MB max RSS
Minco all-vs-all ani -s -5: 0.02 s wall, 10.4 MB max RSS
Mash sketch -s 10000: 0.73 s wall, 70.9 MB max RSS
Mash sketch -s 10000 -k 32: 0.72 s wall, 70.4 MB max RSS
Mash dist -s 10000: see time_mash_dist.txt
Mash dist -s 10000 -k 32: see time_mash_k32_dist.txt
```
