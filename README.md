# minco

`minco` is a pure C fixed-size context-minhash tool for fast genome and read ANI
analysis.

The public sketch keeps the bottom 10,000 context hashes per sample by default,
independent of genome length. Change this with `--sketch-size`. Pairwise reports use shared contexts plus object
differences for ANI-style comparisons, and the `matrix` command provides faster
context-distance output for large screens.

Raw assembled-genome CtxMoE ANI (`ani --raw-output -s3`) uses a minco-trained
MoE calibration against Nayfach ANIm labels for the default 10,000-context
sketch. Default assembled-genome ANI (`ani -s1`, Best) applies a matching HGB
recalibration layer on top of that raw MoE score.

## Build

```bash
make
```

The binary is written to `bin/minco`.

Useful targets:

```bash
make minco
make test
make clean
```

## Quick Start

Sketch assembled genomes:

```bash
bin/minco sketch -p8 --ctxmeta both -o genomes.minco genomes/*.fna.gz
```

Use a larger or smaller fixed sketch:

```bash
bin/minco sketch -p8 --sketch-size 20000 -o genomes.20k.minco genomes/*.fna.gz
```

Sketch from a path list:

```bash
bin/minco sketch -p8 -l genomes.list -o genomes.minco
```

Build the index before large reference or all-vs-all runs:

```bash
bin/minco sketch -i genomes.minco
```

Query ANI against a reference sketch:

```bash
bin/minco ani -p8 -r ref.minco -q query.minco -m0 -o ani.tsv
```

All-vs-all ANI lower triangle:

```bash
bin/minco ani -p8 -q genomes.minco -m2 -s -1 -d -o ani.triangle.tsv
```

Direct FASTA/FASTQ ANI with a non-default sketch size:

```bash
bin/minco ani -S 20000 -f0 -n0 -t0 -o pair.20k.tsv ref.fna query.fna
```

Fast all-vs-all context-distance lower triangle:

```bash
bin/minco matrix --format triangle -q genomes.minco -d -p8 -o ctxdist.triangle.tsv
```

Noisy long-read query ANI:

```bash
bin/minco sketch --conflict --readsQC -p8 -o reads.minco reads.fastq.gz
bin/minco ani -p8 -r ref.minco --qraw reads.minco -m0 -o reads_vs_ref.tsv
```

Inspect a sketch:

```bash
bin/minco sketch --psmp genomes.minco
```

## Help And Manual

```bash
bin/minco --help
bin/minco examples
bin/minco sketch --help
bin/minco ani --help
bin/minco matrix --help
```

The longer user manual is in `docs/USER_MANUAL.md`.

## Metadata

`--ctxmeta` controls the context-cardinality sidecar `minco.ctxmeta.tsv`.
Modes are `preconflict` (default), `postconflict`, `both`, and `none`.

The density-based estimator uses the final bottom-k threshold:

```text
observed_contexts_under_threshold / hash_threshold
```

This sidecar is metadata only. It does not change `comblco` or `comblco.index`.

## Benchmark Helper

```bash
scripts/run_minco_benchmark.sh selected_genomes.tsv outdir 8
```

The helper builds one sketch set, creates the index, and writes indexed
all-vs-all ANI output to `outdir/minco_ani.tsv`.
