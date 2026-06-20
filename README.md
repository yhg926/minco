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
recalibration layer on top of that raw MoE score. `Best` does not silently
switch to containment/Mash-style metrics under low aligned fraction; use
`-s5`/`-s6` explicitly when you want MashD or AafD.

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

Downsample an existing larger bottom-k sketch without re-sketching:

```bash
bin/minco set --downsample -S 1000 -o genomes.1k.minco genomes.20k.minco
bin/minco sketch -i genomes.1k.minco
```

The target sketch size and compact density summary are stored in `minco.stat`.
When `minco ani` uses an existing reference sketch with a direct FASTA/FASTQ
query, the temporary query sketch inherits that reference size automatically.

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

Direct sequence-query ANI at reference density:

```bash
bin/minco ani -r ref.minco -q query.fna.gz --query-density ref -m0 -o ani.tsv
```

Raw FASTQ query ANI can use the same reference-density model without first
building a query sketch:

```bash
bin/minco ani -r ref.minco --qraw reads.fq.gz --query-density ref -m0 -o reads_vs_ref.tsv
```

In this mode minco streams reads, accumulates readwise mutation features, and
uses auxiliary unique context coverage to report AF. Detail output appends
`Reads_with_ctx_match`, `Total_reads`, `Read_match_fraction`,
`Unique_query_ctx`, `Unique_query_ctx_hit`, `Unique_ref_ctx_hit`, and
experimental density-block counters.
With OpenMP builds and `-p > 1`, direct readwise FASTQ density scans process
reads in parallel batches and merge batch-level state as they stream, avoiding a
large end-of-run merge. During long FASTQ runs minco writes progress to stderr:
reads processed, read rate, and input-file percent when the query is a regular
file.
By default minco uses `--density-block-ctx 100`: it accumulates consecutive
reads until at least 100 retained density contexts are available before lookup.
For example, on an `S=1000` reference this tests 100-context pseudo-read blocks
while keeping `Total_reads` as the real read count. Use `--density-block-ctx 0`
to force exact per-read lookup, or set another `N` with
`--density-block-ctx N`.
Add `--abundance-est depth` to also append experimental per-reference
breadth/depth metrics, including `Relative_abundance_depth` and
`Normalized_abundance_depth`, plus `Default_call`. With no custom `-f`, `-n`,
`-t`, or `--top`, minco applies the default readwise abundance report: it uses
a sketch-size-scaled unique-context cutoff (`S/100`), requires ANI >= 0.95,
and prints only `major` and `low_abundance` calls. The normalized column sums
to 1 across printed rows.
Use explicit filters such as `-f0 -n0 -t0` when every candidate comparison must
be reported. This abundance option requires the direct readwise FASTQ density
path, so do not combine it with `--readsQC`, `--abundance`, or
`--save-query-sketch`.

Keep the generated query sketch for debugging:

```bash
bin/minco ani -r ref.minco -q query.fna.gz --query-density ref \
  --save-query-sketch query.debug.minco -m0 -o ani.tsv
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
Release notes are in `CHANGELOG.md`.

## Metadata

New sketches write the `minco.*` / `minco.ctxobj64*` filenames. Readers still
accept earlier KSSD-style filenames for compatibility with existing sketch
directories.

`--ctxmeta` controls context-cardinality metadata. It writes compact binary
per-sample records in `minco.ctxmeta` and stores the operational set-level
density summary in `minco.stat`. Modes are `preconflict` (default),
`postconflict`, `both`, and `none`.

The unique-context estimator uses the final bottom-k threshold:

```text
observed_contexts_under_threshold * hash_space / (hash_threshold + 1)
```

This metadata is optional. It does not change `minco.ctxobj64` or
`minco.ctxobj64.offsets`. The optional search index is
`minco.refindex.ctxgid64obj32`. `minco.ctxmeta` stores per-sample density
estimates as fixed-width records without text keys or headers. Use
`minco sketch --pctxmeta DIR` and `minco sketch --pctxsetmeta DIR` to print TSV
views when needed.
`minco ani --query-density ref` uses this metadata only for FASTA/FASTQ query
inputs: one-sample references use that sample's density threshold, while
combined references use the largest sample density threshold. The query
extraction is internal to `ani`; no persistent query sketch is written unless
`--save-query-sketch DIR` is used for debugging. The save target must not
already exist.
`minco ani` uses `minco.ctxmeta` when available to estimate asymmetric real
aligned fractions from the common bottom-k hash threshold. Detail output appends
`Real_Qry_align_fraction`, `Real_Ref_align_fraction`,
`Real_min_align_fraction`, and `AF_source`; without ctxmeta these fall back to
the fixed-sketch aligned fractions.

Sketch maintenance commands preserve this metadata when possible:
`minco sketch --keep` and `minco sketch --remove` filter `minco.ctxmeta`
alongside the sketch entries and refresh the density summary in `minco.stat`.
Filtered reference sketches can therefore still be used directly with
`minco ani --query-density ref`.

## Benchmark Helper

```bash
scripts/run_minco_benchmark.sh selected_genomes.tsv outdir 8
```

The helper builds one sketch set, creates the index, and writes indexed
all-vs-all ANI output to `outdir/minco_ani.tsv`.
