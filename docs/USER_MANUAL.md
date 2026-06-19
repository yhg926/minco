# minco User Manual

`minco` builds fixed-size context-minhash sketches and estimates genome or read
ANI from those sketches. It is a pure C command line tool.

The public command set is:

```text
minco sketch   create, inspect, index, filter, append, and deduplicate sketches
minco ani      estimate ANI from sketches or direct FASTA/FASTQ inputs
minco matrix   report context-distance matrices, sparse edges, clusters, plans
minco examples print common workflows
minco doctor   check basic build/runtime information
```

Use `minco <subcommand> --help` for the exact option list.

## Core Defaults

The default public sketch is a bottom-k context-minhash sketch:

```text
context/object pattern: coden pattern
sketch size:            10,000 context hashes per sample
ctxmeta:                preconflict
threads:                1 unless -p/--threads is set
```

Change the fixed sketch size with `-S/--sketch-size`:

```bash
minco sketch --sketch-size 20000 -o genomes.20k.minco genomes/*.fna.gz
```

The final sketch keeps the bottom `--sketch-size` context hashes per sample.

## Sketch Directory Layout

A sketch output is a directory. Common files are:

```text
comblco                 sketch entries
lcofiles.stat           sketch parameters and sample table
comblco.index           optional inverted index built by `minco sketch -i`
minco.ctxmeta.tsv       optional density and unique-context estimates
lcofiles.infilemeta     optional input metadata
lcofiles.qc             optional readsQC count ranges
lcofiles.anno           optional FASTA/FASTQ header annotations
comblco.a               optional per-entry abundance/counts
comblco.position        optional source positions
```

Build `comblco.index` before large reference searches, all-vs-all ANI, sparse
matrix reports, clusters, or dedup plans:

```bash
minco sketch -i genomes.minco
```

## Sketch Assembled Genomes

Sketch one or more FASTA files:

```bash
minco sketch -p8 --ctxmeta both -o genomes.minco genomes/*.fna.gz
```

Sketch with a non-default size:

```bash
minco sketch -p8 --sketch-size 20000 -o genomes.20k.minco genomes/*.fna.gz
```

Sketch paths from a list:

```bash
minco sketch -p8 -l genomes.list -o genomes.minco
```

Treat all input files as one final sample:

```bash
minco sketch --asone -o one_sample.minco contigs_part1.fna contigs_part2.fna
```

Split a multi-FASTA so each record becomes one sample:

```bash
minco sketch --splitmfa -o records.minco assembly_records.fna
```

## Sketch Noisy Reads

For Nanopore or other noisy long-read FASTQ, keep conflicting context objects
and apply read-count QC before taking the final bottom-k sketch:

```bash
minco sketch --conflict --readsQC -p8 -o reads.minco reads.fastq.gz
```

For a two-step QC workflow, keep abundance first, then apply stored QC ranges:

```bash
minco sketch -A --conflict --readsQC -p8 -o reads_abundance.minco reads.fastq.gz
minco sketch --sketchQC -o reads_qc.minco reads_abundance.minco
```

Stream input from another tool:

```bash
samtools fastq reads.bam | minco sketch --conflict --readsQC -o reads.minco -
```

Use `--pipecmd` when each path needs conversion:

```bash
minco sketch --pipecmd 'samtools fastq {}' --conflict --readsQC \
  -o reads.minco reads.bam
```

## Inspect Sketches

List samples and sketch sizes:

```bash
minco sketch --psmp genomes.minco
```

Print sketch entries or an index:

```bash
minco sketch --psketch genomes.minco
minco sketch --pindex genomes.minco
```

Write and inspect source positions:

```bash
minco sketch --position -o positioned.minco genomes/*.fna.gz
minco sketch --ppos positioned.minco
```

## ANI From Sketches

Reference/query ANI:

```bash
minco sketch -i ref.minco
minco ani -p8 -r ref.minco -q query.minco -m0 -o ani.tsv
```

All-vs-all ANI lower triangle:

```bash
minco sketch -i genomes.minco
minco ani -p8 -q genomes.minco -m2 -s -1 -d -o ani.triangle.tsv
```

`ani -m2` computes ANI-style comparison fields. For faster large screening when
context distance is enough, use `matrix --format triangle`.

## ANI From Direct FASTA/FASTQ

Direct sequence inputs are auto-sketched into temporary sketch directories.
When any input is already a sketch, direct sequence inputs reuse that sketch's
parameters. Otherwise they use the default public sketch settings.

Pairwise direct input:

```bash
minco ani -f0 -n0 -t0 -o pair.tsv ref.fasta query.fasta
```

Pairwise direct input with a non-default auto-sketch size:

```bash
minco ani -S 20000 -f0 -n0 -t0 -o pair.20k.tsv ref.fasta query.fasta
```

Reference sketch against one direct query:

```bash
minco ani -p8 -r ref.minco -q query.fasta.gz -m0 -o ani.tsv
```

When `-r` is an existing sketch and `-q` or `--qraw` is a FASTA/FASTQ file,
`-S/--sketch-size` controls the temporary query sketch. Use the same size as
the reference sketch if you want symmetric sketch sizes.

Current sketch metadata does not store the requested sketch size separately.
`minco sketch --psmp DIR` shows the actual entry count for each sample.

Path-list mode:

```bash
minco ani --reflist refs.txt --qrylist qrys.txt -p8 -o ani.tsv
```

## Raw-Read ANI

Use `--qraw` for reads or read sketches. For sequence inputs, query conflicts
are kept automatically. For noisy long reads, add `--readsQC`.

Reads against assembled references:

```bash
minco sketch -i ref.minco
minco ani -p8 -r ref.minco --qraw reads.fastq.gz --readsQC -m0 -o reads_vs_ref.tsv
```

Pre-sketched reads against references:

```bash
minco sketch --conflict --readsQC -p8 -o reads.minco reads.fastq.gz
minco ani -p8 -r ref.minco --qraw reads.minco -m0 -o reads_vs_ref.tsv
```

Reads against reads:

```bash
minco sketch --conflict --readsQC -p8 -o ref_reads.minco ref_reads/*.fastq.gz
minco sketch --conflict --readsQC -p8 -o qry_reads.minco qry_reads/*.fastq.gz
minco sketch -i ref_reads.minco
minco ani -p8 -r ref_reads.minco --qraw qry_reads.minco -m0 -o reads_to_reads.tsv
```

## ANI Filters And Metrics

Default detail ANI filters are:

```text
-n/--anicut  0.95
-f/--afcut   0.5, or 0.2 for --qraw
-t/--ctxcut  3 shared contexts
```

Disable these filters when you need every pair:

```bash
minco ani -f0 -n0 -t0 -r ref.minco -q query.minco -o all_pairs.tsv
```

Selected metric values for `-s/--slmetrics`:

```text
1  Best
2  Recalibrated
3  CtxMoE
4  Naive
5  MashD
6  AafD
7  MashD_if_far
8  AafD_if_far
```

In `-m1` and `-m2` formats, positive `-s` values print distance and negative
`-s` values print ANI. Detail output prints both.

For normal assembled-genome ANI, `CtxMoE` uses a minco-trained MoE calibration
against Nayfach ANIm labels for the default 10,000-context sketch. `Best` and
`Recalibrated` apply a matching HGB calibration layer on top of that raw MoE
score. Use `--raw-output -s3` to inspect the raw CtxMoE value.

## Context-Distance Matrix Reports

`matrix` reports context-distance outputs and is intended for fast screening,
edge lists, clustering, and dedup planning.

Self lower triangle:

```bash
minco sketch -i genomes.minco
minco matrix --format triangle -q genomes.minco -d -p8 -o ctxdist.triangle.tsv
```

Rectangular full matrix:

```bash
minco matrix -r ref.minco -q query.minco --format full -p8 -o matrix.tsv
```

Sparse edges:

```bash
minco matrix --sparse --cut 0.05 genomes.minco > sparse_edges.tsv
```

Clusters:

```bash
minco matrix --format clusters --cut 0.05 genomes.minco > clusters.tsv
```

Dedup plan:

```bash
minco matrix --format dedup-plan --cut 0.001 --keep-out keep.txt \
  --remove-out remove.txt --keep-matrix-out keep_matrix.tsv genomes.minco
```

Default metric is `ctx-naive`. Other accepted metrics include `ctx-moe`, `mash`,
`aaf`, and expressions such as `ctx-naive&aaf` or `ctx-moe|mash` for sparse
candidate filtering.

## Sketch Maintenance

Append sketches:

```bash
minco sketch --append -o merged.minco base.minco add.minco
```

Keep or remove listed samples:

```bash
minco sketch --keep keep_names.txt -o kept.minco genomes.minco
minco sketch --remove remove_names.txt -o filtered.minco genomes.minco
```

Deduplicate:

```bash
minco sketch --dedup 0.001 --metric ctx-moe -o dedup.minco genomes.minco
minco sketch --dedup 0.001 --dedup-index -o dedup_indexed.minco genomes.minco
```

Use `--drop-position` with keep/remove/dedup output if the position sidecar is
not needed in the result.

## Context Metadata

`--ctxmeta` writes `minco.ctxmeta.tsv`.

Modes:

```text
none          do not write context metadata
preconflict   estimate before conflict resolution, default
postconflict  estimate after conflict resolution
both          write both estimates
```

The density estimator uses the final bottom-k threshold:

```text
observed_contexts_under_threshold / hash_threshold
```

The sidecar helps estimate total unique contexts but does not change sketch
entries, ANI output, or the index.

## Choosing `ani` Or `matrix`

Use `ani` when the output must contain ANI, selected metric, confidence, context
counts, and alignment fractions.

Use `matrix` when the task is large all-vs-all screening, sparse candidate
edges, clusters, or dedup planning and a context-distance metric is sufficient.

For a large indexed genome set, this is the fast screening pattern:

```bash
minco sketch -p8 -l genomes.list -o genomes.minco
minco sketch -i genomes.minco
minco matrix --format triangle -q genomes.minco -d -p8 -o ctxdist.triangle.tsv
```
