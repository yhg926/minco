# minco User Manual

`minco` builds fixed-size context-minhash sketches and estimates genome or read
ANI from those sketches. It is a pure C command line tool.

The public command set is:

```text
minco sketch   create, inspect, index, filter, append, and deduplicate sketches
minco ani      estimate ANI from sketches or direct FASTA/FASTQ inputs
minco matrix   report context-distance matrices, sparse edges, clusters, plans
minco set      combine, subset, and downsample existing sketches
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
minco.ctxobj64                 sketch entries
minco.ctxobj64.offsets         per-sample offsets into minco.ctxobj64
minco.stat                     sketch parameters, target sketch size, sample table, density summary
minco.ctxmeta                  optional binary density and unique-context records
minco.infilemeta               optional input metadata
minco.qc                       optional readsQC count ranges
minco.anno                     optional FASTA/FASTQ header annotations
minco.ctxobj64.abund           optional per-entry abundance/counts
minco.ctxobj64.position        optional source positions
minco.refindex.ctxgid64obj32   optional reference index built by `minco sketch -i`
```

Build `minco.refindex.ctxgid64obj32` before large reference searches,
all-vs-all ANI, sparse matrix reports, clusters, or dedup plans:

```bash
minco sketch -i genomes.minco
```

Downsample an existing larger bottom-k sketch without re-sketching:

```bash
minco set --downsample -S 1000 -o genomes.1k.minco genomes.20k.minco
minco sketch -i genomes.1k.minco
```

`set --downsample` keeps the first `S` entries from each already sorted
bottom-k sample sketch, rewrites `minco.stat` and context metadata for the new
target size, and preserves sample metadata, annotations, abundance, and
position sidecars when present. It does not copy `minco.refindex.ctxgid64obj32`
because the index must match the resized sketch.

New sketches write only the `minco.*` / `minco.ctxobj64*` names. Readers still
accept earlier KSSD-style sketch filenames so existing sketch directories remain
usable.

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
the temporary query sketch inherits the reference sketch size recorded in
`minco.stat`. Use `-S/--sketch-size` for sequence-vs-sequence auto sketching
or when no existing sketch is available as a template.

Reference-density sequence-query ANI:

```bash
minco ani -p8 -r ref.minco -q query.fasta.gz --query-density ref -m0 -o ani.tsv
```

Default `--query-density fixed` keeps the normal fixed bottom-k query sketch.
`--query-density ref` instead keeps all query context hashes up to a
reference-derived threshold inside the ANI run. The temporary query sketch is
deleted automatically, so users do not need to create or manage it. If the
reference sketch has one sample, minco uses that sample's density threshold
from `minco.stat` or `minco.ctxmeta`. If the reference is a combined sketch,
minco uses the largest sample density threshold from `minco.stat` or scans
`minco.ctxmeta` when needed. Legacy TSV metadata is accepted as a fallback for
old experiment sketches. This option affects only FASTA/FASTQ query inputs;
query sketches are used as-is.

To inspect the generated query sketch for debugging, save it explicitly:

```bash
minco ani -p8 -r ref.minco -q query.fasta.gz --query-density ref \
  --save-query-sketch query.debug.minco -m0 -o ani.tsv
```

`--save-query-sketch` saves the exact generated query sketch used by ANI. For
multiple query sequence inputs, this is the merged query sketch. The target
directory must not already exist.

For direct raw-read FASTQ queries, reference density can be applied without
materializing a global query sketch:

```bash
minco ani -p8 -r ref.minco --qraw reads.fastq.gz --query-density ref -m0 -o reads_vs_ref.tsv
```

When `--qraw` receives FASTQ input, `--query-density ref` streams reads
readwise by default unless `--save-query-sketch`, `--readsQC`, or `--abundance`
requires the materialized sketch path. The readwise path tracks unique
reference/query context coverage for AF and appends `Reads_with_ctx_match`,
`Total_reads`,
`Read_match_fraction`, `Unique_query_ctx`, `Unique_query_ctx_hit`, and
`Unique_ref_ctx_hit` to detail output. It also reports
`Density_block_ctx`, `Total_density_blocks`, `Blocks_with_ctx_match`, and
`Block_match_fraction`; by default `Density_block_ctx` is 100. `AF_source` is
reported as `readwise_coverage`. Long runs report progress to stderr with reads
processed, read rate, and input-file percent when the query is a regular file.
When minco is built with OpenMP and `-p` is greater than 1, the direct readwise
FASTQ density path reads batches of 16,384 reads, processes retained density
contexts in parallel, and merges each batch before reading the next one. This
keeps memory bounded for large streamed FASTQ files while preserving the same
output fields and read counts.
By default minco accumulates consecutive reads into pseudo-read blocks before
lookup. The block is processed when it has at least 100 retained density
contexts, and the final partial block is processed at end of file. Use
`--density-block-ctx N` to choose another block size, or
`--density-block-ctx 0` to force exact per-read lookup. For example,
`--density-block-ctx 100` with an `S=1000` reference tests 100-context blocks
without forcing every individual read to contribute a context.
In this readwise density path, `-t/--ctxcut` is applied to unique context
overlap support, not the summed read-hit `XnY_ctx`. Use a sketch-size-scaled
cutoff for comparable discovery: `-t100` at `S=10000` is roughly comparable to
`-t10` at `S=1000`.

Experimental depth abundance can be requested on the same readwise path:

```bash
minco ani -p8 -r ref.minco --qraw reads.fastq.gz --query-density ref \
  --abundance-est depth -m0 -o reads_vs_ref.abundance.tsv
```

This appends `Ref_breadth`, `Ref_mean_depth`, `Ref_hit_mean_depth`,
`Ref_depth_variance`, `Ref_depth_cv`, `Ref_zero_fraction`,
`Relative_abundance_depth`, `Normalized_abundance_depth`, and `Default_call`.
Depth statistics use the occurrence coverage of each reference context entry.
ANI mutation features (`XnY_ctx`, `N_diff_obj`, `N_diff_obj_section`, and
`N_mut2_ctx`) are reduced once per unique reference context entry using the best
observed context-object difference across reads or density blocks, so repeated
coverage does not inflate the mutation-distance term.
With no custom `-f`, `-n`, `-t`, or `--top`, minco applies the default readwise
abundance report: it sets the unique-context support cutoff to
`min(S, max(100, ceil(S/100)))`, removes the hidden AF cutoff, uses ANI cutoff
0.96, and prints only `major` and `low_abundance` calls.
Default calls use these rules:

```text
support        min(S, max(100, ceil(S/100))); applied to XnY_ctx and Unique_ref_ctx_hit
major          Ref_breadth >= 0.5, Relative_abundance_depth >= 1e-4, support pass, ANI >= 0.96
low_abundance  Ref_breadth >= 0.5, Relative_abundance_depth >= 1e-5, support pass, ANI >= 0.96
weak           below the default call thresholds
```

`Normalized_abundance_depth` sums to 1 across the rows printed after filters and
`--top`. Use explicit filters such as `-f0 -n0 -t0` when every candidate
comparison must be reported. The estimate is a transparent breadth/depth
baseline; it does not yet deconvolve shared contexts with EM or fit NB/ZINB
mixture models. It currently requires direct FASTQ `--qraw --query-density ref`
without `--readsQC`, `--abundance`, or `--save-query-sketch`.

For very large metagenomes, use `--readwise-profile-only` together with
`--abundance-est depth`. This mode keeps the reference coverage/depth
accumulators but does not build the exact global query-context and query-ref
context sets. It is intended for CAMI-style profiling and other large
metagenome abundance runs where exact readwise query AF would otherwise grow
with all unique sample contexts. In profile-only mode `Unique_query_ctx` is 0,
`Unique_query_ctx_hit` is approximated from unique reference-context hits,
query AF is reported as the reference breadth, and naive ANI is still computed
from unique-best context-object differences.

### CAMI Taxonomic Profile Output

`minco ani` can write a CAMI-style taxonomic profile from the same final rows
printed by direct readwise FASTQ abundance mode:

```bash
minco ani -p16 -r ref.minco --qraw reads.fastq.gz --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --cami-taxmap ref.cami_taxmap.tsv \
  --cami-profile sample.profile --cami-sample-id sample_1 \
  -m0 -o reads_vs_ref.tsv
```

`--cami-profile FILE` writes a profile with CAMI header lines:

```text
@SampleID:sample_1
@Version:0.9.1
@Ranks:...
@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE	_CAMI_genomeID	_CAMI_OTU
```

`PERCENTAGE` is `Normalized_abundance_depth * 100` after normal minco filters
and `--top` are applied. If several printed references map to the same
`TAXID/RANK/TAXPATH`, minco sums them into one CAMI row. A taxmap may contain
multiple rows for the same `ref_key`, one per rank, so minco can emit full
lineage profiles when those lineage rows are provided. In a complete lineage
taxmap, percentages are expected to sum to about 100 within each rank, not
across all lineage rows together. References missing from the taxmap are
skipped with a warning.

`--cami-taxmap` is a tab-delimited file with this schema:

```text
ref_key	TAXID	RANK	TAXPATH	TAXPATHSN	[_CAMI_genomeID]	[_CAMI_OTU]
```

Header lines are accepted. `ref_key` may be the stored minco reference path,
the reference basename, an assembly accession such as `GCF_009858895.2`, or a
sequence accession from `minco.anno` such as `NC_045512.2`.
To report a full lineage, repeat the same `ref_key` for each desired rank:

```text
GCF_000000001.1	2	superkingdom	2	Bacteria
GCF_000000001.1	1224	phylum	2|1224	Bacteria|Pseudomonadota
GCF_000000001.1	561	species	2|1224|...|561	Bacteria|Pseudomonadota|...|Escherichia coli
```

For a RefSeq viral sketch, a simple species-level map can be generated from
`virus_assembly_summary.txt`:

```bash
awk -F'\t' '!/^#/ {
  print $1 "\t" $7 "\tspecies\t10239|" $7 "\tViruses|" $8
}' virus_assembly_summary.txt > refseq_virus.cami_taxmap.tsv
```

The CAMI writer is intentionally map-driven. minco can aggregate provided full
lineages, but it does not infer missing lineage ranks from FASTA headers alone.

Current `minco.stat` records the requested target sketch size and the compact
density summary used by `--query-density ref`. `minco sketch --psmp DIR` shows
the actual entry count for each sample, which can be smaller than the target
for very small inputs.

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
-t/--ctxcut  3 contexts; in direct readwise FASTQ density mode, this is unique context overlap
```

When `minco.ctxmeta` is available, `-f/--afcut` uses density-estimated
real aligned fractions rather than the fixed-sketch fractions. Detail output
keeps `Qry_align_fraction` and `Ref_align_fraction` for compatibility and
appends `Real_Qry_align_fraction`, `Real_Ref_align_fraction`,
`Real_min_align_fraction`, and `AF_source`. Without ctxmeta, `AF_source` is
`sketch` and the real-AF columns fall back to the fixed-sketch values.
For direct readwise FASTQ density ANI, AF is computed from auxiliary unique
context coverage and `AF_source` is `readwise_coverage`.

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
score. In a narrow low-AF complete-assembly case, `Best` can mark
`guarded_low_confidence` and use a density-aware exact-context AAF fallback
instead of the HGB value. Use `-s5` or `-s6` explicitly when you want raw MashD
or AafD.
Use `--raw-output -s3` to inspect the raw CtxMoE value.

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

`--keep` and `--remove` filter per-sample sidecars together with the sketch
payload. If the input sketch has `minco.ctxmeta`, the output keeps only the
selected metadata records and refreshes the compact density summary in
`minco.stat`; filtered reference sketches therefore remain compatible with
`minco ani --query-density ref`.

Deduplicate:

```bash
minco sketch --dedup 0.001 --metric ctx-moe -o dedup.minco genomes.minco
minco sketch --dedup 0.001 --dedup-index -o dedup_indexed.minco genomes.minco
```

Use `--drop-position` with keep/remove/dedup output if the position sidecar is
not needed in the result.

## Context Metadata

`--ctxmeta` writes binary per-sample context metadata to `minco.ctxmeta` and
stores the compact sketch-set density summary in `minco.stat`.

Modes:

```text
none          do not write context metadata
preconflict   estimate before conflict resolution, default
postconflict  estimate after conflict resolution
both          write both estimates
```

The unique-context estimator uses the final bottom-k threshold:

```text
observed_contexts_under_threshold * hash_space / (hash_threshold + 1)
```

The per-sample metadata helps estimate total unique contexts but does not
change sketch entries or the index. The compact operational density summary in
`minco.stat` keeps both the minimum and largest sample density; for combined
references, the largest sample density is the query/universal threshold used by
`--query-density ref`. Print text views with `minco sketch --pctxmeta DIR` and
`minco sketch --pctxsetmeta DIR`.

`minco ani --query-density ref` uses this metadata when a query is supplied as
FASTA/FASTQ. A one-sample reference applies that sample's threshold; a combined
reference applies the largest sample density threshold so the auto-sketched
query is extracted internally at a reference-dependent density instead of a
fixed bottom-k size. The temporary query sketch is removed after ANI output is
written unless `--save-query-sketch DIR` is used.

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
