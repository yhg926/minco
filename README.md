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

Build a context-marker reference sketch that keeps only contexts present in
exactly one reference:

```bash
bin/minco sketch -i genomes.minco
bin/minco matrix --format dedup-plan -m aaf --cut 0.001 \
  --keep-out keep.txt --remove-out remove.txt -o dedup_plan.tsv genomes.minco
bin/minco sketch --keep keep.txt -o genomes.dedup.minco genomes.minco
bin/minco set --uniq_union --markerdb-ctx -o genomes.ctxmarker.minco genomes.dedup.minco
bin/minco sketch -i genomes.ctxmarker.minco
```

`matrix --format dedup-plan` predicts the post-dedup context-markerdb size for
kept references and warns for kept references below 500 entries by default.
Change this with `--markerdb-warn-threshold N`, or use `0` to disable the
dedup-stage warning.

Markerdb creation warns when any reference keeps fewer than 500 marker entries;
change this with `--markerdb-warn-threshold N`, or use `0` to disable it.
Markerdb sketches can have different retained marker counts per reference.
Those counts are read from `minco.ctxobj64.offsets`, so `minco sketch --psmp`
reports the real post-markerdb size for each sample. Markerdb generation
preserves the conflict-object mode of the input sketch.

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

In this mode minco streams reads and uses auxiliary unique context coverage to
report AF. With `--abundance-est depth`, depth remains coverage-based while ANI
mutation features are reduced once per unique reference context entry using the
best observed context-object difference. Detail output appends
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
For large combined reference sets, direct readwise profiling uses
`--readwise-assign best-diff-unique` by default. For each retained query
context, minco compares all matching reference context-objects and counts the
context only when exactly one reference has the best object-difference score.
This avoids inflating ANI/support for conserved contexts shared by many
references. Use `--readwise-assign all` to reproduce the legacy behavior,
`best-diff` to keep all tied best matches, or `best-diff-split` to split depth
coverage and breadth across tied best matches while leaving XnY-style support
counts unweighted.
With depth abundance, the selected readwise ANI defaults to
`--readwise-ani zip-aaf`, which estimates latent context AF from
`Ref_breadth` and `Ref_mean_depth` with a zero-inflated Poisson model and then
converts that AF to context AAF ANI. Use `--readwise-ani naive` to report the
object-difference ANI instead.
Experimental reliability modes can be tested with `--readwise-ctx-filter`.
`poisson-diff` hard-filters suspicious contexts used for readwise
naive/object-difference ANI by combining breadth-implied depth surprise with
best object difference. `poisson-depth`, `poisson-product`, and `product-nb`
test Bonferroni-adjusted high-tail depth or `best_diff * depth` outliers.
`product-topfrac` ranks all hit contexts by `best_diff * depth` and drops the
top per-reference fraction, but only positive-product contexts are removed; set
that fraction with `--readwise-fake-threshold` such as `0.25`.
`product-topfrac-median` uses the same top-fraction detector but keeps those
contexts and replaces their `best_diff` and depth with per-reference medians.
`product1-topfrac-median` ranks by `(best_diff + 1) * depth`, allowing high-depth
exact-match contexts to be median-replaced too.
`fake-prob` uses the same score family to estimate `P(fake context)` and
soft-weights nonzero-diff contexts for naive ANI diagnostics. Tune strictness
with `--readwise-fake-threshold`; lower values are more aggressive. These modes
are research options and are not enabled by default.
Add `--abundance-est depth` to also append experimental per-reference
breadth/depth metrics, including `Relative_abundance_depth` and
`Normalized_abundance_depth`, `Ref_zip_af`, `Ref_zip_aaf_ani`, plus
post-filter diagnostic `Reliable_Ref_*` columns and `Default_call`.
With no custom `-f`, `-n`, `-t`, or `--top`, minco applies
the default readwise abundance report: it uses
a sketch-size-scaled support cutoff
`min(S, max(100, ceil(S/100)))`, requires `Ref_breadth >= 0.5` and
`ANI >= 0.96`, and prints only `major` and `low_abundance` calls. The
normalized column sums to 1 across printed rows.
Use explicit filters such as `-f0 -n0 -t0` when every candidate comparison must
be reported. This abundance option requires the direct readwise FASTQ density
path, so do not combine it with `--readsQC`, `--abundance`, or
`--save-query-sketch`.
For very large metagenomes and large reference sets, add
`--readwise-profile-only` with `--abundance-est depth`. This keeps the
coverage/depth abundance path but skips exact global query-context sets, which
prevents memory from growing with every unique metagenome context. In this
mode `Unique_query_ctx` is reported as 0, query AF is approximated from
reference breadth for filtering/reporting, and naive ANI still uses
unique-best matched context-object differences.

Write Kraken-like read tracking as a sidecar table:

```bash
bin/minco ani -p16 -r ref.minco --qraw reads.fq.gz --query-density ref \
  --readwise-track reads.track.tsv \
  --readwise-taxonomy both \
  --gtdb-taxmap ref.gtdb_taxmap.tsv \
  --ncbi-taxmap ref.ncbi_taxmap.tsv \
  -m0 -o reads_vs_ref.tsv
```

`--readwise-track` writes one row for each read with at least one selected
reference hit. It records the read id, read length, retained density-context
count, matched/selected context counts, selected context offsets from the
0-based read start, target reference count, target reference list, and optional
GTDB/NCBI lowest common ancestor labels. Multi-reference read hits are assigned
to the closest common ancestry over the selected best-diff target references.
Tracking forces exact per-read density units (`--density-block-ctx 0`), because
block mode intentionally loses read identity and context offsets. A summary TSV
is written to `<reads.track.tsv>.summary.tsv` unless
`--readwise-track-summary` is set; it includes tracked-read percentage,
density-positive no-hit read percentage, and context-level reference-hit
fractions. `estimated_ref_absent_ctx_pct` is estimated from sampled read
contexts as `1 - total_matched_ctx / total_density_ctx`; it estimates
whole-genome reference absence only when the reference database contains
whole-genome density-sampled contexts. On markerdbs it measures absence from
the retained marker context set.

Write a CAMI taxonomic profile from the printed readwise abundance rows:

```bash
bin/minco ani -p16 -r ref.minco --qraw reads.fq.gz --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --cami-taxmap ref.cami_taxmap.tsv \
  --cami-profile sample.profile --cami-sample-id sample_1 \
  -m0 -o reads_vs_ref.tsv
```

For sensitive CAMI-style species discovery on a large reference set, the
current experimental setting is:

```bash
bin/minco ani -p16 -r ref.minco --qraw reads.fq.gz --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique \
  --readwise-ani zip-aaf \
  -f0.05 -n0.94 -t10 -m0 -o reads_vs_ref.tsv
```

In CAMI II marine samples 0-2, a GTDB-only `S=1000` reference with this recipe
matched `S=10000` accuracy while using much less memory; see
`research/PROGRESS.md` and
`research/experiments/2026-06-20_cami2_marine_s1000_choice/NOTE.md`.

`--cami-taxmap` is a tab-delimited map with at least five columns:
`ref_key`, `TAXID`, `RANK`, `TAXPATH`, and `TAXPATHSN`, followed by optional
`_CAMI_genomeID` and `_CAMI_OTU`. `ref_key` can be the stored minco reference
path, the reference basename, an assembly accession such as `GCF_...`, or a
sequence accession such as `NC_...`. To emit a full lineage profile, include
multiple rows with the same `ref_key`, one per rank; minco aggregates rows with
the same `TAXID/RANK/TAXPATH`. With a complete lineage taxmap, percentages sum
to about 100 within each rank. Minco does not infer missing ranks on its own.

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
