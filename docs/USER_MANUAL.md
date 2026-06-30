# minco User Manual

`minco` builds fixed-size context-minhash sketches and estimates genome or read
ANI from those sketches. It is a pure C command line tool.

The public command set is:

```text
minco sketch   create, inspect, index, filter, append, and deduplicate sketches
minco ani      estimate ANI from sketches or direct FASTA/FASTQ inputs
minco profile  profile raw reads with conservative direct depth defaults
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

Build a marker reference sketch from already sketched references:

```bash
minco sketch -i genomes.minco
minco matrix --format dedup-plan -m aaf --cut 0.001 \
  --keep-out keep.txt --remove-out remove.txt -o dedup_plan.tsv genomes.minco
minco sketch --keep keep.txt -o genomes.dedup.minco genomes.minco
minco set --uniq_union --markerdb-ctx -o genomes.marker.minco genomes.dedup.minco
minco sketch -i genomes.marker.minco
```

`set --uniq_union --markerdb` keeps full `ctxobj64` entries that occur in only
one reference. `set --uniq_union --markerdb-ctx` is stricter at the context
level: it keeps all entries for a context only when that context is present in
exactly one reference, which removes contexts shared across references even when
their object bits differ.

For dedup-first workflows, `matrix --format dedup-plan` can use AAF distance
(`-m aaf`) to write `--keep-out` and `--remove-out` lists before `sketch --keep`
or `sketch --remove`. During that dedup stage it predicts the final
context-markerdb size across the kept references and warns on stderr when a kept
reference is expected to retain fewer than 500 marker entries. Set
`--markerdb-warn-threshold N` to change the cutoff, or set it to `0` to disable
that dedup-stage warning.

Both markerdb modes warn on stderr when a reference retains fewer than 500
marker entries. Set `--markerdb-warn-threshold N` to change the cutoff, or set
it to `0` to disable the warning.

Markerdb output is not fixed-size across samples. The authoritative per-sample
marker count is `minco.ctxobj64.offsets[i + 1] - minco.ctxobj64.offsets[i]`;
`minco sketch --psmp` prints that value. Markerdb creation preserves the input
sketch conflict-object mode: conflict-free input remains conflict-free, while
input sketched with `--conflict` keeps its retained conflicting context objects.

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

For calibrated species profiling, use the default wrapper:

```bash
scripts/minco_profile \
  -r ref.minco --reads reads.fastq.gz -p16 \
  -o calibrated.profile.tsv
```

Before an expensive run, check that a reference bundle is self-contained:

```bash
scripts/minco_profile --check-ref -r ref.minco --reads reads.fastq.gz
```

This is the recommended no-manual-strategy entry point for the current MinCO
species default. `scripts/minco_profile` is a short launcher over
`minco_profile_default.py`, which delegates to
`minco_profile_calibrated.py`, adds `--strategy universal-auto-exact`, and uses
the `candidate` preset unless an advanced caller explicitly chooses another
preset. The candidate preset enables candidate rescue/surface calls and
`normalized-depth-alpha2` candidate abundance as one default path. Existing
scripts can continue to call `scripts/minco_profile_default.py` directly. Use
`--profile-preset current` or `MINCO_PROFILE_PRESET=current` to reproduce the
previous calibrated default without those candidate additions. Packaged species
databases can ship `species_taxmap.tsv` and `joined_feature_training/` beside
`ref.minco`, so the command above has no strategy or calibration-path flags.
For unpackaged layouts, pass `--taxmap` and `--train-features`, or set
`MINCO_PROFILE_TAXMAP` and `MINCO_PROFILE_TRAIN_FEATURES`. Packaged databases
can also ship a fitted RF/HGB model cache such as
`minco_profile_rf_hgb.train12.unfiltered.joblib`; pass `--model-cache` or set
`MINCO_PROFILE_MODEL_CACHE` to skip refitting the calibration models for every
sample. If a packaged database includes `candidate_surface_taxmap.tsv` or
`accession_species_taxmap.tsv` beside `ref.minco`, the launcher passes it as
the candidate-surface label map; use
`MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP` for an explicit deployment path. Set
`MINCO_PROFILE_MINCO` to choose a non-default `minco` binary, or pass `--minco`
explicitly.
`--check-ref` and its alias `--preflight` stop after validating these packaged
inputs and print the resolved sidecars and delegated strategy.
Advanced validation runs can override `--candidate-surface-switch
accession-ani93-xny650-br20` to test the stricter raw-retention surface mode;
`accession-current-or-ani93-xny650-br20` preserves the current candidate
surface while adding that stricter validation surface. The default preset does
not select either mode. `--candidate-surface-max-called-species N` can be used
in validation runs to disable candidate-surface additions on very large
pre-surface call sets; the default `0` leaves this guard off.

When raw reads are supplied, the wrapper auto-selects the initial evidence
scheduler: `-p1` uses same-stream unique sidecar generation, while `-p>=2`
runs unique and split readwise MinCO passes concurrently and splits the thread
budget between them. Use `--sequential-readwise-passes`,
`--same-stream-readwise-passes`, or `--legacy-dual-readwise-passes` only for
debugging or scheduler benchmarks.
Advanced speed experiments can also pass `--same-stream-exact-split` to write
the exact split table during the initial split pass and avoid a later exact
rerun when the auto-exact gate fires. This is not the default because it adds
work when the exact gate does not fire.

For direct species, AMR, virus, gene, or mixed-domain profiling without the
calibrated species gate, use the simpler C `profile` subcommand:

```bash
minco profile -p16 -r ref.minco reads.fastq.gz -o profile.tsv
```

`profile` expands to the conservative direct readwise density path with depth
abundance, profile-only memory bounds, `best-diff-split` assignment,
`--readwise-ani naive`, and `product-topfrac-median` context defaking at
threshold 0.25. It prints the default called rows unless `--report-all` is set.
This keeps `minco ani` available for assembled-genome ANI, matrices, and
benchmark sweeps that need explicit low-level flags.

The calibrated wrapper joins unfiltered `best-diff-unique` and
`best-diff-split` readwise evidence and applies the current RF/HGB gate.

`--train-pool train12` includes strainmadness joined-feature tables when the
training directory contains `test.joined_features.tsv`. Use `--train-pool
train9` when strainmadness itself is the holdout test. The calibrated wrapper
also accepts precomputed `--unique-table` and `--split-table` inputs so
benchmark runs do not need to repeat the two MinCO readwise passes.
`--strategy universal` is the current fixed experimental gate: adaptive+sub95
calls by default, a
high-confidence raw-unique fallback, guarded high-extra/low-unique-AF tail
rescue, guarded low-extra split rescue, and panel zip-corrected depth abundance
normalized over the final retained species. The default `--strategy
universal-auto-exact` evaluates the universal gate on normal density blocks,
then reruns only the split pass with `--density-block-ctx 0` when
`probability_extra_mass_ratio <= 0.10` and the high-uAF/raw-unique guard passes.
The default also skips that exact rerun when the block-mode low-extra split
rescue already added candidates; use `--exact-split-low-extra-mode allow` to
keep the older exact behavior. This is the current best documented MinCO
F1-priority default, not a general Sylph-beating abundance/default claim. The
default launcher also applies candidate rescue/surface calls plus
normalized-depth candidate abundance. Use `--profile-preset current` to
reproduce the previous MinCO default, and use `--strategy probability` to
reproduce the legacy RF/HGB threshold-only wrapper output.

For speed-priority calibrated screening, use `--strategy universal` with a
model cache and leave `--report-all` off. This always skips the exact split
rerun and writes only called rows. With the default low-extra exact-skip rule,
CAMI II Toy Mouse sample6 ran in 1:47.90 with 3.43 GiB peak RSS versus Sylph
sketch+profile at 1:50.98 with 18.81 GiB; the same sample still favored Sylph
for F1 and abundance accuracy, so this is a runtime win rather than a broad
accuracy win.

The calibrated wrapper prints `reported_ani` and `reported_ani_source` near the
front of the output table. `reported_ani` is a diagnostic continuous ANI value:
it uses split `Ref_zip_aaf_ani` by default, or unique `Ref_zip_aaf_ani` when the
raw-unique fallback is active. It is preferred over the raw/emitted readwise ANI
column for reporting because that raw field can saturate in metagenome
profiles, but it is not a separate default call gate.

For retained species, `calibrated_abundance` is normalized from a broad-panel
depth rule: normal rows use the larger of split and unique
`Ref_mean_depth / Ref_zip_af`, while tail-rescue-added rows use split
`Ref_mean_depth` directly. Fixed-call sweeps over Toy Mouse, CAMI3
source-readmap, HMP airskin, and HMP gastrooral found no all-panel abundance
replacement, so this remains the default. `--abundance-genus-xny-blend-alpha`
is an experimental abundance-only diagnostic; the default `0.0` leaves the
selected default abundance unchanged. The experimental
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002` option
applies a guarded within-genus abundance reallocation using hit-depth and
breadth features when the output-derived guard passes; its default is `off`.

This appends `Ref_breadth`, `Ref_mean_depth`, `Ref_hit_mean_depth`,
`Ref_depth_variance`, `Ref_depth_cv`, `Ref_zero_fraction`,
`Relative_abundance_depth`, `Normalized_abundance_depth`, and `Default_call`.
Depth statistics use the occurrence coverage of each reference context entry.
ANI mutation features (`XnY_ctx`, `N_diff_obj`, `N_diff_obj_section`, and
`N_mut2_ctx`) are reduced once per unique reference context entry using the best
observed context-object difference across reads or density blocks, so repeated
coverage does not inflate the mutation-distance term.
In combined reference sketches, retained readwise query contexts may hit many
references. minco now uses `--readwise-assign best-diff-unique` by default:
it first finds the best object-difference score among matching references and
counts the context only if that best hit is unique. This prevents conserved
shared contexts from inflating support and ANI for many references at once.
Use `--readwise-assign all` for legacy behavior, `best-diff` to keep all tied
best hits, or `best-diff-split` to split depth coverage across tied best hits
and compute fractional breadth/depth while keeping XnY-style support counts
unweighted.
When depth abundance is available, selected readwise ANI defaults to
`--readwise-ani zip-aaf`. This fits a zero-inflated Poisson correction from
`Ref_breadth` and `Ref_mean_depth` to estimate latent reference AF, then
converts that AF to context AAF ANI. Use `--readwise-ani naive` to select the
object-difference readwise ANI instead.
`--readwise-ctx-filter` exposes experimental context reliability modes for
readwise naive/object-difference ANI. `poisson-diff` hard-filters suspicious
contexts by comparing per-context depth to the depth expected from reference
breadth and combining that surprise with the best object difference.
`poisson-depth`, `poisson-product`, and `product-nb` test Bonferroni-adjusted
high-tail depth or `best_diff * depth` outliers. `product-topfrac` ranks all
hit contexts by `best_diff * depth` and drops the top per-reference fraction,
but only positive-product contexts are removed; set that fraction with
`--readwise-fake-threshold`, for example `0.25`. `product-topfrac-median` uses
the same detector but keeps those contexts and replaces their `best_diff` and
depth with per-reference medians. `product1-topfrac-median` ranks by
`(best_diff + 1) * depth`, allowing high-depth exact-match contexts to be
median-replaced too.
`fake-prob` uses the same score family to estimate `P(fake context)` and
soft-weights nonzero-diff contexts for diagnostic naive ANI.
`--readwise-fake-threshold` controls strictness; lower values are more
aggressive for probability modes. These modes are off by default and should be
treated as research options until validated across more CAMI samples.
With no custom `-f`, `-n`, `-t`, or `--top`, minco applies the default readwise
abundance report: it sets the unique-context support cutoff to
`min(S, max(100, ceil(S/100)))`, removes the hidden AF cutoff, uses ANI cutoff
0.95 for the default species/prokaryotic profile, and prints only `major` and
`low_abundance` calls.
Default calls use these rules:

```text
support        min(S, max(100, ceil(S/100))); applied to XnY_ctx and Unique_ref_ctx_hit
major          Ref_breadth >= 0.5, Relative_abundance_depth >= 1e-4, support pass, ANI >= 0.95
low_abundance  Ref_breadth >= 0.5, Relative_abundance_depth >= 1e-5, support pass, ANI >= 0.95
weak           below the default call thresholds
```

`Normalized_abundance_depth` sums to 1 across the rows printed after filters and
`--top`. The abundance table also includes `Ref_zip_af`, `Ref_zip_aaf_ani`, and
post-filter diagnostic `Reliable_Ref_*` columns. The `Ref_zip_*` values are the
model values used by `--readwise-ani zip-aaf`.
Use explicit filters such as `-f0 -n0 -t0` when every candidate
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

### AMR And Gene-Panel Reporting

AMR, viral, MGE, and gene-panel databases should not be interpreted with the
same reporting semantics as species profiling. For AMR reads, use the direct
readwise depth table as determinant evidence, then report determinant/family
calls rather than independent exact alleles. Exact allele names from the best
matching reference are useful as representative hints, but close alleles can
share enough contexts that exact-label precision is much lower than
determinant-level precision.

Domain-specific reporting rules are stored in an optional `minco.domain`
sidecar, one byte per reference. Tag an AMR or gene-panel database before
distribution or before merging it with species, viral, or other databases:

```bash
minco set --set-domain amr minco-db-amr-s100.minco
```

Supported profiles are `species`/`default`, `amr`, `virus`, and `gene`.
`minco.domain` is preserved by sketch append, merge, downsample, keep/remove,
and markerdb operations. If a tagged database is merged with an untagged
database, tagged references keep their domain profile and untagged references
use the default species profile. This allows one merged refdb to keep
domain-specific readwise `Default_call` behavior.

The current AMR/gene readwise defaults are:

```text
major      Ref_breadth >= 0.60, support >= max(30, ceil(0.03*S)), ANI >= 0.96
screening  Ref_breadth >= 0.35, support >= max(20, ceil(0.02*S)), ANI >= 0.95
weak       below those AMR/gene thresholds
```

Here `S` is the reference sketch size after any markerdb/downsample operation,
and support requires both `XnY_ctx` and `Unique_ref_ctx_hit` to pass the cutoff.
The default/species profile keeps the species-style `major` and `low_abundance`
rules based on `Ref_breadth >= 0.5`, relative depth, and the
`min(S, max(100, ceil(S/100)))` support cutoff.

Recommended user-facing AMR statuses are:

```text
detected/high_confidence     strong determinant evidence
detected/screening_positive  sensitive screen-positive call, lower confidence
not_detected                 determinant is in the refdb but did not pass evidence gates
cannot_assess/not_in_refdb   benchmark-only: expected by an external truth list but absent from the refdb
```

A normal MinCO AMR run can only report determinants represented in the
installed AMR reference database. It cannot know that `acrF`, `mdtM`, or any
other determinant is absent unless the user supplies an external expected list.
Therefore `not_in_refdb` rows belong in benchmark/comparator reports, not in
ordinary clinical screen output.

Example normal AMR report shape after determinant grouping:

```text
sample	determinant	status	confidence	best_allele_hint	ANI	AF	XnY_ctx	mean_depth
SRR39268261	aph(3'')-Ib	detected	high	aph(3'')-Ib	0.998627	0.970246	750	55.84
SRR39268261	aph(6)-Id	detected	high	aph(6)-Id	1.000000	1.000000	806	49.20
SRR39268261	blaTEM	detected	high	blaTEM-1	1.000000	1.000000	830	56.33
SRR39268261	dfrA	detected	high	dfrA8	1.000000	1.000000	479	93.47
```

Example comparator-only rows when an external NCBI AMR table is supplied:

```text
sample	expected_by	determinant	benchmark_status
SRR39268261	NCBI	acrF	expected_but_absent_from_minco_refdb
SRR39268261	NCBI	mdtM	expected_but_absent_from_minco_refdb
```

### Read Tracking Output

`minco profile` can also write a Kraken-like read tracking sidecar:

```bash
minco profile -p16 -r ref.minco reads.fastq.gz \
  --track reads.track.tsv \
  --track-summary reads.track.summary.tsv \
  --taxonomy both \
  --gtdb-taxmap ref.gtdb_taxmap.tsv \
  --ncbi-taxmap ref.ncbi_taxmap.tsv \
  -o profile.tsv
```

`--track FILE` writes one row per read with at least one selected
reference hit. It uses the same retained density contexts and
assignment rule as the readwise profiler. The sidecar does not change
the main ANI or abundance table.

Read tracking forces `--density-block-ctx 0`, because density blocks merge
consecutive reads and cannot preserve exact read ids or context offsets. The
tracking columns are:

```text
read_id
read_ord
read_len
possible_ctx
density_ctx
matched_ctx
selected_ctx
target_ref_count
target_refs
target_ref_ids
gtdb_rank
gtdb_name
ncbi_rank
ncbi_name
ctx_offsets
ctx_offsets_truncated
assignment_status
```

`ctx_offsets` are 0-based offsets from the read start for selected hit
contexts. Long offset and target-ref lists are capped in the row and marked by
`ctx_offsets_truncated` or a trailing `...` in the reference list. For reads
that select more than one target reference, taxonomy labels are the lowest
common ancestor over the selected target references, not over every raw
candidate. Use `--taxonomy gtdb`, `ncbi`, or `both`; GTDB and NCBI maps
use the same tab-delimited schema as `--cami-taxmap`.

The summary file defaults to `<FILE>.summary.tsv` unless
`--track-summary` is given. It reports `total_reads`,
`reads_with_density_ctx`, `reads_with_ref_hit`, `tracked_read_pct`,
`density_positive_no_ref_hit_read_pct`, `total_density_ctx`,
`total_matched_ctx`, `sampled_ctx_ref_hit_pct`, and, when `minco.ctxmeta` is
available, `estimated_unknown_reads_pct`. This estimates full-context reference
absence from a downsampled full reference sketch by Horvitz-Thompson correction:

```text
capture_probability(ctx) =
  max(reference_sketch_density for candidate refs containing ctx) /
  query_density

sketch_corrected_present_ctx =
  sum(observed matched sampled-read contexts / capture_probability(ctx))

estimated_unknown_reads_pct =
  100 * (1 - sketch_corrected_present_ctx / total_density_ctx)
```

The percentage is clamped to `[0, 100]` for reporting. Use this as a
whole-genome unknown-read estimate only with a full density-sampled reference
sketch such as S2000. On markerdbs, it is not a whole-genome unknown-read
estimate because the reference no longer contains all sampled genome contexts.

### Context Ambiguity Edge Output

For debugging abundance assignment and testing context-level EM, `minco
profile` can write an experimental ambiguity-edge sidecar:

```bash
minco profile -p16 -r ref.minco reads.fastq.gz \
  --edge-out reads.edges.tsv \
  --edge-max 10000000 \
  -o profile.tsv
```

`--edge-out FILE` writes one row per candidate reference edge in each retained
read/context group. The columns are:

```text
edge_id
read_id
unit_id
qctx
edge_rank
ref_begin
gid
diff
best_diff
candidate_refs
selected_refs
selected_by_mode
cov_inc
```

By default only ambiguous context groups are written, groups with more than 64
candidate refs are skipped, and output is capped at 10,000,000 rows. Use
`--edge-all`, `--edge-selected-only`, `--edge-max`, and
`--edge-max-candidates` to change those limits. Edge output forces
`--density-block-ctx 0`, because density blocks deduplicate contexts across
multiple reads and cannot provide exact per-read ambiguity events. This sidecar
can be large and slow; it is a research/debug substrate, not a default
calibrated species-profile output.

### CAMI Taxonomic Profile Output

`minco profile` can write a CAMI-style taxonomic profile from the same final
rows printed by direct readwise FASTQ abundance mode:

```bash
minco profile -p16 -r ref.minco reads.fastq.gz \
  --cami-taxmap ref.cami_taxmap.tsv \
  --cami-profile sample.profile --cami-sample-id sample_1 \
  -o profile.tsv
```

For sensitive species discovery with a large `S=10000` reference database,
start with the current experimental threshold:

```bash
minco ani -p16 -r ref.minco --qraw reads.fastq.gz --query-density ref \
  --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique \
  --readwise-ani zip-aaf \
  -f0.05 -n0.94 -t10 -m0 -o reads_vs_ref.tsv
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
