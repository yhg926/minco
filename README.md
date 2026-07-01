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

Calibrated species profiling uses the default wrapper:

```bash
scripts/minco_profile \
  -r ref.minco --reads reads.fq.gz -p16 \
  -o calibrated.profile.tsv
```

Before an expensive run, check that a reference bundle is self-contained:

```bash
scripts/minco_profile --check-ref -r ref.minco --reads reads.fq.gz
```

This is the recommended no-manual-strategy entry point for the current MinCO
species default. `scripts/minco_profile` is a short launcher over
`minco_profile_default.py`, which delegates to
`minco_profile_calibrated.py`, adds `--strategy universal-auto-exact`, and uses
the `candidate` preset unless the caller explicitly chooses another preset.
Existing scripts can continue to call `scripts/minco_profile_default.py`
directly.
The candidate preset enables strict split-evidence profile rescue, candidate
surface calls, abundance-only ANI/sparse-depth guards, and
`normalized-depth-alpha2` candidate abundance as one default path. Use
`--no-profile-rescue` to keep the candidate preset but disable profile-rescue
additions. Use `--profile-preset current` or
`MINCO_PROFILE_PRESET=current` to reproduce the previous calibrated default
without candidate additions. Packaged species databases can ship
`species_taxmap.tsv` and `joined_feature_training/` beside `ref.minco`, so the
command above has no strategy or calibration-path flags.
For unpackaged layouts, pass `--taxmap` and `--train-features`, or set
`MINCO_PROFILE_TAXMAP` and `MINCO_PROFILE_TRAIN_FEATURES`. A packaged database
can also ship a fitted RF/HGB model cache such as
`minco_profile_rf_hgb.train12.unfiltered.joblib`; pass `--model-cache` or set
`MINCO_PROFILE_MODEL_CACHE` to skip refitting the calibration models for every
sample. If the database ships `candidate_surface_taxmap.tsv` or
`accession_species_taxmap.tsv` beside `ref.minco`, the launcher passes it as
the candidate-surface label map; the same path can also be supplied with
`MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP`. Set `MINCO_PROFILE_MINCO` to choose a
non-default `minco` binary, or pass `--minco` explicitly.
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

Direct raw-read profiling against a species, AMR, virus, gene, or mixed
reference sketch uses the C subcommand:

```bash
bin/minco profile -p16 -r ref.minco reads.fq.gz -o profile.tsv
```

`profile` is the conservative C wrapper for direct readwise abundance
profiling. It uses reference-density extraction, depth abundance, profile-only
memory bounds, best-diff-split shared-context assignment, naive readwise ANI,
and product-topfrac-median context defaking by default. Use `--report-all` when
a benchmark or debugging run needs every reference comparison instead of only
the default called rows.

In the calibrated wrapper, `train12` includes the strainmadness joined-feature
tables when present. Use `train9` for locked strainmadness holdout testing.
`--strategy universal`
applies the current fixed experimental gate: adaptive+sub95 calls, a
high-confidence raw-unique fallback, guarded high-extra/low-unique-AF tail
rescue, guarded low-extra split rescue, and panel zip-corrected depth abundance
over the final retained species. The default `--strategy universal-auto-exact`
first runs that block-mode gate, then reruns the split pass with
`--density-block-ctx 0` when the block-mode probability-rescue mass is low and
the high-uAF/raw-unique guard passes. By default it skips that exact rerun when
the block-mode low-extra split rescue already added candidates; pass
`--exact-split-low-extra-mode allow` to keep the older exact behavior. This is
the current best documented MinCO F1-priority default, at the cost of variable
runtime. The default launcher also applies the candidate preset: strict
split-evidence profile rescue, candidate-surface calls, abundance-only
ANI/sparse-depth guards, plus normalized-depth candidate abundance. This should
not be read as a general Sylph-beating abundance/default claim; Sylph remains
stronger on some held-out panels. Use `--no-profile-rescue` to disable only the
strict rescue, use
`--profile-preset current` to reproduce the previous MinCO default, and use
`--strategy probability` only to reproduce the legacy RF/HGB threshold-only
wrapper output.

For speed-priority calibrated screening, use `--strategy universal` with a
model cache and leave `--report-all` off. This always skips the exact split
rerun and writes only called rows. With the default low-extra exact-skip rule,
CAMI II Toy Mouse sample6 ran in 1:47.90 with 3.43 GiB peak RSS versus Sylph
sketch+profile at 1:50.98 with 18.81 GiB; the same sample still favored Sylph
for F1 and abundance accuracy, so this is a runtime win rather than a broad
accuracy win.

The calibrated wrapper reports `reported_ani` as a diagnostic continuous ANI
field. It uses split `Ref_zip_aaf_ani` by default, or unique `Ref_zip_aaf_ani`
when the raw-unique fallback is active. This avoids the saturated raw/emitted
readwise ANI field in metagenome profiles, but it is not used as a separate
default call gate.

For retained species, `calibrated_abundance` is normalized from a broad-panel
depth rule: normal rows use the larger of split and unique
`Ref_mean_depth / Ref_zip_af`, while tail-rescue-added rows use split
`Ref_mean_depth` directly. Fixed-call sweeps over Toy Mouse, CAMI3
source-readmap, HMP airskin, and HMP gastrooral found no all-panel abundance
replacement, so this remains the default. `--abundance-genus-xny-blend-alpha`
is available only as an experimental abundance-only diagnostic; its default is
`0.0`, which leaves the selected default abundance unchanged. The experimental
`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002` option
applies a guarded within-genus abundance reallocation using hit-depth and
breadth features when the output-derived guard passes; its default is `off`.

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
For conservative direct readwise profiling, prefer the shorter `minco profile`
command above; it selects these readwise abundance settings automatically and
keeps `ani` as the lower-level ANI/debug interface. For the current
F1-priority calibrated species default candidate, use
`scripts/minco_profile` as described above.
With no custom `-f`, `-n`, `-t`, or `--top`, minco applies
the default readwise abundance report: it uses
a sketch-size-scaled support cutoff
`min(S, max(100, ceil(S/100)))`, requires `Ref_breadth >= 0.5` and
`ANI >= 0.95`, and prints only `major` and `low_abundance` calls. The
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

For AMR/gene reference databases, treat the readwise table as a determinant
screen rather than exact allele truth. A practical clinical-style report should
group close allele hits to determinant/family names and emit detected rows such
as `high_confidence` and `screening_positive`; exact allele labels are best
shown as representative hints. Normal MinCO runs cannot know which genes are
missing from the installed reference database, so they should not print absent
determinant names. A `not_in_refdb` status is meaningful only in a benchmark or
comparator report where an external expected list, such as an NCBI AMR table,
is supplied.

Tag AMR or gene-panel sketches before distribution or merge so these rules
travel with the database:

```bash
bin/minco set --set-domain amr minco-db-amr-s100.minco
```

`minco.domain` is one byte per reference and is preserved by sketch append,
merge, downsample, keep/remove, and markerdb operations. Untagged references in
a mixed merged database remain the default species profile. AMR/gene readwise
`Default_call` uses domain-specific gates: `major` requires
`Ref_breadth >= 0.60`, `support >= max(30, ceil(0.03*S))`, and `ANI >= 0.96`;
`screening` requires `Ref_breadth >= 0.35`,
`support >= max(20, ceil(0.02*S))`, and `ANI >= 0.95`.

Write Kraken-like read tracking as a sidecar table:

```bash
bin/minco profile -p16 -r ref.minco reads.fq.gz \
  --track reads.track.tsv \
  --taxonomy both \
  --gtdb-taxmap ref.gtdb_taxmap.tsv \
  --ncbi-taxmap ref.ncbi_taxmap.tsv \
  -o profile.tsv
```

`--track` writes one row for each read with at least one selected
reference hit. It records the read id, read length, retained density-context
count, matched/selected context counts, selected context offsets from the
0-based read start, target reference count, target reference list, and optional
GTDB/NCBI lowest common ancestor labels. Multi-reference read hits are assigned
to the closest common ancestry over the selected best-diff target references.
Tracking forces exact per-read density units (`--density-block-ctx 0`), because
block mode intentionally loses read identity and context offsets. A summary TSV
is written to `<reads.track.tsv>.summary.tsv` unless
`--track-summary` is set; it includes tracked-read percentage,
density-positive no-hit read percentage, and context-level reference-hit
fractions. When `minco.ctxmeta` is available, the summary reports
`estimated_unknown_reads_pct`, a Horvitz-Thompson estimate of the percentage of
read contexts absent from the full reference genome set. It divides each
observed context hit by its sketch capture probability. For shared contexts,
the capture probability uses the largest reference sketch density among
candidate references. This is the efficient S2000-style estimate for
full-context unknown read percentage, but it is a whole-genome estimate only
for full density-sampled reference sketches, not markerdbs.

For abundance-method debugging, `profile` can also dump context-level
ambiguity edges:

```bash
bin/minco profile -p16 -r ref.minco reads.fq.gz \
  --edge-out reads.edges.tsv \
  --edge-max 10000000 \
  -o profile.tsv
```

`--edge-out` writes one row per candidate reference edge in each retained
read/context ambiguity group. Columns include `read_id`, `unit_id`, `qctx`,
candidate `gid`, `diff`, `best_diff`, `candidate_refs`, `selected_refs`,
`selected_by_mode`, and `cov_inc`. By default only ambiguous groups are
written and groups with more than 64 candidate refs are skipped; use
`--edge-all`, `--edge-selected-only`, `--edge-max`, and
`--edge-max-candidates` to change this. Edge output forces
`--density-block-ctx 0` so the context groups are exact per-read events. It is
an experimental diagnostic/EM substrate and is not part of the calibrated
species default.

Write a CAMI taxonomic profile from the printed readwise abundance rows:

```bash
bin/minco profile -p16 -r ref.minco reads.fq.gz \
  --cami-taxmap ref.cami_taxmap.tsv \
  --cami-profile sample.profile --cami-sample-id sample_1 \
  -o profile.tsv
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
