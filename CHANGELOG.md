# Changelog

## 2026-06-24

### Added

- Added `minco ani --readwise-track FILE`, a Kraken-like read tracking sidecar
  for direct `--qraw --query-density ref` readwise mode. It writes per-read
  selected reference hits, selected context offsets from the read start, target
  reference counts/lists, and optional GTDB/NCBI LCA labels.
- Added `--readwise-taxonomy none|gtdb|ncbi|both`, `--gtdb-taxmap`,
  `--ncbi-taxmap`, and `--readwise-track-summary`. Taxonomy maps reuse the
  existing CAMI taxmap schema.
- Read tracking writes a summary with tracked-read percentage,
  density-positive no-hit percentage, and a density/read-length based
  diagnostic estimate of reads from non-reference organisms.

### Changed

- `--readwise-track` forces per-read density units (`--density-block-ctx 0`)
  because block mode cannot preserve exact read ids or context offsets.

### Validation

- `make` and `bash tests/smoke.sh` passed after adding read tracking.

## 2026-06-21

### Added

- Added experimental `minco ani --readwise-ctx-filter fake-prob`, which
  estimates a fake-context probability from the existing breadth/depth plus
  object-difference score and soft-weights nonzero-diff contexts for readwise
  naive ANI diagnostics.
- Direct readwise detail output now appends `Fake_ctx_prob_mean` and
  `Fake_ctx_prob_weighted` when readwise context reliability statistics are
  available.

### Validation

- `make -C minco_core test` passed after adding the fake-context probability
  mode.
- On CAMI Toy Mouse Gut sample0 S=1000, `L. crispatus` target
  `GCF_018987235.1` improved from raw readwise naive ANI `0.934031` and
  hard-filtered ANI `0.962932` to probability-weighted ANI `0.966331` at
  `--readwise-fake-threshold 3.0`, with similar runtime (`1:51.57`) and peak
  RSS (`3.63 GB`).

## 2026-06-20

### Important Progress

- Established the first replicated CAMI II marine profiling win for minco
  readwise metagenome profiling. On marine short-read samples 0-2, the
  `best-diff-unique` shared-context assignment plus ZIP AAF ANI recipe beat
  local Sylph species-level mean F1 (`0.860` vs `0.832`) mainly by reducing
  mean false positives (`41.3` vs `60.7`) at nearly equal recall.
- Found that the practical GTDB-only `S=1000` reference sketch matches or
  slightly exceeds the `S=10000` reference on CAMI marine samples 0-2:
  mean F1 `0.862` vs `0.860`, while reducing mean runtime from `169.8s` to
  `124.5s` and peak RSS from `33.65GB` to `3.70GB`. This makes `S=1000` the
  current benchmark target for readwise profiling, pending validation on a
  different dataset.

### Added

- Added CAMI taxonomic profile output for direct readwise FASTQ abundance ANI:
  `minco ani --abundance-est depth --cami-taxmap MAP.tsv --cami-profile OUT.profile`.
  The profile is generated from the final printed/filtered rows, uses
  `Normalized_abundance_depth * 100` as `PERCENTAGE`, and aggregates rows that
  map to the same `TAXID/RANK/TAXPATH`.
- CAMI taxmaps can now contain multiple rows for the same reference key, one
  per rank, allowing full lineage profile output when the lineage rows are
  supplied.
- Added `minco ani --readwise-profile-only` for large direct FASTQ abundance
  profiling. It skips exact global query-context/query-ref-context sets while
  preserving reference breadth/depth coverage and CAMI profile generation.

### Changed

- Accelerated direct readwise FASTQ density ANI:
  `minco ani -r REF --qraw READS --query-density ref` now processes FASTQ
  reads in parallel 16,384-read batches when OpenMP threads are available.
  Per-thread readwise accumulators, query context sets, reference hit sets, and
  coverage hits are merged after each batch, keeping memory bounded during large
  streamed FASTQ runs.
- Fixed threaded readwise abundance coverage merging so per-reference context
  coverage hits are counted once per hit instead of twice.
- Readwise abundance ANI now stores coverage and best observed object
  difference in one packed per-reference-entry array. Coverage/depth keeps the
  occurrence count, while ANI mutation features are computed once per unique
  reference context entry from the best observed object difference.
- Retuned the automatic readwise abundance report for fixed-size sketches:
  support is now `min(S, max(100, ceil(S/100)))` instead of stale fixed
  `XnY` thresholds, and the default report uses `ANI >= 0.96` with
  `Ref_breadth >= 0.5`.
- `minco sketch --keep` and `minco sketch --remove` now preserve filtered
  `minco.ctxmeta` records and refresh the compact density summary in
  `minco.stat`. Filtered/subset reference sketches therefore remain usable with
  `minco ani --query-density ref`.

### Documentation

- Updated README, user manual, and command-line help text for readwise FASTQ
  density mode, `--density-block-ctx`, progress reporting, and the filtered
  sketch context-metadata behavior.
- Added CAMI II clinical pathogen detection benchmark notes under
  `research/experiments/2026-06-20_cami2_clinical_pathogen_minco/`.

### Validation

- `make test` passed after the readwise batching and ctxmeta filtering changes.
- Large streamed FASTQ analyses completed with bounded memory on the strict
  non-phage viral S1000 reference database.
