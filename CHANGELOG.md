# Changelog

## 2026-06-20

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
