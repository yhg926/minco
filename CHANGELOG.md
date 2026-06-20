# Changelog

## 2026-06-20

### Changed

- Accelerated direct readwise FASTQ density ANI:
  `minco ani -r REF --qraw READS --query-density ref` now processes FASTQ
  reads in parallel 16,384-read batches when OpenMP threads are available.
  Per-thread readwise accumulators, query context sets, reference hit sets, and
  coverage hits are merged after each batch, keeping memory bounded during large
  streamed FASTQ runs.
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
