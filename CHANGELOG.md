# Changelog

## 2026-07-01

### Added

- Added hidden C readwise density-cache replay support:
  `--readwise-density-cache-out` writes retained per-read density vectors and
  `--readwise-density-cache-in` replays them for exact split evaluation without
  rereading raw input. The Python default wrapper uses this for lazy exact
  split abundance/call replay.
- Added `--no-profile-rescue` to the user-facing `scripts/minco_profile`
  launcher. It keeps the candidate preset but disables profile-rescue
  additions for audit/replay runs.

### Changed

- Promoted the calibrated `scripts/minco_profile` candidate preset to the
  strict split-evidence rescue switch
  `split-p002-x300-ani095-af06-b025-d1-top1`.
- `universal-auto-exact` now also has an abundance-only exact split trigger at
  `--exact-split-abundance-trigger 0.30`, selected to avoid the raw validation
  regression observed at `0.10`.
- The candidate preset now injects the validated abundance-only reliability
  guards: ANI floor `0.90`, sparse-depth cap `poisson-breadth`, sparse breadth
  max `0.15`, and sparse depth-ratio min `200`.
- Strict rescued rows keep native profile raw abundance mass. Normalized-depth
  candidate mass is still used for older loose candidate additions and
  candidate-surface rows.

### Validation

- `make minco`, `bash tests/smoke.sh`, `bash tests/full_cli.sh`, and the
  Python raw default smoke passed after adding density-cache replay.
- `scripts/minco_profile --check-ref` passed on the packaged GTDB reference
  sidecars, confirming default scope, species taxmap, model cache metadata,
  candidate-surface taxmap, strict rescue switch, and abundance guards are
  discovered without expert options.

## 2026-06-26

### Added

- Added opt-in `--candidate-surface-switch accession-ani93-xny650-br20` to
  `scripts/minco_profile_calibrated.py`. It implements a stricter raw-retention
  candidate-surface mode for validation runs and does not change the selected
  default preset.
- Added validation-only `--candidate-surface-switch
  accession-current-or-ani93-xny650-br20`, which preserves the current
  candidate surface and adds the stricter raw-retention surface for replay
  audits. It is not selected by the default preset.
- Added opt-in `--candidate-surface-max-called-species`, an output-derived
  guard that skips candidate-surface additions when the pre-surface called
  species count exceeds a supplied threshold. The default `0` disables it.
- Added `--strategy {probability,adaptive-sub95,universal,universal-auto-exact}` to
  `scripts/minco_profile_calibrated.py`. The experimental `universal` strategy
  uses adaptive+sub95 calls, a high-confidence raw-unique fallback, and panel
  zip-corrected depth abundance over the final retained species.
- Added experimental `universal-auto-exact` wrapper mode. It evaluates the
  normal block-mode universal gate first, then reruns the split pass with
  `--density-block-ctx 0` when `probability_extra_mass_ratio <= 0.10` by
  default.
- Added `--abundance-genus-xny-blend-alpha` as an experimental abundance-only
  diagnostic for calibrated universal strategies. The default is `0.0`, so the
  selected calibrated abundance rule is unchanged; alpha 0.25 was rejected as
  the default after four-panel fixed-call validation.
- Added `--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002`
  as an experimental abundance-only candidate. It is off by default, preserves
  candidate-added row mass, and requires independent holdout validation before
  any default promotion.

### Changed

- Clarified the public profiling entrypoints: `scripts/minco_profile_default.py`
  is the recommended no-manual-strategy launcher for calibrated species
  profiling with `universal-auto-exact`, while `minco profile` remains the
  conservative direct path for species/AMR/virus/gene/mixed-domain profiling
  and read tracking.
- Promoted `scripts/minco_profile_calibrated.py` default `--strategy` from
  legacy `probability` to guarded `universal-auto-exact`. Use
  `--strategy probability` to reproduce the legacy RF/HGB threshold-only gate.
- Fixed calibrated-profile output so `profile_strategy` reports the actual
  selected/default strategy (`universal-auto-exact`) instead of the internal
  block-mode sub-strategy (`universal`) after auto-exact evaluation.

### Validation

- Added a regression check that the implicit no-`--strategy` calibrated wrapper
  output matches explicit `--strategy universal-auto-exact`.

## 2026-06-25

### Added

- Added `minco profile`, a simplified read-profiling subcommand over the
  direct readwise engine. It sets reference-density extraction, depth
  abundance, profile-only mode, `best-diff-split` assignment, naive readwise
  ANI, and `product-topfrac-median` context defaking by default.
- `minco profile` supports CAMI profile output, read tracking, GTDB/NCBI
  tracking taxmaps, `--report-all` for benchmark/debug tables, and
  `--dual-evidence` for optional Marker_* columns from full-index references.

### Validation

- `make minco`, `bash tests/smoke.sh`, and `bash tests/full_cli.sh` passed
  after adding `profile`.

## 2026-06-24

### Added

- Added `minco ani --readwise-track FILE`, a Kraken-like read tracking sidecar
  for direct `--qraw --query-density ref` readwise mode. It writes per-read
  selected reference hits, selected context offsets from the read start, target
  reference counts/lists, and optional GTDB/NCBI LCA labels.
- Added `minco set --set-domain species|amr|virus|gene|default`, which writes
  a `minco.domain` sidecar so domain-specific readwise reporting rules survive
  downsample, append/merge, keep/remove, and markerdb operations.
- Documented AMR/gene-panel reporting semantics: normal MinCO AMR runs should
  report detected determinants from the installed refdb, while `not_in_refdb`
  rows are benchmark/comparator-only and require an external truth list.
- Added `--readwise-taxonomy none|gtdb|ncbi|both`, `--gtdb-taxmap`,
  `--ncbi-taxmap`, and `--readwise-track-summary`. Taxonomy maps reuse the
  existing CAMI taxmap schema.
- Read tracking writes a summary with tracked-read percentage,
  density-positive no-hit read percentage, and context-level reference-hit /
  reference-absent percentages from sampled read contexts.
- Read tracking summaries now include `estimated_unknown_reads_pct` when
  `minco.ctxmeta` is available. This uses per-reference sketch densities to
  Horvitz-Thompson-correct S2000-style sampled hits toward a full-context
  reference-absence estimate, instead of reporting the raw sampled-sketch absent
  percentage as the user-facing unknown fraction.

### Changed

- `--readwise-track` forces per-read density units (`--density-block-ctx 0`)
  because block mode cannot preserve exact read ids or context offsets.
- Renamed the previous read-level no-hit diagnostic in the tracking summary so
  it is not confused with whole-genome non-reference context absence.
- Direct readwise FASTQ density mode now keeps completed output when all
  records are parsed but `gzclose` fails during the final gzip integrity check;
  MinCO emits a warning instead of discarding the finished result. Malformed
  FASTQ records still fail with the underlying `kseq_read` error.
- Direct readwise `Default_call` now uses AMR/gene-specific `major` and
  `screening` gates for references tagged in `minco.domain`; untagged,
  species, and virus references retain the existing species-style rules.

### Validation

- `make`, `bash tests/smoke.sh`, and `bash tests/full_cli.sh` passed after
  adding read tracking and domain-specific reporting profiles.

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
- Changed the default species/prokaryotic readwise abundance report ANI cutoff
  from `0.96` to `0.95`, matching the standard prokaryotic species boundary.
  AMR/gene-specific `major` and `screening` thresholds are unchanged.
- Retuned the automatic readwise abundance report for fixed-size sketches:
  support is now `min(S, max(100, ceil(S/100)))` instead of stale fixed
  `XnY` thresholds, and the default report uses `ANI >= 0.95` with
  `Ref_breadth >= 0.5`.
- `minco sketch --keep` and `minco sketch --remove` now preserve filtered
  `minco.ctxmeta` records and refresh the compact density summary in
  `minco.stat`. Filtered/subset reference sketches therefore remain usable with
  `minco ani --query-density ref`.
- Added experimental exact-split sidecar generation:
  `minco ani --readwise-exact-split-out FILE` can write exact per-read
  best-diff-split evidence while the main split pass uses block mode, and
  `scripts/minco_profile_calibrated.py --same-stream-exact-split` can reuse
  that sidecar when `universal-auto-exact` fires.
- Added `scripts/minco_profile_calibrated.py --exact-split-low-extra-mode`.
  The default `skip` mode avoids the expensive exact split rerun when the
  block-mode low-extra split rescue already added candidates; `allow` restores
  the older exact behavior.

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
