# KSSDmini FASTQ density estimator accuracy

Date: 2026-06-18
Author/agent: ubuntu
Project: KSSD3mini
Code commit: NA
KSSD3A binary/tool version: kssd3mini_stage3_native + kssd3mini_exact20m

## Question

How accurately does the KSSDmini bottom-k density estimator predict unique
context count on FASTQ read input?

## Hypothesis

For FASTQ read sketches using `--conflict`, the estimator
`observed_ctx_under_threshold * 2^44 / (threshold + 1)` should remain accurate
for distinct context hashes, but conflict/context multiplicity means
`sketch_entries` can be 10,000 while distinct selected contexts are lower.

## Dataset

- Input data: `/mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq`.
- Input sketch or matrix: none; sketches generated from FASTQ.
- Sample count: one FASTQ source, evaluated as first 10k, 50k, 100k, and 200k read subsets.
- Selection criteria: first N reads of the R1 FASTQ.
- Storage location: source FASTQ under `/mnt/new3T/kssd3test/patmg_CAMI2`.

## Methods

Key parameters:

```text
Normal estimator:
  bin/kssd3mini_stage3_native sketch -p 4 --conflict --ctxmeta both
  fixed sketch size = 10000 context-object entries

Exact validation:
  bin/kssd3mini_exact20m sketch -p 2 --conflict --ctxmeta both
  KSSD3MINI_SKETCH_SIZE=20000000
  All exact subset runs had exact_ctxobj_entries < 20,000,000, so no cap hit.

Full FASTQ check:
  bin/kssd3mini_stage3_native sketch -p 8 --conflict --ctxmeta both
  Exact count not computed for the full 846 MB FASTQ because the estimate was
  about 132M unique contexts, above the 20M debug cap.
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
Subset accuracy summary:
  10k reads:  exact 526,144 unique contexts, estimate 520,874, error -1.002%
  50k reads:  exact 2,473,699 unique contexts, estimate 2,465,496, error -0.332%
  100k reads: exact 4,809,687 unique contexts, estimate 4,801,715, error -0.166%
  200k reads: exact 9,315,040 unique contexts, estimate 9,380,419, error +0.702%

Across four subsets:
  MAPE = 0.550%
  mean relative bias = -0.199%
  relative RMSE = 0.639%

Full FASTQ normal sketch:
  input = 3,444,570 reads, 846 MB
  wall time = 5.64 s
  max RSS = 872,728 KB
  selected_estimated_unique_ctx = 132,218,033
```

## Validation

- Checks run: exact large-sketch sidecar compared with normal 10k density sidecar.
- Expected behavior: exact subset `sketch_entries` should be below 20M.
- Observed behavior: max exact subset `exact_ctxobj_entries` was 10,668,319 at 200k reads.
- Failure modes checked: `--conflict` used consistently for normal and exact sketches.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, matrices, trees, figures, and sketches.

## Conclusion

On these FASTQ read subsets, the density estimator predicted exact distinct
context hashes with sub-1% average error. This supports using the ctxmeta
density estimate as a practical unique-context/read-complexity estimate for
FASTQ sketches.

## Paper-Relevant Claim

Exploratory: KSSDmini's 10,000-entry bottom-k density estimator estimates FASTQ
distinct context cardinality with about 0.55% MAPE on four subsets from one real
R1 read file.

## Caveats

- Tested on one R1 FASTQ source and prefix subsets, not independent FASTQ datasets.
- Exact validation was only run up to 200k reads; full 846 MB exact count was not computed.
- The estimate is for distinct 44-bit context hashes, not un-hashed context strings.
- `--conflict` keeps multiple objects for the same context, so `sketch_entries` is not the same as distinct observed contexts.

## Next Experiment

Implement a counter-only exact context cardinality mode for FASTQ, then validate
full read files and multiple datasets without writing huge exact sketches.
