# minco ONT yield density estimator accuracy

Date: 2026-06-18
Author/agent: ubuntu
Project: minco
Code commit: NA
minco binary/tool version: minco_stage3_native + exact20m/exact100m

## Question

Does the minco bottom-k density estimator remain accurate for error-prone
Nanopore FASTQ reads across ONT yield downsample levels from 5m to 250m?

## Hypothesis

Even though ONT read errors inflate unique contexts as yield increases, the
bottom-k threshold estimator should remain unbiased for distinct context-hash
cardinality, provided the exact validation sketch does not cap.

## Dataset

- Input data: ONT yield-curve FASTQ files under `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield`.
- Input sketch or matrix: none; sketches generated in this experiment.
- Sample count: one isolate across seven yield levels.
- Selection criteria: sample `2009K-1545__SRR19787631`, present in all yield subdirectories.
- Storage location: source FASTQ under `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/{5m,10m,20m,50m,75m,100m,250m}`.

## Methods

Key parameters:

```text
Normal minco:
  bin/minco_stage3_native sketch -p 4 --conflict --ctxmeta both
  Fixed sketch size = 10,000 context-object entries.

Exact validation:
  5m,10m,20m:
    bin/minco_exact20m sketch -p 2 --conflict --ctxmeta both
    MINCO_SKETCH_SIZE=20000000
  50m,75m,100m,250m:
    bin/minco_exact100m sketch -p 2 --conflict --ctxmeta both
    MINCO_SKETCH_SIZE=100000000

All final exact validation rows had exact_capped=False.
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
Across seven yields:
  MAPE = 0.521%
  mean relative bias = -0.171%
  RMSE relative error = 0.758%
  max absolute relative error = 1.558%

Per-yield relative error:
  5m:   +0.172%
  10m:  +0.926%
  20m:  -0.030%
  50m:  +0.042%
  75m:  -1.558%
  100m: -0.833%
  250m: +0.088%

Estimated unique contexts increased from 4.11M at 5m to 83.17M at 250m,
showing the expected error-driven growth in distinct read contexts with ONT
yield.
```

## Validation

- Checks run: exact-debug sketches with caps 20M or 100M were compared with normal 10k density estimates.
- Expected behavior: final exact validation should not hit the cap.
- Observed behavior: all final rows had exact_capped=False; exact_ctxobj_entries at 250m was 88,729,615 under the 100M cap.
- Failure modes checked: the smaller exact20m build did cap at 50m and above, so those rows were superseded by exact100m.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, matrices, trees, figures, and sketches.

## Conclusion

For this ONT Salmonella yield series, density estimation remains accurate
despite error-prone Nanopore reads. Across 5m-250m, the 10,000-entry sketch
estimated exact distinct context hashes with 0.52% MAPE.

## Paper-Relevant Claim

Exploratory: minco's fixed 10,000-entry bottom-k density estimator can
estimate ONT read distinct context cardinality over a 5m-250m yield curve with
about 0.5% MAPE for this Salmonella sample.

## Caveats

- Tested one isolate only; needs multiple samples from the 50-sample yield folders.
- Exact validation is for distinct 44-bit context hashes, not un-hashed context strings.
- `--conflict` was used to match raw-read sketching; sketch entries can exceed distinct contexts because multiple objects may share one context.
- Disk was tight, so large exact sketches are in `/tmp` and should not be treated as durable.

## Next Experiment

Run normal density sketches for all 50 samples at all seven yields, then exact
validate a stratified subset or implement counter-only exact cardinality to
avoid writing large exact sketches.
