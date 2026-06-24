# threshold_combo_search_all_metrics

Date: 2026-06-23
Author/agent: ubuntu
Project: KSSD3mini

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Commit: `NA`
- Ref/describe: `NA`
- Working tree: `NA`
- Binary/script/tool: `minco-dev`
- Provenance files: `provenance/code_status.txt`

## Command and Parameter Provenance

- Working directory: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Exact commands: `commands.sh`
- Parameters, arguments, environment settings, seeds, thresholds, and config files: `parameters.tsv`
- Rule: commands must be complete enough to rerun without consulting chat history.

## Question

Can threshold search over the current MinCO S2000 rows find a universal rule
that improves presence F1 while also accounting for abundance L1 and available
ANI accuracy, across CAMI II toy mouse GTDB truth and CAMI3 ToyGut NCBI truth?

## Hypothesis

The previous active gate was conservative for CAMI3. A lower support threshold
and AF floor, plus robust effective depth and controlled ctx+obj add-back, may
raise CAMI3 recall enough to improve the overall multi-metric result without
destroying mouse abundance.

## Dataset

- Input data: cached MinCO readwise result TSVs and Sylph profiles from earlier
  S2000 experiments.
- Mouse samples: CAMI II toy mouse samples 0-2, scored against GTDB species truth.
- CAMI3 samples: ToyGut samples 0-2, scored against bacterial NCBI species truth.
- Ref modes:
  - coden11 S2000 ctx-marker plus full-sketch fallback for marker-poor refs.
  - coden11 S2000 ctx+obj marker rows.
- ANI truth: mouse sample0 source-genome to GTDB representative ANIm table,
  available for matched positive GTDB species only.
- Availability status: cached local files available on 2026-06-23; some large
  inputs live in `/tmp` and `/mnt/new3T`.

## Methods

Parameter summary:

```text
Primary search: 5000 focused combos for ctx-marker/full and 5000 focused combos
for ctx+obj.
Focused grid included ANI threshold 0.945-0.97, XnY threshold 5-20, fixed AF
0.25-0.35, formula AF from effective context lengths 22-30, delta gate VMR
10-100, delta max 0.015-0.04, robust-depth median cutoff 10-30, abundance AF
exponent 1.0-1.05, and strict/moderate/off rescue.

Post-search ensemble: ctx-marker/full default plus ctx+obj add-back with
absolute/relative add-back thresholds and abundance scaling 0.05-1.0.
```

Commands are recorded in `commands.sh`; detailed parameters are recorded in `parameters.tsv`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
method                                      all_F1    all_L1    mouse_F1  mouse_L1  cami3_F1  cami3_L1
Sylph                                       0.783546  15.038388 0.963606  2.287032  0.603485  27.789743
ctx-marker/full best F1                     0.777458  18.493025 0.939060  2.892808  0.615856  34.093242
ctx-marker/full best L1 with F1>=0.775      0.775173  18.061240 0.938907  2.947439  0.611439  33.175042
ctx+obj best F1                             0.775598  18.268764 0.914134  3.747801  0.637062  32.789727
ensemble scaled add-back best F1            0.782475  18.034939 0.918887  2.964606  0.646062  33.105272
ensemble scaled add-back best L1 F1>=0.78   0.780616  18.024015 0.934958  2.953692  0.626275  33.094338
```

## Validation

- Checks run: `python3 -m py_compile` for the scorer; smoke sweep with 1000
  combos; focused ctx-marker/full and ctx+obj sweeps; ensemble add-back grids.
- Expected behavior: previous known rules should be reproduced in the summary.
- Observed behavior: known current and hybrid rules are present; focused search
  improved MinCO overall F1 but did not beat Sylph when L1 is included.
- Failure modes checked: high-CAMI3 ctx+obj add-back was tested separately; it
  did not close the remaining overall F1 gap.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, matrices, trees, figures, and sketches.

Per-experiment outputs should go in `results/` when they are small enough to keep here. Large outputs may remain external, but must be listed in `artifacts.md` with availability and retention status.

## Conclusion

The best searched MinCO rule is the scaled ensemble add-back:
ctx-marker/full marker-L1 gate plus ctx+obj add-back with all ctx+obj new calls
included at abundance scale 0.2. It reaches all-sample F1 0.782475, close to
Sylph 0.783546, and beats Sylph on CAMI3 F1 (0.646062 vs 0.603485). However it
still loses badly on abundance L1 (18.034939 vs Sylph 15.038388) and mouse F1
(0.918887 vs Sylph 0.963606). Therefore this is not yet the universal default.

For a conservative all-metric default, the prior ctx-marker/full rule remains
better for mouse abundance, while ctx+obj is a useful CAMI3 recall/ANI diagnostic
source rather than a default replacement.

## Paper-Relevant Claim

Exploratory: ctx+obj add-back can recover CAMI3 presence recall beyond Sylph,
but current MinCO abundance normalization and mouse precision are not yet good
enough for a universal all-metric default.

## Caveats

- The search used cached result rows, not a full MinCO rerun for every strategy.
- ANI was evaluated only on available mouse sample0 source-to-representative
  ANIm matches; CAMI3 true source ANI is only partially available.
- The ensemble add-back is a post-hoc scorer prototype, not yet implemented in
  MinCO.
- CAMI3 truth here is NCBI species bacterial truth; mouse truth is GTDB species.

## Next Experiment

Target abundance normalization for ensemble calls. The main remaining gap is L1,
not F1: ctx+obj add-back recovers presence calls, but abundance mass assignment
still over- or mis-weights added species.
