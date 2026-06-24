# Full S2000 Marked-Unique Full Benchmark

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory used: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: `main`, dirty working tree
- Working-tree status and diffstat are recorded in `provenance/code_status.txt` and `provenance/code_diffstat.txt`.
- Main script: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/assemble_full_benchmark.py`

## Question

Using the best known MinCO strategies where possible, how does the requested
full S2000 marked-unique strategy compare with the current best MinCO rows and
Sylph for:

1. species presence/absence F1;
2. abundance L1 error;
3. read-based ANI accuracy against source-genome truth?

## Scope And Important Caveat

This benchmark assembles the best available cached results. It does not run a
new integrated MinCO binary that computes marker and full counters inside the
same full S2000 inverted-index scan.

The full S2000 marked-unique row is therefore a measured proxy:

- marker branch: physical ctx-markerdb results;
- full branch: cached full S2000 marker-poor fallback results;
- justification: the prior virtual-index benchmark showed `0` mismatched
  marker-size refs between the full S2000 inverted-index calculation and the
  physical ctx-markerdb for all `200,527` refs.

That makes the marker branch semantically valid as a proxy, but integrated
runtime and implementation behavior still need a dedicated MinCO run after the
dual-counter code exists.

## Datasets

F1 and abundance use the same six-sample mixed panel as the current abundance
calibration note:

```text
mouse_gtdb: samples 0,1,2; GTDB species truth
cami3_ncbi: samples 0,1,2; NCBI species-taxid truth
```

ANI accuracy uses two existing source-aware panels:

```text
CAMI3 ToyGut selected calls:
  read-based ANI versus skani source-genome-to-called-reference ANI

Toy Mouse sample0 source-positive GTDB representatives:
  read-based ANI versus source-to-GTDB-representative ANIm
```

## Compared Strategies

F1 and abundance:

```text
Sylph baseline

MinCO current mixed F1-priority:
  marker_l1_ctxobj_f1_blend + raw_value_sum

MinCO current mixed L1-priority:
  marker_l1_ctxobj_cami_blend + raw_value_sum

MinCO old ctx-marker robust rescue:
  coden11 S2000 physical ctx-markerdb

MinCO full S2000 marked-unique proxy:
  old ctx-marker robust rescue plus marker-poor full fallback
  fallback: ctx_marker_size < 200, XnY_ctx >= 700
  fallback abundance: 0.4 * robust_effective_depth
```

ANI:

```text
MinCO coden15 Ref_zip_aaf_ani
Sylph Adjusted_ANI
MinCO old ctx-marker Ref_zip_aaf_ani
Sylph Naive_ANI
MinCO coden15 product1 p90 ANI_median2_else_max
```

## Results: F1 And Abundance

Six-sample mouse+CAMI3 comparison:

```text
method                                all_F1      all_L1_pp   mouse_F1    mouse_L1    cami3_F1   cami3_L1
Sylph baseline                        0.783546    15.038388   0.963606    2.287032    0.603485   27.789743
MinCO current mixed F1-priority       0.782475    18.014617   0.918887    3.267440    0.646062   32.761794
MinCO current mixed L1-priority       0.770982    17.823364   0.884563    2.989748    0.657402   32.656979
MinCO old ctx-marker robust rescue    0.747536    20.513677   0.939759    1.794503    0.555312   39.232850
MinCO full S2000 marked-unique proxy  0.750078    20.478490   0.944844    1.724130    0.555312   39.232850
```

Winner decisions:

```text
F1:           Sylph = 0.783546, MinCO best = 0.782475
Abundance L1: Sylph = 15.038388 pp, MinCO best = 17.823364 pp
```

The full S2000 marked-unique proxy is a small MinCO-only gain over the old
ctx-marker row:

```text
F1 delta: +0.002543
L1 delta: -0.035187 percentage points
```

The gain is driven by Toy Mouse. CAMI3 is neutral relative to the old ctx-marker
row and remains behind Sylph.

## Results: ANI Accuracy

Primary CAMI3 source/ref ANI comparison using `read_major_source`:

```text
method                                      n    Pearson    Spearman   MAE
MinCO coden15 Ref_zip_aaf_ani               164  0.478258   0.560073   0.007627
Sylph Adjusted_ANI                          165  0.462588   0.526159   0.008038
MinCO old ctx-marker Ref_zip_aaf_ani        146  0.254005   0.487373   0.010709
Sylph Naive_ANI                             165  0.081465   0.206605   0.029642
```

Toy Mouse sample0 source-positive ANIm sweep:

```text
method                                      n   Pearson    MAE       RMSE      >0.02 err  >0.03 err
MinCO coden15 product1 p90 median2/max      75  0.469986   0.006585  0.010137  5          1
MinCO coden15 ZIP-AAF baseline              75  0.421777   0.008147  0.014246  11         6
```

Interpretation: the best ANI reporting signal remains independent from the
abundance/presence gate. Use coden15 ZIP-AAF or the coden15 product1 p90
median2/max reporter for ANI reporting. Do not use the old emitted/naive ANI as
a continuous accuracy signal because previous checks showed it saturates at
`1.0`.

## Decision

The requested full S2000 marked-unique direction is still useful as an
implementation architecture because it can avoid the extra physical markerdb
copy and exactly reproduce marker sizes from the full inverted index.

It is not yet the overall best accuracy strategy:

```text
F1:       Sylph still slightly wins the six-sample mixed panel.
Abund:    Sylph clearly wins the six-sample mixed L1 panel.
ANI:      MinCO coden15 still slightly wins the CAMI3 source/ref ANI subset.
MinCO-only: full S2000 marked-unique proxy slightly improves old ctx-marker, but not enough.
```

Current practical recommendation:

```text
For mixed F1:        use current MinCO mixed F1-priority row only as a near-Sylph candidate.
For mixed abundance: Sylph remains the baseline to beat.
For MinCO abundance experiments: keep ctx-marker robust depth/rescue as first Toy Mouse baseline.
For ANI reporting:  use coden15 Ref_zip_aaf_ani or coden15 product1 p90 median2/max.
For implementation: build full S2000 marked-unique dual counters next, but do not expect it alone to beat Sylph.
```

## Validation

- `python3 -m py_compile` passed for `assemble_full_benchmark.py`.
- Assembler completed with exit status `0`.
- Runtime: `0:00.03` wall clock.
- Peak RSS: `12,032` KB.
- Generated tables:
  - `results/f1_abundance_comparison.tsv`
  - `results/ani_accuracy_comparison.tsv`
  - `results/benchmark_decision.tsv`
  - `summary.tsv`

## Caveats

- The full S2000 marked-unique row is a proxy until MinCO has integrated
  marker/full dual counters in `ani`.
- F1 and abundance use mixed truth namespaces: mouse GTDB species and CAMI3
  NCBI species taxids.
- The six-sample panel is small and partly reused for model/rule selection.
- CAMI3 source/ref ANI uses skani source/ref truth, not full ANIm.
- Toy Mouse ANIm ANI has no Sylph row in the cached sweep, so the MinCO versus
  Sylph ANI comparison relies on the CAMI3 source/ref panel.

## Next Experiment

Implement full S2000 marked-unique dual counters in MinCO `ani`, emitting both:

```text
full_XnY/full_depth/full_breadth
marker_XnY/marker_depth/marker_breadth/marker_size
```

Then rerun the same six-sample F1/abundance panel and the source/ref ANI
benchmark from one integrated binary output.
