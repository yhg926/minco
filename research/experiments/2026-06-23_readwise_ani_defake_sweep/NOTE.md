# Readwise ANI Defake Sweep

## Purpose

Test whether reporting readwise naive ANI after median replacement of suspicious
contexts improves source-genome to GTDB-representative ANI accuracy.

The sweep is focused on top-fraction replacement using product scores:

- product0: `best_diff * depth`, implemented as `product-topfrac-median`
- product1: `(best_diff + 1) * depth`, implemented as `product1-topfrac-median`

The scoring target is the existing toy mouse gut sample 0 source-to-GTDB
representative ANIm table.

## Implementation Note

The normal MinCO code path kept raw readwise ANI in the report, even when
`--readwise-ctx-filter product-topfrac-median` produced filtered support
features. For this experiment only, binaries were rebuilt with
`MINCO_REPORT_FILTERED_READWISE_ANI=1`, so the reported `ANI` field uses the
filtered/replaced feature accumulator.

Normal builds still default to raw reported readwise ANI unless compiled with
that switch.

## Scoring Rules

For each run, the scorer reports:

- `ANI_adjusted_naive`: reported MinCO `ANI` from the adjusted-ANI binary
- `ANI_zip_aaf`: MinCO `Ref_zip_aaf_ani`
- `ANI_median2_else_zip`: use adjusted naive only when reliable median hit
  depth is at least 2; otherwise use ZIP-AAF ANI
- `ANI_median2_else_max`: same median-depth gate, but use max(adjusted naive,
  ZIP-AAF) for low-depth references

The primary error metrics are MAE and RMSE against source-to-representative
ANIm ANI.

## Result

Best observed combo on toy mouse gut sample 0 source-positive GTDB
representatives:

`coden15 + product1-topfrac-median + p90 + ANI_median2_else_max`

This means:

- use coden15 ctxobj96 markerdb
- rank suspicious contexts by `(best_diff + 1) * depth`
- median-replace the top 10% per reference
- for ANI point reporting, use adjusted naive when reliable median hit depth is
  at least 2; otherwise use `max(adjusted naive, ZIP-AAF)`

Top results:

| coden | filter | kept percentile | estimator | MAE | RMSE | Pearson | >0.02 err | >0.03 err |
|---|---|---:|---|---:|---:|---:|---:|---:|
| c15 | product1 | 0.90 | ANI_median2_else_max | 0.006585 | 0.010137 | 0.469986 | 5 | 1 |
| c11 | product1 | 0.90 | ANI_median2_else_max | 0.007539 | 0.012113 | 0.319042 | 8 | 3 |
| c11 | product1 | 0.90 | ANI_adjusted_naive | 0.007735 | 0.012232 | 0.300543 | 8 | 3 |
| c15 | product1 | 0.85 | ANI_median2_else_max | 0.007751 | 0.011410 | 0.314894 | 7 | 1 |
| c15 | product1 | 0.85 | ANI_adjusted_naive | 0.007901 | 0.011399 | 0.333105 | 6 | 1 |
| c15 | ZIP-AAF baseline | n/a | ANI_zip_aaf | 0.008147 | 0.014246 | 0.421777 | 11 | 6 |

The strict rule `median_depth < 2 -> ZIP-AAF` was not best for source-positive
ANI accuracy. It avoids some low-depth saturation, but ZIP-AAF underestimates
several true positives. The `max(adjusted naive, ZIP-AAF)` low-depth fallback
gave the best ANIm agreement in this test. This should be treated as an ANI
reporting rule, not as a species-call gate, because using the maximum can be too
permissive for false-positive gating.

For the Lactobacillus crispatus check:

| run | ANIm | adjusted naive | ZIP-AAF | final max rule | median depth | XnY | breadth | replaced ctx |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| c11 product0 p80 | 0.981000 | 1.000000 | 0.901884 | 1.000000 | 1.0 | 100 | 0.108696 | 21 |
| c11 product1 p90 | 0.981000 | 0.981035 | 0.901884 | 0.981035 | 1.0 | 100 | 0.108696 | 10 |
| c15 product1 p90 | 0.981000 | 0.957539 | 0.980053 | 0.980053 | 1.0 | 73 | 0.068933 | 9 |

Product0 does not solve the saturation. p85 improves over product0 in some
coden15 rows, but p90 is better than p85 for both c11 and c15 in this sweep.
