# Universal MinCO Strategy Search

Date: 2026-06-23

## Goal

Find a universal strategy ranked by:

1. presence/absence F1;
2. abundance L1;
3. ANI reporting accuracy.

ANI is treated as a reporting layer and can use a separate estimator from the
call/abundance gate.

## Initial Data

Cached joined unique/split S1000 feature tables were used for the first F1
search:

- marine0-2, CAMI3 toy gut0-2, plant0-2 training panel:
  `/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/train.joined_features.tsv`
- strain0-2 external panel:
  `/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/test.joined_features.tsv`

These feature tables do not contain the newer markerdb robust effective-depth
or adjusted-ANI fields, so this first search is for universal call strategy.

## Candidate Scope

The search uses deployable signals available in both unique and split MinCO
feature tables:

- direct and relaxed rule flags;
- ANI;
- XnY;
- real minimum align fraction;
- breadth;
- ZIP-AAF ANI and ZIP AF;
- depth CV;
- normalized abundance depth.

Dataset labels are not used in candidate rules.

## Result: Active Universal Strategy

The best current all-metrics strategy is layered:

1. Presence gate:
   - Use the previous no-leak RF/HGB average probability model on joined
     unique/split S1000 MinCO candidate features.
   - Use one global probability threshold: `P >= 0.35`.
   - Do not union direct calls into this gate; `direct_or_prob` increased
     false positives and reduced mean F1.
2. Abundance on the S1000 joined-feature panel:
   - Keep the same selected taxa from the `P >= 0.35` gate.
   - Use `s_Normalized_abundance_depth_max` as raw abundance.
   - Renormalize over selected taxa within each sample.
3. Markerdb abundance module:
   - Keep the existing ctx-marker robust effective-depth logic for markerdb
     runs where reliable depth fields are available:

     ```text
     if Reliable_Ref_hit_median_depth >= 20:
         Effective_abundance_depth = Reliable_Ref_hit_median_depth
     else:
         Effective_abundance_depth = Reliable_Ref_mean_depth / Reliable_Ref_zip_af^1.05
     ```

   - Keep the conservative intra-genus winner rescue as an optional markerdb
     rescue rule. It improved CAMI II toy mouse abundance, but did not trigger
     on CAMI3 ToyGut.
4. ANI reporting, independent from call gating:
   - Use coden15 ctxobj96 reporting with `product1-topfrac-median`.
   - Rank suspicious contexts by `(best_diff + 1) * depth`.
   - Median-replace the top 10 percent per reference (`p90` kept).
   - Report adjusted naive ANI when reliable median hit depth is at least 2;
     otherwise report `max(adjusted naive, ZIP-AAF ANI)`.

## Presence Benchmarks

First, a deployable rule sweep over raw joined features did not beat the
previous RF/HGB model. The best simple rule was `u_or_s_direct`:

| method | mean F1 | mean FP | mean FN | marine F1 | toy gut F1 | plant F1 | strain F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| previous RF/HGB avg | 0.639618 | 24.83 | 27.67 | 0.851772 | 0.626735 | 0.614334 | 0.465632 |
| Sylph default | 0.639078 | 26.25 | 29.75 | 0.832332 | 0.616695 | 0.579912 | 0.527374 |
| `u_or_s_direct` | 0.584417 | 24.25 | 34.25 | 0.845092 | 0.506027 | 0.510036 | 0.476514 |

Second, logreg/HGB sample-held-out and dataset-held-out model variants were
tested with base and enriched feature sets. None beat the previous RF/HGB avg
mean F1. The best of this staged pass was:

| method | mean F1 | mean FP | mean FN | min dataset F1 |
|---|---:|---:|---:|---:|
| `loso_enriched_model_hgb` | 0.614120 | 35.50 | 29.25 | 0.526825 |

Finally, retuning the saved previous RF/HGB probabilities with one global
threshold improved the overall F1:

| method | threshold | mean F1 | precision | recall | mean FP | mean FN | marine F1 | toy gut F1 | plant F1 | strain F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RF/HGB avg, global prob | 0.35 | 0.646593 | 0.706213 | 0.624539 | 22.17 | 28.17 | 0.858970 | 0.624104 | 0.629183 | 0.474116 |
| previous RF/HGB avg | native | 0.639618 | 0.685842 | 0.629220 | 24.83 | 27.67 | 0.851772 | 0.626735 | 0.614334 | 0.465632 |
| Sylph default | default | 0.639078 | 0.676871 | 0.627612 | 26.25 | 29.75 | 0.832332 | 0.616695 | 0.579912 | 0.527374 |

Interpretation: `P >= 0.35` is the new mean-F1 winner on the available
12-sample panel. It improves mean F1 over both the previous MinCO strategy and
Sylph. It does not solve the strain-domain weakness: Sylph still has higher
strain F1.

## Abundance Benchmarks

On the same 12-sample S1000 joined-feature panel, using the F1-winning gate
and split normalized depth gave the best tested MinCO abundance L1:

| method | mean L1 | Pearson | Spearman | marine L1 | toy gut L1 | plant L1 | strain L1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `rf_hgb_t035_s_norm` | 65.371844 | 0.785375 | 0.344144 | 29.951621 | 71.007567 | 101.799489 | 58.728697 |
| `rf_hgb_native_threshold_s_norm` | 65.429243 | 0.785379 | 0.336266 | 29.918486 | 71.489201 | 101.507444 | 58.801842 |
| `u_or_s_direct_s_norm` | 68.535351 | 0.750363 | 0.267323 | 30.383358 | 74.392640 | 111.459537 | 57.905869 |

This is a small abundance-L1 improvement over the native RF/HGB thresholds, and
it preserves the better F1. Exact Sylph abundance L1 on this 12-sample panel
was not rescored because the previous accession-to-taxid taxmap cache was lost
from `/tmp`; MinCO abundance could still be scored because the joined features
already contain taxids. The CAMI marine truth profile was restored from CAMI
for this run.

Existing markerdb abundance results remain important:

| dataset | method | mean F1 | mean L1 | note |
|---|---|---:|---:|---|
| CAMI II toy mouse GTDB 0-2 | ctx-marker robust depth + intra-genus rescue | 0.937927 | 1.7945 | prior best MinCO markerdb abundance, better L1 than Sylph 2.287 |
| CAMI3 ToyGut 0-2 | ctx+obj robust depth | 0.601000 | 32.6323 | below Sylph L1 27.7897 |

## ANI Reporting

The ANI reporting layer uses the independent result from
`2026-06-23_readwise_ani_defake_sweep`:

| ANI reporter | n | MAE | RMSE | Pearson | >0.02 err | >0.03 err |
|---|---:|---:|---:|---:|---:|---:|
| coden15 product1 p90, `ANI_median2_else_max` | 75 | 0.006585 | 0.010137 | 0.469986 | 5 | 1 |

This reporter should not be used as the call gate, because the low-depth
`max(adjusted naive, ZIP-AAF)` fallback is intentionally permissive for ANI
reporting.

## Current Conclusion

The current active strategy is:

```text
Call gate:       RF/HGB average probability >= 0.35
Abundance:       selected taxa, raw = split Normalized_abundance_depth, renorm per sample
Markerdb depth:  robust effective depth + optional conservative intra-genus rescue
ANI report:      coden15 product1 p90 adjusted-ANI reporter
```

This is better than the previous MinCO strategy in mean F1 on the available
12-sample panel and slightly better in L1 under the same S1000 abundance
scoring. The unresolved weakness is external-domain transfer, especially
strain-like samples, where Sylph still has better F1.
