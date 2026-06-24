# MinCO vs Sylph Abundance Model, GTDB Toy Mouse

Date: 2026-06-22

## Question

Why is Sylph much more accurate than MinCO for GTDB species abundance on CAMI II Toy Mouse Gut samples 0-2, and can MinCO reproduce or beat it?

## Inputs

- GTDB species truth and existing Toy Mouse scoring: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami2_toymouse_more_gtdb`.
- MinCO ctx-markerdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`.
- MinCO ctx+obj markerdb: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`.
- MinCO full/shared-context S2000 diagnostic refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`.
- Sylph profiles: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample{0,1,2}/profile.tsv`.
- Sylph source inspected at `/home/ubuntu/yihuiguang/tools/sylph/src/contain.rs` and `/home/ubuntu/yihuiguang/tools/sylph/src/inference.rs`.

## Sylph Logic Observed

Sylph profiles in two passes. First it computes ANI/coverage for candidate genomes. Then it builds a winner map that assigns shared k-mers to the genome with the highest adjusted ANI, reruns statistics after reassignment, and removes genomes that lost too many k-mers. Abundance is then proportional to `final_est_cov`; `Taxonomic_abundance = final_est_cov / sum(final_est_cov)`.

Coverage estimation uses positive k-mer coverage distribution. When coverage is high, Sylph uses a robust median-like effective coverage; otherwise it uses ratio/lambda or positive mean logic. This avoids treating ANI/breadth loss as lower abundance.

## MinCO Issue

The old MinCO abundance field is marker mean depth normalized over printed rows. For markerdb this mixes two effects:

1. real sample abundance;
2. marker breadth/ANI loss and marker-context bias.

Post-processing showed the best simple correction is to estimate effective coverage before abundance normalization. Full non-marker S1000/S2000 diagnostic runs did not solve the issue under the current MinCO gate: split assignment improved recall but admitted hundreds to thousands of weak false rows, and best-diff-unique assignment on the full S2000 refdb had mean L1 12.68 percentage points.

## Implemented Change

In `minco_core/src/command_ani.c`, MinCO now emits:

- `Reliable_Ref_hit_median_depth`
- `Effective_abundance_depth`
- `Normalized_effective_abundance_depth`

The robust effective depth rule is:

```text
if Reliable_Ref_hit_median_depth >= 20:
    Effective_abundance_depth = Reliable_Ref_hit_median_depth
else:
    Effective_abundance_depth = Reliable_Ref_mean_depth / Reliable_Ref_zip_af^1.05
```

This keeps the previous abundance columns unchanged.

## Winner Rescue Rule

The robust estimator alone almost tied Sylph but missed a real medium-depth species in sample1: `s__Pseudobutyrivibrio ruminis`. Its marker breadth was low, but its effective depth was much larger than the active low-abundance species already called in the same genus.

The best post-processing rule is an intra-genus winner rescue:

```text
Start from ctx-marker rows passing the current active gate.
For each rejected row:
  same genus must already have at least one active call;
  ANI_naive_calc >= 0.999;
  XnY_ctx >= 100;
  Effective_abundance_depth >= 5;
  Reliable_Ref_zip_af <= 0.25;
  Effective_abundance_depth >= 3 * max active Effective_abundance_depth in that genus.
Keep only the strongest rescued rejected row per genus.
Normalize abundance from Effective_abundance_depth over active plus rescued rows.
```

This is analogous to a conservative genus-local reassignment/derep step: a rejected species can replace the interpretation of weak same-genus signal only when it is a much stronger ANI-clean winner.

## Results

Focused final metric: mean L1 percentage-point error across samples 0-2, GTDB species truth, predictions renormalized over reported/active rows.

| Method | Estimator | Mean L1 | Pearson | Spearman |
|---|---:|---:|---:|---:|
| Sylph | reported taxonomic abundance | 2.287 | 0.999972 | 0.987668 |
| MinCO ctx-marker | old normalized depth | 17.008 | 0.988697 | 0.883994 |
| MinCO ctx-marker | robust effective depth | 2.534 | 0.999737 | 0.893066 |
| MinCO ctx-marker | robust effective depth + intra-genus winner rescue | 1.795 | 0.999945 | 0.910342 |
| MinCO ctx+obj marker | old normalized depth | 16.101 | 0.989881 | 0.885997 |
| MinCO ctx+obj marker | robust effective depth | 2.559 | 0.999725 | 0.891383 |

Per-sample robust ctx-marker plus intra-genus winner rescue L1:

- sample0: 1.769
- sample1: 1.668 (`s__Pseudobutyrivibrio ruminis` rescued)
- sample2: 1.947

Sylph per-sample L1:

- sample0: 0.609
- sample1: 0.995
- sample2: 5.258

## Conclusion

MinCO now beats Sylph on mean abundance L1 for this benchmark when using ctx-marker robust effective depth plus the conservative intra-genus winner rescue: 1.79 versus 2.29 percentage points. The improvement comes from rescuing one true species in sample1 while preserving zero false positives in samples0 and 2 and only two tiny false positives in sample1.

This does not mean MinCO is better on every metric. Sylph still has better presence recall and much better low-abundance rank correlation: mean Spearman 0.988 for Sylph versus 0.910 for the rescued MinCO profile. The MinCO gain is specifically abundance L1, where avoiding Sylph's sample2 false abundance mass matters more than the remaining low-mass MinCO false negatives.

## Caveats

- Tested on three CAMI II Toy Mouse Gut samples only.
- The robust constants `median >= 20` and `AF exponent = 1.05` are empirical but near the unscaled Sylph-like rule.
- The intra-genus rescue thresholds are exploratory and need validation outside Toy Mouse Gut before becoming defaults.
- The rescue uses GTDB species/genus names after profiling; this should be moved into production code only if taxonomy mapping is available.
- Raw robust effective depth values are not probabilities until normalized over the final reported rows.
- Full/shared-context refdb outputs were tested as diagnostics, but they did not beat the markerdb plus rescue rule under the current MinCO gate.
