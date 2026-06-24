# CAMI3 ToyGut Current MinCO vs Sylph

Date: 2026-06-22

## Question

Does the current best MinCO S2000 global context markerdb strategy remain competitive with Sylph on the newer CAMI3 ToyGut samples?

## Inputs

- MinCO reference: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- MinCO strategy: product0 top-25% median replacement, `best-diff-split`, naive ANI, active gate from the Toy Mouse sample0 benchmark.
- Sylph database: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`
- Samples: CAMI3 ToyGut sample0, sample1, sample2.
- Truth: bacterial species rows from CAMI taxonomic profiles `taxonomic_profile_0.txt`, `taxonomic_profile_1.txt`, and `taxonomic_profile_2.txt`.

Only bacterial species were scored because both MinCO and the selected Sylph database are GTDB bacterial/archaeal references, while the CAMI3 ToyGut profiles contain substantial viral, plasmid, and eukaryotic abundance.

## Method

MinCO was run on all three read sets with the current best S2000 global context markerdb command:

```text
bin/minco ani -p16 -r <S2000_global_ctx_markerdb> \
  --qraw <reads> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0
```

The offline active gate was:

```text
XnY_ctx >= 15
ANI_naive_calc > 0.95
Reliable_ztp_af >= 0.40
if Reliable_Ref_hit_mean_depth > 3 and Reliable_Ref_hit_depth_variance / Reliable_Ref_hit_mean_depth > 50:
    require ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03
```

Sylph sample0 used the existing profile. Sylph sample1 and sample2 were newly sketched and profiled against the same local GTDB r226 Sylph database.

Predictions from both tools were mapped from GTDB accession to NCBI species taxid using GTDB r232 metadata. CAMI truth was parsed at the bacterial species rank and scored by NCBI species taxid. If multiple references mapped to the same NCBI species taxid, predictions were collapsed to one species.

## Results

Presence/absence against bacterial NCBI species taxids:

```text
sample  method          truth  pred  TP  FP  FN  precision  recall  F1
0       minco_current   118    61    49  12  69  0.803279   0.415254 0.547486
0       sylph           118    70    57  13  61  0.814286   0.483051 0.606383
1       minco_current   125    70    55  15  70  0.785714   0.440000 0.564103
1       sylph           125    77    63  14  62  0.818182   0.504000 0.623762
2       minco_current   127    57    51   6  76  0.894737   0.401575 0.554348
2       sylph           127    66    56  10  71  0.848485   0.440945 0.580311
```

Mean over the three samples:

```text
method          precision  recall    F1        mean TP  mean FP  mean FN
minco_current  0.827910   0.418943  0.555312  51.67    11.00    71.67
sylph          0.826984   0.475999  0.603485  58.67    12.33    64.67
```

Totals over the three samples:

```text
method          TP   FP  FN   predicted_species  truth_species
minco_current  155  33  215  188                370
sylph          176  37  194  213                370
```

Runtime and memory:

```text
MinCO sample0: 2:11.14, peak RSS 4,156,672 KB
MinCO sample1: 2:10.12, peak RSS 4,198,256 KB
MinCO sample2: 2:10.23, peak RSS 4,132,944 KB

Sylph sample1 sketch: 0:54.33, peak RSS 500,804 KB
Sylph sample1 profile: 2:28.92, peak RSS 19,317,948 KB
Sylph sample2 sketch: 0:56.26, peak RSS 500,668 KB
Sylph sample2 profile: 2:27.55, peak RSS 19,216,020 KB
```

Abundance correlation also favored Sylph on these samples. Using bacterial-truth-normalized abundance and renormalized prediction vectors:

```text
sample  method          Pearson   Spearman  MAE percentage points  L1 percentage points
0       minco_current   0.916350  0.533859  0.496323               58.566162
0       sylph           0.954895  0.654911  0.230905               27.246783
1       minco_current   0.836393  0.552008  0.393600               49.199957
1       sylph           0.903143  0.701399  0.192412               24.051467
2       minco_current   0.899290  0.551258  0.361774               45.945296
2       sylph           0.912813  0.648452  0.252527               32.070980
```

## Conclusion

On these CAMI3 ToyGut samples, the current best MinCO S2000 global context markerdb strategy does not beat Sylph. Sylph has nearly the same mean precision but higher recall, giving higher F1 on all three samples:

```text
mean F1: MinCO 0.555312, Sylph 0.603485
```

This weakens the sample0-only conclusion from the CAMI2 Toy Mouse benchmark. The current MinCO gate remains high precision, especially on sample2, but is too conservative for these CAMI3 ToyGut profiles and misses more bacterial species than Sylph.

## Caveats

- This benchmark scores bacterial NCBI species taxids from CAMI profiles, not GTDB species truth. It is appropriate for direct comparison to the CAMI3 profiles, but it is not the same scoring axis as the earlier CAMI2 GTDB-species benchmark.
- The CAMI3 ToyGut profiles contain roughly half non-bacterial abundance. Viral, plasmid, and eukaryotic truth rows were excluded because the MinCO and Sylph references tested here are GTDB references.
- Sylph sample0 was reused from an existing profile; sample1 and sample2 were newly profiled in this run.
- GTDB r232 metadata was used for accession-to-NCBI mapping for both tools, while Sylph used a GTDB r226 database. Assembly-core fallback in the metadata loader handles common GCA/GCF swaps, but version/database differences can still affect mapping.
- The current MinCO active gate was tuned on CAMI2 Toy Mouse sample0. These results suggest it needs cross-sample recalibration before being treated as a default.
