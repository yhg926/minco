# S2000 Ctx+Obj Markerdb on CAMI II Toy Mouse Sample0

Date: 2026-06-22

## Question

Test whether MinCO markerdb built on full context+object entries, matching the original KSSD-style uniqueness logic more closely, performs better than the current ctx-only markerdb on the Toy Mouse sample0 read-wise GTDB species benchmark.

## Markerdb Mode

No MinCO core change was needed for this test. The existing mode already supports the requested behavior:

- `minco set --uniq_union --markerdb` keeps full `ctxobj64` entries that occur in exactly one reference.
- `minco set --uniq_union --markerdb-ctx` is stricter. It keeps entries only when the context is present in exactly one reference, so a context shared by multiple references is removed even if object bits differ.

This experiment tests `--markerdb` against the current best S2000 ctx-only global markerdb baseline.

## Inputs

- Workspace: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Binary: `bin/minco`, version `minco 0.1`
- Binary sha256: `40af778fb7b9a122e2db513a27d652ef52caeffb7196bba4ea1f936258c65f2b`
- Deduplicated S2000 sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Toy Mouse sample0 reads: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- GTDB truth species profile: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/mouse0_gtdb_species_profile.tsv`

## Commands

The runnable command record is in:

`/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse/commands.sh`

Core markerdb build command:

```bash
bin/minco set --uniq_union --markerdb -p 8 \
  -o /tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker \
  /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
```

Toy Mouse profile command:

```bash
bin/minco ani -p16 \
  -r /tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker \
  --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0 \
  -o /tmp/gtdb232_s2000_ctxobj_marker_20260622/toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.tsv
```

## Gate Used For Primary Comparison

The primary comparison reuses the current active product0 strategy:

- `XnY_ctx >= 15`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.40`
- Conditional delta rule: if `Reliable_Ref_hit_mean_depth > 3` and `Reliable_Ref_hit_depth_variance / Reliable_Ref_hit_mean_depth > 50`, require `ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03`

`ANI_from_Reliable_ztp_af = 1 + log(Reliable_ztp_af) / 24`.

## Marker Retention

Ctx+obj markerdb retained many more markers than ctx-only:

| markerdb | refs below 500 | total refs | zero-size refs | min | median | mean | max |
|---|---:|---:|---:|---:|---:|---:|---:|
| ctx-only global markerdb | 9936 | 200527 | not recorded here | not recorded here | not recorded here | 1192.52 | not recorded here |
| ctx+obj markerdb | 2547 | 200527 | 0 | 22 | 1842 | 1724.48 | 2000 |

The ctx+obj build warning reported `2547/200527` references below the 500-marker threshold.

## Results

Primary GTDB species score on Toy Mouse sample0:

| method | predicted taxa | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ctx-only current best recheck | 58 | 58 | 0 | 7 | 1.000000 | 0.892308 | 0.943089 |
| ctx+obj markerdb, active gate | 58 | 56 | 2 | 9 | 0.965517 | 0.861538 | 0.910569 |

The ctx+obj markerdb improved marker retention but degraded species-call F1.

## Error Pattern

Ctx+obj active-gate false positives:

| GTDB species | accession | ANI_naive_calc | XnY_ctx | Reliable_Ref_breadth | Reliable_ztp_af |
|---|---|---:|---:|---:|---:|
| `s__Bacteroides caecimuris_A` | `GCA_948742525.1` | 0.973198 | 20 | 0.022528 | 1.000000 |
| `s__Parabacteroides distasonis_A` | `GCF_004793765.1` | 0.992981 | 21 | 0.015895 | 1.000000 |

Both false positives are low-breadth, low-depth rows where the positive-depth ZTP breadth adjustment saturates to 1.0.

Ctx+obj active-gate false negatives:

- `s__Bacteroides fragilis`
- `s__Bacteroides ovatus`
- `s__Blautia_A wexlerae`
- `s__Lactobacillus crispatus`
- `s__Lactobacillus johnsonii`
- `s__Lactobacillus sp910589675`
- `s__Limosilactobacillus kinnaridis`
- `s__Limosilactobacillus reuteri_L`
- `s__Prevotella sp900115435`

## Small Gate Retune

A small sweep over `XnY_ctx`, `Reliable_ztp_af`, `Ref_zip_af`, and `Reliable_Ref_breadth` found the best simple retune:

- `XnY_ctx >= 15`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.35`
- `Reliable_Ref_breadth >= 0.03`
- Same conditional delta rule as the active gate

This gives `TP=57`, `FP=0`, `FN=8`, `F1=0.934426`. It removes the two false positives and recovers one true species, but still does not beat the ctx-only current best (`F1=0.943089`).

## Abundance

On the selected truth species, ctx+obj had slightly better renormalized abundance error but worse unnormalized mass recovery:

| method | renormalized | Pearson | Spearman | MAE percentage points | L1 percentage points | predicted sum on truth |
|---|---:|---:|---:|---:|---:|---:|
| ctx-only current best | false | 0.998262 | 0.888085 | 0.461136 | 29.973822 | 0.705454 |
| ctx-only current best | true | 0.998262 | 0.888085 | 0.155189 | 10.087294 | 1.000000 |
| ctx+obj active gate | false | 0.998378 | 0.889160 | 0.580205 | 37.713333 | 0.628046 |
| ctx+obj active gate | true | 0.998378 | 0.889160 | 0.152421 | 9.907339 | 0.999936 |

The main species-call objective is worse for ctx+obj despite this small renormalized abundance improvement.

## Runtime

| step | elapsed | max RSS KB | exit |
|---|---:|---:|---:|
| build ctx+obj markerdb | 1:06.22 | 6319360 | 0 |
| index ctx+obj markerdb | 0:07.80 | 6808064 | 0 |
| Toy Mouse readwise profile | 2:04.33 | 6002996 | 0 |

The indexed ctx+obj markerdb directory used about `6.6G` on disk.

## Conclusion

Ctx+obj markerdb is better for marker retention, but it is not better for Toy Mouse sample0 GTDB species detection under the current gate. The current best MinCO S2000 markerdb remains the ctx-only global markerdb with active product0 gate: `F1=0.943089`, `TP=58`, `FP=0`, `FN=7`.

## Caveats

- This is one CAMI II Toy Mouse sample. It should not be generalized to all samples without additional replication.
- The ctx+obj result may need a different gate because the marker population and depth distribution changed.
- The active ZTP adjustment can overcorrect very low-depth, low-breadth rows; the two ctx+obj false positives are examples.
- Large generated artifacts are under `/tmp` and may need to be regenerated if `/tmp` is cleaned.
