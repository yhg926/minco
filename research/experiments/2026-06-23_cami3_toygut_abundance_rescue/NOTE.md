# CAMI3 ToyGut Abundance Rescue Validation

Date: 2026-06-23

## Question

Does the CAMI2 Toy Mouse abundance strategy, robust effective depth plus conservative intra-genus winner rescue, transfer to CAMI3 Toy Human Gut samples 0-2 and beat Sylph?

## Inputs

- Workspace: `/home/ubuntu/yihuiguang/tools/KSSD3mini`.
- MinCO binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/bin/minco`.
- Ctx-marker refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`.
- Ctx+obj marker refdb: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`.
- CAMI3 ToyGut reads:
  - sample0: `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`
  - sample1: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz`
  - sample2: `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz`
- CAMI3 bacterial species truth: `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_{0,1,2}.txt`.
- Sylph profiles from the previous CAMI3 comparison:
  - sample0: `/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv`
  - sample1: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv`
  - sample2: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample2/profile.tsv`

This benchmark uses CAMI bacterial NCBI species-taxid scoring, matching the previous CAMI3 ToyGut comparison. It is not the same axis as the GTDB-species Toy Mouse abundance benchmark.

## Method

Fresh MinCO profiles were generated with the current binary so the outputs contain:

- `Effective_abundance_depth`
- `Normalized_effective_abundance_depth`
- `Reliable_Ref_hit_median_depth`

The MinCO command family was unchanged from the previous CAMI3 comparison except for the binary path:

```text
minco_core/bin/minco ani -p16 -r <markerdb> --qraw <reads> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0
```

The transferred Toy Mouse rescue rule was:

```text
Start from rows passing the current active gate.
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

## Results

Mean abundance metrics over CAMI3 ToyGut samples 0-2, with predictions renormalized to sum to 1 over reported bacterial species:

| method | Pearson | Spearman | L1 percentage points |
|---|---:|---:|---:|
| Sylph | 0.923617 | 0.668254 | 27.789743 |
| MinCO ctx+obj active, robust depth | 0.917313 | 0.565025 | 32.632308 |
| MinCO ctx+obj active, robust depth + intra-genus rescue | 0.917313 | 0.565025 | 32.632308 |
| MinCO ctx+obj active, old normalized depth | 0.894380 | 0.571113 | 46.613987 |
| MinCO ctx-marker active, robust depth | 0.898686 | 0.546870 | 39.232850 |
| MinCO ctx-marker active, robust depth + intra-genus rescue | 0.898686 | 0.546870 | 39.232850 |
| MinCO ctx-marker active, old normalized depth | 0.884011 | 0.545709 | 51.237138 |

Mean presence/absence metrics:

| method | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| Sylph | 58.667 | 12.333 | 64.667 | 0.826984 | 0.475999 | 0.603485 |
| MinCO ctx+obj active, robust depth | 59.000 | 13.667 | 64.333 | 0.812318 | 0.477906 | 0.601000 |
| MinCO ctx-marker active, robust depth | 51.667 | 11.000 | 71.667 | 0.827910 | 0.418943 | 0.555312 |

The intra-genus rescue added zero rows on CAMI3 for both markerdb variants. `rescued_species.tsv` contains only the header.

## Rescue Blocker Check

Rejected same-genus rows existed, but the Toy Mouse pattern did not recur. The rule was blocked mainly by low `XnY_ctx` for high-ratio rows, and by the final 3x same-genus dominance test for rows with enough XnY:

| ref_kind | sample | rejected with active genus and eff>0 | after ANI>=0.999 | after XnY>=100 | after eff>=5 | after zip<=0.25 | after eff>=3x active genus max |
|---|---:|---:|---:|---:|---:|---:|---:|
| ctxmarker | 0 | 1317 | 563 | 9 | 7 | 4 | 0 |
| ctxmarker | 1 | 1509 | 612 | 13 | 9 | 5 | 0 |
| ctxmarker | 2 | 1276 | 515 | 9 | 8 | 5 | 0 |
| ctxobj | 0 | 2696 | 224 | 8 | 4 | 2 | 0 |
| ctxobj | 1 | 3147 | 266 | 8 | 6 | 2 | 0 |
| ctxobj | 2 | 2612 | 234 | 9 | 9 | 6 | 0 |

## Runtime

Fresh MinCO profile runtimes:

| sample | markerdb | elapsed | peak RSS KB |
|---:|---|---:|---:|
| 0 | ctx-marker | 2:12.37 | 4155816 |
| 0 | ctx+obj | 2:02.82 | 5946760 |
| 1 | ctx-marker | 2:11.37 | 4197860 |
| 1 | ctx+obj | 2:03.47 | 6006348 |
| 2 | ctx-marker | 2:11.97 | 4132368 |
| 2 | ctx+obj | 2:04.07 | 5883684 |

## Conclusion

The Toy Mouse abundance rescue does not transfer to CAMI3 ToyGut under the same rule. It makes no calls because CAMI3 rejected same-genus candidates do not satisfy the safe high-support winner pattern.

Robust effective depth still helps MinCO abundance substantially:

- ctx-marker L1 improves from 51.24 to 39.23 percentage points.
- ctx+obj L1 improves from 46.61 to 32.63 percentage points.

But Sylph remains better on CAMI3 abundance:

```text
Sylph L1                 = 27.79 percentage points
best MinCO CAMI3 L1      = 32.63 percentage points
best MinCO CAMI3 method  = ctx+obj active robust effective depth
```

## Caveats

- This is CAMI bacterial NCBI species-taxid scoring, not GTDB species scoring.
- CAMI3 ToyGut contains large non-bacterial truth mass; only bacterial species were scored because the tested MinCO and Sylph references are GTDB bacterial/archaeal references.
- The strict rescue rule was transferred unchanged from Toy Mouse. A different CAMI3-tuned rescue may improve metrics, but that would be a new calibration experiment rather than validation of the Toy Mouse rule.
- An exploratory broad threshold sweep was interrupted because it was too slow; it was not used for the reported metrics.
