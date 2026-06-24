# CAMI II Toy Mouse More Samples With GTDB Truth

## Question

Does the current Toy Mouse sample0 GTDB-species result generalize to more CAMI
II Toy Mouse Gut samples after transferring CAMISIM source genomes to GTDB r232
species?

## Data

Downloaded two additional official CAMI II Toy Mouse Gut short-read archives from
the PUBLISSO CAMI data host:

- `2017.12.29_11.37.26_sample_1_reads.tar`
- `2017.12.29_11.37.26_sample_2_reads.tar`

Each archive contains:

- `anonymous_reads.fq.gz`
- `reads_mapping.tsv.gz`
- `taxonomic_profile_N.txt`
- `distributions/distribution_N.txt`

Local extracted read paths:

- `/mnt/new3T/minco_cami2_toymouse_20260621/sample_1/2017.12.29_11.37.26_sample_1/reads/anonymous_reads.fq.gz`
- `/mnt/new3T/minco_cami2_toymouse_20260621/sample_2/2017.12.29_11.37.26_sample_2/reads/anonymous_reads.fq.gz`

## GTDB Truth

Truth was built from CAMISIM distributions plus `setup/internal/genome_locations.tsv`,
mapping each source genome accession to GTDB r232 metadata. Species abundance is
renormalized over source genomes that map to GTDB species.

| sample | positive genomes | mapped | unmapped | GTDB species |
| --- | ---: | ---: | ---: | ---: |
| 0 | 75 | 75 | 0 | 65 |
| 1 | 90 | 89 | 1 | 80 |
| 2 | 87 | 86 | 1 | 73 |

Unmapped source caveat:

- sample1: `GCF_900113595.1`, NCBI `Lachnospiraceae bacterium NLAE-zl-G231`, raw abundance `1/4567`.
- sample2: `GCF_000245055.1`, NCBI `Desulfovibrio sp. U5L`, raw abundance `219/4042`.

## Methods

Scored three profilers against GTDB species truth:

- `ctx_only_current_best`: old global ctx-only S2000 markerdb with the current active product0 gate.
- `ctxobj_markerdb_product0_active`: ctx+obj S2000 markerdb with the same gate.
- `sylph_gtdb_profile`: Sylph GTDB r226 c200 profile mapped to GTDB species.

MinCO active gate:

- `XnY_ctx >= 15`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.40`
- if `Reliable_Ref_hit_mean_depth > 3` and `Reliable_depth_vmr > 50`, require `ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03`

## Results

Per-sample GTDB-species presence/absence:

| sample | method | TP | FP | FN | F1 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0 | ctx-only | 58 | 0 | 7 | 0.943089 |
| 0 | ctx+obj | 56 | 2 | 9 | 0.910569 |
| 0 | Sylph | 59 | 2 | 6 | 0.936508 |
| 1 | ctx-only | 70 | 2 | 10 | 0.921053 |
| 1 | ctx+obj | 69 | 6 | 11 | 0.890323 |
| 1 | Sylph | 78 | 2 | 2 | 0.975000 |
| 2 | ctx-only | 66 | 0 | 7 | 0.949640 |
| 2 | ctx+obj | 67 | 2 | 6 | 0.943662 |
| 2 | Sylph | 71 | 1 | 2 | 0.979310 |

Mean across samples:

| method | mean precision | mean recall | mean F1 |
| --- | ---: | ---: | ---: |
| ctx-only | 0.990741 | 0.890472 | 0.937927 |
| ctx+obj | 0.952177 | 0.880616 | 0.914851 |
| Sylph | 0.976108 | 0.951765 | 0.963606 |

Pooled totals across samples:

| method | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| ctx-only | 194 | 2 | 24 | 0.937198 |
| ctx+obj | 192 | 10 | 26 | 0.914286 |
| Sylph | 208 | 5 | 10 | 0.965197 |

## Abundance

MinCO abundance remains under-calibrated in raw mass. For ctx-only, predicted
mass on truth species was `0.705`, `0.565`, and `0.668` for samples 0-2. Sylph
was much closer: `0.9995`, `0.9992`, and `0.9481`.

After renormalization, MinCO abundance correlations are still weaker than Sylph:

- ctx-only sample1 renormalized MAE: `0.2667` percentage points.
- ctx-only sample2 renormalized MAE: `0.2685` percentage points.
- Sylph sample1 renormalized MAE: `0.0124` percentage points.
- Sylph sample2 renormalized MAE: `0.0720` percentage points.

## Interpretation

The sample0 conclusion does not generalize as a win over Sylph. The old global
ctx-only markerdb is still the best MinCO variant on these Toy Mouse samples,
but Sylph has much higher recall on samples 1 and 2 and better abundance
calibration.

The ctx+obj markerdb again does not improve Toy Mouse performance. It slightly
rescues one true species on sample2, but adds enough false positives and sample1
misses that its mean F1 is lower than ctx-only.

The dominant MinCO failure mode is recall from the reliable adjusted breadth
gate: most MinCO false negatives fail `Reliable_ztp_af >= 0.40`. The conditional
delta rule explains only a small number of misses, including `s__Lactobacillus
amylovorus` in sample2. This suggests the current breadth adjustment is still
too strict or not abundance-aware enough for low-to-mid abundance species.

Runtime and memory were stable:

- ctx-only sample1/sample2: about 2.1-2.3 minutes, about 4.5 GB RSS.
- ctx+obj sample1/sample2: about 2.1 minutes, about 6.3 GB RSS.
- Sylph profile peak RSS was about 19.5 GB.
