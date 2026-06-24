# CAMI3 ToyGut Coden15 Formula-AF Transfer Test

Date: 2026-06-23

## Question

Does the Toy Mouse relative abundance ranking reproduce on CAMI3 ToyGut?

Toy Mouse reference result:

| method | L1 pp | Pearson |
|---|---:|---:|
| old ctx-marker robust rescue | 1.7945 | 0.999945 |
| coden15 formula-AF gate | 2.0091 | 0.999936 |
| Sylph | 2.2870 | 0.999972 |
| coden15 original gate | 2.4022 | 0.999917 |

## Inputs

- CAMI3 ToyGut samples 0-2 reads and taxonomic profiles from `/mnt/new3T/minco_cami3_toygut_20260620` plus extra sample read paths from `/mnt/new3T/minco_cami3_toygut_extra_20260621`.
- Coden15 ctx-marker refdb: `/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker`.
- Existing CAMI3 old ctx-marker and Sylph baseline: `research/experiments/2026-06-23_cami3_toygut_abundance_rescue/summary.tsv`.
- CAMI3 scoring axis: bacterial NCBI species taxid, mapped from GTDB metadata. This is not GTDB species-name truth.

## Coden15 Gates

Both coden15 methods used:

- `--readwise-assign best-diff-split`
- `--readwise-ani naive`
- `--readwise-ctx-filter product-topfrac-median`
- `--readwise-fake-threshold 0.25`
- `Reliable_ztp_af` as depth-adjusted AF
- the same intra-genus winner rescue template as the Toy Mouse abundance model

Original gate:

- `XnY_ctx >= 15`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.40`
- if mean depth `>3` and VMR `>50`, require `ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03`
- robust depth: median if median `>=20`, otherwise `mean / Reliable_Ref_zip_af^1.05`

Formula-AF gate:

- `XnY_ctx >= 10`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.95^24 = 0.291989`
- if mean depth `>3` and VMR `>50`, require `ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.02`
- robust depth: median if median `>=10`, otherwise `mean / Reliable_Ref_zip_af`

No intra-genus rescue rows were added for these CAMI3 coden15 runs.

## Mean Results

| method | L1 pp | Pearson | Spearman | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sylph | 27.7897 | 0.923617 | 0.668254 | 0.603485 | 58.67 | 12.33 | 64.67 |
| coden15 formula-AF gate robust rescue | 33.7127 | 0.909277 | 0.570194 | 0.597670 | 58.33 | 13.33 | 65.00 |
| old ctx-marker robust rescue | 39.2328 | 0.898686 | 0.546870 | 0.555312 | 51.67 | 11.00 | 71.67 |
| coden15 original gate robust rescue | 45.9801 | 0.896053 | 0.521866 | 0.508266 | 45.33 | 9.67 | 78.00 |

Per-sample formula-AF L1 pp: sample0 `34.5184`, sample1 `31.9916`, sample2 `34.6280`.

Per-sample original-gate L1 pp: sample0 `48.7655`, sample1 `34.3245`, sample2 `54.8503`.

## Interpretation

The Toy Mouse relative ranking does not reproduce on CAMI3. On CAMI3, Sylph remains best by L1, Pearson, Spearman, and mean F1. The coden15 formula-AF gate is still a real improvement over coden15 original and over old ctx-marker robust-rescue on this CAMI3 benchmark, but it does not beat Sylph.

This supports treating formula-AF as a useful coden15 diagnostic and retuning direction, not as the current abundance default. The old ctx-marker robust-rescue rule remains the first MinCO abundance baseline for Toy Mouse, but CAMI3 suggests the coden15 formula-AF gate transfers better than old ctx-marker on this specific NCBI species-taxid benchmark.

## Runtime

Each CAMI3 coden15 profile used 16 threads and finished successfully:

| sample | wall time | max RSS |
|---|---:|---:|
| 0 | 7:26.00 | 6.42 GB |
| 1 | 7:18.11 | 6.57 GB |
| 2 | 7:17.19 | 6.46 GB |

## Artifacts

- Profile outputs: `cami3_sample*_coden15_ctxmarker_split_naive_product.tsv`
- Scorer: `score_cami3_coden15_formula_af.py`
- Combined result: `summary.tsv`
- Per-sample presence: `coden15_cami3_sample_presence.tsv`
- Per-sample abundance: `coden15_cami3_sample_abundance.tsv`
- Selected species: `coden15_cami3_selected_species.tsv`
- Rescue audit: `coden15_cami3_rescued_species.tsv`
