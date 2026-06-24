# CAMI3 Read ANI Versus Source-Genome-to-Reference ANI

Date: 2026-06-23

## Question

For CAMI3 ToyGut selected calls, compare the read-based ANI to the actual source-genome-to-reference ANI. Report correlation and MAE for Sylph and MinCO.

## Method

For each selected call, I paired the called reference with CAMI3 source genomes from the same NCBI species taxid in that sample's read mapping file. Source-to-reference ANI was computed with `skani dist --min-af 0`.

Two source choices are reported:

- `read_major_source`: source genome with the most reads for that species taxid in the sample.
- `best_ani_source`: source genome with the highest skani ANI to the called reference. This is optimistic when multiple source strains exist.

Methods included:

- Sylph `Adjusted_ANI`
- Sylph `Naive_ANI`
- MinCO old ctx-marker `Ref_zip_aaf_ani`
- MinCO coden15 formula-AF `Ref_zip_aaf_ani`
- MinCO emitted/naive ANI columns as a saturation check

Diagnostics:

| metric | value |
|---|---:|
| selected calls considered | 616 |
| selected calls with reference path | 616 |
| source genomes needed | 216 |
| source genomes resolved | 216 |
| source/ref pairs requested | 218 |
| source/ref pairs scored by skani | 218 |
| calls with any scored source | 475 |

## Main Result

Using `read_major_source`:

| read ANI estimator | n | Pearson | Spearman | MAE | MAE pp | mean error |
|---|---:|---:|---:|---:|---:|---:|
| MinCO coden15 `Ref_zip_aaf_ani` | 164 | 0.4783 | 0.5601 | 0.00763 | 0.763 | -0.00221 |
| Sylph `Adjusted_ANI` | 165 | 0.4626 | 0.5262 | 0.00804 | 0.804 | -0.00441 |
| MinCO old ctx-marker `Ref_zip_aaf_ani` | 146 | 0.2540 | 0.4874 | 0.01071 | 1.071 | -0.00398 |
| Sylph `Naive_ANI` | 165 | 0.0815 | 0.2066 | 0.02964 | 2.964 | -0.02697 |

Using `best_ani_source`:

| read ANI estimator | n | Pearson | Spearman | MAE | MAE pp | mean error |
|---|---:|---:|---:|---:|---:|---:|
| MinCO coden15 `Ref_zip_aaf_ani` | 164 | 0.4833 | 0.5563 | 0.00781 | 0.781 | -0.00315 |
| Sylph `Adjusted_ANI` | 165 | 0.4767 | 0.5610 | 0.00833 | 0.833 | -0.00570 |
| MinCO old ctx-marker `Ref_zip_aaf_ani` | 146 | 0.2531 | 0.4528 | 0.01113 | 1.113 | -0.00502 |
| Sylph `Naive_ANI` | 165 | 0.1251 | 0.2629 | 0.03017 | 3.017 | -0.02827 |

## Saturation Check

MinCO emitted `ANI` and raw `ANI_naive_calc` for these selected rows are saturated at `1.0`, so Pearson and Spearman are undefined. Their MAE is worse than MinCO `Ref_zip_aaf_ani`:

| estimator | source choice | n | MAE |
|---|---|---:|---:|
| MinCO old emitted/naive ANI | read_major_source | 146 | 0.01293 |
| MinCO coden15 emitted/naive ANI | read_major_source | 164 | 0.01380 |
| MinCO old emitted/naive ANI | best_ani_source | 146 | 0.01189 |
| MinCO coden15 emitted/naive ANI | best_ani_source | 164 | 0.01286 |

## Interpretation

For source-genome-to-reference ANI accuracy, MinCO coden15 `Ref_zip_aaf_ani` is slightly better than Sylph `Adjusted_ANI` by Pearson and MAE on this paired-source subset. Sylph adjusted ANI is close and has nearly identical or slightly better Spearman under the optimistic `best_ani_source` pairing.

The important MinCO finding is that the selected/emitted ANI used in the gate is not useful for continuous ANI accuracy here because it saturates at `1.0`. The MinCO ANI signal that tracks source/ref ANI is `Ref_zip_aaf_ani`.

## Caveats

- Source pairing uses NCBI species taxid from the read mapping and the called reference metadata. This is source-genome truth, but still ambiguous when a species taxid contains multiple source strains.
- `best_ani_source` is optimistic; `read_major_source` is closer to the dominant read signal.
- Source/ref ANI is skani ANI, not full ANIm. It is fast and consistent across all scored pairs, but not an exact alignment-based truth.
- This evaluates selected calls with matching source species taxids, not all false-positive calls.

## Artifacts

- `score_source_ref_ani.py`: builds source/ref pairings, resolves source FASTAs, runs skani.
- `score_alternative_read_ani.py`: compares alternate read ANI columns against cached source/ref ANI.
- `summary.tsv`: selected-call ANI summary using the selected `pred_ani` column.
- `alternative_read_ani_summary.tsv`: main table with Sylph adjusted/naive and MinCO ZIP/emitted/naive ANI.
- `read_vs_source_ref_ani_detail.tsv`: per-call selected source pairing and errors.
- `call_source_candidates_with_skani.tsv`: all source candidates per call with skani ANI.
- `source_ref_skani_pairs.tsv`: cached skani source/ref ANI table.
- `diagnostics.tsv`: counts and resolution diagnostics.
