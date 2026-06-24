# Toy Mouse sample0 source genome to GTDB representative ANI

Date: 2026-06-23

## Question

For one Toy Mouse gut sample, use the known positive source genomes as ground truth and compare each source genome to its GTDB representative with ANIm. Then list the corresponding MinCO ANI values and Sylph ANI values for the same GTDB representative.

## Dataset

- Sample: Toy Mouse gut sample0.
- Positive source genomes: 75 rows from `research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/mouse0_gtdb_source_species_ground_truth.tsv`.
- GTDB metadata: `/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz` and `/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz`.
- GTDB representative FASTA paths: existing resolved path lists in `research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse/`.
- MinCO profile: `research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_ctxmarker_rawani_report.tsv`.
- MinCO reference: coden11, ctx-marker, ctxobj64 storage. The stderr run log reports `coden_len=11; ctx_bits=44; obj_bits=20; storage=ctxobj64`.
- Sylph profile: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/profile.tsv`.

## Method

1. Map each positive source genome accession to its GTDB metadata row.
2. Extract `gtdb_genome_representative` and resolve the representative FASTA path.
3. Compute source-genome-to-representative ANIm with `nucmer --mum`, `delta-filter -1`, and identity from the filtered delta records:
   `ANIm_ANI = 1 - similarity_errors / aligned_bp`.
4. Join MinCO by representative accession and report:
   `minco_naive_ANI_calc`, `minco_aafANI_Ref_zip_aaf_ani`, and the profile `ANI` column.
5. Join Sylph by representative accession and report:
   `sylph_Adjusted_ANI` and `sylph_Naive_ANI`.

The MinCO readwise reporting code was fixed for this run. Before the fix, product-topfrac-median context filtering supplied the feature vector used for reported naive ANI, so many representative rows had `N_diff_obj=0` and `ANI=1.0`. After the fix, filtered/defaked features still drive support and abundance, but reported naive ANI uses the raw readwise diff features. This removes the accidental saturation.

ANIm intermediate files are cached under `/tmp/toymouse_source_rep_anim_20260623/`.

## Results

Validation counts:

- Positive source genomes: 75.
- GTDB representative paths found: 75/75.
- ANIm truth values computed: 75/75.
- MinCO representative rows found: 75/75.
- Sylph representative rows with ANI values: 69/75.
- Source genome is already the GTDB representative: 19/75.
- ANIm ANI range: 0.9612852911 to 1.0.
- Fixed MinCO coden length in joined table: 11.
- Fixed MinCO naive ANI values across the 75 joined rows: 64 unique values, 2 rows at 1.0.

### Coden11 Summary

Summary against source-to-representative ANIm, using the fixed coden11 ctx-marker profile:

| comparison_set | estimator | n | Pearson | Spearman | MAE |
|---|---:|---:|---:|---:|---:|
| all_available | sylph_Adjusted_ANI | 69 | 0.802689 | 0.775086 | 0.003862 |
| all_available | minco_aafANI_Ref_zip_aaf_ani | 75 | 0.564614 | 0.694266 | 0.011812 |
| all_available | minco_ANI_col | 75 | 0.422997 | 0.359613 | 0.016520 |
| all_available | minco_naive_ANI_calc | 75 | 0.422995 | 0.359613 | 0.016520 |
| all_available | sylph_Naive_ANI | 69 | 0.271463 | 0.215225 | 0.032675 |
| common_minco_sylph | sylph_Adjusted_ANI | 69 | 0.802689 | 0.775086 | 0.003862 |
| common_minco_sylph | minco_aafANI_Ref_zip_aaf_ani | 69 | 0.649036 | 0.736683 | 0.008967 |
| common_minco_sylph | minco_ANI_col | 69 | 0.393297 | 0.300733 | 0.015330 |
| common_minco_sylph | minco_naive_ANI_calc | 69 | 0.393294 | 0.300733 | 0.015331 |
| common_minco_sylph | sylph_Naive_ANI | 69 | 0.271463 | 0.215225 | 0.032675 |

### Coden15 Trial

I reran Toy Mouse sample0 with the rebuilt coden15 binary against `/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker`.

The coden15 stderr log reports:

- `coden_len=15`
- `ctx_bits=60`
- `obj_bits=32`
- `storage=ctxobj96`

Runtime and output:

- Wall time: 7:08.82.
- Peak RSS: 6,362,188 KB.
- Reported profile rows: 3,071.
- Joined source-positive representatives: 75/75.
- Fixed coden15 MinCO naive ANI values across the 75 joined rows: 64 unique values, 2 rows at 1.0.

Common 69-row comparison:

| method | estimator | Pearson | Spearman | MAE |
|---|---|---:|---:|---:|
| coden11 | minco_aafANI_Ref_zip_aaf_ani | 0.649036 | 0.736683 | 0.008967 |
| coden15 | minco_aafANI_Ref_zip_aaf_ani | 0.492383 | 0.569539 | 0.006771 |
| coden11 | minco_naive_ANI_calc | 0.393294 | 0.300733 | 0.015331 |
| coden15 | minco_naive_ANI_calc | 0.423006 | 0.435913 | 0.025913 |
| Sylph | Adjusted_ANI | 0.802689 | 0.775086 | 0.003862 |

Coden15 improves MinCO AAF ANI absolute error on this source-to-representative ANIm benchmark, but the raw naive ANI is worse than coden11 and strongly underestimates ANIm.

### Assembly Source-to-Representative Naive ANI

To separate intrinsic MinCO source/ref sketch-model bias from readwise/sequencing bias, I also compared each positive source genome assembly directly to its GTDB representative assembly. This uses explicit source and representative sketches and then joins only each source's true representative target.

Summary against source-to-representative ANIm:

| estimator | n | Pearson | Spearman | MAE | Mean error |
|---|---:|---:|---:|---:|---:|
| asm_minco_naive_coden11 | 75 | 0.987161 | 0.981968 | 0.001006 | 0.000353 |
| asm_minco_naive_coden15 | 75 | 0.981530 | 0.974275 | 0.003485 | -0.003379 |
| sylph_Adjusted_ANI | 69 | 0.802689 | 0.775086 | 0.003862 | 0.000034 |
| readwise_minco_naive_coden11 | 75 | 0.422995 | 0.359613 | 0.016520 | -0.015315 |
| readwise_minco_naive_coden15 | 75 | 0.419866 | 0.438069 | 0.026172 | -0.025145 |
| sylph_Naive_ANI | 69 | 0.271463 | 0.215225 | 0.032675 | -0.031335 |

Interpretation: the large negative MinCO naive ANI bias seen in readwise mode is not intrinsic to the source genome vs representative genome relationship. Assembly-to-assembly MinCO naive coden11 is nearly unbiased and tracks ANIm very tightly. The readwise bias therefore comes from sequencing/readwise context assignment, context depth/filtering, markerdb reduction, or their interaction.

Caveat: in this current local build, the generated coden11 assembly source sketches report variable sketch-entry counts despite `-S 2000`, while coden15 reports exactly 2000 entries. The existing GTDB coden11 S2000 reference still reports 2000 entries. So the coden11 assembly result should be read as a direct source/ref assembly sanity check, not a perfectly matched S2000 markerdb readwise simulation.

### Readwise Bias vs Abundance and Depth

I stratified source-positive pairs by source abundance and joined readwise support/depth metrics from the MinCO profiles.

Coden11 readwise naive ANI error by source abundance:

| abundance bin | n | c11 MAE | c11 mean error | mean breadth | mean depth | mean hit depth |
|---|---:|---:|---:|---:|---:|---:|
| <=1 | 31 | 0.022788 | -0.022467 | 0.161429 | 0.394347 | 1.527851 |
| 1-5 | 22 | 0.017477 | -0.016808 | 0.340905 | 7.376476 | 9.546196 |
| 5-10 | 6 | 0.013026 | -0.006033 | 0.594693 | 1.958750 | 3.015773 |
| 10-20 | 6 | 0.009208 | -0.009208 | 0.834184 | 2.138997 | 2.573511 |
| 20-50 | 3 | 0.002092 | 0.000479 | 0.653405 | 3.581669 | 5.055784 |
| 50-100 | 4 | 0.001376 | 0.001299 | 0.908950 | 9.482950 | 10.387921 |
| 100-500 | 1 | 0.000000 | 0.000000 | 0.997863 | 66.264245 | 66.406138 |
| >1000 | 2 | 0.001433 | 0.001206 | 0.833175 | 246.867250 | 297.696302 |

High-abundance examples:

- `s__Lacticaseibacillus paracasei`: abundance 3057, ANIm 0.989315, c11 readwise naive 0.991954, c11 breadth 0.823102, mean depth 340.31.
- `s__Lactobacillus helveticus`: abundance 1264.08, ANIm 0.992292, c11 readwise naive 0.992065, c11 breadth 0.843248, mean depth 153.42.
- `s__Secundilactobacillus pentosiphilus`: abundance 427, ANIm 1.0, c11 readwise naive 1.0, c11 breadth 0.997863, mean depth 66.26.

Interpretation: the main negative readwise naive ANI bias is strongly associated with low source abundance and low marker breadth/depth. For coden11, Spearman correlation between absolute readwise error and source abundance is -0.596; with hit depth it is -0.458. High abundance plus high breadth/depth gives much better readwise ANI.

Abundance alone is not sufficient: `s__Lactobacillus sp910589675` has source abundance 28 but only c11 breadth 0.130435 and mean depth 0.5, so it still has weak readwise support. There are also representative-level mixture effects: two `s__Staphylococcus aureus` source genomes map to the same GTDB representative; the readwise representative row is aggregate, so it matches the high-ANI source well but can look inflated for the lower-ANI source. This is not simple low coverage; it is source/strain mixture at the representative row.

Main table:

- `toymouse_sample0_source_rep_ani.tsv`
- `toymouse_sample0_source_rep_ani_coden15.tsv`
- `summary_coden11_vs_coden15.tsv`
- `source_ref_assembly_minco_naive.tsv`
- `source_ref_assembly_minco_naive_summary.tsv`
- `readwise_bias_by_abundance_bin.tsv`
- `high_abundance_readwise_ani_check.tsv`

The table contains one row per positive source genome and includes source accession, GTDB representative accession, GTDB species, source abundance, ANIm truth, MinCO naive ANI, MinCO AAF ANI, Sylph adjusted ANI, and Sylph naive ANI.

## Interpretation

Sylph adjusted ANI is closest to direct ANIm truth when using readwise profiles. The MinCO naive ANI saturation bug is fixed, but raw MinCO readwise naive ANI underestimates source-to-representative ANIm mostly at low abundance/low support. Assembly-to-assembly MinCO naive does not show this bias, so the remaining problem is in the readwise/sequencing path rather than the core source/ref naive ANI formula. High abundance with high marker breadth/depth gives much better readwise ANI; representative-level source mixtures can still distort per-source comparisons.

## Caveats

- This is one Toy Mouse gut sample only.
- Sylph has six missing representative ANI values in this profile, so fair direct comparison uses the 69-row common subset.
- ANIm here is implemented from filtered MUMmer delta records rather than a pyANI wrapper; the identity formula matches the ANIm-style filtered-alignment calculation used in earlier Toy Mouse notes.
- The MinCO profile used here is a newly rerun coden11 ctx-marker readwise profile after the raw-ANI reporting fix.
- The same source change is compiled into `minco_core/bin_coden15/minco`; coden15-specific outputs are kept as separate `*_coden15.tsv` files.
