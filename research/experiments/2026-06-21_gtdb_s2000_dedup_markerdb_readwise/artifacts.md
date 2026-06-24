# Artifacts

Date: 2026-06-21

Large generated data are stored under `/tmp/gtdb232_s2000_dedup_marker.qKJofv`. `/tmp/gtdb232_s2000_dedup_marker.latest` points to the same work directory. These paths are temporary and should be copied to durable storage if needed.

## Source Data

- S10000 GTDB MinCO sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`
- Toy Mouse sample0 reads: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- Toy Mouse gold profile: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/taxonomic_profile_0.txt`
- Toy Mouse positive source genomes and abundances: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/positive_source_genomes.tsv`
- GTDB r232 bacterial metadata: `/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz`
- GTDB r232 archaeal metadata: `/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz`
- Source-aware taxmap: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/mouse0_sourceaware_taxmap_ani95_minaf50.tsv`
- Prior S1000 unique readwise output: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_s1000_gtdb_unique_zip_unfiltered.tsv`
- Prior S1000 split readwise output: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_s1000_gtdb_split_zip_unfiltered.tsv`

## Generated Sketches

- S2000 sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_anno`
- S2000 deduped sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- S2000 deduped context markerdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`

## Dedup Outputs

- Dedup plan: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/dedup_plan.aaf003.tsv`
- Kept references: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/keep.aaf003.txt`
- Removed references: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/remove.aaf003.txt`
- Dedup log: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/dedup_aaf003.log`
- Dedup timing: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/dedup_aaf003.time.log`

## Markerdb Outputs

- Markerdb log: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.log`
- Markerdb psmp summary input: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv`
- Markerdb timing: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.time.log`
- Markerdb index timing: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/index_ctxmarker_s2000_dedup.time.log`

## Readwise Outputs

- S2000 markerdb unique output: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip_unfiltered.tsv`
- S2000 markerdb unique log: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip.log`
- S2000 markerdb unique timing: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_unique_zip.time.log`
- S2000 markerdb split output: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_split_zip_unfiltered.tsv`
- S2000 markerdb split log: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_split_zip.log`
- S2000 markerdb split timing: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/toymouse_sample0_s2000_dedup_ctxmarker_split_zip.time.log`
- S2000 markerdb split naive unfiltered output: `toymouse_sample0_s2000_marker_split_naive_unfiltered.tsv`
- S2000 markerdb split naive unfiltered timing: `toymouse_sample0_s2000_marker_split_naive_unfiltered.time.log`
- S2000 markerdb split zip-aaf with Poisson-depth p<=0.05 output: `toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.tsv`
- S2000 markerdb split zip-aaf with Poisson-depth p<=0.05 timing: `toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.time.log`
- S2000 markerdb split naive with Poisson-depth p<=0.05 output: `toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.tsv`
- S2000 markerdb split naive with Poisson-depth p<=0.05 timing: `toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.time.log`
- S2000 markerdb split naive with Poisson-product p<=0.05 output: `toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.tsv`
- S2000 markerdb split naive with Poisson-product p<=0.05 timing: `toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.time.log`
- S2000 markerdb split naive with product-NB p<=0.05 output: `toymouse_sample0_s2000_marker_split_naive_product_nb_p005.tsv`
- S2000 markerdb split naive with product-NB p<=0.05 timing: `toymouse_sample0_s2000_marker_split_naive_product_nb_p005.time.log`
- S2000 markerdb split naive with product-NB p<=0.05 and post-defake reliable abundance columns: `toymouse_sample0_s2000_marker_split_naive_product_nb_reliable_abundance_p005.tsv`
- S2000 markerdb split naive with product-NB p<=0.05 and post-defake reliable abundance timing: `toymouse_sample0_s2000_marker_split_naive_product_nb_reliable_abundance_p005.time.log`
- S2000 markerdb split naive with product top-25%-fraction filter output: `toymouse_sample0_s2000_marker_split_naive_product_topfrac025.tsv`
- S2000 markerdb split naive with product top-25%-fraction filter timing: `toymouse_sample0_s2000_marker_split_naive_product_topfrac025.time.log`
- S2000 markerdb split naive with product top-25%-fraction median replacement output: `toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv`
- S2000 markerdb split naive with product top-25%-fraction median replacement timing: `toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.time.log`
- S2000 markerdb split naive with `(best_diff+1)*depth` top-25%-fraction median replacement output: `toymouse_sample0_s2000_marker_split_naive_product1_topfrac_median025.tsv`
- S2000 markerdb split naive with `(best_diff+1)*depth` top-25%-fraction median replacement timing: `toymouse_sample0_s2000_marker_split_naive_product1_topfrac_median025.time.log`
- S2000 markerdb split per-read trace output: `toymouse_sample0_s2000_marker_split_zip_unfiltered_perread_trace_lcri.tsv`
- S2000 markerdb split per-read trace timing: `toymouse_sample0_s2000_marker_split_perread_trace_lcri.time.log`
- S2000 markerdb `GCF_018987235.1` assignment trace: `l_crispatus_s2000_marker_split_readwise_assignment_trace.tsv`
- S2000 markerdb assignment trace joined to CAMISIM truth: `l_crispatus_s2000_marker_split_readwise_assignment_trace_with_truth.tsv`
- S2000 markerdb origin summary by species: `l_crispatus_s2000_marker_split_readwise_assignment_trace_origin_by_species.tsv`
- S2000 markerdb origin summary by genome: `l_crispatus_s2000_marker_split_readwise_assignment_trace_origin_by_genome.tsv`
- S2000 markerdb final-context origin summary: `l_crispatus_s2000_marker_split_readwise_assignment_trace_final_context_origin_summary.tsv`

## Scoring Outputs In This Note Folder

- GTDB truth and score helper: `build_and_score_gtdb_ground_truth.py`
- Source-aware score helper: `score_sourceaware_direct.py`
- Reliable truncated-depth score helper: `score_reliable_truncated_depth_gate.py`
- Read-origin trace join helper: `join_l_crispatus_trace_to_truth.py`
- GTDB per-source genome truth map: `mouse0_gtdb_source_species_ground_truth.tsv`
- GTDB aggregated species profile: `mouse0_gtdb_species_profile.tsv`
- GTDB truth summary: `mouse0_gtdb_truth_summary.tsv`
- Direct GTDB-species scores: `gtdb_direct_scores.tsv`
- Direct GTDB false positives: `gtdb_direct_score_false_positives.tsv`
- Direct GTDB false negatives: `gtdb_direct_score_false_negatives.tsv`
- Direct GTDB-species scores using recalculated naive ANI for every method: `gtdb_direct_scores_naive_ani.tsv`
- Direct GTDB naive-ANI false positives: `gtdb_direct_score_naive_ani_false_positives.tsv`
- Direct GTDB naive-ANI false negatives: `gtdb_direct_score_naive_ani_false_negatives.tsv`
- Direct GTDB-species scores using recalculated naive ANI with strict `ANI > 0.95`: `gtdb_direct_scores_naive_ani_gt095.tsv`
- Direct GTDB strict naive-ANI false positives: `gtdb_direct_score_naive_ani_gt095_false_positives.tsv`
- Direct GTDB strict naive-ANI false negatives: `gtdb_direct_score_naive_ani_gt095_false_negatives.tsv`
- Direct GTDB strict naive-ANI scores with breadth-depth filter: `gtdb_direct_scores_naive_ani_gt095_breadth_depth.tsv`
- Direct GTDB strict naive-ANI breadth-depth false positives: `gtdb_direct_score_naive_ani_gt095_breadth_depth_false_positives.tsv`
- Direct GTDB strict naive-ANI breadth-depth false negatives: `gtdb_direct_score_naive_ani_gt095_breadth_depth_false_negatives.tsv`
- Direct GTDB strict naive-ANI scores with adjusted breadth `Ref_zip_af >= 0.5`: `gtdb_direct_scores_naive_ani_gt095_zipaf05.tsv`
- Direct GTDB strict naive-ANI adjusted-breadth false positives: `gtdb_direct_score_naive_ani_gt095_zipaf05_false_positives.tsv`
- Direct GTDB strict naive-ANI adjusted-breadth false negatives: `gtdb_direct_score_naive_ani_gt095_zipaf05_false_negatives.tsv`
- Direct GTDB reliable truncated-depth scores: `gtdb_direct_scores_reliable_truncated_depth.tsv`
- Direct GTDB reliable truncated-depth FP/FN details: `gtdb_direct_score_reliable_truncated_depth_details.tsv`
- Direct GTDB reliable truncated-depth threshold sweep: `gtdb_direct_score_reliable_truncated_depth_sweep.tsv`
- Direct GTDB reliable ZTP AF-to-ANI consistency scores: `gtdb_direct_score_ztp_ani_consistency.tsv`
- Direct GTDB product top-25%-fraction scores: `gtdb_direct_scores_product_topfrac025.tsv`
- Direct GTDB product top-25%-fraction FP/FN details: `gtdb_direct_score_product_topfrac025_details.tsv`
- Direct GTDB product top-25%-fraction median replacement scores: `gtdb_direct_scores_product_topfrac_median025.tsv`
- Direct GTDB product top-25%-fraction median replacement FP/FN details: `gtdb_direct_score_product_topfrac_median025_details.tsv`
- Direct GTDB product top-25%-fraction median replacement comparison against product-NB: `gtdb_direct_score_product_topfrac_median025_vs_product_nb.tsv`
- Direct GTDB product top-25%-fraction median replacement context-support sweep: `gtdb_direct_score_product_topfrac_median025_ctxcut_sweep.tsv`
- Direct GTDB product top-25%-fraction median replacement best-gate FP/FN details: `gtdb_direct_score_product_topfrac_median025_ctxcut15_floor040_details.tsv`
- Direct GTDB `(best_diff+1)*depth` top-25%-fraction median replacement scores: `gtdb_direct_scores_product1_topfrac_median025.tsv`
- Direct GTDB `(best_diff+1)*depth` top-25%-fraction median replacement context-support sweep: `gtdb_direct_score_product1_topfrac_median025_ctxcut_sweep.tsv`
- Direct GTDB `(best_diff+1)*depth` top-25%-fraction median replacement best-gate FP/FN details: `gtdb_direct_score_product1_topfrac_median025_ctxcut15_floor040_details.tsv`
- `L. mulieris_A` `(best_diff+1)*depth` median-replacement metrics: `l_mulieris_a_product1_topfrac_median025_metrics.tsv`
- Current product0-vs-Sylph GTDB abundance correlation summary: `mouse0_gtdb_abundance_correlation_current_product0_vs_sylph.tsv`
- Current product0-vs-Sylph GTDB abundance by-species detail: `mouse0_gtdb_abundance_current_product0_vs_sylph_by_species.tsv`
- Current product0-vs-S1000 non-markerdb GTDB abundance debug summary: `mouse0_gtdb_abundance_debug_s1000_vs_current_product0.tsv`
- Current product0-vs-S1000 non-markerdb dominant-species abundance detail: `mouse0_gtdb_abundance_debug_main_errors_s1000_vs_current.tsv`
- Current product0 GTDB ZIP-lambda abundance test summary: `mouse0_gtdb_abundance_zip_lambda_test.tsv`
- Current product0 GTDB ZIP-lambda abundance by-species detail: `mouse0_gtdb_abundance_zip_lambda_test_by_species.tsv`
- Pairwise context markerdb run directory: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c`
- Pairwise context markerdb sketch: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/sketch_T_S2000_aaf003_dedup_pairctx_af005`
- Pairwise context markerdb marker-size summary: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/markerdb_pairctx_af005.psmp.tsv`
- Pairwise context markerdb build log: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/markerdb_pairctx_af005.log`
- Pairwise context markerdb Toy Mouse ZIP direct output: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/toymouse_sample0_pairctx_af005_split_zip_unfiltered.tsv`
- Pairwise context markerdb Toy Mouse product-NB output: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/toymouse_sample0_pairctx_af005_split_naive_product_nb_p005.tsv`
- Pairwise context markerdb Toy Mouse product0 top-fraction median output: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/toymouse_sample0_pairctx_af005_split_naive_product_topfrac_median025.tsv`
- Pairwise context markerdb GTDB scores: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/pairctx_af005_gtdb_scores.tsv`
- Pairwise context markerdb GTDB FP/FN details: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/pairctx_af005_gtdb_score_details.tsv`
- Pairwise context markerdb active-gate raw-support sweep: `/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/pairctx_active_rawbreadth_sweep.tsv`
- Direct source-aware scores: `sourceaware_direct_scores.tsv`
- `L. crispatus` target comparison: `l_crispatus_markerdb_comparison.tsv`
- Main metric table: `summary.tsv`

## Failed Or Non-Authoritative Logs

- `/tmp/gtdb232_s2000_dedup_marker.qKJofv/score_exact.log` failed because `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv` was missing; this is superseded by `build_and_score_gtdb_ground_truth.py`.
- `/tmp/gtdb232_s2000_dedup_marker.qKJofv/score_sourceaware.log` came from a broader wrapper attempt and failed on `/tmp/gs_marine_short.profile`; use `sourceaware_direct_scores.tsv` for the valid comparison.
