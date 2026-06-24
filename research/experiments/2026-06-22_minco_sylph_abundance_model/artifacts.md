# Artifacts

- `effective_abundance_mean_metrics.tsv`: focused final mean comparison against GTDB species abundance truth.
- `effective_abundance_sample_metrics.tsv`: per-sample focused comparison.
- `intragenus_winner_rescue_mean_metrics.tsv`: final mean comparison for robust ctx-marker abundance plus conservative intra-genus winner rescue.
- `intragenus_winner_rescue_sample_metrics.tsv`: per-sample metrics for the intra-genus winner rescue; sample1 rescues `s__Pseudobutyrivibrio ruminis`.
- `score_effective_abundance.py`: focused scoring script for old MinCO abundance, robust MinCO abundance, and Sylph.
- `score_intragenus_winner_rescue.py`: scoring script for the final winning post-processing rule.
- `evaluate_abundance_estimators.py`: broader exploratory estimator sweep.
- `abundance_estimator_mean.tsv` and `abundance_estimator_metrics.tsv`: broad sweep outputs.
- `toymouse_sample{0,1,2}_ctxmarker_median_col.tsv`: rebuilt ctx-marker MinCO outputs with `Reliable_Ref_hit_median_depth`.
- `toymouse_sample{0,1,2}_ctxobjmarker_median_col.tsv`: rebuilt ctx+obj marker MinCO outputs with `Reliable_Ref_hit_median_depth`.
- `toymouse_sample0_ctxmarker_effective_abundance.tsv`: final binary smoke-test output with `Effective_abundance_depth` and `Normalized_effective_abundance_depth`.
- `toymouse_sample{0,1,2}_s1000_full_split_naive_product_topfrac_median025.tsv`: full S1000 non-marker diagnostic outputs.
- `toymouse_sample{0,1,2}_s2000_dedup_full_split_naive_product_topfrac_median025.tsv`: S2000 dedup non-marker diagnostic outputs.
- `toymouse_sample{0,1,2}_s2000_dedup_full_unique_naive_product_topfrac_median025.tsv`: S2000 full/shared-context best-diff-unique diagnostic outputs; mean L1 remained poor at 12.68.
- External truth and prior benchmark folder: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_cami2_toymouse_more_gtdb`.
- Sylph profiles: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample{0,1,2}/profile.tsv`.
