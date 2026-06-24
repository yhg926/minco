# Artifacts

| Kind | Path | Description | Availability | Preserve? | Notes |
| --- | --- | --- | --- | --- | --- |
| script | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py` | Offline threshold scorer for ctx-marker/full, ctx+obj, and ANI summaries | available | yes | Added in this experiment |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/results/combo_search_summary.tsv` | focused ctx-marker/full combo summary | available | yes | 5003 combos |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/results/ctxobj/combo_search_summary.tsv` | focused ctx+obj combo summary | available | yes | 5003 combos |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/results/ensemble_scaled_grid.tsv` | scaled ctx+obj add-back grid | available | yes | best F1 result came from this file |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/results/ensemble_scaled_highcami_grid.tsv` | high-CAMI3 ctx+obj add-back grid | available | yes | did not improve over scaled grid |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/final_comparison.tsv` | compact final comparison table | available | yes | used in NOTE.md |
| log | `/tmp/minco_threshold_combo_focused5000.stderr` | ctx-marker/full focused sweep timing/progress | temporary | maybe | peak RSS 3.90 GB, wall 8:51.96 |
| log | `/tmp/minco_threshold_combo_ctxobj_focused5000.stderr` | ctx+obj focused sweep timing/progress | temporary | maybe | peak RSS 3.83 GB, wall 7:39.50 |
| input | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_minco_sylph_abundance_model/toymouse_sample{0,1,2}_ctxmarker_median_col.tsv` | mouse ctx-marker MinCO rows | available | yes | cached prior experiment |
| input | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_minco_sylph_abundance_model/toymouse_sample{0,1,2}_s2000_dedup_full_split_naive_product_topfrac_median025.tsv` | mouse full S2000 MinCO rows | available | yes | cached prior experiment |
| input | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-22_minco_sylph_abundance_model/toymouse_sample{0,1,2}_ctxobjmarker_median_col.tsv` | mouse ctx+obj MinCO rows | available | yes | cached prior experiment |
| input | `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s{0,1,2}_ctxmarker_current_binary.tsv` | CAMI3 ctx-marker MinCO rows | temporary | yes | should be copied if `/tmp` is cleaned |
| input | `/tmp/cami3_toygut_hybrid_full_20260623/minco_s{0,1,2}_full_product0.tsv` | CAMI3 full S2000 MinCO rows | temporary | yes | should be copied if `/tmp` is cleaned |
| input | `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s{0,1,2}_ctxobj_current_binary.tsv` | CAMI3 ctx+obj MinCO rows | temporary | yes | should be copied if `/tmp` is cleaned |
| input | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-23_toymouse_source_rep_ani/concrete_ani_check_coden11_coden15_vs_sylph.tsv` | mouse source-to-representative ANIm and ANI estimator table | available | yes | used for matched ANI MAE |
