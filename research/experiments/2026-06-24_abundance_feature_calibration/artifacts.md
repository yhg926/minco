# Artifacts

| Kind | Path or Link | Source Path or URI | Description | Availability | Preserve? | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| script | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/search_feature_calibration.py` | local workspace | Main feature-calibration experiment script | available | yes | Imports previous abundance-normalization scorer |
| input-script | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py` | local workspace | Previous callset and abundance-normalization logic | available | yes | Used to build selected-call rows |
| command-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/commands.sh` | local workspace | Exact commands and working directory | available | yes | Rerunnable |
| parameter-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/parameters.tsv` | local workspace | Parameters, model grid, and runtime summary | available | yes |  |
| code-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/provenance/code_status.txt` | local workspace | Commit, branch, dirty status, and diffstat | available | yes | Full diff not captured here |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/selected_call_features.tsv` | generated | Feature rows for selected calls | available | yes | 4,319 lines including header |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/truth_profiles.tsv` | generated | Full per-sample truth abundance rows for correct FN/L1 scoring | available | yes | 589 lines including header |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/calibration_rule_summary.tsv` | generated | All raw and LOSO model summaries | available | yes | 1,117 lines including header |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/calibration_sample_metrics.tsv` | generated | Per-sample metrics for raw and LOSO models | available | yes | 6,697 lines including header |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/calibration_top200_by_l1.tsv` | generated | Top 200 rows sorted by all_L1 then all_F1 | available | yes |  |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/calibration_top200_by_f1.tsv` | generated | Top 200 rows sorted by all_F1 then all_L1 | available | yes |  |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/results/final_comparison.tsv` | generated | Compact comparison against Sylph | available | yes |  |
| temporary-log | `/tmp/minco_feature_calibration_features.stdout` | generated | Feature/truth table generation stdout | temporary | no | Regenerable |
| temporary-log | `/tmp/minco_feature_calibration_features.stderr` | generated | Feature/truth table generation timing and warnings | temporary | no | Regenerable |
| temporary-log | `/tmp/minco_feature_calibration_120.stdout` | generated | Broad 120-model run stdout | temporary | no | Regenerable |
| temporary-log | `/tmp/minco_feature_calibration_120.stderr` | generated | Broad 120-model run timing and warnings | temporary | no | Regenerable |
