# Artifacts

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| input | `/tmp/cami_marine_sample0_reads.fq.gz` | CAMI II marine sample0 reads cache | maybe |
| input | `/tmp/cami_marine_sample1_reads.fq.gz` | CAMI II marine sample1 reads cache | maybe |
| input | `/tmp/cami_marine_sample2_reads.fq.gz` | CAMI II marine sample2 reads cache | maybe |
| input | `/tmp/gs_marine_short.profile` | CAMI II marine gold profile cache | maybe |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz` | CAMI III Toy Human Gut sample0 reads | yes |
| input | `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz` | Downloaded CAMI III Toy Human Gut sample1 reads | yes |
| input | `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz` | Downloaded CAMI III Toy Human Gut sample2 reads | yes |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_*.txt` | Toy Human Gut gold profiles | yes |
| input | `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno` | S1000 GTDB minco reference sketch | yes |
| input | `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb` | Sylph GTDB database | yes |
| output | `/tmp/minco_multisample_calibration_20260621/model_marine_toy0_2` | Six-sample model outputs | maybe |
| output | `/tmp/minco_multisample_calibration_20260621/threshold_marine_toy0_2` | Six-sample threshold-grid outputs | maybe |
| output | `/tmp/minco_multisample_calibration_20260621/sylph_toy_sample1/profile.tsv` | Sylph toy1 profile | maybe |
| output | `/tmp/minco_multisample_calibration_20260621/sylph_toy_sample2/profile.tsv` | Sylph toy2 profile | maybe |
| output | `summary_marine_toy0_2.tsv` | Per-sample baseline/model metrics copied into repo | yes |
| output | `mean_summary_marine_toy0_2.tsv` | Mean six-sample baseline/model metrics copied into repo | yes |
| output | `threshold_loso_summary_marine_toy0_2.tsv` | Per-sample threshold-grid metrics copied into repo | yes |
| output | `threshold_loso_mean_marine_toy0_2.tsv` | Mean threshold-grid metrics copied into repo | yes |
| output | `rf_feature_importance_marine_toy0_2.tsv` | Final RF feature importances copied into repo | yes |
| script | `scripts/calibrate_multisample_calls.py` | Multi-sample model calibration script | yes |
| script | `scripts/grid_loso_thresholds.py` | Interpretable threshold grid script | yes |
