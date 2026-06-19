# Artifacts

| Kind | Path | Size | Description | Preserve? |
| --- | --- | ---: | --- | --- |
| input | `/mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv` | 120,544,386 bytes | Source pair table with ANIm labels and historical KSSD3 T10 fields. | yes |
| input | `/mnt/new3T/skani_data/Nayfach_data/fna` | 52,515 `.fna` files | Nayfach genome FASTA directory used to run direct minco pairs. | yes |
| script | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/minco_nayfach_calibration_pilot.py` | repo file | Stratified sampler, minco feature extractor, model trainer, C HGB exporter. | yes |
| script | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/optimize_moe_model.py` | repo file | Optimizes MoE denominator bases and linear coefficients for `model_ani.h`. | yes |
| script | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/train_fixedp_moe.py` | repo file | Earlier fixed-denominator MoE coefficient refit; kept as a candidate comparison. | yes |
| script | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/scripts/validate_recalibrated_binary.py` | repo file | Runs compiled `minco ani` on held-out pairs and compares to Python predictions/ANIm. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/sample_pairs.tsv` | 2.4M | 12,000 sampled pair list, 2,000 per ANI bin. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/sample_bin_counts.tsv` | 4.0K | Available source-table pair counts per ANI bin. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/minco_features.tsv` | 5.3M | Initial raw minco feature table before MoE replacement. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/model_metrics.tsv` | 4.0K | Initial train/test/all metrics before MoE replacement. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/model_metrics_by_bin.tsv` | 16K | Metrics split by ANI band. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/model_predictions.tsv` | 1.4M | Per-pair train/test predictions for all evaluated models. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/ridge_coefficients.tsv` | 4.0K | Coefficients for the ridge/MoE-style full-feature correction. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/fixedp_moe_metrics.tsv` | 4.0K | Metrics for the fixed-denominator MoE coefficient refit. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/fixedp_moe_coefficients.tsv` | 4.0K | Raw-basis candidate coefficients for `model_ani.h`; not installed. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/moe_opt_long/optimized_moe_metrics.tsv` | 4.0K | Optimized MoE train/test/all metrics before source installation. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/moe_opt_long/optimized_moe_params.tsv` | 4.0K | Optimized MoE denominator parameters and coefficients. | yes |
| generated C | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/moe_opt_long/optimized_moe_model_ani_arrays.hfrag` | 4.0K | C array fragment used to patch `model_ani.h`. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/minco_features.tsv` | 5.3M | Feature table regenerated after installing optimized MoE. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_metrics.tsv` | 4.0K | Final metrics after MoE replacement and HGB retraining. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_metrics_by_bin.tsv` | 16K | Final metrics split by ANI band. | yes |
| output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv` | 1.4M | Final per-pair predictions after MoE replacement. | yes |
| generated C | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_refaf_hgb.minco_generated.h` | 524K | Final generated HGB header trained on recalibrated MoE raw score. | yes |
| installed C | `minco_core/src/model_ani.h` | 17K | Active MoE model header with optimized `NUM_CODENS == 11` constants. | yes |
| installed C | `minco_core/src/model_refaf_hgb.h` | 522K | Active HGB header regenerated after MoE replacement. | yes |
| backup | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_ani.before_moe_recalibration.h` | 17K | Previous `model_ani.h` before optimized MoE installation. | yes |
| backup | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_refaf_hgb.before_minco_recalibration.h` | 514K | Previous generated HGB header before Nayfach minco recalibration. | yes |
| backup | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_refaf_hgb.before_moe_recalibrated_hgb.h` | 521K | HGB header before retraining on optimized MoE raw score. | yes |
| validation | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s3_raw_alltest.tsv` | 912K | Compiled raw CtxMoE output on all 3,600 held-out pairs. | yes |
| validation | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s3_raw_alltest.metrics.tsv` | 4.0K | Binary-level raw CtxMoE held-out metrics. | yes |
| validation | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s1_best_alltest.tsv` | 920K | Compiled default Best output on all 3,600 held-out pairs. | yes |
| validation | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/binary_validation_s1_best_alltest.metrics.tsv` | 4.0K | Binary-level default Best held-out metrics. | yes |
| smoke output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin5` | directory | 30-pair end-to-end smoke run. | maybe |
| pilot output | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin500` | directory | 3,000-pair pilot used before final 12,000-pair run. | maybe |
| temporary | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/tmp` | 4.0K | Empty temp directory for per-pair minco output. | no |
| temporary | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000/tmp_binary_validation` | 4.0K | Empty temp directory for binary validation output. | no |
