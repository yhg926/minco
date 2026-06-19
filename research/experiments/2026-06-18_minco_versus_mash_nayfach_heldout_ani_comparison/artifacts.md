# Artifacts

| Kind | Path | Size | Description | Preserve? |
| --- | --- | ---: | --- | --- |
| input | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv` | 2.4M | Source sampled pair table with FASTA paths and ANIm labels. | yes |
| input | `research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv` | 1.4M | Final minco raw MoE and Best/HGB predictions. | yes |
| script | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/scripts/compare_mash_minco_anim.py` | repo file | Samples held-out pairs, caches Mash sketches, computes Mash distances, writes metrics. | yes |
| output | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/comparison_pairs.selected.tsv` | 112K | Selected 600-pair balanced held-out set. | yes |
| output | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/mash_sketch_status.tsv` | 504K | Per-genome Mash sketch status for 1,109 unique genomes. | yes |
| output | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/mash_minco_anim_pairs.tsv` | 136K | Per-pair ANIm, Mash distance/ANI, and minco ANI values. | yes |
| output | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/mash_minco_anim_metrics.tsv` | 4.0K | Aggregate and by-band metrics. | yes |
| sketch cache | `research/experiments/2026-06-18_minco_versus_mash_nayfach_heldout_ani_comparison/run_perbin100/mash_sketches_s10000` | 87M | Cached Mash `-s 10000` sketches for selected genomes. | maybe |
