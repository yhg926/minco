# Artifacts

| Kind | Path or Link | Source Path or URI | Description | Availability | Preserve? | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| input-ref | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker` | local temporary ref build | coden11 S2000 ctx-marker reference | available | yes | old current-best MinCO reference |
| input-ref | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup` | local temporary ref build | coden11 S2000 full dedup reference | available | yes | used for full fallback |
| input-map | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv` | local temporary ref build | ctx-marker size per reference | available | yes | accession -> unique marker size |
| input-mouse | `research/experiments/2026-06-22_minco_sylph_abundance_model/toymouse_sample*_ctxmarker_median_col.tsv` | previous experiment | mouse ctx-marker readwise outputs | available | yes | samples 0-2 |
| input-mouse | `research/experiments/2026-06-22_minco_sylph_abundance_model/toymouse_sample*_s2000_dedup_full_split_naive_product_topfrac_median025.tsv` | previous experiment | mouse full S2000 readwise outputs | available | yes | samples 0-2 |
| input-cami3 | `/tmp/cami3_toygut_abundance_rescue_20260623/minco_s*_ctxmarker_current_binary.tsv` | previous experiment | CAMI3 ctx-marker readwise outputs | available | yes | samples 0-2 |
| generated-cami3 | `/tmp/cami3_toygut_hybrid_full_20260623/minco_s*_full_product0.tsv` | this experiment | CAMI3 full S2000 readwise outputs | available | yes | generated with `minco_core/bin/minco`, `-p8` |
| script | `score_hybrid_mouse.py` | this experiment | exploratory mouse cutoff/XnY scorer | available | yes | includes fast metadata mapping and vectorized rescue |
| script | `score_hybrid_mouse_baseline_plus.py` | this experiment | final mouse baseline-plus fallback scorer | available | yes | selected rule comparison |
| script | `score_hybrid_cami3.py` | this experiment | CAMI3 baseline-plus fallback scorer | available | yes | selected rule comparison |
| output | `hybrid_mouse_mean_metrics.tsv` | this experiment | naive size-cutoff sweep on mouse | available | yes | showed simple fallback creates many FP |
| output | `hybrid_mouse_xnygrid_mean_metrics.tsv` | this experiment | mouse cutoff and XnY grid | available | yes | found XnY filtering |
| output | `hybrid_mouse_abundscale_mean_metrics.tsv` | this experiment | mouse fallback abundance scale grid | available | yes | found scale 0.4 |
| output | `hybrid_mouse_baseline_plus_mean_metrics.tsv` | this experiment | final mouse selected-strategy metrics | available | yes | main mouse result |
| output | `hybrid_mouse_baseline_plus_sample_metrics.tsv` | this experiment | final mouse per-sample metrics | available | yes | main mouse result |
| output | `hybrid_cami3_summary.tsv` | this experiment | CAMI3 selected-strategy summary | available | yes | main CAMI3 result |
| output | `hybrid_cami3_sample_presence.tsv` | this experiment | CAMI3 per-sample presence metrics | available | yes | main CAMI3 result |
| output | `hybrid_cami3_sample_abundance.tsv` | this experiment | CAMI3 per-sample abundance metrics | available | yes | main CAMI3 result |
| output | `results/score_hybrid_mouse_baseline_plus.stdout.tsv` | this experiment | captured final mouse scorer stdout | available | yes | timing stderr beside it |
| output | `results/score_hybrid_cami3.stdout.tsv` | this experiment | captured final CAMI3 scorer stdout | available | yes | timing stderr beside it |
| command-provenance | `commands.sh` | this experiment | exact commands and command blocks | available | yes | rerun starting point |
| parameter-provenance | `parameters.tsv` | this experiment | selected thresholds and inputs | available | yes | rerun metadata |
| code-provenance | `provenance/code_status.txt` | init script | repository status | available | yes | git metadata was not available |
