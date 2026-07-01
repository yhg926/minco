# Artifacts

## Repo-Tracked Small Outputs

| path | status | description |
| --- | --- | --- |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/sweep_exact_hit_abundance_trigger.py` | tracked candidate | Offline scorer for cached exact-hit abundance trigger sweep. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/score_trigger_validation_reruns.py` | tracked candidate | Raw-rerun validation scorer for samples that triggered threshold 0.10. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_scores.tsv` | tracked candidate | Per-sample offline cached scores for current and trigger thresholds. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_panel_summary.tsv` | tracked candidate | Offline cached panel summaries. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_panel_summary_with_sylph.tsv` | tracked candidate | Offline cached panel summaries plus external Sylph baseline rows. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_overall.tsv` | tracked candidate | Offline cached cross-panel threshold summary. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_trace.tsv` | tracked candidate | Per-sample offline trigger decisions and exact-table availability. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_mapping.tsv` | tracked candidate | Offline profile and exact-table paths used by sample. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_rerun_sample_scores.tsv` | tracked candidate | Scores for four raw rerun samples. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_rerun_sample_deltas.tsv` | tracked candidate | Baseline-vs-trigger deltas for four raw rerun samples. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_validated_scores.tsv` | tracked candidate | Combined validated per-sample scores. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_validated_panel_summary.tsv` | tracked candidate | Validated panel summaries after raw rerun replacement. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_validated_panel_summary_with_sylph.tsv` | tracked candidate | Validated panel summaries plus external Sylph baseline rows. |
| `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/results/exact_hit_abundance_trigger_validated_overall.tsv` | tracked candidate | Validated cross-panel summary used for the 0.30 default decision. |

## External / Temporary Inputs

| path | status | description |
| --- | --- | --- |
| `/tmp/minco_current_code_toymouse_refresh_20260627/sample{5,6,7}_current.tsv` | temporary | Current toy mouse MinCO profiles used by offline sweep. |
| `/tmp/cami2_toymouse_current_default_20260626/run/sample{5,6,7}_default/minco.best_diff_split.exact.unfiltered.tsv` | temporary | Cached toy mouse exact split tables. |
| `/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin*/minco_sample*_current_default.tsv` | external/local | HMP airskin current profiles for most samples. |
| `/tmp/minco_current_code_hmp_airskin*/work/minco.best_diff_split.exact.unfiltered.tsv` | temporary | Cached HMP airskin exact split tables where present. |
| `/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample{0,6}_current_default.tsv` | temporary | HMP gastrooral current profiles. |
| `/tmp/minco_exactsplit_universal_20260626/hmp_gastrooral6_s1000_gtdb_split_exact_unfiltered.tsv` | temporary | Cached HMP gastrooral sample 6 exact split table. |
| `/tmp/cami3_toy_human_gut_20260626/run/sample{0,1,2}_universal_autoexact*.tsv` | temporary | CAMI3 current cached profiles. |
| `/tmp/cami3_toy_human_gut_20260626/run/sample2_exact_forced/minco.best_diff_split.exact.unfiltered.tsv` | temporary | Cached CAMI3 sample 2 exact split table. |

## Raw Rerun Validation Artifacts

| path | status | description |
| --- | --- | --- |
| `/tmp/minco_exact_trigger_validation_20260630/reads` | temporary, 18 GB | Restored/regenerated FASTQ inputs for CAMI3 samples 0-2 and HMP airskin sample 13. |
| `/tmp/minco_exact_trigger_validation_20260630/run` | temporary, 522 MB | Raw rerun workdirs, profile outputs, baseline table outputs, and exact split sidecars. |
| `/tmp/minco_exact_trigger_validation_20260630/logs` | temporary | `/usr/bin/time -v` logs for read restoration, trigger-0.10 raw reruns, and baseline table replays. |
| `/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample0_trigger010.tsv` | temporary | CAMI3 sample 0 raw rerun with `--exact-split-abundance-trigger 0.10`. |
| `/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample1_trigger010.tsv` | temporary | CAMI3 sample 1 raw rerun with `--exact-split-abundance-trigger 0.10`. |
| `/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample2_trigger010.tsv` | temporary | CAMI3 sample 2 raw rerun with `--exact-split-abundance-trigger 0.10`; this is the regression case. |
| `/tmp/minco_exact_trigger_validation_20260630/run/hmp_airskin13_trigger010.tsv` | temporary | HMP airskin sample 13 raw rerun with `--exact-split-abundance-trigger 0.10`. |
