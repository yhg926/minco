# Artifacts

Repo-tracked:

- Experiment note: `research/experiments/2026-06-20_minco_readwise_correction_research/NOTE.md`
- Analysis script: `research/experiments/2026-06-20_minco_readwise_correction_research/scripts/analyze_readwise_corrections.py`
- Summary table: `research/experiments/2026-06-20_minco_readwise_correction_research/summary.tsv`
- Commands: `research/experiments/2026-06-20_minco_readwise_correction_research/commands.sh`

External inputs:

- CAMI marine sample0 reads: `/tmp/cami_marine_sample0_reads.fq.gz`
- CAMI gold profile: `/tmp/gs_marine_short.profile`
- minco S10000 unfiltered detail table: `/tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_unfiltered.tsv`
- Sylph profile: `/tmp/sylph_marine_sample0/profile.tsv`
- Full-lineage taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`

Generated `/tmp` outputs:

- Output directory: `/tmp/minco_readwise_correction_20260620`
- Joined minco/Sylph/gold feature table: `/tmp/minco_readwise_correction_20260620/joined_s10000_rows.tsv`
- Strategy summary: `/tmp/minco_readwise_correction_20260620/strategy_summary.tsv`
- Strategy grid: `/tmp/minco_readwise_correction_20260620/strategy_grid.tsv`
- Model summary: `/tmp/minco_readwise_correction_20260620/model_summary.tsv`
- Feature strata summary: `/tmp/minco_readwise_correction_20260620/feature_strata_summary.tsv`
- Sylph-only TP correction features: `/tmp/minco_readwise_correction_20260620/sylph_only_tp_minco_correction_features.tsv`
- CV predictions for Sylph-matched distance calibration: `/tmp/minco_readwise_correction_20260620/sylph_matched_distance_cv_predictions.tsv`
- Same-sample HGB Sylph-calibrated ANI predictions: `/tmp/minco_readwise_correction_20260620/hgb_sylphcal_ani_allrows.tsv`
- Same-sample row classifier probabilities: `/tmp/minco_readwise_correction_20260620/row_gold_classifier_allrows_probability.tsv`
- Gold species table: `/tmp/minco_readwise_correction_20260620/gold_species_sample0.tsv`
- Sylph profile with mapped taxids: `/tmp/minco_readwise_correction_20260620/sylph_with_taxids.tsv`
- Manifest: `/tmp/minco_readwise_correction_20260620/manifest.tsv`
