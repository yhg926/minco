# Artifacts

## Downloaded CAMI Data

- Setup archive: `/mnt/new3T/minco_cami2_toymouse_20260621/CAMISIM_setup.tar.gz`
- Sample0 reads tar: `/mnt/new3T/minco_cami2_toymouse_20260621/2017.12.29_11.37.26_sample_0_reads.tar`
- FASTQ: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- Official gold profile: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/taxonomic_profile_0.txt`
- Setup distribution: `/mnt/new3T/minco_cami2_toymouse_20260621/setup/distributions/distribution_0.txt`
- Setup metadata: `/mnt/new3T/minco_cami2_toymouse_20260621/setup/internal/meta_data.tsv`
- Setup genome locations: `/mnt/new3T/minco_cami2_toymouse_20260621/setup/internal/genome_locations.tsv`

## Tool Outputs

- minco unique unfiltered: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_s1000_gtdb_unique_zip_unfiltered.tsv`
- minco split unfiltered: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_s1000_gtdb_split_zip_unfiltered.tsv`
- minco unique time log: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_unique.time.log`
- minco split time log: `/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/toymouse_sample0_split.time.log`
- Sylph sketch: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/anonymous_reads.fq.gz.sylsp`
- Sylph profile: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/profile.tsv`
- Sylph sketch log: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/sketch.log`
- Sylph profile log: `/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/profile.log`

## Repo Results

- Manifest: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/manifest_mouse0.tsv`
- Train9/test mouse score directory: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/train9_test_mouse0/`
- Method metrics: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_method_metrics.tsv`
- Abundance metrics: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_abundance_metrics.tsv`
- Call details: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_call_details.tsv`
- Gold species absent from taxmap: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_gold_not_in_taxmap.tsv`
- Reconstructed profile from CAMISIM setup: `/mnt/new3T/minco_cami2_toymouse_20260621/taxonomic_profile_mouse0_species.txt`
- Source-aware crosswalk working directory: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/`
- Positive source genome metadata: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/positive_source_genomes.tsv`
- Selected called representatives: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/selected_called_refs.tsv`
- skani representative-vs-source table: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/skani_called_vs_sources.tsv`
- Source-aware overrides: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/mouse0_sourceaware_overrides_ani95_minaf50.tsv`
- Source-aware full taxmap: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/mouse0_sourceaware_taxmap_ani95_minaf50.tsv`
- Source-aware score directory: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/train9_test_mouse0_sourceaware_ani95_minaf50/`
- Original versus source-aware comparison: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_comparison.tsv`
- Gold call diagnostics: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_gold_call_diagnostics.tsv`
- Sylph-only true positives versus best minco model: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_sylph_only_vs_train9_rf_hgb_avg.tsv`
- Concise Sylph advantage summary: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_sylph_advantage_summary.tsv`
- minco RF/HGB false positives: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_minco_rfhgb_false_positives.tsv`
- RF/HGB threshold/guard sweep: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_rfhgb_threshold_guard_sweep.tsv`
- ANIm source-aware pairs: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_sylph_advantage_anim.tsv`
- ANIm versus minco/Sylph ANI comparison: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_sylph_advantage_anim_compare.tsv`
- ANIm/minco/Sylph comparison with coverage: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/mouse0_sourceaware_sylph_advantage_anim_compare_with_coverage.tsv`
- `L. crispatus` per-read trace minco output: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_zip_unfiltered_perread_trace_lcri.tsv`
- `L. crispatus` selected assignment trace: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_gcf018987235_readwise_assignment_trace.tsv`
- `L. crispatus` selected assignment trace joined to CAMISIM truth: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_gcf018987235_readwise_assignment_trace_with_truth.tsv`
- `L. crispatus` assignment origin summary by genome: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_gcf018987235_trace_origin_by_genome.tsv`
- `L. crispatus` final context origin summary: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_gcf018987235_final_context_origin_summary.tsv`
- `L. crispatus` context fake-model feature prototype: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_gcf018987235_context_fake_model_features.tsv`
- Filtered naive per-read output, threshold 3.5: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_poisson_diff_t35.tsv`
- Filtered naive per-read time log, threshold 3.5: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_poisson_diff_t35.time.log`
- Filtered naive source-aware score directory: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/train9_test_mouse0_sourceaware_poisson_diff_t35/`
- Fake-probability weighted naive per-read output, threshold 3.0: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_fakeprob_t30.tsv`
- Fake-probability weighted naive per-read output, threshold 3.5: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/toymouse_sample0_s1000_gtdb_split_naive_perread_fakeprob_t35.tsv`
- `L. crispatus` fake-probability target summary: `research/experiments/2026-06-21_minco_toy_mouse_gut_sample0/l_crispatus_fakeprob_summary.tsv`
