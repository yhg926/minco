# Artifacts

## Code

- MinCO source patched for edge dump: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/src/command_ani.c`
- Built binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/bin/minco`
- EM scorer: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/score_mouse0_edge_em.py`
- Filtered edge-EM search: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/search_mouse0_filtered_edge_em.py`
- Best-filter beta refinement: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/refine_mouse0_filtered_edge_em_beta.py`
- Three-sample validation scorer: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/score_mouse012_filtered_edge_em.py`
- Run script: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/commands.sh`

## Inputs

- Full GTDB S2000 dedup refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Toy Mouse sample0 reads: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- GTDB metadata loaded by existing mouse scorer: `/mnt/new3T/gtdbr220/bac120_metadata_r226.tsv.gz`
- Truth profile loaded by existing code path: `evaluate_abundance_estimators.load_truth_profile(0)`

## Temporary Large Outputs

Stored under `/tmp/minco_context_em_20260624`; these can be regenerated with `commands.sh`.

- `/tmp/minco_context_em_20260624/mouse_s0_edges.tsv`, 356 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s0_profile.tsv`, 103 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s0_edge.time.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s0_edge.stderr.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s0_edge.stdout.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s1_edges.tsv`, 462 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s1_profile.tsv`, 132 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s1_edge.time.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s1_edge.stderr.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s1_edge.stdout.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s2_edges.tsv`, 487 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s2_profile.tsv`, 134 MB, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s2_edge.time.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s2_edge.stderr.log`, available now, temporary
- `/tmp/minco_context_em_20260624/mouse_s2_edge.stdout.log`, available now, temporary

Smoke-test temporary files:

- `/tmp/minco_edge_smoke_1k.fq.gz`
- `/tmp/minco_edge_smoke_edges.tsv`
- `/tmp/minco_edge_smoke_profile.tsv`
- `/tmp/minco_edge_smoke.time.log`

## Experiment Outputs

- Summary table: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/summary.tsv`
- Filtered follow-up summary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/filtered_summary.tsv`
- Focused EM grid: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_edge_em_focused_grid.tsv`
- Filtered edge-EM grid: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_filtered_edge_em_grid.tsv`
- Filtered edge-EM top rows: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_filtered_edge_em_top50.tsv`
- Filtered edge-EM best broad-grid row: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_filtered_edge_em_best.tsv`
- Filtered edge diagnostics: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_edge_filter_diagnostics.tsv`
- Best-filter beta refinement: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_filtered_edge_em_beta_refine.tsv`
- Mouse0-2 filtered validation sample metrics: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse012_filtered_edge_em_sample_metrics.tsv`
- Mouse0-2 filtered validation summary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse012_filtered_edge_em_summary.tsv`
- Mouse0-2 filtered validation diagnostics: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse012_filtered_edge_em_diagnostics.tsv`
- Mouse0-2 filtered validation top summary copy: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/mouse012_filtered_summary.tsv`
- Edge extraction stats: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_edge_stats.tsv`
- Species-level marker-filtered edge cache: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse0_marker_candidate_edge_species.tsv`
- Sample1 edge extraction stats: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse1_edge_stats.tsv`
- Sample1 species-level marker-filtered edge cache: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse1_marker_candidate_edge_species.tsv`
- Sample2 edge extraction stats: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse2_edge_stats.tsv`
- Sample2 species-level marker-filtered edge cache: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/mouse2_marker_candidate_edge_species.tsv`
- Scorer stdout copy: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse0_edge_em.stdout.tsv`
- Scorer stderr: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse0_edge_em.stderr.txt`
- Scorer time log: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse0_edge_em.time.log`
- Filtered search stdout copy: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/search_mouse0_filtered_edge_em.stdout.tsv`
- Filtered search stderr: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/search_mouse0_filtered_edge_em.stderr.txt`
- Beta refinement stdout copy: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/refine_mouse0_filtered_edge_em_beta.stdout.tsv`
- Beta refinement stderr: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/refine_mouse0_filtered_edge_em_beta.stderr.txt`
- Mouse0-2 validation stdout copy: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse012_filtered_edge_em.stdout.tsv`
- Mouse0-2 validation stderr: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse012_filtered_edge_em.stderr.txt`
- Mouse0-2 validation time log: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/results/score_mouse012_filtered_edge_em.time.log`

## Provenance

- Code status: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/provenance/code_status.txt`
- Code diff: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/provenance/code_diff.patch`
- Code diff stat: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_context_level_ambiguity_em_prototype/provenance/code_diffstat.txt`
