#!/usr/bin/env bash
set -euo pipefail

cd /home/ubuntu/yihuiguang/tools/KSSD3mini

# Code provenance.
git --git-dir=.minco_git --work-tree=. rev-parse HEAD
git --git-dir=.minco_git --work-tree=. describe --always --dirty --tags
git --git-dir=.minco_git --work-tree=. status --short

# Inspect implementation of the selected default.
sed -n '1,280p' scripts/minco_profile_calibrated.py
sed -n '280,760p' scripts/minco_profile_calibrated.py
sed -n '1,220p' scripts/minco_profile_default.py
python3 -B scripts/minco_profile_default.py --help
sed -n '1,180p' minco_core/src/command_profile_wrapper.c

# Inspect source evidence records used in this decision note.
sed -n '1,220p' research/experiments/2026-06-26_cami2_hmp_unseen_transfer/results/default_strategy_readiness.tsv
sed -n '1,220p' research/experiments/2026-06-26_cami2_hmp_unseen_transfer/summary.tsv
sed -n '1,220p' research/experiments/2026-06-26_cami2_toymouse_current_default/summary.tsv
sed -n '1,220p' research/experiments/2026-06-26_cami2_toymouse_current_default/NOTE.md
sed -n '1,120p' research/experiments/2026-06-25_cami2_marine_samples3_5_minco_vs_sylph/summary_samples0_5.tsv
sed -n '1,120p' research/experiments/2026-06-25_cami2_plant_samples3_5_minco_vs_sylph/summary.tsv
sed -n '1,120p' research/experiments/2026-06-25_cami2_hmp_pilot_default_strategy/summary.tsv
sed -n '1,120p' research/experiments/2026-06-25_cami3_toygut_more_minco_vs_sylph/summary.tsv

# Rebuild the holdout-bundle audit and consolidated stage-decision tables from
# source metrics.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_cami3_gtdb_taxid_transfer.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_cami3_gtdb_source_readmap.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_marine_gtdb_taxid_transfer.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py --help
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_toymouse_current_refresh.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_plant_holdout_exactsplit.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_plant_gtdb_transfer_feasibility.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_strain_holdout_exactsplit.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_strain_gtdb_transfer_feasibility.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_edge_em_cross_domain.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cami3_source_readmap_abundance.py
python3 -B research/experiments/2026-06-26_cami2_hmp_unseen_transfer/search_fixed_call_abundance.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/cross_validate_cami3_loose_split_rescue.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/preflight_gtdb232_sylph_db.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_holdout_resources.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_block_exact_split_semantics.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_restricted_exact_feasibility.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/summary.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/decision_checks.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/objective_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/OBJECTIVE_COMPLETION_AUDIT.md
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_bundle_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_bundle_summary.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_db_preflight.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_db_build_command.sh
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_build_command.sh
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_hmp_profile_commands.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_hmp_profile_commands.sh
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_hmp_profile_commands.tsv
sed -n '1,100p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_hmp_profile_commands.sh
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_hmp_rescore_command.sh
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_hmp_rescore_command.sh
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_resource_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_resource_audit_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/block_exact_split_semantics.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/block_exact_split_semantics_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_restricted_exact_feasibility.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_restricted_exact_feasibility_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/runtime_strategy_cases.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/runtime_readiness.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/NEXT_RELEASE_GRADE_ACTIONS.md
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_gtdb_taxid_transfer_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_gtdb_taxid_transfer_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_gtdb_source_readmap_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_gtdb_source_readmap_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_gtdb_taxid_transfer_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_gtdb_taxid_transfer_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gtdb_source_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gtdb_source_abundance_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_scores.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_abundance_sweep_summary.tsv
sed -n '1,80p' research/experiments/2026-06-26_cami2_hmp_unseen_transfer/results/fixed_call_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/plant_holdout_exactsplit_current_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/plant_gtdb_transfer_feasibility_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/strain_holdout_exactsplit_current_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/strain_gtdb_transfer_feasibility_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/edge_em_cross_domain_methods.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/edge_em_cross_domain_policy_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_loose_split_rescue_crossval_summary.tsv
sed -n '1,80p' research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/alternative_read_ani_summary.tsv
sed -n '1,80p' research/experiments/2026-06-25_minco_vs_mashscreen_toymouse/ani_accuracy_summary.tsv

# Documentation and claim checks.
rg -n "universal-auto-exact|current_best_minco_default|LOW_EXTRA|low-extra|reported_ani|Ref_zip_aaf_ani|Sylph|HMP|marine|plant|profile default|conservative direct|recommended no-manual-strategy|sidecar beside --ref|joined_feature_training" \
  minco_core/src/command_profile_wrapper.c \
  minco_core/src/global_wrapper.c \
  scripts/minco_profile_calibrated.py \
  scripts/minco_profile_default.py \
  README.md \
  docs/USER_MANUAL.md \
  research/claims.tsv \
  research/experiments/2026-06-26_cami2_toymouse_current_default \
  research/experiments/2026-06-27_universal_strategy_decision

# Focused validation.
make minco
bin/minco profile --help
python3 -B -m py_compile \
  scripts/minco_profile_calibrated.py \
  scripts/minco_profile_default.py \
  tests/test_minco_profile_calibrated_auto_exact.py \
  tests/test_cami3_extension_after_recovery.py \
  research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py \
  research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py \
  research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py \
  research/experiments/2026-06-27_universal_strategy_decision/preflight_gtdb232_sylph_db.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_holdout_resources.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_block_exact_split_semantics.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_restricted_exact_feasibility.py \
  research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py \
  research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py \
  research/experiments/2026-06-27_universal_strategy_decision/write_current_default_strategy_manifest.py \
  research/experiments/2026-06-27_universal_strategy_decision/write_abundance_allocation_gap_manifest.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend_guard.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_alpha_sweep.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_cami3_gtdb_taxid_transfer.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_cami3_gtdb_source_readmap.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_marine_gtdb_taxid_transfer.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_toymouse_current_refresh.py \
  research/experiments/2026-06-27_universal_strategy_decision/summarize_plant_holdout_exactsplit.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_plant_gtdb_transfer_feasibility.py \
  research/experiments/2026-06-27_universal_strategy_decision/summarize_strain_holdout_exactsplit.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_strain_gtdb_transfer_feasibility.py \
  research/experiments/2026-06-27_universal_strategy_decision/summarize_edge_em_cross_domain.py \
  research/experiments/2026-06-27_universal_strategy_decision/sweep_cami3_source_readmap_abundance.py \
  research/experiments/2026-06-26_cami2_hmp_unseen_transfer/search_fixed_call_abundance.py \
  research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_call_filters.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_switch.py \
  research/experiments/2026-06-27_universal_strategy_decision/score_adaptive_call_filter_wrapper_validation.py \
  research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_external_exactsplit.py \
  research/experiments/2026-06-27_universal_strategy_decision/cross_validate_cami3_loose_split_rescue.py
python3 -B -m pytest -q tests/test_minco_profile_calibrated_auto_exact.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_cami3_extension_after_recovery.py -q
bash tests/smoke.sh

# Runtime scheduling benchmark for the calibrated default wrapper.
mkdir -p /tmp/minco_parallel_default_bench_20260627
/usr/bin/time -v -o /tmp/minco_parallel_default_bench_20260627/sample6_parallel.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 \
  --train-pool train12 \
  --scope bacteria \
  --workdir /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work \
  -p 16 \
  --report-all \
  -o /tmp/minco_parallel_default_bench_20260627/sample6_parallel.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel.time.log \
  /tmp/cami2_toymouse_current_default_20260626/logs/sample6_default.time.log
md5sum \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_unique.unfiltered.tsv \
  /tmp/cami2_toymouse_current_default_20260626/run/sample6_default/minco.best_diff_unique.unfiltered.tsv \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.unfiltered.tsv \
  /tmp/cami2_toymouse_current_default_20260626/run/sample6_default/minco.best_diff_split.unfiltered.tsv \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.exact.unfiltered.tsv \
  /tmp/cami2_toymouse_current_default_20260626/run/sample6_default/minco.best_diff_split.exact.unfiltered.tsv

# Experimental fused same-stream unique sidecar p16 full-sample benchmark.
# This produced identical raw unique/split/exact table MD5s to the concurrent
# scheduler but was slower, so p>=2 auto scheduling remains concurrent
# dual-pass.
mkdir -p /tmp/minco_same_stream_fused_bench_20260627
/usr/bin/time -v -o /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 \
  --train-pool train12 \
  --scope bacteria \
  --workdir /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused_work \
  -p 16 \
  --same-stream-readwise-passes \
  --report-all \
  -o /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused.time.log \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel.time.log \
  /tmp/cami2_toymouse_current_default_20260626/logs/sample6_default.time.log
md5sum \
  /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused_work/minco.best_diff_unique.unfiltered.tsv \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_unique.unfiltered.tsv \
  /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused_work/minco.best_diff_split.unfiltered.tsv \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.unfiltered.tsv \
  /tmp/minco_same_stream_fused_bench_20260627/sample6_same_stream_fused_work/minco.best_diff_split.exact.unfiltered.tsv \
  /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.exact.unfiltered.tsv

# Low-thread scheduler benchmark on a bounded first50k-read subset. This
# supports the auto -p1 sidecar scheduler.
mkdir -p /tmp/minco_scheduler_lowthread_20260627
set +o pipefail
gzip -cd /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz \
  | head -n 200000 \
  > /tmp/minco_scheduler_lowthread_20260627/sample6.first50k.fq
set -o pipefail
/usr/bin/time -v -o /tmp/minco_scheduler_lowthread_20260627/unique_p1.time.log bin/minco ani \
  -p 1 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /tmp/minco_scheduler_lowthread_20260627/sample6.first50k.fq \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-unique \
  --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0 \
  -o /tmp/minco_scheduler_lowthread_20260627/unique_p1.tsv
/usr/bin/time -v -o /tmp/minco_scheduler_lowthread_20260627/split_p1.time.log bin/minco ani \
  -p 1 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /tmp/minco_scheduler_lowthread_20260627/sample6.first50k.fq \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0 \
  -o /tmp/minco_scheduler_lowthread_20260627/split_p1.tsv
/usr/bin/time -v -o /tmp/minco_scheduler_lowthread_20260627/same_stream_p1.time.log bin/minco ani \
  -p 1 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /tmp/minco_scheduler_lowthread_20260627/sample6.first50k.fq \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0 \
  -o /tmp/minco_scheduler_lowthread_20260627/same_stream_split_p1.tsv \
  --readwise-unique-out /tmp/minco_scheduler_lowthread_20260627/same_stream_unique_p1.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/minco_scheduler_lowthread_20260627/unique_p1.time.log \
  /tmp/minco_scheduler_lowthread_20260627/split_p1.time.log \
  /tmp/minco_scheduler_lowthread_20260627/same_stream_p1.time.log
md5sum \
  /tmp/minco_scheduler_lowthread_20260627/unique_p1.tsv \
  /tmp/minco_scheduler_lowthread_20260627/same_stream_unique_p1.tsv \
  /tmp/minco_scheduler_lowthread_20260627/split_p1.tsv \
  /tmp/minco_scheduler_lowthread_20260627/same_stream_split_p1.tsv
column -ts $'\t' research/experiments/2026-06-27_universal_strategy_decision/results/same_stream_sidecar_benchmark.tsv

# Current-code HMP airskin refresh used by the decision table. Sample6 FASTQ
# was restored by streaming the selected BAMs from the CAMI II archive into
# samtools fastq, then sample6 and sample11 both used exact split tables.
mkdir -p /tmp/minco_current_code_hmp_refresh_20260627
/usr/bin/time -v -o /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample6_stream_bam_to_fastq_20260627.time.log \
  bash -c 'set -euo pipefail; curl -fL "$1" | tar -xzf - -T "$2" --to-command="samtools fastq -" | gzip -1 > "$3"' _ \
  https://frl.publisso.de/data/frl:6425518/airskinurogenital/sample_6.tar.gz \
  /tmp/cami2_hmp_unseen_transfer_20260626/nonzero_bam_paths_sample6.txt \
  /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample6.nonzero.fastq.gz
gzip -t /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample6.nonzero.fastq.gz
/usr/bin/time -v -o /tmp/minco_current_code_hmp_refresh_20260627/sample6_exact_split_rerun_20260627.time.log bin/minco ani \
  -p8 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample6.nonzero.fastq.gz \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani zip-aaf \
  -m0 \
  -f0 \
  -n0 \
  -t0 \
  --density-block-ctx 0 \
  -o /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_split_exact_unfiltered.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_refresh_20260627/sample11_exact_split.time.log bin/minco ani \
  -p8 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample11.nonzero.fastq.gz \
  --query-density ref \
  --abundance-est depth \
  --readwise-profile-only \
  --readwise-assign best-diff-split \
  --readwise-ani zip-aaf \
  -m0 \
  -f0 \
  -n0 \
  -t0 \
  --density-block-ctx 0 \
  -o /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_split_exact_unfiltered.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_refresh_20260627/sample6_table_refresh_exact_rerun_20260627.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  --unique-table /tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample6_unique_zip_unfiltered.tsv \
  --split-table /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_split_exact_unfiltered.tsv \
  --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv \
  --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 \
  --train-pool train12 \
  --scope bacteria \
  --report-all \
  -o /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_refresh_20260627/sample11_table_refresh_exact.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  --unique-table /tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample11_unique_zip_unfiltered.tsv \
  --split-table /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_split_exact_unfiltered.tsv \
  --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv \
  --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 \
  --train-pool train12 \
  --scope bacteria \
  --report-all \
  -o /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-method minco_current_code_hmp_refresh_exact6 \
  --output-prefix hmp_current_refresh_source_abundance
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_source_abundance_scores.tsv

# Same-release GTDB r232 Sylph build/profile/rescore provenance. The chunked
# build and HMP chunked profile commands were run once as long external jobs;
# do not rerun them casually because they write under /mnt/new3T and /tmp.
LOG=/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/chunked_build_driver.log
echo "START $(date -Is)" > "$LOG"
bash research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_build_command.sh >> "$LOG" 2>&1
echo "END $(date -Is) status=$?" >> "$LOG"
bash research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_hmp_profile_commands.sh
bash research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_hmp_rescore_command.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-method minco_current_code_hmp_refresh_exact6 \
  --output-prefix hmp_current_refresh_r232_source_abundance
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_quality.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_scores.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/gtdb232_sylph_chunked_build_runtime.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_r232_chunked_sylph_profile_runtime.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/diagnose_hmp_abundance_gap.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_r232_abundance_gap_variant_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_r232_abundance_gap_top_errors.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_r232_abundance_gap_mapping.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_read_context_ambiguity_outputs.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/read_context_ambiguity_output_audit.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_one_scan_feature_feasibility.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/one_scan_feature_feasibility.tsv

# Keep repo experiment records lightweight.
find research/experiments/2026-06-26_cami2_toymouse_current_default -type f -size +1M -print
find research/experiments/2026-06-27_universal_strategy_decision -type f -size +1M -print

# Model-cache and speed-priority calibrated mode benchmark on CAMI II Toy Mouse
# sample6. The cache file is temporary and should be packaged beside a released
# refdb rather than committed to the code repository.
mkdir -p /tmp/minco_model_cache_bench_20260627
/usr/bin/time -v -o /tmp/minco_model_cache_bench_20260627/write_cache.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py \
  --unique-table /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_unique.unfiltered.tsv \
  --split-table /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.unfiltered.tsv \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 \
  --train-pool train12 \
  --scope bacteria \
  --report-all \
  --write-model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  -o /tmp/minco_model_cache_bench_20260627/table_fit.tsv
/usr/bin/time -v -o /tmp/minco_model_cache_bench_20260627/table_cache.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py \
  --unique-table /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_unique.unfiltered.tsv \
  --split-table /tmp/minco_parallel_default_bench_20260627/sample6_parallel_work/minco.best_diff_split.unfiltered.tsv \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  --report-all \
  -o /tmp/minco_model_cache_bench_20260627/table_cache.tsv
/usr/bin/time -v -o /tmp/minco_model_cache_bench_20260627/sample6_fast_calls_cache.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  --strategy universal \
  --workdir /tmp/minco_model_cache_bench_20260627/sample6_fast_calls_cache_work \
  -p 16 \
  -o /tmp/minco_model_cache_bench_20260627/sample6_fast_calls_cache.tsv
/usr/bin/time -v -o /tmp/minco_model_cache_bench_20260627/sample6_fast_universal_cache.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz \
  --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  --strategy universal \
  --workdir /tmp/minco_model_cache_bench_20260627/sample6_fast_universal_cache_work \
  -p 16 \
  --report-all \
  -o /tmp/minco_model_cache_bench_20260627/sample6_fast_universal_cache.tsv
python3 -B - <<'PY'
from pathlib import Path
import importlib.util
import pandas as pd

ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
scorer_path = ROOT / "research/experiments/2026-06-26_cami2_toymouse_current_default/score_sample7_current_default.py"
spec = importlib.util.spec_from_file_location("toy_scorer", scorer_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
score = mod.score
by_accession, by_core = score.truth.load_gtdb_metadata(score.truth.GTDB_METADATA)
genome_taxid, genome_location = mod.load_meta()
score.EXP_DIR = Path("/tmp/minco_model_cache_bench_20260627")
score.BASE = mod.DATA_BASE
score.distribution_path = mod.distribution_path
sample = 6
source_rows = score.build_gtdb_truth(sample, by_accession, by_core, genome_taxid, genome_location)
profile = score.write_truth_files(sample, source_rows)
gold = {str(row["gtdb_species"]) for row in profile}
methods = [
    ("minco_fast_universal_cache_calls", Path("/tmp/minco_model_cache_bench_20260627/sample6_fast_calls_cache.tsv"), "minco"),
    ("minco_fast_universal_cache_report_all", Path("/tmp/minco_model_cache_bench_20260627/sample6_fast_universal_cache.tsv"), "minco"),
    ("minco_current_autoexact_parallel", Path("/tmp/minco_parallel_default_bench_20260627/sample6_parallel.tsv"), "minco"),
    ("sylph_gtdb_profile", Path("/tmp/cami2_toymouse_samples5_7_20260625/run/sylph_sample6/profile.tsv"), "sylph"),
]
rows = []
abund = []
for method, path, kind in methods:
    data = score.load_sylph_rows(path, by_accession, by_core) if kind == "sylph" else mod.score_calibrated_minco(path)
    selected = data.loc[data["active_gate_pass"]].copy()
    tp, fp, fn, precision, recall, f1 = score.score_sets(selected["gtdb_species"], gold)
    rows.append({"sample": sample, "method": method, "path": str(path), "gold_taxa": len(gold), "pred_taxa": len(tp | fp), "TP": len(tp), "FP": len(fp), "FN": len(fn), "precision": precision, "recall": recall, "F1": f1, "selected_rows": len(selected)})
    for rec in score.abundance_metrics(selected, profile):
        abund.append({"sample": sample, "method": method, **rec})
score_df = pd.DataFrame(rows)
abund_df = pd.DataFrame(abund)
merged = score_df.merge(abund_df.loc[abund_df["renorm_pred"].astype(str).str.lower().eq("true")], on=["sample", "method"], how="left")
merged.to_csv("/tmp/minco_model_cache_bench_20260627/sample6_fast_vs_sylph_score.tsv", sep="\t", index=False)
print(merged[["method", "TP", "FP", "FN", "precision", "recall", "F1", "l1_pct_points", "pearson", "selected_rows"]].to_string(index=False))
PY

# Exact-skip audit for CAMI II Toy Mouse samples5-7. This compares current-code
# cached table-mode calls from normal block split versus exact split, with the
# same unique table and same fitted RF/HGB model. The large generated profile
# tables stay under /tmp; only summary TSVs are copied into this note.
mkdir -p /tmp/minco_exact_skip_audit_20260627
for sample in 5 6 7; do
  for mode in block exact; do
    split=/tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_split.unfiltered.tsv
    if [ "$mode" = exact ]; then
      split=/tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_split.exact.unfiltered.tsv
    fi
    /usr/bin/time -v -o /tmp/minco_exact_skip_audit_20260627/mouse${sample}_${mode}.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py \
      --unique-table /tmp/cami2_toymouse_current_default_20260626/run/sample${sample}_default/minco.best_diff_unique.unfiltered.tsv \
      --split-table "$split" \
      --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv \
      --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
      --train-pool train12 \
      --scope bacteria \
      --strategy universal \
      -o /tmp/minco_exact_skip_audit_20260627/mouse${sample}_${mode}_universal.tsv
  done
done
python3 -B - <<'PY'
from pathlib import Path
import importlib.util
import pandas as pd

ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
scorer_path = ROOT / "research/experiments/2026-06-26_cami2_toymouse_current_default/score_sample7_current_default.py"
spec = importlib.util.spec_from_file_location("toy_scorer", scorer_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
score = mod.score
by_accession, by_core = score.truth.load_gtdb_metadata(score.truth.GTDB_METADATA)
genome_taxid, genome_location = mod.load_meta()
score.EXP_DIR = Path("/tmp/minco_exact_skip_audit_20260627")
score.BASE = mod.DATA_BASE
score.distribution_path = mod.distribution_path
rows = []
abund = []
diffs = []
for sample in [5, 6, 7]:
    source_rows = score.build_gtdb_truth(sample, by_accession, by_core, genome_taxid, genome_location)
    profile = score.write_truth_files(sample, source_rows)
    gold = {str(row["gtdb_species"]) for row in profile}
    method_paths = [
        ("minco_block_universal", Path(f"/tmp/minco_exact_skip_audit_20260627/mouse{sample}_block_universal.tsv"), "minco"),
        ("minco_exact_split_universal", Path(f"/tmp/minco_exact_skip_audit_20260627/mouse{sample}_exact_universal.tsv"), "minco"),
        ("sylph_gtdb_profile", Path(f"/tmp/cami2_toymouse_samples5_7_20260625/run/sylph_sample{sample}/profile.tsv"), "sylph"),
    ]
    callsets = {}
    loaded = {}
    for method, path, kind in method_paths:
        data = score.load_sylph_rows(path, by_accession, by_core) if kind == "sylph" else mod.score_calibrated_minco(path)
        loaded[method] = data
        selected = data.loc[data["active_gate_pass"]].copy()
        callsets[method] = set(selected["gtdb_species"].astype(str))
        tp, fp, fn, precision, recall, f1 = score.score_sets(selected["gtdb_species"], gold)
        rows.append({"sample": sample, "method": method, "path": str(path), "gold_taxa": len(gold), "pred_taxa": len(tp | fp), "TP": len(tp), "FP": len(fp), "FN": len(fn), "precision": precision, "recall": recall, "F1": f1, "selected_rows": len(selected)})
        for rec in score.abundance_metrics(selected, profile):
            abund.append({"sample": sample, "method": method, **rec})
    for kind, vals in [("block_only", sorted(callsets["minco_block_universal"] - callsets["minco_exact_split_universal"])), ("exact_only", sorted(callsets["minco_exact_split_universal"] - callsets["minco_block_universal"]))]:
        for sp in vals:
            data = loaded["minco_block_universal"] if kind == "block_only" else loaded["minco_exact_split_universal"]
            hit = data.loc[data["gtdb_species"].astype(str).eq(sp)]
            rec = {"sample": sample, "kind": kind, "gtdb_species": sp, "in_gold": sp in gold}
            if not hit.empty:
                r = hit.iloc[0]
                for col in ["calibrated_probability", "calibrated_abundance", "s_XnY_ctx_max", "s_ANI_max", "s_Ref_breadth_max", "s_Ref_mean_depth_max", "s_Ref_zip_af_max", "u_XnY_ctx_max", "u_Real_min_align_fraction_max"]:
                    if col in hit.columns:
                        rec[col] = r[col]
            diffs.append(rec)
score_df = pd.DataFrame(rows)
abund_df = pd.DataFrame(abund)
merged = score_df.merge(abund_df.loc[abund_df["renorm_pred"].astype(str).str.lower().eq("true")], on=["sample", "method"], how="left")
summary = []
for method, sub in merged.groupby("method"):
    tp = int(sub["TP"].sum())
    fp = int(sub["FP"].sum())
    fn = int(sub["FN"].sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    summary.append({"method": method, "samples": ",".join(map(str, sorted(sub["sample"].unique()))), "mean_F1": sub["F1"].mean(), "mean_L1_pp": sub["l1_pct_points"].mean(), "mean_Pearson": sub["pearson"].mean(), "pooled_TP": tp, "pooled_FP": fp, "pooled_FN": fn, "pooled_F1": f1})
outdir = Path("/tmp/minco_exact_skip_audit_20260627")
merged.to_csv(outdir / "toymouse_exact_skip_compare.tsv", sep="\t", index=False)
pd.DataFrame(summary).sort_values("mean_F1", ascending=False).to_csv(outdir / "toymouse_exact_skip_summary.tsv", sep="\t", index=False)
pd.DataFrame(diffs).to_csv(outdir / "toymouse_exact_skip_call_diffs.tsv", sep="\t", index=False)
PY
cp /tmp/minco_exact_skip_audit_20260627/toymouse_exact_skip_compare.tsv research/experiments/2026-06-27_universal_strategy_decision/results/toymouse_exact_skip_compare.tsv
cp /tmp/minco_exact_skip_audit_20260627/toymouse_exact_skip_summary.tsv research/experiments/2026-06-27_universal_strategy_decision/results/toymouse_exact_skip_summary.tsv
cp /tmp/minco_exact_skip_audit_20260627/toymouse_exact_skip_call_diffs.tsv research/experiments/2026-06-27_universal_strategy_decision/results/toymouse_exact_skip_call_diffs.tsv

# 2026-06-27 exact-split sidecar speed check.
make -j2
mkdir -p /tmp/minco_exact_split_sidecar_bench_20260627
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_bench_20260627/sidecar_50k.time.log bin/minco ani -p1 -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --qraw /tmp/toymouse_sample0_50000_reads.fq.gz --query-density ref --abundance-est depth --readwise-profile-only --readwise-assign best-diff-split --readwise-ani zip-aaf --density-block-ctx 100 -m0 -f0 -n0 -t0 --readwise-exact-split-out /tmp/minco_exact_split_sidecar_bench_20260627/exact_sidecar_50k.tsv -o /tmp/minco_exact_split_sidecar_bench_20260627/block_with_sidecar_50k.tsv
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_bench_20260627/exact_standalone_50k.time.log bin/minco ani -p1 -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --qraw /tmp/toymouse_sample0_50000_reads.fq.gz --query-density ref --abundance-est depth --readwise-profile-only --readwise-assign best-diff-split --readwise-ani zip-aaf --density-block-ctx 0 -m0 -f0 -n0 -t0 -o /tmp/minco_exact_split_sidecar_bench_20260627/exact_standalone_50k.tsv
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_bench_20260627/block_standalone_50k.time.log bin/minco ani -p1 -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --qraw /tmp/toymouse_sample0_50000_reads.fq.gz --query-density ref --abundance-est depth --readwise-profile-only --readwise-assign best-diff-split --readwise-ani zip-aaf --density-block-ctx 100 -m0 -f0 -n0 -t0 -o /tmp/minco_exact_split_sidecar_bench_20260627/block_standalone_50k.tsv
md5sum /tmp/minco_exact_split_sidecar_bench_20260627/exact_sidecar_50k.tsv /tmp/minco_exact_split_sidecar_bench_20260627/exact_standalone_50k.tsv
md5sum /tmp/minco_exact_split_sidecar_bench_20260627/block_with_sidecar_50k.tsv /tmp/minco_exact_split_sidecar_bench_20260627/block_standalone_50k.tsv
mkdir -p /tmp/minco_exact_split_sidecar_wrapper_20260627
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_wrapper_20260627/legacy_p1_50k.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --exact-split-trigger 999 --scope bacteria --workdir /tmp/minco_exact_split_sidecar_wrapper_20260627/legacy_work -p1 --report-all -o /tmp/minco_exact_split_sidecar_wrapper_20260627/legacy_p1_50k.tsv
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sidecar_p1_50k.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --exact-split-trigger 999 --scope bacteria --same-stream-exact-split --workdir /tmp/minco_exact_split_sidecar_wrapper_20260627/sidecar_work -p1 --report-all -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sidecar_p1_50k.tsv
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --scope bacteria --same-stream-exact-split --workdir /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work -p16 -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16.tsv
python3 -B - <<'PY'
from pathlib import Path
import importlib.util
import pandas as pd

ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
scorer_path = ROOT / "research/experiments/2026-06-26_cami2_toymouse_current_default/score_sample7_current_default.py"
spec = importlib.util.spec_from_file_location("toy_scorer", scorer_path)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
score = mod.score
sample = 6
minco_profile = Path("/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16.tsv")
sylph_profile = Path("/tmp/cami2_toymouse_samples5_7_20260625/run/sylph_sample6/profile.tsv")
by_accession, by_core = score.truth.load_gtdb_metadata(score.truth.GTDB_METADATA)
genome_taxid, genome_location = mod.load_meta()
score.EXP_DIR = ROOT / "research/experiments/2026-06-27_universal_strategy_decision/results"
score.BASE = mod.DATA_BASE
score.distribution_path = mod.distribution_path
source_rows = score.build_gtdb_truth(sample, by_accession, by_core, genome_taxid, genome_location)
profile = score.write_truth_files(sample, source_rows)
gold = {str(row["gtdb_species"]) for row in profile}
rows = []
for method, path, kind in [
    ("minco_sidecar_autoexact_sample6", minco_profile, "minco"),
    ("sylph_gtdb_profile_sample6", sylph_profile, "sylph"),
]:
    calls = score.load_sylph_rows(path, by_accession, by_core) if kind == "sylph" else mod.score_calibrated_minco(path)
    selected = calls.loc[calls["active_gate_pass"]].copy()
    tp, fp, fn, precision, recall, f1 = score.score_sets(selected["gtdb_species"], gold)
    renorm = [r for r in score.abundance_metrics(selected, profile) if str(r.get("renorm_pred")).lower() == "true"][0]
    rows.append({"sample": sample, "method": method, "path": str(path), "gold_taxa": len(gold), "pred_taxa": len(tp | fp), "TP": len(tp), "FP": len(fp), "FN": len(fn), "precision": precision, "recall": recall, "F1": f1, "l1_pct_points": renorm["l1_pct_points"], "pearson": renorm["pearson"]})
pd.DataFrame(rows).to_csv(ROOT / "research/experiments/2026-06-27_universal_strategy_decision/results/exact_split_sidecar_sample6_score.tsv", sep="\t", index=False)
PY

# 2026-06-27 low-extra exact-skip speed check.
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_onestream_p16.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --scope bacteria --same-stream-readwise-passes --same-stream-exact-split --workdir /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_onestream_p16_work -p16 -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_onestream_p16.tsv
/usr/bin/time -v -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_lowextra_skip_p16.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz --taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --scope bacteria --workdir /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_lowextra_skip_p16_work -p16 -o /tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_lowextra_skip_p16.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_exact_lowextra_skip.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_split_lowextra_skip_runtime.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_split_lowextra_skip_sample6_score.tsv

# 2026-06-27 abundance error decomposition from cached official GTDB truth/profile tables.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_error_decomposition_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_error_decomposition_delta.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_error_decomposition_validation.tsv

# 2026-06-27 cross-panel fixed-call abundance formula sweep.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_call_filters.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_calibrator.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_sample28_holdout.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_external_exactsplit.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_cached_exactsplit_diagnostics.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_abundance_variant_validation.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_abundance_variant_overall.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_variant_safety_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_abundance_switch_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_call_filter_validation.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_call_filter_safety_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_call_filter_overall.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_switch_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_switch_lopo.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/supervised_abundance_calibrator_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/supervised_abundance_sample28_holdout_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/supervised_abundance_external_exactsplit_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cached_exactsplit_diagnostic_audit.tsv

# 2026-06-27 wrapper-realistic experimental abundance blend validation.
# The wrapper default remains --abundance-genus-xny-blend-alpha 0.0; alpha 0.25
# is tested here only as an explicit experiment and is not promoted.
bash research/experiments/2026-06-27_universal_strategy_decision/run_abundance_blend_wrapper_validation.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_abundance_blend_wrapper_validation.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_blend_wrapper_validation_panel_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_blend_wrapper_validation_vs_offline.tsv

# 2026-06-29 wrapper-realistic experimental adaptive call-filter validation.
# The wrapper default remains --adaptive-call-filter-switch off; lopo-min-xny25
# is tested here only as an explicit experiment and is not promoted.
bash research/experiments/2026-06-27_universal_strategy_decision/run_adaptive_call_filter_wrapper_validation.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_adaptive_call_filter_wrapper_validation.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_wrapper_validation_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_wrapper_validation_panel_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_wrapper_validation_vs_offline.tsv

# 2026-06-29 external cached exact-split stress test for the same adaptive
# call-filter candidate. This is diagnostic/non-release but is outside the
# four GTDB panels that selected the rule.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_external_exactsplit.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_external_exactsplit_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_external_exactsplit_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_external_exactsplit_validation.tsv

# 2026-06-27 current HMP gastrooral table-mode universal check.
python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/cami2_hmp_pilot_20260625/run/minco_sample6_unique_zip_unfiltered.tsv --split-table /tmp/cami2_hmp_pilot_20260625/run/minco_sample6_split_zip_unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --train-features /mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2 --train-pool train12 --scope bacteria --strategy universal-auto-exact --out /tmp/cami2_hmp_pilot_20260625/run/minco_sample6_universal_autoexact_tablemode_current.tsv --report-all
python3 -B research/experiments/2026-06-25_cami2_hmp_pilot_default_strategy/score_hmp_pilot.py --run /tmp/cami2_hmp_pilot_20260625/run --truth /tmp/cami2_hmp_pilot_20260625/truth/taxonomic_profile_0.txt --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --sample-id 0 --sample-label gastrooral0 --calibrated-table /tmp/cami2_hmp_pilot_20260625/run/minco_sample0_universal_autoexact_tablemode_check.tsv --out /tmp/cami2_hmp_pilot_20260625/run/summary_sample0_current_universal_tablemode.tsv
python3 -B research/experiments/2026-06-25_cami2_hmp_pilot_default_strategy/score_hmp_pilot.py --run /tmp/cami2_hmp_pilot_20260625/run --truth /tmp/cami2_hmp_pilot_20260625/truth/taxonomic_profile_6.txt --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --sample-id 6 --sample-label gastrooral6 --calibrated-table /tmp/cami2_hmp_pilot_20260625/run/minco_sample6_universal_autoexact_tablemode_current.tsv --out /tmp/cami2_hmp_pilot_20260625/run/summary_sample6_current_universal_tablemode.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_hmp_gastrooral_current.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_current_universal_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_current_universal_mapping.tsv

# 2026-06-27 HMP gastrooral same-release chunked GTDB r232 Sylph profiles and source-abundance scorer.
mkdir -p /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample0 /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample6 /tmp/minco_hmp_gastrooral_r232_profile_inputs
sed '$a /tmp/cami2_hmp_pilot_20260625/run/sylph_sample0/gastrooral_sample0.nonzero.fastq.gz.sylsp' /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /tmp/minco_hmp_gastrooral_r232_profile_inputs/sample0.profile_inputs.list
sed '$a /tmp/cami2_hmp_pilot_20260625/run/sylph_sample6/gastrooral_sample6.nonzero.fastq.gz.sylsp' /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /tmp/minco_hmp_gastrooral_r232_profile_inputs/sample6.profile_inputs.list
/usr/bin/time -v -o /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample0/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /tmp/minco_hmp_gastrooral_r232_profile_inputs/sample0.profile_inputs.list -o /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample0/profile.chunked.tsv
/usr/bin/time -v -o /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample6/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /tmp/minco_hmp_gastrooral_r232_profile_inputs/sample6.profile_inputs.list -o /tmp/cami2_hmp_pilot_20260625/run_sylph_r232/sylph_sample6/profile.chunked.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_r232_source_abundance_quality.tsv

# 2026-06-27 HMP gastrooral raw-read current default rerun and source-abundance scorer.
mkdir -p /tmp/minco_current_code_hmp_gastrooral_20260627
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample0.nonzero.fastq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_work -p16 -o /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_20260627/sample6_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample6.nonzero.fastq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /tmp/minco_current_code_hmp_gastrooral_20260627/sample6_work -p16 -o /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --minco-sample0 /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv --minco-sample6 /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv --minco-method minco_current_raw_default_gastrooral_source_abundance --output-prefix hmp_gastrooral_raw_default_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_hmp_gastrooral_raw_default.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_raw_default_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_raw_default_runtime.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_raw_default_audit.tsv

# 2026-06-27 HMP gastrooral exact-needed speed substrate check.
# Replays sample0 without exact from cached block/unique tables, then reruns
# sample0 with --same-stream-exact-split to test whether sidecar output matches
# the exact rerun while avoiding the third FASTQ pass.
python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_work/minco.best_diff_split.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --report-all -o /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_block_universal_replay.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --minco-sample0 /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_block_universal_replay.tsv --minco-sample6 /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv --minco-method minco_hmp_gastrooral_sample0_block_universal_replay_sample6_current --output-prefix hmp_gastrooral_block_replay_r232_source_abundance
mkdir -p /tmp/minco_hmp_gastrooral_sidecar_speed_20260627
/usr/bin/time -v -o /tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_sidecar.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample0.nonzero.fastq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --scope bacteria --train-pool train12 --workdir /tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work -p16 --same-stream-exact-split --report-all -o /tmp/minco_hmp_gastrooral_sidecar_speed_20260627/minco_sample0_sidecar.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --minco-sample0 /tmp/minco_hmp_gastrooral_sidecar_speed_20260627/minco_sample0_sidecar.tsv --minco-sample6 /tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv --minco-method minco_hmp_gastrooral_sample0_sidecar_sample6_current --output-prefix hmp_gastrooral_sidecar_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_hmp_gastrooral_sidecar_speed.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_exact_sidecar_policy.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_exact_sidecar_speed.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_exact_sidecar_speed_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_sidecar_policy_instances.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_sidecar_policy_audit.tsv

# 2026-06-27 cheap prefix-read exact preflight audit.
mkdir -p /tmp/minco_exact_preflight_20260627
bash -lc 'gzip -dc /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample0.nonzero.fastq.gz | head -n 200000 > /tmp/minco_exact_preflight_20260627/gastrooral_sample0.first50k.fq'
bash -lc 'gzip -dc /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample6.nonzero.fastq.gz | head -n 200000 > /tmp/minco_exact_preflight_20260627/gastrooral_sample6.first50k.fq'
bash -lc 'gzip -dc /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample0.nonzero.fastq.gz | head -n 800000 > /tmp/minco_exact_preflight_20260627/gastrooral_sample0.first200k.fq'
bash -lc 'gzip -dc /tmp/cami2_hmp_pilot_20260625/reads/gastrooral_sample6.nonzero.fastq.gz | head -n 800000 > /tmp/minco_exact_preflight_20260627/gastrooral_sample6.first200k.fq'
/usr/bin/time -v -o /tmp/minco_exact_preflight_20260627/sample0_preflight_universal.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/minco_exact_preflight_20260627/gastrooral_sample0.first50k.fq --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --workdir /tmp/minco_exact_preflight_20260627/sample0_preflight_work -p4 --report-all -o /tmp/minco_exact_preflight_20260627/sample0_preflight_universal.tsv
/usr/bin/time -v -o /tmp/minco_exact_preflight_20260627/sample6_preflight_universal.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/minco_exact_preflight_20260627/gastrooral_sample6.first50k.fq --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --workdir /tmp/minco_exact_preflight_20260627/sample6_preflight_work -p4 --report-all -o /tmp/minco_exact_preflight_20260627/sample6_preflight_universal.tsv
/usr/bin/time -v -o /tmp/minco_exact_preflight_20260627/sample0_preflight200k_universal.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/minco_exact_preflight_20260627/gastrooral_sample0.first200k.fq --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --workdir /tmp/minco_exact_preflight_20260627/sample0_preflight200k_work -p4 --report-all -o /tmp/minco_exact_preflight_20260627/sample0_preflight200k_universal.tsv
/usr/bin/time -v -o /tmp/minco_exact_preflight_20260627/sample6_preflight200k_universal.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/minco_exact_preflight_20260627/gastrooral_sample6.first200k.fq --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --workdir /tmp/minco_exact_preflight_20260627/sample6_preflight200k_work -p4 --report-all -o /tmp/minco_exact_preflight_20260627/sample6_preflight200k_universal.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_exact_preflight_policy.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_preflight_policy_instances.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/exact_preflight_policy_audit.tsv

# 2026-06-27 HMP airskin source-truth candidate audit for the next same-release r232 sample.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_source_truth_candidate_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_source_truth_candidate_audit.tsv

# 2026-06-27 HMP airskin sample28 restore, current MinCO default, r232 Sylph, and source-abundance scoring.
/usr/bin/time -v -o /tmp/cami2_hmp_airskin_20260625/logs/sample28_stream_bam_to_fastq_20260627.time.log bash -c 'set -euo pipefail; mkdir -p "$4"; curl -fL "$1" | tar -xzf - -T "$2" --to-command="samtools fastq -" | gzip -1 > "$3"' _ https://frl.publisso.de/data/frl:6425518/airskinurogenital/sample_28.tar.gz /tmp/cami2_hmp_airskin_20260625/nonzero_bam_paths_sample28.txt /tmp/cami2_hmp_airskin_20260625/reads/airskinurogenital_sample28.nonzero.fastq.gz /tmp/cami2_hmp_airskin_20260625/logs
gzip -t /tmp/cami2_hmp_airskin_20260625/reads/airskinurogenital_sample28.nonzero.fastq.gz
/usr/bin/time -v -o /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_hmp_airskin_20260625/reads/airskinurogenital_sample28.nonzero.fastq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal-auto-exact --scope bacteria --train-pool train12 --workdir /tmp/minco_current_code_hmp_airskin28_20260627/work -p16 --report-all -o /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv
sed '$a /tmp/cami2_hmp_airskin_20260625/run/sylph_sample28/airskinurogenital_sample28.nonzero.fastq.gz.sylsp' /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /tmp/minco_hmp_airskin28_r232_profile_inputs/sample28.profile_inputs.list
/usr/bin/time -v -o /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /tmp/minco_hmp_airskin28_r232_profile_inputs/sample28.profile_inputs.list -o /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py --samples 28 --sylph-source r232 --r232-run /tmp/cami2_hmp_airskin_20260625/run_sylph_r232 --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv --minco-method minco_current_default_gtdb_source_abundance_sample28 --output-prefix hmp_airskin28_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py --samples 6,11,28 --sylph-source r232 --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv --minco-method minco_current_default_gtdb_source_abundance_exact6_sample28 --output-prefix hmp_current_refresh_r232_source_abundance_3sample
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin28_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_3sample_summary.tsv

# 2026-06-28 HMP airskin sample22 restore, current MinCO default, r232 Sylph,
# and source-abundance scoring. This fixed script preserves the exact sample22
# command after the moving next-action script advanced beyond sample22.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample22_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin22_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample22_6_11_28_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,90p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample22_stream_bam_to_fastq_20260627.time.log \
  /tmp/minco_current_code_hmp_airskin22_20260627/minco_sample22_current_default.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample22/sketch.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample22/profile.time.log

# 2026-06-28 HMP airskin sample5 restore, current MinCO default, r232 Sylph,
# source-abundance scoring, and regenerated decision/readiness tables. This
# fixed script preserves sample5 after the moving next-action script advanced
# to sample0.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample5_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 5,6,11,22,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample5_6_11_22_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample5_6_11_22_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_calibrator.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_sample28_holdout.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_supervised_abundance_external_exactsplit.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin5_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample5_6_11_22_28_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,90p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample5_stream_bam_to_fastq_20260627.time.log \
  /tmp/minco_current_code_hmp_airskin5_20260627/minco_sample5_current_default.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample5/sketch.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample5/profile.time.log

# 2026-06-28 HMP airskin sample0 restore, current MinCO default, r232 Sylph,
# source-abundance scoring, corrected aggregate rescoring, and regenerated
# decision/readiness tables. This fixed script preserves sample0 after the
# moving next-action script advanced to sample21.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample0_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 6,11 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_exact6 \
  --output-prefix hmp_current_refresh_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 6,11,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_exact6_sample28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_3sample
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 22,6,11,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample22_6_11_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample22_6_11_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 5,6,11,22,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample5_6_11_22_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample5_6_11_22_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 0,5,6,11,22,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample0_5_6_11_22_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample0_5_6_11_22_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin0_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample0_5_6_11_22_28_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,90p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample0_stream_bam_to_fastq_20260627.time.log \
  /tmp/minco_current_code_hmp_airskin0_20260627/minco_sample0_current_default.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample0/sketch.time.log \
  /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample0/profile.time.log

# 2026-06-28 non-destructive /tmp headroom audit before running the then-next
# release-grade sample21 command.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
sed -n '1,30p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_cleanup_candidates.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv

# 2026-06-28 completed alternate-root sample21 run. The moving next-action
# script later advanced to sample18, so sample21 is preserved in a fixed script.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample21_release_grade_commands.sh

# The first sample21 run completed restore/profile but hit an older scorer path
# assumption for raw unique/split tables under /tmp. After adding generic
# --minco-unique-sample and --minco-split-sample overrides, scoring was rerun
# from the completed alternate-root artifacts without rerunning the heavy
# restore/profile steps.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 21 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv \
  --minco-unique-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 21=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample21 \
  --output-prefix hmp_airskin21_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 21,0,5,6,11,22,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv \
  --minco-unique-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 21=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.tsv \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample21_0_5_6_11_22_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample21_0_5_6_11_22_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin21_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample21_0_5_6_11_22_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample21_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.time.log

# 2026-06-28 completed alternate-root sample18 run. The moving next-action
# script later advanced to sample13, so sample18 is preserved in a fixed script.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
cp research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh \
  research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample18_release_grade_commands.sh
chmod 755 research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample18_release_grade_commands.sh

# Corrected aggregate including both sample18 and sample21.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gtdb_source_abundance.py \
  --samples 18,21,0,5,6,11,22,28 \
  --sylph-source r232 \
  --r232-run /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232 \
  --minco-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/minco_sample18_current_default.tsv \
  --minco-unique-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 18=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 18=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample18/profile.tsv \
  --minco-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv \
  --minco-unique-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_unique.unfiltered.tsv \
  --minco-split-sample 21=/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_split.unfiltered.tsv \
  --sylph-sample 21=/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.tsv \
  --sylph-sample28 /tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv \
  --minco-sample6 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv \
  --minco-sample11 /tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv \
  --minco-sample28 /tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv \
  --minco-method minco_current_default_gtdb_source_abundance_sample18_21_0_5_6_11_22_28 \
  --output-prefix hmp_current_refresh_r232_source_abundance_sample18_21_0_5_6_11_22_28
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin18_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample18_21_0_5_6_11_22_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample18_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/minco_sample18_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample18/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample18/profile.time.log

# 2026-06-28 completed alternate-root sample13 run. The moving next-action
# script later advanced to sample25, so sample13 is preserved in a fixed script.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample13_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin13_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample13_0_5_6_11_18_21_22_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample13_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/minco_sample13_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample13/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample13/profile.time.log

# 2026-06-28 completed alternate-root sample25 run. The moving next-action
# script later advanced to sample3, so sample25 is preserved in a fixed script.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample25_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin25_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample25_0_5_6_11_13_18_21_22_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample25_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/minco_sample25_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample25/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample25/profile.time.log

# 2026-06-28 completed alternate-root sample3 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample1.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample3_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_call_filters.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin3_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample3_0_5_6_11_13_18_21_22_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample3_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/minco_sample3_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample3/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample3/profile.time.log

# 2026-06-28 completed alternate-root sample1 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample17.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample1_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin1_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample1_0_3_5_6_11_13_18_21_22_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample1_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/minco_sample1_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample1/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample1/profile.time.log

# 2026-06-28 completed alternate-root sample17 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample24.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample17_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_call_filters.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin17_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample17_0_1_3_5_6_11_13_18_21_22_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample17_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/minco_sample17_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample17/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample17/profile.time.log

# 2026-06-28 completed alternate-root sample24 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample16.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample24_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin24_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample24_0_1_3_5_6_11_13_17_18_21_22_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample24_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/minco_sample24_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample24/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample24/profile.time.log

# 2026-06-28 completed alternate-root sample16 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample15.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample16_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin16_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample16_0_1_3_5_6_11_13_17_18_21_22_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample16_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/minco_sample16_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample16/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample16/profile.time.log

# 2026-06-28 completed alternate-root sample15 run. The original run used the
# then-current moving next-action script; reruns should use the fixed script
# because the moving script later advanced to sample4.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample15_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin15_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample15_0_1_3_5_6_11_13_16_17_18_21_22_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample15_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin15_20260627/minco_sample15_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample15/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample15/profile.time.log

# 2026-06-28 completed alternate-root sample4 run. The fixed script preserves
# the sample4 command after the moving next-action script advanced to sample23.
cp research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh \
  research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample4_release_grade_commands.sh
chmod 755 research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample4_release_grade_commands.sh
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample4_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin4_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample4_0_1_3_5_6_11_13_15_16_17_18_21_22_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample4_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin4_20260627/minco_sample4_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample4/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample4/profile.time.log

# 2026-06-28 completed alternate-root sample23 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample20.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample23_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin23_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample23_0_1_3_4_5_6_11_13_15_16_17_18_21_22_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample23_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin23_20260627/minco_sample23_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample23/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample23/profile.time.log

# 2026-06-28 completed alternate-root sample20 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample7.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample20_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin20_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample20_0_1_3_4_5_6_11_13_15_16_17_18_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample20_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin20_20260627/minco_sample20_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample20/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample20/profile.time.log

# 2026-06-28 completed alternate-root sample7 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample14.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample7_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin7_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample7_0_1_3_4_5_6_11_13_15_16_17_18_20_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample7_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin7_20260627/minco_sample7_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample7/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample7/profile.time.log

# 2026-06-28 completed alternate-root sample14 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample9.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample14_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin14_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample14_0_1_3_4_5_6_7_11_13_15_16_17_18_20_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample14_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin14_20260627/minco_sample14_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample14/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample14/profile.time.log

# 2026-06-28 completed alternate-root sample9 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample10.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample9_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin9_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample9_0_1_3_4_5_6_7_11_13_14_15_16_17_18_20_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample9_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin9_20260627/minco_sample9_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample9/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample9/profile.time.log

# 2026-06-28 completed alternate-root sample10 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample19.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample10_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin10_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample10_0_1_3_4_5_6_7_9_11_13_14_15_16_17_18_20_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample10_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin10_20260627/minco_sample10_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample10/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample10/profile.time.log

# 2026-06-29 completed alternate-root sample19 run. Reruns should use the
# fixed script because the moving next-action script later advanced to sample26.
MINCO_HMP_WORK_ROOT=/mnt/new3T/minco_release_holdouts_20260628 \
  bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_sample19_release_grade_commands.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/decompose_abundance_errors.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_abundance_variants.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_abundance_variant_safety.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_abundance_switch.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_runtime_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin19_r232_source_abundance_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28_summary.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
rg -n "Elapsed|Maximum resident|User time|System time|Command being timed" \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/logs/sample19_stream_bam_to_fastq_20260627.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin19_20260627/minco_sample19_current_default.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample19/sketch.time.log \
  /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample19/profile.time.log

# 2026-06-29 fixed-call truth-aware abundance oracle diagnostic.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_abundance_oracle_bounds.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_oracle_bounds_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_oracle_bounds_panel_delta.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv

# 2026-06-29 emitted-profile missed-truth candidate-surface diagnostic.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_missed_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/missed_truth_candidate_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/missed_truth_candidate_panel_summary.tsv
sed -n '1,30p' research/experiments/2026-06-27_universal_strategy_decision/results/missed_truth_candidate_top.tsv

# 2026-06-29 HMP raw-table follow-up for profile-absent high-abundance misses.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_missed_truth_raw_tables.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_missed_truth_raw_table_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_missed_truth_raw_table_panel_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_missed_truth_raw_table_detail.tsv

# 2026-06-29 HMP raw-side-channel rescue threshold sweep.
# Full per-sample score grids are intentionally written to /tmp.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_hmp_raw_candidate_rescue.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_raw_candidate_rescue_audit.tsv
sed -n '1,12p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_raw_candidate_rescue_top100.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_raw_candidate_rescue_runtime.tsv

# 2026-06-29 cross-panel zero-mass candidate rescue sweep.
# Full per-sample score grids are intentionally written to /tmp.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_candidate_rescue.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_audit.tsv
sed -n '1,14p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_top100.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_validation.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_runtime.tsv

# 2026-06-29 detail audit for the best cross-panel candidate rescue rule.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cross_panel_candidate_rescue_hits.py
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_best_hits_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_best_hits_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_best_hits.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_rescue_best_hits_runtime.tsv

# 2026-06-29 emitted-profile candidate rescue wrapper replay.
/usr/bin/time -v bash research/experiments/2026-06-27_universal_strategy_decision/run_candidate_rescue_wrapper_validation.sh
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_candidate_rescue_wrapper_validation.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_rescue_wrapper_validation_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_rescue_wrapper_validation_panel_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_rescue_wrapper_validation_runtime.tsv

# 2026-06-29 raw-side rescue candidate taxid/accession collapse diagnostic.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_raw_rescue_taxid_collapse.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_rescue_taxid_collapse_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_rescue_taxid_collapse_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_rescue_taxid_collapse_runtime.tsv

# 2026-06-29 focused accession-level raw-side candidate surface replay.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_raw_side_candidate_surface_gastrooral.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_gastrooral_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_gastrooral_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_gastrooral_vs_offline.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_gastrooral_runtime.tsv

# 2026-06-29 all-HMP accession-level raw-side candidate surface replay.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_raw_side_candidate_surface_hmp.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_hmp_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_hmp_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_hmp_vs_offline.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_side_candidate_surface_hmp_runtime.tsv

# 2026-06-29 focused integrated wrapper candidate-surface replay.
mkdir -p /tmp/minco_candidate_surface_integrated_gastrooral_20260629
/usr/bin/time -v -o /tmp/minco_candidate_surface_integrated_gastrooral_20260629/sample0.gtdb_surface.time.log python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample0_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --candidate-surface-switch accession-ani90-xny100-br01-af70 --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv -o /tmp/minco_candidate_surface_integrated_gastrooral_20260629/minco_sample0_candidate_surface_gtdb.tsv
/usr/bin/time -v -o /tmp/minco_candidate_surface_integrated_gastrooral_20260629/sample6.gtdb_surface.time.log python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample6_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_20260627/sample6_work/minco.best_diff_split.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --strategy universal --scope bacteria --train-pool train12 --candidate-surface-switch accession-ani90-xny100-br01-af70 --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv -o /tmp/minco_candidate_surface_integrated_gastrooral_20260629/minco_sample6_candidate_surface_gtdb.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --minco-sample0 /tmp/minco_candidate_surface_integrated_gastrooral_20260629/minco_sample0_candidate_surface_gtdb.tsv --minco-sample6 /tmp/minco_candidate_surface_integrated_gastrooral_20260629/minco_sample6_candidate_surface_gtdb.tsv --minco-method minco_integrated_candidate_surface_gtdb_gastrooral --output-prefix candidate_surface_integrated_gtdb_gastrooral
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_gtdb_gastrooral_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_gtdb_gastrooral_summary.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_gtdb_gastrooral_runtime.tsv

# 2026-06-29 all-HMP integrated wrapper candidate-surface replay.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_integrated_candidate_surface_hmp.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_vs_postprocessor.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_vs_offline.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_outer_runtime.tsv

# 2026-06-29 HMP-only abundance sweep for integrated candidate-surface rows.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_candidate_surface_abundance_policy.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_abundance_policy_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_abundance_policy_overall.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_abundance_policy_runtime.tsv

# 2026-06-29 cross-panel abundance sweep for the selected candidate call set.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_candidate_abundance_policy.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_policy_audit.tsv
sed -n '1,25p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_policy_overall.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_policy_runtime.tsv

# 2026-06-29 HMP integrated wrapper validation for normalized-depth candidate abundance.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_candidate_surface_abundance_policy.py
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_integrated_candidate_surface_hmp_abundance_policy.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_abundance_policy_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_abundance_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_abundance_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_abundance_vs_policy_sweep.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_surface_integrated_hmp_abundance_outer_runtime.tsv

# 2026-06-29 cross-panel wrapper validation for normalized-depth candidate abundance.
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_cross_panel_candidate_abundance_wrapper.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_wrapper_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_wrapper_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_wrapper_vs_zero.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_wrapper_vs_policy_sweep.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_candidate_abundance_wrapper_outer_runtime.tsv

# 2026-06-29 candidate default comparison against current MinCO and Sylph.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_candidate_abundance_default_decision.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_default_decision_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_default_vs_current.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_default_vs_sylph.tsv

# 2026-06-29 selected candidate preset launcher implementation validation.
python3 -B -m py_compile scripts/minco_profile_default.py scripts/minco_profile_calibrated.py tests/test_cami3_extension_after_recovery.py research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py research/experiments/2026-06-27_universal_strategy_decision/write_current_default_strategy_manifest.py research/experiments/2026-06-27_universal_strategy_decision/write_abundance_allocation_gap_manifest.py research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_callset_oracle_feasibility.py research/experiments/2026-06-27_universal_strategy_decision/sweep_selected_call_mass_transforms.py research/experiments/2026-06-27_universal_strategy_decision/audit_selected_call_mass_transform_guard.py research/experiments/2026-06-27_universal_strategy_decision/sweep_selected_call_feature_allocators.py research/experiments/2026-06-27_universal_strategy_decision/audit_selected_call_feature_allocator_guard.py research/experiments/2026-06-27_universal_strategy_decision/validate_feature_allocator_wrapper_parity.py research/experiments/2026-06-27_universal_strategy_decision/validate_refined_feature_allocator_score_replay.py research/experiments/2026-06-27_universal_strategy_decision/validate_refined_feature_allocator_marine_diagnostic.py research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend.py research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend_guard.py research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_alpha_sweep.py
python3 -B -m pytest -q tests/test_minco_profile_calibrated_auto_exact.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_cami3_extension_after_recovery.py -q
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/decision_checks.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv

# 2026-06-29 compact selected-default strategy manifest.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_current_default_strategy_manifest.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/current_default_strategy_manifest.tsv

# 2026-06-29 one-flag candidate preset cached replay.
# A first report-all replay filled /tmp on HMP airskin sample25; the validator
# now uses normal called-row output, which is sufficient for scoring.
rm -rf /tmp/minco_candidate_preset_replay_20260629
/usr/bin/time -v python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_candidate_preset_replay.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/summarize_candidate_abundance_default_decision.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_replay_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_replay_summary.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_default_decision_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/release_readiness.tsv

# 2026-06-29 selected-default abundance blend rejection audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_blend_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_blend_panel_delta.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_blend_guard.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_blend_guard_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_blend_guard_lopo.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_preset_genus_xny_alpha_sweep.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_alpha_sweep_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_preset_genus_xny_alpha_sweep_overall.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_allocation_gap_manifest.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_allocation_gap_manifest.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_ALLOCATION_GAP.md
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_strategy_decision_matrix.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 selected candidate call-set abundance oracle feasibility.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_candidate_callset_oracle_feasibility.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_callset_oracle_feasibility_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/candidate_callset_oracle_feasibility_summary.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_NEXT_TARGET.md

# 2026-06-29 selected call-set base-row mass-transform sweep.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_selected_call_mass_transforms.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_mass_transform_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_mass_transform_overall.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_MASS_TRANSFORM_SWEEP.md

# 2026-06-29 selected call-set guarded mass-transform audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_selected_call_mass_transform_guard.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_mass_transform_guard_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_mass_transform_guard_lopo.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_MASS_TRANSFORM_GUARD.md
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 selected call-set feature-derived allocator sweep.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_selected_call_feature_allocators.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_feature_allocator_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_feature_allocator_overall.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_SWEEP.md
sed -n '1,180p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 selected call-set feature-allocator guard audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_selected_call_feature_allocator_guard.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_feature_allocator_guard_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/selected_call_feature_allocator_guard_lopo.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_GUARD.md
sed -n '1,200p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 implemented guarded feature-allocator wrapper parity.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_feature_allocator_wrapper_parity.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_wrapper_parity_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_wrapper_parity_summary.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_WRAPPER_PARITY.md
sed -n '1,220p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 guarded feature-allocator external exact-split diagnostic.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_feature_allocator_external_exactsplit.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_external_exactsplit_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_external_exactsplit_summary.tsv
sed -n '1,180p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_EXTERNAL_EXACTSPLIT.md
sed -n '1,240p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_STRATEGY_DECISION.md

# 2026-06-29 consolidated current-default release snapshot.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_default_release_snapshot.py
sed -n '1,100p' research/experiments/2026-06-27_universal_strategy_decision/results/default_release_snapshot_audit.tsv
sed -n '1,100p' research/experiments/2026-06-27_universal_strategy_decision/results/default_release_snapshot_panel_metrics.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/DEFAULT_RELEASE_SNAPSHOT.md

# 2026-06-29 universal strategy release-gate validator.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_release_blocker_audit.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_RELEASE_BLOCKER.md
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_RELEASE_GATE.md

# 2026-06-29 holdout manifest refresh and gap action plan.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_default_release_snapshot.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_bundle_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_gap_action_plan.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/HOLDOUT_GAP_ACTION_PLAN.md

# 2026-06-29 CAMI3 source-readmap scope audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_scope.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_scope_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_scope_summary.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_SCOPE_AUDIT.md

# 2026-06-29 CAMI3 source-readmap resolver-candidate audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_resolver_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_resolver_candidate_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_resolver_candidate_summary.tsv
sed -n '1,180p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_RESOLVER_CANDIDATES.md

# 2026-06-29 CAMI3 source-readmap exact-binomial fallback rescore.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/rescore_cami3_source_readmap_binomial_fallback.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_binomial_fallback_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_binomial_fallback_summary.tsv
sed -n '1,220p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_BINOMIAL_FALLBACK_RESCORING.md

# 2026-06-29 CAMI3 selected-default binomial-fallback rescore and truth-policy audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/rescore_cami3_binomial_fallback_candidate_default.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_binomial_fallback_truth_policy.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_default_release_snapshot.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_binomial_fallback_candidate_default_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_binomial_fallback_truth_policy_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/holdout_bundle_summary.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/results/default_release_snapshot_audit.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_RELEASE_GATE.md

# 2026-06-29 CAMI3 source-readmap extension cache audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_extension_cache.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_extension_cache_audit.tsv
sed -n '150,175p' research/experiments/2026-06-27_universal_strategy_decision/results/decision_checks.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/OBJECTIVE_COMPLETION_AUDIT.md
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_EXTENSION_CACHE_AUDIT.md
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/HOLDOUT_GAP_ACTION_PLAN.md

# 2026-06-29 marine exact-binomial truth-transfer feasibility audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_binomial_truth_transfer.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_binomial_transfer_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_binomial_transfer_quality.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/MARINE_BINOMIAL_TRANSFER_AUDIT.md

# 2026-06-29 guarded feature-allocator release-candidate audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_feature_allocator_release_candidate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_release_candidate_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/abundance_release_blocker_audit.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_RELEASE_CANDIDATE.md
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_RELEASE_BLOCKER.md
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_RELEASE_GATE.md

# 2026-06-29 refined guarded feature-allocator candidate.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_feature_allocator_combined_guard_refinement.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_refined_feature_allocator_guard_parity.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_refined_feature_allocator_score_replay.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_refined_feature_allocator_marine_diagnostic.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_combined_guard_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_guard_parity_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_score_replay_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_marine_diagnostic_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_independent_holdout_inventory_audit.tsv
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_COMBINED_GUARD_REFINEMENT.md
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_REFINED_GUARD_PARITY.md
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_REFINED_SCORE_REPLAY.md
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_REFINED_MARINE_DIAGNOSTIC.md
sed -n '1,140p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_REFINED_HOLDOUT_INVENTORY.md
sed -n '1,150p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_RELEASE_BLOCKER.md

# 2026-06-29 refined allocator release-holdout recovery plan.
/usr/bin/time -v -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample3.candidate_preset.time.log python3 -B scripts/minco_profile_default.py \
  --profile-preset candidate \
  --unique-table /tmp/cami3_toy_human_gut_20260626/run/sample3_autoexact/minco.best_diff_unique.unfiltered.tsv \
  --split-table /tmp/cami3_toy_human_gut_20260626/run/sample3_autoexact/minco.best_diff_split.unfiltered.tsv \
  --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample3.candidate_preset.tsv
/usr/bin/time -v -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample4.candidate_preset.time.log python3 -B scripts/minco_profile_default.py \
  --profile-preset candidate \
  --unique-table /tmp/cami3_toy_human_gut_20260626/run/sample4_autoexact/minco.best_diff_unique.unfiltered.tsv \
  --split-table /tmp/cami3_toy_human_gut_20260626/run/sample4_autoexact/minco.best_diff_split.unfiltered.tsv \
  --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample4.candidate_preset.tsv
/usr/bin/time -v -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample5.candidate_preset.time.log python3 -B scripts/minco_profile_default.py \
  --profile-preset candidate \
  --unique-table /tmp/cami3_toy_human_gut_20260626/run/sample5_autoexact/minco.best_diff_unique.unfiltered.tsv \
  --split-table /tmp/cami3_toy_human_gut_20260626/run/sample5_autoexact/minco.best_diff_split.unfiltered.tsv \
  --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv \
  --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib \
  --train-pool train12 \
  --scope bacteria \
  -o /tmp/minco_candidate_preset_replay_20260629/cami3_toy_human_gut_gtdb_source_readmap.sample5.candidate_preset.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_extension_cache.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_refined_allocator_release_holdout_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_holdout_recovery_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_holdout_recovery_runtime.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/ABUNDANCE_FEATURE_ALLOCATOR_REFINED_HOLDOUT_RECOVERY_PLAN.md

# 2026-06-29 CAMI3 samples3-5 local source-profile fallback diagnostic.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_cami3_source_profile_binomial_extension.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_refined_allocator_release_holdout_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_profile_binomial_extension_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_profile_binomial_extension_quality.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_profile_binomial_extension_summary.tsv
sed -n '1,180p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_PROFILE_BINOMIAL_EXTENSION.md

# 2026-06-29 CAMI3 samples3-5 source-readmap recovery input audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_recovery_inputs.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_cami3_source_readmap_extension_after_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_refined_allocator_release_holdout_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_abundance_strategy_decision_matrix.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_recovery_inputs_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_recovery_inputs.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_readmap_only_recovery_commands.sh
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_RECOVERY_INPUTS.md
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_extension_after_recovery_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_extension_after_recovery_inputs.tsv
sed -n '1,180p' research/experiments/2026-06-27_universal_strategy_decision/CAMI3_SOURCE_READMAP_EXTENSION_AFTER_RECOVERY.md

# 2026-06-29 CAMI3 samples3-5 readmap-only recovery and post-recovery scoring.
# The recovery script streams remote archives and extracts only
# reads_mapping.tsv.gz for each sample; it does not store full archives.
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_readmap_only_recovery_commands.sh
bash research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_readmap_only_recovery_commands.sh
gzip -t /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_3_reads/reads_mapping.tsv.gz
gzip -t /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_4_reads/reads_mapping.tsv.gz
gzip -t /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_5_reads/reads_mapping.tsv.gz
du -ch \
  /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_3_reads/reads_mapping.tsv.gz \
  /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_4_reads/reads_mapping.tsv.gz \
  /mnt/new3T/minco_cami3_toygut_extra_20260621/sample_5_reads/reads_mapping.tsv.gz
df -h /mnt/new3T /tmp
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_cami3_source_readmap_recovery_inputs.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_cami3_source_readmap_extension_after_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_refined_allocator_release_holdout_recovery.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_extension_after_recovery_summary.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/cami3_source_readmap_extension_after_recovery_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/feature_allocator_refined_holdout_recovery_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-30 wrapper-output replay for combined candidate surface mode.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_raw_retention_surface_wrapper_cross_panel.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_wrapper_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_wrapper_overall.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_wrapper_panel_summary.tsv

# 2026-06-30 max-called-species guard for the combined candidate surface.
for sample in 3 4 5; do
  /usr/bin/time -v -o /tmp/minco_raw_retention_surface_maxcalls_guard_20260630.marine${sample}.time.log python3 -B scripts/minco_profile_calibrated.py --unique-table /mnt/new3T/minco_marine_selected_default_20260630/sample${sample}_work/minco.best_diff_unique.unfiltered.tsv --split-table /mnt/new3T/minco_marine_selected_default_20260630/sample${sample}_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --train-pool train12 --strategy universal-auto-exact --candidate-rescue-switch emitted-ani90-xny100-br01-af70 --candidate-surface-switch accession-current-or-ani93-xny650-br20 --candidate-surface-max-called-species 250 --candidate-abundance-policy normalized-depth-alpha2 -o /tmp/minco_raw_retention_surface_maxcalls_guard_20260630.marine${sample}.tsv
done
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_raw_retention_surface_maxcalls_guard.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_maxcalls_guard_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_maxcalls_guard_sweep.tsv

# 2026-06-30 full wrapper-output replay for max-called-species guard.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_raw_retention_surface_maxcalls_guard_wrapper_cross_panel.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_maxcalls_guard_wrapper_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_maxcalls_guard_wrapper_overall.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_retention_surface_maxcalls_guard_wrapper_panel_summary.tsv

# 2026-06-30 opt-in raw-retention candidate-surface implementation smoke.
python3 -B -m py_compile scripts/minco_profile_calibrated.py scripts/minco_profile_default.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_minco_profile_calibrated_auto_exact.py -q
python3 scripts/minco_profile_calibrated.py --help | rg -n "candidate-surface-switch|accession-ani93"
mkdir -p /tmp/minco_raw_retention_surface_smoke_20260630
/usr/bin/time -v python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/cami2_toymouse_current_default_20260626/run/sample5_default/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/cami2_toymouse_current_default_20260626/run/sample5_default/minco.best_diff_split.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --candidate-surface-taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --strategy universal-auto-exact --candidate-surface-switch accession-ani93-xny650-br20 --candidate-abundance-policy zero -o /tmp/minco_raw_retention_surface_smoke_20260630/toymouse_sample5_ani93_xny650_br20_zero.tsv
python3 -B - <<'PY'
import pandas as pd
p='/tmp/minco_raw_retention_surface_smoke_20260630/toymouse_sample5_ani93_xny650_br20_zero.tsv'
df=pd.read_csv(p, sep='\t')
added=df.get('candidate_surface_added', pd.Series(False, index=df.index)).astype(str).str.lower().eq('true')
called=df.get('calibrated_call', pd.Series(False, index=df.index)).astype(str).str.lower().eq('true')
print('rows', len(df))
print('called', int(called.sum()))
print('surface_added', int(added.sum()))
print('surface_added_called', int((added & called).sum()))
print('surface_added_abundance_raw_sum', float(pd.to_numeric(df.loc[added, 'calibrated_abundance_raw'], errors='coerce').fillna(0).sum()) if added.any() else 0.0)
for col in ['candidate_surface_switch','candidate_surface_ani_min','candidate_surface_xny_min','candidate_surface_breadth_min','candidate_surface_real_af_min']:
    if col in df.columns:
        print(col, sorted(df[col].dropna().astype(str).unique())[:5])
PY

# 2026-06-30 broad cached raw-candidate retention replay across current panels.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_raw_candidate_retention_cross_panel.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_candidate_retention_cross_panel_audit.tsv
sed -n '1,25p' research/experiments/2026-06-27_universal_strategy_decision/results/raw_candidate_retention_cross_panel_overall.tsv

# 2026-06-30 marine setup+assembly truth profile rescore and route closure.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_marine_setup_truth_profile_rescore.py --samples 3,4,5 --profile-set-label cached_existing --minco-method minco_cached_s1000_unique_zipaaf_setup_truth --sylph-method sylph_cached_marine_setup_truth
mkdir -p /mnt/new3T/minco_marine_selected_default_20260630/sample3_work /mnt/new3T/minco_marine_selected_default_20260630/logs && /usr/bin/time -v -o /mnt/new3T/minco_marine_selected_default_20260630/logs/sample3_minco_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py --profile-preset candidate -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_marine_samples3_5_20260625/data/cami_marine_sample3_reads.fq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /mnt/new3T/minco_marine_selected_default_20260630/sample3_work -p16 -o /mnt/new3T/minco_marine_selected_default_20260630/minco_sample3_current_default.tsv
mkdir -p /mnt/new3T/minco_marine_selected_default_20260630/sample4_work /mnt/new3T/minco_marine_selected_default_20260630/logs && /usr/bin/time -v -o /mnt/new3T/minco_marine_selected_default_20260630/logs/sample4_minco_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py --profile-preset candidate -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_marine_samples3_5_20260625/data/cami_marine_sample4_reads.fq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /mnt/new3T/minco_marine_selected_default_20260630/sample4_work -p16 -o /mnt/new3T/minco_marine_selected_default_20260630/minco_sample4_current_default.tsv
mkdir -p /mnt/new3T/minco_marine_selected_default_20260630/sample5_work /mnt/new3T/minco_marine_selected_default_20260630/logs && /usr/bin/time -v -o /mnt/new3T/minco_marine_selected_default_20260630/logs/sample5_minco_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py --profile-preset candidate -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_marine_samples3_5_20260625/data/cami_marine_sample5_reads.fq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /mnt/new3T/minco_marine_selected_default_20260630/sample5_work -p16 -o /mnt/new3T/minco_marine_selected_default_20260630/minco_sample5_current_default.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_marine_setup_truth_profile_rescore.py --samples 3,4,5 --profile-set-label selected_default_same_namespace --minco-method minco_selected_default_setup_truth --sylph-method sylph_cached_marine_setup_truth --minco-sample 3=/mnt/new3T/minco_marine_selected_default_20260630/minco_sample3_current_default.tsv --minco-sample 4=/mnt/new3T/minco_marine_selected_default_20260630/minco_sample4_current_default.tsv --minco-sample 5=/mnt/new3T/minco_marine_selected_default_20260630/minco_sample5_current_default.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_remaining_holdout_route_options.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_universal_strategy_next_evidence_routes.py tests/test_cami3_extension_after_recovery.py tests/test_minco_profile_calibrated_auto_exact.py -q

# 2026-06-30 completed-negative route failure-mode audit and raw-side diagnostic sweep.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_negative_route_failure_modes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_raw_candidate_retention_negative_routes.py
python3 -B -m py_compile research/experiments/2026-06-27_universal_strategy_decision/audit_negative_route_failure_modes.py research/experiments/2026-06-27_universal_strategy_decision/sweep_raw_candidate_retention_negative_routes.py

# 2026-06-30 marine setup metadata truth-upgrade audit.
mkdir -p /tmp/cami2_marine_samples3_5_20260625/setup_metadata
curl -fL https://frl.publisso.de/data/frl%3A6425521/marine/short_read/marmgCAMI2_setup.tar.gz | tar -xzf - -C /tmp/cami2_marine_samples3_5_20260625/setup_metadata simulation_short_read/genome_to_id.tsv simulation_short_read/metadata.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_setup_metadata_truth_upgrade.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_source_mapping_local_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_remaining_holdout_route_options.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_setup_metadata_truth_upgrade_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_source_mapping_local_inventory_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv

# 2026-06-30 route-state correction after setup+assembly truth upgrade.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_remaining_holdout_route_options.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/remaining_holdout_route_options_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-29 marine local source-mapping inventory after readmap feasibility.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_source_mapping_local_inventory.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_source_mapping_local_inventory_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_source_mapping_local_inventory.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv

# 2026-06-29 remaining local holdout route audit after marine inventory failure.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_remaining_holdout_route_options.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/remaining_holdout_route_options_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/remaining_holdout_route_options.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-29 marine truth-upgrade candidate-rule audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_truth_upgrade_candidate_rules.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_truth_upgrade_candidate_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_truth_upgrade_candidate_rules.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_truth_upgrade_remaining_blockers.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-29 next-evidence route audit for the active universal-strategy goal.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_universal_strategy_next_evidence_routes.py -q
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_NEXT_EVIDENCE_ROUTES.md

# 2026-06-29 active-goal completion audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv
sed -n '1,160p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_GOAL_COMPLETION_AUDIT.md

# 2026-06-29 refreshed headroom-aware next-route audit after CAMI3 recovery.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_cleanup_candidates.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_universal_strategy_next_evidence_routes.py tests/test_cami3_extension_after_recovery.py tests/test_minco_profile_calibrated_auto_exact.py -q

# 2026-06-29 HMP airskin sample26 recovery, scoring, and post-run audit.
# Executed with explicit approval as a one-time destructive scratch cleanup:
# rm -rf /tmp/cami_marine_sample1_reads.fq.gz /tmp/cami_marine_sample2_reads.fq.gz /tmp/minco_context_em_external_20260624/cami_marine_sample0_reads.fq.gz
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
# At execution time this generated script targeted sample26. It has since been
# regenerated by plan_next_release_grade_actions.py to target sample2.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
df -h /tmp /mnt/new3T
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample26_stream_bam_to_fastq_20260627.time.log
tail -40 /tmp/minco_current_code_hmp_airskin26_20260627/minco_sample26_current_default.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample26/sketch.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample26/profile.time.log
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample26_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_28_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample26_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_28_scores.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_universal_strategy_next_evidence_routes.py tests/test_cami3_extension_after_recovery.py tests/test_minco_profile_calibrated_auto_exact.py -q

# 2026-06-29 HMP airskin sample2 recovery, scoring, and post-run audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
# At execution time this generated script targeted sample2. It has since been
# regenerated by plan_next_release_grade_actions.py to target sample8.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
df -h /tmp /mnt/new3T
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample2_stream_bam_to_fastq_20260627.time.log
tail -40 /tmp/minco_current_code_hmp_airskin2_20260627/minco_sample2_current_default.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample2/sketch.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample2/profile.time.log
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin2_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample2_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_28_summary.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py

# 2026-06-29 HMP airskin sample8 recovery, scoring, and post-run audit.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
# Executed with explicit approval as a one-time destructive scratch cleanup:
# rm -f /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample0.nonzero.fastq.gz /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample5.nonzero.fastq.gz
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
# At execution time this generated script targeted sample8. It has since been
# regenerated by plan_next_release_grade_actions.py to target sample0.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
df -h /tmp /mnt/new3T
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample8_stream_bam_to_fastq_20260627.time.log
tail -40 /tmp/minco_current_code_hmp_airskin8_20260627/minco_sample8_current_default.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample8/sketch.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample8/profile.time.log
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin8_r232_source_abundance_scores.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample8_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_28_summary.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_abundance_release_blocker.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_holdout_bundle.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_strategy_decision_summary.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/build_release_readiness.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_NEXT_EVIDENCE_ROUTES.md

# 2026-06-29 FASTQ cleanup, HMP omitted samples27/12 recovery, and route-state fix.
# Executed with explicit user request to remove old generated FASTQ caches.
rm -f /tmp/cami2_hmp_airskin_20260625/reads/airskinurogenital_sample28.nonzero.fastq.gz /tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample22.nonzero.fastq.gz
df -h /tmp /mnt/new3T
find /mnt/new3T/minco_release_holdouts_20260628 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) -print
rm -f /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample1.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample3.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample4.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample7.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample9.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample10.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample13.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample14.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample15.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample16.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample17.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample18.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample19.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample20.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample21.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample23.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample24.nonzero.fastq.gz /mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample25.nonzero.fastq.gz
df -h /mnt/new3T
find /mnt/new3T/minco_release_holdouts_20260628 -type f \( -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.fastq' -o -name '*.fq' \) -print
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
# At execution time this generated script targeted sample27. It was later
# regenerated by plan_next_release_grade_actions.py to target sample12.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample27_stream_bam_to_fastq_20260627.time.log
tail -40 /tmp/minco_current_code_hmp_airskin27_20260627/minco_sample27_current_default.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample27/sketch.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample27/profile.time.log
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin27_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample27_0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_26_28_summary.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
# At execution time this generated script targeted sample12.
bash research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin_next_release_sample_commands.sh
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/logs/sample12_stream_bam_to_fastq_20260627.time.log
tail -40 /tmp/minco_current_code_hmp_airskin12_20260627/minco_sample12_current_default.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample12/sketch.time.log
tail -40 /tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample12/profile.time.log
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin12_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_current_refresh_r232_source_abundance_sample12_0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_summary.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-29 HMP airskin sample8 extra-call diagnostic.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_sample8_extra_calls.py
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin8_minco_call_status_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin8_minco_filter_diagnostic.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_airskin8_minco_missing_truth.tsv

# 2026-06-29 opt-in cross-panel call-filter audit with HMP omitted samples2/8/26.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/sweep_cross_panel_call_filters.py --include-hmp-omitted
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_call_filter_with_hmp_omitted_safety_audit.tsv
sed -n '1,40p' research/experiments/2026-06-27_universal_strategy_decision/results/cross_panel_call_filter_with_hmp_omitted_overall.tsv

# 2026-06-29 opt-in adaptive call-filter switch audit with HMP omitted samples2/8/26.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_adaptive_call_filter_switch.py --include-hmp-omitted
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_switch_with_hmp_omitted_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/adaptive_call_filter_switch_with_hmp_omitted_lopo.tsv

# 2026-06-29 next-action planner correction after HMP omitted samples2/8/26.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_hmp_airskin_source_truth_candidates.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_tmp_headroom_for_next_action.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/plan_next_release_grade_actions.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_refined_allocator_independent_holdout_inventory.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/next_release_grade_actions.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/results/tmp_headroom_audit.tsv
sed -n '1,80p' research/experiments/2026-06-27_universal_strategy_decision/UNIVERSAL_STRATEGY_NEXT_EVIDENCE_ROUTES.md

# 2026-06-29 marine source-readmap feasibility after HMP omitted route completion.
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_marine_source_readmap_feasibility.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_source_readmap_feasibility_audit.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/marine_source_readmap_feasibility_summary.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_next_evidence_routes.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_holdout_gap_action_plan.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_next_evidence_routes_audit.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_release_gate.tsv
sed -n '1,120p' research/experiments/2026-06-27_universal_strategy_decision/results/universal_strategy_goal_completion_audit.tsv

# 2026-06-30 HMP gastrooral sample1 independent same-namespace holdout for the guarded raw-retention candidate.
mkdir -p /tmp/cami2_hmp_gastrooral_sample1_20260630/reads /tmp/cami2_hmp_gastrooral_sample1_20260630/logs /tmp/cami2_hmp_gastrooral_sample1_20260630/run /tmp/minco_current_code_hmp_gastrooral_sample1_20260630 /tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1 /tmp/cami2_hmp_gastrooral_sample1_20260630/run/sylph_sample1 /tmp/minco_hmp_gastrooral_sample1_r232_profile_inputs
awk '$2+0>0{print "2017.12.04_18.45.54_sample_1/bam/" $1 ".bam"}' /tmp/cami2_hmp_pilot_20260625/truth/abundance1.tsv > /tmp/cami2_hmp_gastrooral_sample1_20260630/nonzero_bam_paths_sample1.txt
/usr/bin/time -v -o /tmp/cami2_hmp_gastrooral_sample1_20260630/logs/sample1_stream_bam_to_fastq_20260630.time.log bash -c 'set -euo pipefail; curl -fL https://frl.publisso.de/data/frl:6425518/gastrooral/sample_1.tar.gz | tar -xzf - -T /tmp/cami2_hmp_gastrooral_sample1_20260630/nonzero_bam_paths_sample1.txt --to-command="samtools fastq -" | gzip -1 > /tmp/cami2_hmp_gastrooral_sample1_20260630/reads/gastrooral_sample1.nonzero.fastq.gz'
gzip -t /tmp/cami2_hmp_gastrooral_sample1_20260630/reads/gastrooral_sample1.nonzero.fastq.gz
python3 -B -m py_compile research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_current_default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/cami2_hmp_gastrooral_sample1_20260630/reads/gastrooral_sample1.nonzero.fastq.gz --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --workdir /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work -p16 -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_current_default.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_current_surface_accession_taxmap.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --train-pool train12 --strategy universal-auto-exact --candidate-rescue-switch emitted-ani90-xny100-br01-af70 --candidate-surface-switch accession-ani90-xny100-br01-af70 --candidate-abundance-policy normalized-depth-alpha2 -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_current_surface_accession_taxmap.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_guarded_candidate.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --train-pool train12 --strategy universal-auto-exact --candidate-rescue-switch emitted-ani90-xny100-br01-af70 --candidate-surface-switch accession-current-or-ani93-xny650-br20 --candidate-surface-max-called-species 250 --candidate-abundance-policy normalized-depth-alpha2 -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_guarded_candidate.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_guarded_candidate_zero.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_calibrated.py --unique-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --candidate-surface-taxmap /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria --train-pool train12 --strategy universal-auto-exact --candidate-rescue-switch emitted-ani90-xny100-br01-af70 --candidate-surface-switch accession-current-or-ani93-xny650-br20 --candidate-surface-max-called-species 250 --candidate-abundance-policy zero -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_guarded_candidate_zero.tsv
/usr/bin/time -v -o /tmp/cami2_hmp_gastrooral_sample1_20260630/run/sylph_sample1/sketch.time.log /home/ubuntu/yihuiguang/bin/sylph sketch -t 16 -r /tmp/cami2_hmp_gastrooral_sample1_20260630/reads/gastrooral_sample1.nonzero.fastq.gz -d /tmp/cami2_hmp_gastrooral_sample1_20260630/run/sylph_sample1
sed '$a /tmp/cami2_hmp_gastrooral_sample1_20260630/run/sylph_sample1/gastrooral_sample1.nonzero.fastq.gz.sylsp' /mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list > /tmp/minco_hmp_gastrooral_sample1_r232_profile_inputs/sample1.profile_inputs.list
/usr/bin/time -v -o /tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.time.log /home/ubuntu/yihuiguang/bin/sylph profile -t 16 -l /tmp/minco_hmp_gastrooral_sample1_r232_profile_inputs/sample1.profile_inputs.list -o /tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.tsv
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --sample 1,gastrooral1,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_current_default.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv,/tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.tsv --minco-method minco_current_raw_default_gastrooral_sample1_source_abundance --output-prefix hmp_gastrooral_sample1_default_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --sample 1,gastrooral1,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_current_surface_accession_taxmap.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv,/tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.tsv --minco-method minco_current_surface_accession_taxmap_gastrooral_sample1_source_abundance --output-prefix hmp_gastrooral_sample1_current_surface_accession_taxmap_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --sample 1,gastrooral1,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_guarded_candidate.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv,/tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.tsv --minco-method minco_guarded_candidate_gastrooral_sample1_source_abundance --output-prefix hmp_gastrooral_sample1_guarded_candidate_r232_source_abundance
python3 -B research/experiments/2026-06-27_universal_strategy_decision/score_hmp_gastrooral_gtdb_source_abundance.py --sample 1,gastrooral1,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_guarded_candidate_zero.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv,/tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv,/tmp/cami2_hmp_gastrooral_sample1_20260630/run_sylph_r232/sylph_sample1/profile.chunked.tsv --minco-method minco_guarded_candidate_zero_gastrooral_sample1_source_abundance --output-prefix hmp_gastrooral_sample1_guarded_candidate_zero_r232_source_abundance
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_sample1_default_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_sample1_current_surface_accession_taxmap_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_sample1_guarded_candidate_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_sample1_guarded_candidate_zero_r232_source_abundance_summary.tsv
sed -n '1,20p' research/experiments/2026-06-27_universal_strategy_decision/results/hmp_gastrooral_sample1_holdout_comparison.tsv

# 2026-06-30 local S1000 reference candidate-surface sidecar installation and discovery check.
cp /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/candidate_surface_taxmap.tsv
sha256sum /tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/candidate_surface_taxmap.tsv
/usr/bin/time -v -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_default_sidecar_discovery.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile_default.py -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --unique-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_unique.unfiltered.tsv --split-table /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/sample1_work/minco.best_diff_split.exact.unfiltered.tsv --taxmap /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv --model-cache /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib --scope bacteria -o /tmp/minco_current_code_hmp_gastrooral_sample1_20260630/minco_sample1_default_sidecar_discovery.tsv

# 2026-06-30 default launcher usability hardening.
chmod +x scripts/minco_profile
python3 -B research/experiments/2026-06-27_universal_strategy_decision/write_current_default_strategy_manifest.py
python3 -B -m py_compile scripts/minco_profile scripts/minco_profile_default.py scripts/minco_profile_calibrated.py
env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile --help
PYTHONDONTWRITEBYTECODE=1 python3 -B -m pytest tests/test_universal_strategy_next_evidence_routes.py tests/test_cami3_extension_after_recovery.py tests/test_minco_profile_calibrated_auto_exact.py -q
bash tests/smoke.sh
python3 -B research/experiments/2026-06-27_universal_strategy_decision/validate_universal_strategy_release_gate.py
python3 -B research/experiments/2026-06-27_universal_strategy_decision/audit_universal_strategy_goal_completion.py
git --git-dir=.minco_git --work-tree=. diff --check

# 2026-06-30 packaged default preflight and no-extra-flags smoke.
install -m 644 /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/species_taxmap.tsv
install -m 644 /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/minco_profile_rf_hgb.train12.unfiltered.joblib
sha256sum /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/species_taxmap.tsv /tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/minco_profile_rf_hgb.train12.unfiltered.joblib /tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib
# Expected failure before minco_profile_defaults.tsv: model cache metadata has scope=bacteria but parser default is scope=all.
if [ ! -e /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/minco_profile_defaults.tsv ]; then
  ! env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile --check-ref -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz
  ! /usr/bin/time -v -o /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz -p16 -o /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.tsv
fi
printf 'option\tvalue\nscope\tbacteria\n' > /tmp/minco_profile_defaults.gtdb_s1000.tsv
install -m 644 /tmp/minco_profile_defaults.gtdb_s1000.tsv /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/minco_profile_defaults.tsv
env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile --check-ref -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz
/usr/bin/time -v -o /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.time.log env PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/minco_profile -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno --reads /tmp/toymouse_sample0_50000_reads.fq.gz -p16 -o /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.tsv
grep -E 'Elapsed \(wall clock\)|Maximum resident set size|Exit status|User time|System time' /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.time.log
wc -l /tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.tsv
python3 -c 'import csv; p="/tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.tsv"; rows=list(csv.DictReader(open(p, newline=""), delimiter="\t")); print("rows", len(rows)); print("called", sum(r.get("calibrated_call","")=="True" for r in rows)); print("candidate_surface_added", sum(r.get("candidate_surface_added","")=="True" for r in rows)); print("calibrated_abundance_sum", sum(float(r.get("calibrated_abundance") or 0) for r in rows)); print("scope_values", sorted(set(r.get("scope","") for r in rows))); print("strategy_values", sorted(set(r.get("profile_strategy","") for r in rows)))'
