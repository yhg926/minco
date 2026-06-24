# Artifacts

## Code

- MinCO source: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/src/command_ani.c`
- ANI options header: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/src/command_ani.h`
- ANI wrapper: `/home/ubuntu/yihuiguang/tools/KSSD3mini/minco_core/src/command_ani_wrapper.c`
- Run script: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/commands.sh`
- Scorer: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/score_integrated_dual_evidence.py`

## External Inputs

- Full S2000 dedup refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Physical markerdb size map for validation: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv`
- Toy Mouse reads:
  - `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_toymouse_20260621/sample_1/2017.12.29_11.37.26_sample_1/reads/anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_toymouse_20260621/sample_2/2017.12.29_11.37.26_sample_2/reads/anonymous_reads.fq.gz`
- CAMI3 ToyGut reads:
  - `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/anonymous_reads.fq.gz`

## Generated Large Outputs

Stored in `/tmp/minco_dual_s2000_20260624`:

- `mouse_s0_dual.tsv`, 103 MB
- `mouse_s1_dual.tsv`, 133 MB
- `mouse_s2_dual.tsv`, 134 MB
- `cami3_s0_dual.tsv`, 86 MB
- `cami3_s1_dual.tsv`, 94 MB
- `cami3_s2_dual.tsv`, 86 MB
- `*_dual.time.log`, `*_dual.stderr.log`, `*_dual.stdout.log`

These are temporary external artifacts under `/tmp`; rerun `commands.sh` if the directory is removed.

## Result Tables

- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_f1_abundance_comparison.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_addback_summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_addback_sample_metrics.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_em_proxy_top50.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_mouse_mean_metrics.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/integrated_cami3_summary.tsv`
- `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_integrated_dual_evidence_benchmark/results/marker_size_validation.tsv`

## Provenance

- Code status: `provenance/code_status.txt`
- Diff stat: `provenance/code_diffstat.txt`
