# Artifacts

Date: 2026-06-17

## Input

- Genome list: `/tmp/kssdmini_1000_speed/genomes.list`
- Source list: `/mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list`
- Genome count: 1000

## Code

- KSSDmini repo: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Copied KSSD3A stage0 source: `/home/ubuntu/yihuiguang/tools/KSSD3mini/kssd3a_stage0`
- Stage0 binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/kssd3mini_stage0`
- Stage1 binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/kssd3mini_stage1`
- Stage2 binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/kssd3mini_stage2`
- Selected Stage3 binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/kssd3mini_stage3_native`
- Native KSSD3A baseline binary: `/home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a`

## Large Temporary Outputs

- Stage0 sketch directory: `/tmp/kssd3a_based_kssdmini_stages/stage0_kssd3a_f8`
- Stage1 sketch directory: `/tmp/kssd3a_based_kssdmini_stages/stage1_f8_hash_bottomk_v2`
- Stage2 exact sketch directory: `/tmp/kssd3a_based_kssdmini_stages/stage2_no_filter_hash_bottomk_1000`
- Selected Stage3 sketch directory: `/tmp/kssd3a_based_kssdmini_stages/stage3_native_default_m2_1000`
- Native KSSD3A baseline sketch directory: `/tmp/kssdmini_1000_speed/baselines/kssd3a_sketch`
- Current KSSDmini compact sketch: `/tmp/kssdmini_1000_speed/compact_bmi2_vec_radix.kmini`

## Tracked Result Files

- `summary.tsv`
- `stage0_kssd3a_f8_sketch.time.log`
- `stage0_kssd3a_f8.files.tsv`
- `stage0_kssd3a_f8.du.tsv`
- `stage0_vs_native.cmp.tsv`
- `stage1_f8_hash_bottomk_sketch.time.log`
- `stage2_exact_postcollect_sketch.time.log`
- `stage3_native_stream_m2_sketch.time.log`
- `stage3_native_stream_m2.files.tsv`
- `stage3_exact_checks.tsv`
- `stage_speed_tuning.tsv`

## Stage0 Byte Check

The stage0 raw payload files `comblco`, `comblco.index`, and
`lcofiles.infilemeta` are byte-identical to the native KSSD3A baseline. The
`lcofiles.stat` file differs.

## Stage3 Byte Check

The selected Stage3 `comblco` payload is byte-identical to the Stage2 exact
post-collection payload and to the current KSSDmini compact sketch after parsing
the `.kmini` entries.
