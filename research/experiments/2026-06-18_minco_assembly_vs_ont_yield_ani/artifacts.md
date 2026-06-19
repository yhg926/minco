# Artifacts

| Type | Path | Notes |
| --- | --- | --- |
| binary | `/home/ubuntu/yihuiguang/tools/minco/bin/minco_stage3_native` | minco stage3 native build. |
| reference assembly | `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz` | ATB/Shovill-SPAdes assembly for matched Illumina BioSample `SAMN29253066`. |
| ONT read inputs | `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/{5m,10m,20m,50m,75m,100m,250m}/2009K-1545__SRR19787631__*.fastq.gz` | Downsampled ONT FASTQs for run `SRR19787631`. |
| assembly sketch | `/tmp/minco_asm_2009K1545_0618` | Temporary minco assembly sketch. |
| no-QC read sketch | `/tmp/minco_ont_yield_2009K1545_noqc_0618` | Temporary minco read sketches without read QC. |
| readsQC read sketch | `/tmp/minco_ont_yield_2009K1545_readsqc_0618` | Temporary minco read sketches with `--readsQC`. |
| no-conflict no-QC read sketch | `/tmp/minco_ont_yield_2009K1545_noconflict_noqc_0618` | Temporary minco read sketches without `--conflict` or read QC. |
| no-conflict readsQC read sketch | `/tmp/minco_ont_yield_2009K1545_noconflict_readsqc_0618` | Temporary minco read sketches with `--readsQC` and without `--conflict`. |
| raw no-QC ANI | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/noqc_ani_detail.tsv` | Copied detail ANI output. |
| raw readsQC ANI | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/readsQC_ani_detail.tsv` | Copied detail ANI output. |
| raw no-conflict no-QC ANI | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/noconflict_noqc_ani_detail.tsv` | Copied detail ANI output. |
| raw no-conflict readsQC ANI | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/noconflict_readsQC_ani_detail.tsv` | Copied detail ANI output. |
| patched readsQC ANI | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/readsQC_kssdqcsample_counts_ani_detail.tsv` | KSSD3A-sampled QC range plus count-preserving bottom-k result. |
| summary | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/summary.tsv` | Compact yield-by-mode table. |
| no-conflict summary | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/summary_noconflict.tsv` | Compact table for no-conflict modes. |
| all-mode summary | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/summary_all_modes.tsv` | Combined conflict-kept and no-conflict table. |
| patched readsQC summary | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/summary_readsqc_kssdqcsample_counts.tsv` | Compact patched readsQC table. |
| patched readsQC debug log | `research/experiments/2026-06-18_minco_assembly_vs_ont_yield_ani/readsQC_kssdqcsample_counts_debug.log` | Saved range-inference log. |
