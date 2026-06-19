# Artifacts

Record paths, sizes, and whether files are temporary or should be preserved.

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| input | `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/{5m,10m,20m,50m,75m,100m,250m}/2009K-1545__SRR19787631__*.fastq.gz` | ONT downsample yield series for one Salmonella isolate. | yes |
| binary | `bin/minco_stage3_native` | Normal 10,000-context minco density sketch. | yes |
| binary | `bin/minco_exact20m` | Exact-debug binary for 5m-20m. | maybe |
| binary | `bin/minco_exact100m` | Exact-debug binary for 50m-250m. | maybe |
| output | `research/experiments/2026-06-18_minco_ont_yield_density_estimator_accuracy/summary.tsv` | Final per-yield accuracy table. | yes |
| output | `research/experiments/2026-06-18_minco_ont_yield_density_estimator_accuracy/per_yield.tsv` | Same final table, intended for plotting. | yes |
| output | `/tmp/minco_ont_density_*_std_0618` | Normal 10k sketches and ctxmeta sidecars. | maybe |
| output | `/tmp/minco_ont_density_*_exact20m_0618` | Exact20m validation sketches; 50m+ are capped and superseded by exact100m. | no |
| output | `/tmp/minco_ont_density_*_exact100m_0618` | Exact100m validation sketches for 50m-250m. | no |
