# Artifacts

Record paths, sizes, and whether files are temporary or should be preserved.

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| input | `/tmp/kssd3a_based_kssdmini_stages/genomes_50.list` | Source GTDB genome path list; first 50 entries used. | yes |
| input | `/tmp/kssdmini_density_test50.list` | Exact 50-genome list used for this run. | maybe |
| binary | `bin/kssd3mini_stage3_native` | Normal 10,000-context KSSDmini stage3 native binary. | yes |
| binary | `bin/kssd3mini_exact20m` | Debug exact-count binary compiled with `KSSD3MINI_SKETCH_SIZE=20000000`. | maybe |
| output | `research/experiments/2026-06-18_kssdmini_density_estimator_predicts_unique_contexts/summary.tsv` | Correlation and error metrics. | yes |
| output | `research/experiments/2026-06-18_kssdmini_density_estimator_predicts_unique_contexts/per_sample.tsv` | Per-sample genome size, exact unique contexts, estimates, and relative errors. | yes |
| output | `/tmp/kssdmini_density50_std_0618` | Normal 10k sketch and ctxmeta sidecar. | maybe |
| output | `/tmp/kssdmini_density50_exact_0618` | Exact-count large sketch; about 1.3 GB. | no |
