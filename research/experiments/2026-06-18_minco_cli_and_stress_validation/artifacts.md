# Artifacts

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| source | `/home/ubuntu/yihuiguang/tools/minco` | minco pure C repo under test. | yes |
| test | `tests/smoke.sh` | Fast public-command smoke test. | yes |
| test | `tests/full_cli.sh` | Expanded CLI/function coverage test. | yes |
| input | `/mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list` | Source genome list; first 1000 used. | yes |
| input | `/tmp/minco_validation_20260618/genomes1000/genomes1000.list` | Exact 1000-genome path list used. | maybe |
| output | `/tmp/minco_validation_20260618/genomes1000/sketch` | 1000-genome sketch directory, 192M. | maybe |
| output | `/tmp/minco_validation_20260618/genomes1000/ani_triangle.tsv` | 1000-genome ANI triangle output, 1000 rows. | maybe |
| output | `/tmp/minco_validation_20260618/genomes1000/matrix_triangle.tsv` | 1000-genome matrix triangle output, 1000 rows. | maybe |
| input | `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/250m` | ONT 250m FASTQ directory; first 4 sorted files used. | yes |
| input | `/tmp/minco_validation_20260618/fastq250m/reads250m_4.list` | Exact FASTQ path list used. | maybe |
| output | `/tmp/minco_validation_20260618/fastq250m/reads_sketch` | Four-read-set conflict readsQC sketch, 340K. | maybe |
| output | `/tmp/minco_validation_20260618/fastq250m/asm_sketch` | Four matching assembly sketches, 808K. | maybe |
| output | `/tmp/minco_validation_20260618/fastq250m/qraw_ani.tsv` | Four read sketches by four assemblies, 16 comparisons. | maybe |
| log | `/tmp/minco_validation_20260618/**/*.time.log` | `/usr/bin/time -v` runtime and memory logs. | maybe |
