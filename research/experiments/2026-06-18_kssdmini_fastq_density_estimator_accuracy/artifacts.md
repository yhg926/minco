# Artifacts

Record paths, sizes, and whether files are temporary or should be preserved.

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| input | `/mnt/new3T/kssd3test/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq` | Real FASTQ used for subsets and full-sketch check; 3,444,570 reads, 846 MB. | yes |
| input | `/tmp/kssdmini_fastq_R1_10k.fq` | First 10,000 reads. | maybe |
| input | `/tmp/kssdmini_fastq_R1_50k.fq` | First 50,000 reads. | maybe |
| input | `/tmp/kssdmini_fastq_R1_100k.fq` | First 100,000 reads. | maybe |
| input | `/tmp/kssdmini_fastq_R1_200k.fq` | First 200,000 reads. | maybe |
| binary | `bin/kssd3mini_stage3_native` | Normal 10,000-context KSSDmini stage3 native binary. | yes |
| binary | `bin/kssd3mini_exact20m` | Debug exact-count binary compiled with `KSSD3MINI_SKETCH_SIZE=20000000`. | maybe |
| output | `research/experiments/2026-06-18_kssdmini_fastq_density_estimator_accuracy/summary.tsv` | FASTQ density accuracy table by read subset. | yes |
| output | `/tmp/kssdmini_fastq_density_*_std_0618` | Normal 10k FASTQ subset sketches. | maybe |
| output | `/tmp/kssdmini_fastq_density_*_exact_0618` | Exact large FASTQ subset sketches. | no |
| output | `/tmp/kssdmini_fastq_full_R1_conflict_0618` | Full R1 normal sketch; exact count not computed. | maybe |
