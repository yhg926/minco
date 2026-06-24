# Artifacts

Repo-tracked:

- Experiment note: `research/experiments/2026-06-20_minco_readwise_assignment_fix/NOTE.md`
- Summary: `research/experiments/2026-06-20_minco_readwise_assignment_fix/summary.tsv`
- Commands: `research/experiments/2026-06-20_minco_readwise_assignment_fix/commands.sh`
- Threshold scorer: `research/experiments/2026-06-20_minco_readwise_assignment_fix/scripts/grid_profile_thresholds.py`

Input data:

- CAMI reads: `/tmp/cami_marine_sample0_reads.fq.gz`
- Gold profile: `/tmp/gs_marine_short.profile`
- Taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Sylph profile: `/tmp/sylph_marine_sample0/profile.tsv`
- S10000 reference sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`

Generated outputs:

- Output directory: `/tmp/minco_readwise_assign_20260620`
- Best-diff unfiltered TSV: `/tmp/minco_readwise_assign_20260620/s10000_bestdiff_unfiltered.tsv`
- Best-diff-split unfiltered TSV: `/tmp/minco_readwise_assign_20260620/s10000_bestdiffsplit_unfiltered.tsv`
- Best-diff-unique naive unfiltered TSV: `/tmp/minco_readwise_assign_20260620/s10000_bestdiffunique_unfiltered.tsv`
- C ZIP AAF unfiltered TSV: `/tmp/minco_readwise_assign_20260620/s10000_unique_zipaaf_unfiltered_c.tsv`
- Final direct CLI TSV: `/tmp/minco_readwise_assign_20260620/s10000_unique_zipaaf_f0.05_n0.94_t10.tsv`
- Final direct CLI threshold table: `/tmp/minco_readwise_assign_20260620/s10000_unique_zipaaf_c.threshold_best.tsv`
- Baseline grid: `/tmp/minco_readwise_assign_20260620/s10000_all_baseline.threshold_best.tsv`
- Best-diff grid: `/tmp/minco_readwise_assign_20260620/s10000_bestdiff.threshold_best.tsv`
- Best-diff-split grid: `/tmp/minco_readwise_assign_20260620/s10000_bestdiffsplit.threshold_best.tsv`
- Best-diff-unique grid: `/tmp/minco_readwise_assign_20260620/s10000_bestdiffunique.threshold_best.tsv`
