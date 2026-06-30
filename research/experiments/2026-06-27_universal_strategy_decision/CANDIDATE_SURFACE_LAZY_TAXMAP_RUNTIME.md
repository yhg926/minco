# Candidate Surface Lazy Taxmap Runtime

Date: 2026-06-30

## Question

Can the current candidate preset reduce Python wrapper overhead without changing MinCO calls, abundances, or candidate-surface output?

## Change

`scripts/minco_profile_calibrated.py` now applies the candidate-surface numeric gate before accession-level taxmap mapping. When an explicit `--candidate-surface-taxmap` is supplied, the file is parsed only if at least one raw unique/split row can pass the candidate-surface numeric thresholds.

This does not change the raw-read MinCO scan count or the exact-split rerun policy.

## Benchmark

Working directory: `/home/ubuntu/yihuiguang/tools/KSSD3mini`

Baseline code: `40d6a6c63a35a7ddfcb51849720a999016a8a705`

Patched code: dirty worktree with lazy candidate-surface numeric prefilter.

Input tables:

- Unique: `/tmp/minco_current_code_hmp_gastrooral_20260627/sample6_work/minco.best_diff_unique.unfiltered.tsv`
- Split: `/tmp/minco_current_code_hmp_gastrooral_20260627/sample6_work/minco.best_diff_split.unfiltered.tsv`
- Taxmap: `/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv`
- Candidate-surface taxmap: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno/candidate_surface_taxmap.tsv`

Results are in `results/candidate_surface_lazy_taxmap_runtime.tsv`.

## Result

On this cached real table-mode profile, the patched path preserved output:

- Rows and columns: identical, `29438 x 177`
- Called rows: identical
- Candidate-surface rows: identical
- Maximum numeric difference: `2.220446049250313e-16`
- Text/boolean columns: identical

Runtime improved from `16.92s` to `14.87s`; peak RSS improved from `2.9340 GiB` to `2.6421 GiB`.

## Caveat

This is a wrapper overhead optimization. It does not solve the larger first-run raw-read speed gap caused by exact-split reruns. The existing runtime boundary remains: unconditional exact sidecar is not promoted because it speeds exact-needed samples but slows samples where exact is skipped.
