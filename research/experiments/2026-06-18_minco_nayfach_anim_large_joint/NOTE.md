# Minco Large Nayfach ANIm Calibration

Date: 2026-06-18
Author/agent: Codex
Project: minco
Code commit: NA; this workspace is not a git repository
Tool version: `bin/minco 0.1`

## Question

Does a larger 60k-120k Nayfach ANIm calibration run improve Minco assembled-genome ANI, and should the MoE and HGB layers be retrained together?

## Hypothesis

A 120k stratified pair set should give a more stable HGB calibration than the 12k pilot. Re-optimizing the raw CtxMoE model on the same large set may improve the HGB input distribution, but it should only be installed if the held-out test split improves.

## Dataset

- Label table: `/mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv`
- FASTA root: `/mnt/new3T/skani_data/Nayfach_data/fna`
- Source rows: 1,075,493 pairwise ANIm rows.
- Sample: 120,000 pairs, 20,000 per ANI bin.
- Bins: `<0.90`, `0.90-0.94`, `0.94-0.95`, `0.95-0.97`, `0.97-0.99`, `>=0.99`.
- Successful feature rows: 119,999 / 120,000.
- Split: stratified 70/30 train/test, giving 83,999 train and 36,000 held-out test pairs.

One feature extraction failed:

```text
P000080777  094_095  exception  minco produced header but no comparison row
```

## Methods

Key parameters:

```text
Sketch size: 10000
Feature extraction: minco ani -S 10000 -p1 -f0 -n0 -t0 --raw-output -s3
Extractor concurrency: 16 direct-pair jobs
MoE optimizer: maxiter 35, popsize 8, polish-iter 160
HGB candidates: current raw MoE HGB, optimized-MoE HGB 160 trees, optimized-MoE HGB 240 trees
Install rule: install only the held-out winner among source-compatible candidates
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

Feature extraction:

```text
120000 direct pair attempts
119999 ok rows
36:06.08 wall time
1,178,304 KB peak RSS
56.22 final pairs/s
```

Held-out raw MoE comparison:

```text
current raw CtxMoE:     MAE 0.0029997561, RMSE 0.0045030119, cutoff errors 1543
large optimized MoE:    MAE 0.0030169095, RMSE 0.0044792826, cutoff errors 1535
```

The large optimized MoE improved RMSE and 0.95 cutoff errors slightly, but it worsened MAE, so it was not installed.

Held-out HGB comparison:

```text
large HGB on current raw MoE:  MAE 0.0025661341, RMSE 0.0038749039, cutoff errors 987
joint optimized-MoE HGB:       MAE 0.0025665830, RMSE 0.0038741155, cutoff errors 999
```

The HGB trained on the current raw MoE won, so only `model_refaf_hgb.h` was replaced.

Compiled binary comparison on the same 36,000 held-out pairs:

```text
previous compiled Best:  MAE 0.0025804308, RMSE 0.0040106008, cutoff errors 1032
installed large HGB:     MAE 0.0025661569, RMSE 0.0038751974, cutoff errors 987
```

This is a small but real improvement: 0.553% lower MAE, 3.376% lower RMSE, and 45 fewer 0.95 cutoff errors.

## Validation

- `make test` passed after installing the large HGB.
- Current compiled Best was evaluated on all 36,000 held-out pairs before replacement.
- Installed compiled Best was evaluated on all 36,000 held-out pairs after replacement.
- All installed validation rows used `BestDist`.
- C-vs-Python exported large HGB agreement was MAE `5.4009e-7` ANI.

## Code Change

Installed:

```text
minco_core/src/model_refaf_hgb.h
```

Preserved previous header:

```text
research/experiments/2026-06-18_minco_nayfach_anim_large_joint/model_refaf_hgb.before_large_hgb.h
```

Not installed:

```text
large_10k_perbin20000/moe_opt_long/optimized_moe_model_ani_arrays.hfrag
large_10k_perbin20000/joint_hgb/model_refaf_hgb.joint_generated.h
```

The raw MoE constants in `model_ani.h` were intentionally left unchanged.

## Important Artifacts

See `artifacts.md`.

## Conclusion

The large-scale run supports installing a larger Nayfach-trained HGB layer, but not replacing the raw MoE constants. Joint MoE/HGB retraining was tested; the optimized MoE branch did not improve held-out MAE after HGB, so the conservative install is HGB-only.

## Paper-Relevant Claim

Tentative: a larger stratified Nayfach ANIm calibration set marginally improves Minco default assembled-genome ANI and 0.95 threshold behavior, but the current raw MoE remains competitive enough that the large retraining gain comes from the HGB layer.

## Caveats

- The gain over the previous 12k-installed default is modest on this 36k held-out split.
- The calibration is still for assembled genomes with default coden context/object pattern and sketch size 10,000.
- Read/`--qraw` modes are not covered.
- One sampled pair failed direct Minco feature extraction and was dropped from training.

## Next Experiment

Run a species/genus-disjoint holdout, not just a random stratified pair holdout from the Nayfach table, to test whether the large HGB improves extrapolation outside the sampled taxonomic mix.
