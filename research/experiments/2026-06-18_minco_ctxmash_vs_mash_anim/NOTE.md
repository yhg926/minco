# Minco Context-Mash Versus Mash ANI On Vibrio ANIm Pairs

Date: 2026-06-18
Author/agent: Codex
Project: Minco/KSSD3mini
Code commit: no git repository present in working directory
KSSD3A binary/tool version: minco 0.1; Mash 2.3

## Question

Does Minco context-Mash ANI have the same accuracy as ordinary Mash ANI when both use sketch size 10,000?

## Hypothesis

For equal-size Minco sketches, Minco `MashD` and Minco `AafD` should be identical because the Mash correction reduces to the same shared-context fraction. Ordinary Mash may differ because it sketches k-mers instead of context-object hashes.

## Dataset

- Input data: `/mnt/new3T/gtdbr220/eval_runs/gtdb_random_genus_vibrio_20260609/vibrio_pair_predictions_with_anim.tsv`
- Ground truth: existing bidirectional ANIm mean in the source table.
- Sample count: 120 selected pairs from 32 unique Vibrio genomes.
- Selection criteria: deterministic ANIm-spread sample, sorted by `anim_identity_mean` and taking 120 evenly spaced pairs.
- Storage location: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-18_minco_ctxmash_vs_mash_anim`

## Methods

Key parameters:

```text
Minco sketch: -S 10000, default coden context/object pattern, -p8
Minco default/best ANI: minco ani -s -1 -f0 -n0 -t0 -m0
Minco context-Mash ANI: minco ani -s -5 -f0 -n0 -t0 -m0
Minco AAF sanity check: minco ani -s -6 -f0 -n0 -t0 -m0
Mash default: mash sketch -s 10000, default k=21; ANI proxy = 1 - Mash distance
Mash k32: mash sketch -s 10000 -k 32; ANI proxy = 1 - Mash distance
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv` and `metrics.with_mash_k32.tsv`.

```text
All 120 pairs:
  Minco ctx-Mash MAE = 0.012506, RMSE = 0.016621
  Mash k21 MAE       = 0.064936, RMSE = 0.071186
  Mash k32 MAE       = 0.025499, RMSE = 0.029501

Same-species/high-ANI 12 pairs:
  Minco ctx-Mash MAE = 0.002874
  Mash k21 MAE       = 0.002432
  Mash k32 MAE       = 0.003042
  Minco best MAE     = 0.001866

Cross-species 108 pairs:
  Minco ctx-Mash MAE = 0.013577
  Mash k21 MAE       = 0.071881
  Mash k32 MAE       = 0.027994
```

## Validation

- Checks run: joined all 120 selected truth pairs without missing Minco or Mash predictions.
- Expected behavior: all Minco assembly sketches retain exactly 10,000 contexts.
- Observed behavior: `qaf == raf` for all selected pairs.
- Sanity check: Minco `-s -5` context-Mash and `-s -6` AAF were identical for all 1024 all-vs-all comparisons.

## Important Artifacts

See `artifacts.md` for generated files and timing notes.

## Conclusion

Minco context-Mash ANI is not the same accuracy profile as ordinary Mash ANI. On this mixed Vibrio spread, Minco context-Mash had lower overall MAE than Mash k21 and Mash k32, mainly because Mash underestimates the lower-ANI cross-species pairs. On the high-ANI same-species subset, Mash k21 was slightly more accurate than Minco context-Mash, while Minco default/best was best.

## Paper-Relevant Claim

Tentative: context-minhash features can give a stronger Mash-style ANI estimate than ordinary k-mer Mash for mixed within-genus comparisons, but this needs replication across more genera and with a larger same-species set.

## Caveats

- This is only one Vibrio-derived pair set.
- The same-species/high-ANI subset has only 12 pairs.
- Mash ANI is reported as `1 - Mash distance`, the standard quick proxy, not a calibrated Mash-specific ANIm model.
- Cross-species Pearson is low for all methods because the ANIm range is narrow, roughly 0.828 to 0.858.

## Next Experiment

Repeat this with a larger same-species-only validation set across multiple species, and include Minco recalibrated/default, Minco ctx-Mash, Mash k21, Mash k32, skani, and kssd3a side by side.
