# minco versus Mash Nayfach heldout ANI comparison

Date: 2026-06-18
Author/agent: Codex
Project: minco
Tool versions: `minco 0.1`; `mash 2.3`

## Question

How does Mash with `-s 10000` compare to minco on the same Nayfach ANIm
held-out pairs?

## Dataset

Input pair source:

```text
research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/sample_pairs.tsv
```

Minco prediction source:

```text
research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/pilot_10k_perbin2000_moe_recalibrated/model_predictions.tsv
```

Sample:

- 600 held-out test pairs.
- 100 pairs from each ANIm band:
  `090_094`, `094_095`, `095_097`, `097_099`, `gte099`, `lt090`.
- 1,109 unique genomes were sketched by Mash.

## Methods

Mash sketches were generated once per unique genome:

```bash
mash sketch -s 10000 -o PREFIX genome.fna
```

Pair distances were computed with:

```bash
mash dist ref.msh qry.msh
```

For direct comparison to ANIm, Mash ANI was computed as:

```text
Mash_ANI = 1 - Mash_distance
```

Compared metrics:

- `mash_ani_1_minus_dist`
- `minco_raw_ctxmoe`
- `minco_best_hgb`

## Results

All 600 Mash pair distances completed successfully.

Aggregate metrics against ANIm:

```text
model                  n    MAE       RMSE      Pearson_r  0.95 errors
mash_ani_1_minus_dist  600  0.007854  0.010525  0.964156   76
minco_raw_ctxmoe       600  0.002822  0.004263  0.993703   26
minco_best_hgb         600  0.002318  0.003510  0.995725   17
```

The 0.95 cutoff errors were:

```text
Mash:      60 FP, 16 FN
minco raw: 21 FP, 5 FN
minco Best: 6 FP, 11 FN
```

By-band MAE:

```text
band      Mash      minco_raw  minco_Best
lt090     0.008117  0.004053   0.003156
090_094   0.009708  0.004455   0.003772
094_095   0.007792  0.003390   0.002038
095_097   0.005622  0.002700   0.002661
097_099   0.007048  0.001515   0.001476
gte099    0.008837  0.000821   0.000806
```

## Validation

- `mash --version` returned `2.3`.
- The comparison script compiled with `python3 -m py_compile`.
- `mash_sketch_status.tsv` shows all 1,109 sketches succeeded or were cached.
- `mash_minco_anim_pairs.tsv` shows all 600 pair distances had status `ok`.

## Conclusion

On this balanced 600-pair Nayfach held-out subset, minco is substantially more
accurate than Mash `-s 10000` when both are compared to ANIm. Minco Best/HGB
has about 3.4-fold lower MAE than Mash ANI (`0.002318` vs `0.007854`) and many
fewer 0.95-threshold mistakes (`17` vs `76`).

## Caveats

- This is a 600-pair balanced subset, not the full 3,600 held-out set.
- Mash used default k-mer length with `-s 10000`; only sketch size was matched
  to minco.
- Mash distance was converted to ANI as `1 - distance`; this is the standard
  simple ANI-style comparison, but Mash itself reports distance.

## Next Experiment

Run the same script with `--per-bin 600` for the full 3,600 held-out split if a
full-table comparison is needed for publication figures.
