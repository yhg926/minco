# minco Nayfach ANIm recalibration pilot

Date: 2026-06-18
Author/agent: Codex
Project: minco
Code commit: NA; workspace `.git` was not readable by `git status`
Tool version: `bin/minco 0.1`

## Question

Can the default `minco ani` assembled-genome ANI calibration be improved by
training on the Nayfach ANIm label table, instead of reusing the stale KSSD3
calibration?

## Hypothesis

The raw minco CtxMoE distance can be improved by refitting the MoE denominator
bases and linear coefficients on minco context-minhash features. A compact HGB
calibration using the recalibrated raw ANI, refAF, shared context count, and
object-difference rates should further reduce ANI error.

## Dataset

- Label table:
  `/mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv`
- FASTA root:
  `/mnt/new3T/skani_data/Nayfach_data/fna`
- Label rows available: 1,075,493 pairwise ANIm rows.
- FASTA files available: 52,515 `.fna` files.
- Stratification bins:
  `<0.90`, `0.90-0.94`, `0.94-0.95`, `0.95-0.97`, `0.97-0.99`, `>=0.99`.
- Final pilot sample: 12,000 pairs, 2,000 per bin.
- Split: deterministic 70/30 train/test stratified by ANI bin:
  8,400 train pairs and 3,600 held-out test pairs.

Available pair counts in the source label table:

```text
lt090       75782
090_094     64049
094_095     40595
095_097     97786
097_099     600933
gte099      196348
total_rows  1075493
```

## Methods

Feature extraction command family:

```bash
minco ani -S 10000 -p1 -f0 -n0 -t0 --raw-output -s3 -o pair.tsv REF.fna QRY.fna
```

Final extraction was run through
`scripts/minco_nayfach_calibration_pilot.py` with:

```text
--sketch-size 10000
--per-bin 2000
--jobs 16
--minco-threads 1
--seed 20260618
```

Models evaluated:

- `minco_raw_ctxmoe`: original CtxMoE ANI from minco raw output.
- `kssd3_t10_raw_reference`: old KSSD3 T10 ANI column from the label table,
  kept as a reference only.
- `moe_fixedp_twostage`: refit the 3-way MoE linear coefficients with the
  current denominator parameters, including an ANI>=0.95 branch.
- `moe_opt_twostage`: optimized MoE denominator bases plus linear coefficients;
  this was installed in `model_ani.h`.
- `moe_style_ridge_full`: linear/ridge correction using the full extracted
  minco diagnostic feature set.
- `hgb_full`: HGB using the full extracted diagnostic feature set.
- `hgb_refaf11_dropin`: compact 11-feature HGB matching the existing C
  `refaf_hgb_predict_ani()` signature.

The drop-in HGB feature order is:

```text
0 raw_ani
1 ref_af
2 log1p_xny
3 diff_obj_rate
4 diff_section_rate
5 mut2_rate
6 ref_aaf_ani
7 raw_minus_ref_aaf
8 ref_af_lt_0_2
9 ref_af_0_2_to_0_5
10 ref_af_ge_0_5
```

## Results

Initial 12,000-pair run before MoE replacement:

- Wall time: 3:11.32.
- Peak RSS: 286,328 KB.
- Feature extraction throughput: ended at 70.50 pair/s with 16 concurrent
  direct-pair minco jobs.
- Successful feature rows: 12,000 / 12,000.

After installing the optimized MoE model, the same 12,000 sampled pairs were
rerun and the HGB layer was retrained from the new raw MoE score:

- Wall time: 1:46.58.
- Peak RSS: 281,476 KB.
- Feature extraction throughput: ended at 115.57 pair/s with 16 concurrent
  direct-pair minco jobs.
- Successful feature rows: 12,000 / 12,000.

Held-out test metrics for the final installed model stack:

```text
model                 n     bias       MAE       RMSE      max_abs   Pearson_r
old_raw_ctxmoe        3600  -0.000352  0.003283  0.005179  0.043127  0.990837
moe_fixedp_twostage   3600  -0.000052  0.003059  0.004597  0.032649  0.992734
moe_opt_twostage      3600   0.000037  0.002916  0.004329  0.028982  0.993561
hgb_refaf11_dropin    3600  -0.000069  0.002503  0.003822  0.030442  0.994980
```

At the 0.95 ANI cutoff, old raw CtxMoE made 152 held-out mistakes
(118 false positives, 34 false negatives). The optimized MoE also made 152
mistakes, but with lower MAE/RMSE and a different error balance
(128 false positives, 24 false negatives). The final drop-in HGB layer made
92 mistakes (45 false positives, 47 false negatives).

By-band held-out MAE, raw vs drop-in HGB:

```text
band      old_raw_MAE  opt_moe_MAE  final_HGB_MAE
lt090     0.005404     0.004192     0.003632
090_094   0.005225     0.004692     0.004186
094_095   0.003795     0.003551     0.002115
095_097   0.002611     0.002523     0.002597
097_099   0.001689     0.001700     0.001673
gte099    0.000972     0.000842     0.000814
```

## Code Change

Replaced the `NUM_CODENS == 11` MoE constants in `minco_core/src/model_ani.h`
with the optimized minco/Nayfach parameters. The previous header is saved as:

```text
research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_ani.before_moe_recalibration.h
```

Replaced `minco_core/src/model_refaf_hgb.h` with the generated 11-feature HGB
model trained from the regenerated new-MoE feature table. Previous HGB headers
are saved as:

```text
research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_refaf_hgb.before_minco_recalibration.h
research/experiments/2026-06-18_minco_nayfach_anim_recalibration_pilot/model_refaf_hgb.before_moe_recalibrated_hgb.h
```

Also updated:

- `README.md`
- `docs/USER_MANUAL.md`
- `minco_core/src/command_ani_wrapper.c` CLI help text

The fixed-p MoE coefficients remain saved in
`pilot_10k_perbin2000/fixedp_moe_coefficients.tsv` as an earlier candidate.
The installed MoE came from
`pilot_10k_perbin2000/moe_opt_long/optimized_moe_model_ani_arrays.hfrag`.

## Validation

Commands/checks:

- `make`
- `make test`
- `bin/minco ani --help | grep -A4 'CtxMoE is'`
- Binary validation of compiled raw CtxMoE (`-s3`) on all 3,600 held-out pairs.
- Binary validation of compiled default Best (`-s1`) on all 3,600 held-out
  pairs using the HGB retrained after MoE replacement.

Binary-level validation over all held-out pairs:

```text
metric                  n     bias       MAE       RMSE      max_abs
compiled_CtxMoE_vs_ANIm 3600   0.000037  0.002916  0.004329  0.028982
compiled_Best_vs_ANIm   3600  -0.000069  0.002503  0.003822  0.030442
```

The compiled Best-vs-Python drop-in MAE was `2.54e-7` ANI, i.e. only
printed-output rounding. All 3,600 default Best outputs used `BestDist`; all
3,600 raw outputs used `CtxMoE`; all validation rows completed with status `ok`.

## Important Artifacts

See `artifacts.md` for paths and sizes.

## Conclusion

The Nayfach-trained MoE replacement improves raw CtxMoE substantially for the
sampled assembly-to-assembly path: held-out MAE drops from `0.003283` to
`0.002916`. The matching HGB layer improves the default Best ANI further to
`0.002503` MAE. Both compiled paths were validated in the binary.

## Paper-Relevant Claim

For fixed-size context-minhash sketches at 10,000 contexts, a minco-trained
MoE distance model improves raw ANIm agreement, and a compact context/object
HGB layer improves the default assembled-genome ANI further while remaining
small enough to embed directly in pure C.

## Caveats

- This is still a stratified 12,000-pair pilot, not full 1.07M-pair training.
- The model is calibrated for default coden context/object sketches and
  `--sketch-size 10000`; larger or smaller sketch sizes should be checked
  separately.
- The old KSSD3 T10 column in the label table is a reference feature only and
  is not a fair direct comparison to current minco unless recomputed under the
  same sampling/design.
- The drop-in HGB intentionally uses only the feature set already available to
  the C calibration function. A freer HGB was slightly more accurate in some
  runs but would require a wider C feature interface.
- The calibration is assembly-to-assembly. Read/`--qraw` modes still need a
  separate calibration path.

## Next Experiment

Run a larger training set, e.g. 60,000 to 120,000 stratified pairs, then validate
against a completely separate species/genus holdout. Also repeat with
`--sketch-size 20000` and with read-derived sketches before changing read-mode
defaults.
