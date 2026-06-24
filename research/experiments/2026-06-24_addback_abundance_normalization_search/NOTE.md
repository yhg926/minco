# addback_abundance_normalization_search

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory used: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: `main`, `b55b4d3-dirty`
- Working tree: dirty; status and diffstat are recorded in `provenance/code_status.txt`.
- Binary/script/tool: `python3 research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py`
- Previous cached scorer imported by this experiment: `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py`

## Command and Parameter Provenance

- Working directory: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Exact commands: `commands.sh`
- Parameters, arguments, environment settings, seeds, thresholds, and config files: `parameters.tsv`

## Question

Can MinCO close or beat Sylph on abundance L1 while retaining the strong ctx-marker/ctx+obj presence-call behavior by applying post-hoc abundance normalization to cached MinCO calls?

## Hypothesis

If the remaining abundance error is mostly caused by simple scale bias between marker-only calls and ctx+obj add-back calls, then global transforms of MinCO effective depth should reduce L1 without changing presence F1. If L1 remains poor despite good callset oracle bounds, then the bottleneck is a stronger per-species abundance estimator, not more global threshold tuning.

## Dataset

- Samples: 6 cached benchmark samples.
- Mouse toy gut: `mouse_gtdb` sample IDs `0`, `1`, `2`, scored against GTDB species truth and GTDB-normalized abundance.
- CAMI3 toy gut: `cami3_ncbi` sample IDs `0`, `1`, `2`, scored against NCBI species taxid truth and bacterial-normalized abundance.
- MinCO cached input rows are loaded indirectly through the previous scorer:
  `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py`
- Sylph baseline is copied from:
  `/home/ubuntu/yihuiguang/tools/KSSD3mini/2026-06-23_threshold_combo_search_all_metrics/results/sylph_baseline.tsv`
- Availability status: available in the current workspace on 2026-06-24.
- Retention risk: cached analysis rows are workspace-local and should be preserved with the project.

## Methods

The script reuses prepared MinCO rows and truth loaders from the previous all-metric threshold search. It fixes a small number of representative callsets, then searches abundance normalization rules over those fixed callsets without rerunning MinCO.

Callsets tested:

```text
marker_l1_only
marker_f1_only
ctxobj_f1_only
marker_l1_ctxobj_f1_all
marker_l1_ctxobj_l1_all
marker_l1_ctxobj_f1_rel0002
marker_f1_ctxobj_f1_all
marker_l1_ctxobj_cami_all
marker_l1_ctxobj_f1_blend
marker_l1_ctxobj_l1_blend
marker_l1_ctxobj_cami_blend
marker_f1_ctxobj_f1_blend
```

Abundance rule families:

```text
linear add-back: marker value plus ctx+obj value only for non-marker taxa
linear_sum: marker and ctx+obj values summed when both exist
group_mass: marker-normalized and addback-normalized groups blended by fixed mass fraction
optional value transforms: power, median cap, ctx+obj scale, uniform smoothing, median-derived floor
```

The search used deterministic downsampling to evaluate 12,000 rules from the larger grid. The script also reports a corrected oracle L1 lower bound for each callset, where the current benchmark L1 is calculated only over truth taxa; the lower bound is therefore `(1 - selected_truth_mass) * 100`.

## Results

Main comparison:

```text
method          all_F1      all_L1      mouse_F1    mouse_L1    cami3_F1    cami3_L1
Sylph           0.783546    15.038388   0.963606    2.287032    0.603485    27.789743
best MinCO L1   0.781261    18.052614   0.916461    3.007268    0.646062    33.097960
best MinCO F1   0.782475    18.119338   0.918887    3.072346    0.646062    33.166330
```

Best-L1 rule:

```text
marker_l1_ctxobj_l1_blend_linear_mp1.0_ap1.2_mc0.0_ac0.0_as0.2_u0.0_f0.01
```

Best-F1 rule:

```text
marker_l1_ctxobj_f1_all_linear_mp1.0_ap0.8_mc0.0_ac0.0_as0.5_u0.0_f0.0
```

Corrected oracle bounds show remaining potential if abundance assignment improves. For example, `marker_l1_ctxobj_cami_blend` has CAMI3 F1 `0.657402`, actual CAMI3 L1 `32.678916`, and corrected oracle CAMI3 L1 `22.602903`, with selected truth mass `0.773971`. This means the presence calls can cover enough truth mass to beat Sylph's CAMI3 L1 in principle, but the current depth-to-abundance mapping does not reach that bound.

## Validation

- `python3 -m py_compile` passed for `search_abundance_normalization.py`.
- A 2,000-rule smoke run completed.
- The final 12,000-rule run completed with exit status 0.
- Runtime for the final run: 5:15.64 wall clock.
- Peak RSS for the final run: 2,651,940 KB.
- Output row counts:
  - `abundance_rule_summary.tsv`: 12,001 lines including header.
  - `abundance_rule_top200_by_l1.tsv`: 201 lines including header.
  - `abundance_rule_top200_by_f1.tsv`: 201 lines including header.
- A SciPy runtime warning about `MessageStream size changed` was printed, but the scoring completed and produced all expected output files.

## Important Artifacts

See `artifacts.md` for paths to scripts, result tables, and temporary logs.

## Conclusion

This search did not beat Sylph on combined abundance L1 or combined F1. It did reproduce the earlier pattern: MinCO remains stronger on CAMI3 presence F1, but abundance L1 is worse, especially on CAMI3.

The negative result is informative. Simple global transforms of marker depth and ctx+obj add-back depth are not enough. The next useful step is a stronger abundance estimator that uses per-species features or read-level reassignment, not another sweep over only global scale/power/cap constants.

## Paper-Relevant Claim

No positive paper claim from this experiment. This is a negative benchmark that narrows the abundance bottleneck: MinCO call coverage has enough theoretical room, but current abundance assignment is the limiting step.

## Caveats

- Only six cached samples were used.
- The abundance rules were tuned and evaluated on the same six samples; this is exploratory, not an independent validation.
- CAMI3 truth is NCBI-taxid based here, while mouse truth is GTDB species based.
- This experiment does not rerun MinCO or Sylph; it recombines cached result rows.
- The current L1 scorer evaluates truth taxa only, so false-positive abundance affects L1 through normalization mass rather than as separate FP rows.
- ANI was not reoptimized in this specific abundance-normalization run.

## Next Experiment

Export a per-selected-species feature table and test a calibrated abundance estimator with leave-one-sample-out evaluation. Features should include marker effective depth, ctx+obj effective depth, XnY, reliable breadth, hit mean/median depth, depth variance/VMR, ZIP/ztp AF, marker size, source mode, and sample-level totals. If feature calibration still cannot approach the oracle bound, move to read-level assignment or EM-style mass redistribution.
