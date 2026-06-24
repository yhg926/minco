# abundance_feature_calibration

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory used: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: `main`, `b55b4d3-dirty`
- Working tree: dirty; status and diffstat are recorded in `provenance/code_status.txt`.
- Main script: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_abundance_feature_calibration/search_feature_calibration.py`
- Input experiment: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_addback_abundance_normalization_search`

## Question

Can a learned abundance estimator using per-selected-species MinCO features improve abundance L1 beyond raw MinCO effective-depth scores and close the gap to Sylph on the six-sample mouse plus CAMI3 benchmark?

## Hypothesis

If MinCO abundance error is recoverable from row-level evidence, a leave-one-sample-out model using depth, breadth, XnY, AF, VMR, marker size, and marker/ctx+obj source features should improve held-out L1 over raw value normalization. If learned models fail to beat the raw baseline, the remaining bottleneck is likely read-level mass assignment or an EM-style abundance solve rather than row-level calibration.

## Dataset

- Samples: 6 cached benchmark samples.
- Mouse toy gut: `mouse_gtdb` sample IDs `0`, `1`, `2`, GTDB species truth.
- CAMI3 toy gut: `cami3_ncbi` sample IDs `0`, `1`, `2`, NCBI species-taxid truth.
- Selected-call feature rows: 4,318 rows plus header in `results/selected_call_features.tsv`.
- Truth profile rows: 588 rows plus header in `results/truth_profiles.tsv`.
- Availability status: available in current workspace on 2026-06-24.

## Methods

The script imports the previous abundance-normalization scorer and builds selected-call feature rows for nine fixed callsets:

```text
marker_l1_only
marker_f1_only
ctxobj_f1_only
marker_l1_ctxobj_f1_all
marker_l1_ctxobj_l1_blend
marker_l1_ctxobj_f1_blend
marker_l1_ctxobj_cami_blend
marker_l1_ctxobj_f1_rel0002
marker_f1_ctxobj_f1_all
```

Presence calls are fixed. The model only changes the positive abundance score assigned to each selected taxon. F1 is therefore a callset property, and L1 measures abundance assignment quality for that fixed callset.

Scoring uses the full per-sample truth profile, including truth taxa absent from selected-call rows. An earlier smoke test incorrectly counted only selected truth taxa; that was fixed before the recorded broad run.

Raw baselines:

```text
value_marker_first
value_sum
value_max
sum_norm
```

Leave-one-sample-out model sweep:

```text
models: Ridge, Ridge(alpha=10), Huber, Tweedie(log), HistGradientBoosting, RandomForest, ExtraTrees
target transforms: raw, sqrt, log100, log1000 where compatible
feature sets: value-only, row features without dataset labels, row features with dataset labels
sample weights: none, positive5, abundance-weighted
sampled configs: 120 model specifications across the full model grid
held-out folds: each of the 6 samples
```

## Results

Main comparison:

```text
method              callset                    all_F1      all_L1      mouse_F1    mouse_L1    cami3_F1    cami3_L1
Sylph               baseline                   0.783546    15.038388   0.963606    2.287032    0.603485    27.789743
best LOSO L1        marker_l1_ctxobj_cami_blend 0.770982   18.302036   0.884563    4.068283    0.657402    32.535789
best LOSO F1->L1    marker_l1_ctxobj_f1_blend   0.782475   18.455096   0.918887    4.214073    0.646062    32.696119
best raw L1         marker_l1_ctxobj_cami_blend 0.770982   17.823364   0.884563    2.989748    0.657402    32.656979
```

Best held-out model by L1:

```text
huber_value_raw_positive5
```

Best raw L1 rule:

```text
raw_value_sum on marker_l1_ctxobj_cami_blend
```

Best F1-priority raw rule:

```text
raw_value_sum on marker_l1_ctxobj_f1_blend
all_F1=0.782475, all_L1=18.014617
```

The learned row-level models did not beat the raw `value_sum` baseline. The new raw `value_sum` result is slightly better than the previous sampled abundance-transform search, but it still does not beat Sylph on combined L1 or combined F1.

## Validation

- `python3 -m py_compile` passed.
- Feature/truth table generation run completed with exit status 0.
- Cached 120-model leave-one-sample-out run completed with exit status 0.
- Broad cached run runtime: 26:45.83 wall clock.
- Broad cached run peak RSS: 231,668 KB.
- Feature generation peak RSS: 2,712,612 KB.
- Output row counts:
  - `calibration_rule_summary.tsv`: 1,117 lines including header.
  - `calibration_sample_metrics.tsv`: 6,697 lines including header.
  - `selected_call_features.tsv`: 4,319 lines including header.
  - `truth_profiles.tsv`: 589 lines including header.

## Conclusion

Row-level learned abundance calibration did not solve the abundance gap. The best held-out model is worse than raw `value_sum`, and raw `value_sum` remains worse than Sylph on combined L1.

This argues that the remaining abundance error is not just a smooth per-row correction problem. The next serious direction should be read-level mass redistribution, EM-style assignment over candidate species, or a constrained abundance solve that uses shared evidence jointly instead of scoring each selected species independently.

## Paper-Relevant Claim

Negative result: On this six-sample mouse plus CAMI3 panel, row-level leave-one-sample-out abundance calibration did not improve over raw MinCO value normalization and did not close the abundance gap to Sylph.

## Caveats

- Only six samples were used.
- Model selection is still over this same six-sample panel; leave-one-sample-out reduces but does not eliminate benchmark overfitting.
- Mouse truth is GTDB species, while CAMI3 truth is NCBI species taxid.
- Feature calibration does not change ANI estimates; ANI was not reevaluated in this experiment.
- Presence F1 is fixed by callset, so this experiment only evaluates abundance scoring conditional on selected taxa.

## Next Experiment

Prototype a read-level or context-level EM abundance estimator for the same fixed callsets. The estimator should distribute ambiguous evidence across selected candidate species using marker uniqueness, best-diff evidence, and depth consistency, then compare L1 against raw `value_sum` and Sylph under the same full-truth scoring.
