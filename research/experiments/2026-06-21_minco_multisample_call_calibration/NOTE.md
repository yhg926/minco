# Minco multi-sample call calibration

Date: 2026-06-21
Author/agent: Codex
Project: minco
Code commit: b55b4d3 plus dirty worktree changes
Tool version: minco 0.1, Sylph `/home/ubuntu/yihuiguang/bin/sylph`

## Question

Can minco improve readwise species calls beyond the fixed direct threshold by
training a multi-sample call model on unique and split features, while keeping
direct CLI scoring separate from diagnostic model scoring?

## Hypothesis

The fixed direct rule is strong for CAMI marine but too strict for CAMI III Toy
Human Gut. A domain-aware model using both unique and split support, breadth,
depth, abundance, and ZIP features should recover gut recall without losing too
much marine precision.

## Dataset

Inputs:

- CAMI II marine samples 0-2, cached reads under `/tmp/cami_marine_sample*.fq.gz`
  and gold profile `/tmp/gs_marine_short.profile`.
- CAMI III Toy Human Gut short-read samples 0-2. Sample0 was already local.
  Samples 1 and 2 were downloaded from the CAMI Toy Human Gut download page to
  `/mnt/new3T/minco_cami3_toygut_extra_20260621`.
- Toy Human Gut gold profiles came from
  `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles.tar.gz`.
- Reference sketch:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`.
- Tax map:
  `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`.
- Sylph database:
  `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`.

External source checked: CAMI states Toy Human Gut has 20 short-read samples,
100 Gbp total, 2x150 bp reads, and provides sample-specific read tarballs plus
gold standards at `https://cami-challenge.org/datasets/toy-human-gut/#download`.

## Methods

Minco candidate generation:

```text
minco ani -p16 -r REF --qraw READS \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0

minco ani -p16 -r REF --qraw READS \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani zip-aaf \
  -m0 -f0 -n0 -t0
```

Baselines:

- `minco_unique_direct`: direct row-level rule matching the proven CLI recipe:
  `XnY_ctx>=10`, `ANI>=0.94`, `Real_min_align_fraction>=0.05`.
- `minco_split_direct`: same direct rule on split candidates.
- `minco_unique_relaxed`: diagnostic relaxed rule:
  `XnY_ctx>=1`, `ANI>=0.90`, `Ref_breadth>=0.005`.
- `sylph_default`: Sylph profile output mapped to the same species taxids.

Experimental model:

- Taxid-level joined unique/split features.
- Whole-sample holdout, never row-random train/test split.
- Models: logistic regression, random forest, HistGradientBoosting.
- Per-fold probability threshold tuned only on training samples.
- Toy scoring is restricted to bacterial gold taxa because the local GTDB/Sylph
  reference is not a full viral/eukaryotic database.

Scripts:

- `scripts/calibrate_multisample_calls.py`
- `scripts/grid_loso_thresholds.py`

## Results

Six-sample leave-one-sample-out mean metrics are in
`mean_summary_marine_toy0_2.tsv`.

```text
method                 samples  mean precision  mean recall  mean F1  mean FP  mean FN
loso_rf                6        0.855           0.684        0.745    25.0     46.8
loso_hgb               6        0.769           0.720        0.740    38.3     42.3
sylph_default          6        0.810           0.682        0.725    36.5     46.8
minco_split_direct     6        0.801           0.614        0.665    32.7     54.3
minco_unique_direct    6        0.828           0.582        0.646    24.5     58.5
minco_unique_relaxed   6        0.618           0.739        0.646    136.0    38.3
```

Interpretable threshold grid:

```text
method          samples  mean precision  mean recall  mean F1  mean FP  mean FN
loso_threshold  6        0.789           0.690        0.718    42.2     44.0
fixed_default   6        0.828           0.582        0.646    24.5     58.5
```

Marine-only leave-one-sample-out showed only a small RF gain:

```text
loso_rf              mean F1 0.863
minco_unique_direct  mean F1 0.861
sylph_default        mean F1 0.832
```

Top RF feature importances from the final all-sample fit:

```text
s_Ref_breadth_max                 0.175
s_Ref_depth_cv_min                0.148
u_Ref_breadth_max                 0.122
u_XnY_ctx_max                     0.094
u_Real_min_align_fraction_max     0.077
s_Real_min_align_fraction_max     0.065
s_XnY_ctx_max                     0.056
s_Ref_mean_depth_max              0.051
```

## Runtime Notes

New minco candidate generation:

```text
marine1 unique  2:10 wall, 3.71 GB RSS
marine1 split   1:55 wall, 4.15 GB RSS
marine2 unique  2:11 wall, 3.71 GB RSS
marine2 split   1:55 wall, 4.15 GB RSS
toy1 unique     2:09 wall, 3.33 GB RSS
toy1 split      1:54 wall, 3.66 GB RSS
toy2 unique     1:52 wall, 3.29 GB RSS
toy2 split      1:54 wall, 3.59 GB RSS
```

New Sylph toy baselines:

```text
toy1 sketch     0:55 wall, 0.50 GB RSS
toy1 profile    2:29 wall, 19.3 GB RSS
toy2 sketch     0:57 wall, 0.50 GB RSS
toy2 profile    1:05 wall, 19.2 GB RSS
```

## Validation

- `python3 -m py_compile` passed for both scripts.
- `calibrate_multisample_calls.py` was rerun after adding RF feature-importance
  export and reproduced the six-sample summary.
- Direct minco row-level scoring is kept separate from model scoring.
- Whole samples, not rows, are held out for model evaluation.

## Conclusion

This is the first mixed-domain calibration result that improves over both the
fixed minco direct threshold and the local Sylph baseline on the tested samples.
The improvement comes from using both split and unique features, especially
split breadth and split depth-CV, while preserving enough unique precision.

The RF model should not be installed as a default yet. It is promising, but it
is still trained and tested on only six local samples from two CAMI domains.
The next step is to add at least one more different CAMI domain, preferably
plant-associated or clinical/pathogen data if a species-level truth profile can
be mapped cleanly.

## Paper-Relevant Claim

Preliminary claim: on CAMI II marine samples 0-2 plus CAMI III Toy Human Gut
samples 0-2, a random-forest minco call model using joined unique/split readwise
features improved local species-level mean F1 to `0.745`, compared with `0.725`
for local Sylph and `0.646` for the fixed minco direct threshold.

## Caveats

- Local species-taxid scoring, not official OPAL web evaluation.
- Minco and Sylph use different GTDB releases/databases.
- Toy scoring is bacteria-only.
- RF is not yet exported to pure C or installed in minco.
- Only marine and toy gut domains are represented in the six-sample model.

## Next Experiment

Add a third distinct domain. Candidate CAMI datasets from the CAMI datasets
page include plant-associated, clinical pathogen detection, strain-madness, toy
mouse gut, and Toy Human Microbiome Project. For release work, either train a
larger model that can be exported to C or distill the RF into domain-specific
threshold rules.
