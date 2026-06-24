# Research Progress

## Current Best Quick Recall

The current stage summary is pinned at
`research/STAGE_SUMMARY_2026-06-24.md`. In short: use S1000 unique ZIP-AAF for
marine F1, ctx-marker robust-depth rescue for Toy Mouse abundance, coden15
`Ref_zip_aaf_ani` for ANI reporting, and full S2000 dual evidence for future
refdb development. Edge-EM remains diagnostic.

For abundance quantification, first try the pinned current best in
`research/CURRENT_BEST.md`: ctx-marker MinCO with robust effective depth plus
conservative intra-genus winner rescue. On GTDB Toy Mouse samples 0-2 it beat
Sylph on mean abundance L1, 1.7945 versus 2.2870 percentage points, although
Sylph still had slightly higher mean Pearson and much higher Spearman.

For the mixed six-sample mouse+CAMI3 abundance benchmark, do not use the
row-feature calibration as a default. The best raw MinCO rule found on
2026-06-24 was `marker_l1_ctxobj_cami_blend` with raw `value_sum`
(`all_L1=17.8234`), still worse than Sylph (`all_L1=15.0384`). The best
F1-priority raw MinCO rule was `marker_l1_ctxobj_f1_blend` with
`all_F1=0.7825`, still just below Sylph `0.7835`. Evidence:
`research/experiments/2026-06-24_abundance_feature_calibration/NOTE.md`.

The full S2000 marked-unique refdb design has now passed the index-equivalence
test: virtual marker sizes computed from the full S2000 inverted index exactly
matched the physical ctx-markerdb sizes for all 200,527 refs. This means the
next implementation can keep only the full S2000 payload and compute
marker/shared status during the full-refdb scan. Evidence:
`research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/NOTE.md`.

The follow-up full benchmark shows that this architecture is not yet an overall
accuracy win. On the six-sample mouse+CAMI3 panel, Sylph still leads both F1
and abundance L1 (`F1=0.7835`, `L1=15.0384`), while the best MinCO mixed rows
are `F1=0.7825` and `L1=17.8234`. The full S2000 marked-unique proxy improves
old ctx-marker only slightly (`F1=0.7501`, `L1=20.4785`) and is neutral on
CAMI3. MinCO coden15 still slightly beats Sylph adjusted ANI on the CAMI3
source/ref ANI subset (`MAE=0.007627` versus `0.008038`). Evidence:
`research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/NOTE.md`.

Integrated dual-counter update: `minco ani --readwise-dual-evidence` now emits
full-refdb evidence plus marker-only evidence from the same full S2000 inverted
index. The integrated run reproduced the proxy exactly on the six-sample panel:
best integrated fallback `F1=0.750078`, `L1=20.478490`; integrated marker-only
`F1=0.747536`, `L1=20.513677`. Marker sizes matched the physical markerdb map
for all comparable emitted rows (`0` comparable mismatches). Evidence:
`research/experiments/2026-06-24_integrated_dual_evidence_benchmark/NOTE.md`.
This removes a storage/runtime blocker but confirms that the next improvement
must change the abundance mechanism, not just the refdb representation.

Context-level edge-EM update: the branch now has external spot scores on CAMI3,
marine, and strain. A guarded no-truth trigger improved L1 versus S2000
edge-marker beta 0 on marine (`-0.3083` pp) and strain (`-0.0631` pp), while
leaving CAMI3 unchanged by using beta 0. It does not improve F1 and is still
only a diagnostic abundance module. Evidence:
`research/experiments/2026-06-24_cross_domain_spot_comparison/NOTE.md`.

## 2026-06-20: Replicated CAMI Marine Readwise Profiling Gain

Minco reached an important milestone for metagenome readwise profiling.

On CAMI II marine short-read samples 0-2, the current recipe:

```bash
minco ani -p16 -r REF --qraw reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique \
  --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10
```

beat local Sylph in species-level F1:

```text
method  mean F1  mean precision  mean recall  mean FP
minco   0.860    0.853           0.867        41.3
sylph   0.832    0.798           0.870        60.7
```

The main improvement is false-positive control from `best-diff-unique`
shared-context assignment. ZIP AAF then keeps low-coverage taxa from being
penalized too strongly.

The practical `S=1000` GTDB-only sketch also matched `S=10000` accuracy while
being much cheaper:

```text
method   mean F1  mean time  peak RSS
S1000    0.862    124.5s     3.70 GB
S10000   0.860    169.8s     33.65 GB
```

Conclusion: for this benchmark, `S=1000` is the current practical target. The
next validation should use a different dataset, because CAMI marine samples are
related replicates.

Evidence:

- `research/experiments/2026-06-20_minco_readwise_assignment_fix/NOTE.md`
- `research/experiments/2026-06-20_cami2_marine_extra_samples/NOTE.md`
- `research/experiments/2026-06-20_cami2_marine_s1000_choice/NOTE.md`
- Claims `C003` and `C004` in `research/claims.tsv`

Next actions are tracked in `research/TODO.md`, including the Chinese patent
draft and non-marine benchmark plan.

## 2026-06-21: Hybrid Unique/Split Rule Did Not Generalize

Tested whether `best-diff-split` can rescue recall while `best-diff-unique`
anchors precision on CAMI3 Toy Human Gut sample0 and CAMI2 Marine sample0.

This was a diagnostic hybrid-rule experiment, not a replacement for the earlier
direct CLI marine benchmark where `S=1000 best-diff-unique + ZIP AAF` beat
local Sylph on samples 0-2.

The split mode recovered extra true taxa, but it also introduced shared-context
false positives. A simple unique-primary/split-rescue threshold family did not
beat the best unique baseline on toy and did not beat local Sylph on either
dataset under the joint parameters:

```text
dataset             method             TP   FP   FN   F1
toy GTDB bacteria   unique_relaxed      57   18   51   0.623
toy GTDB bacteria   hybrid_joint_best   47   12   61   0.563
toy GTDB bacteria   Sylph default       55   13   53   0.625
marine GTDB         unique_current      197  34   59   0.809
marine GTDB         split_current       213  51   43   0.819
marine GTDB         hybrid_joint_best   221  64   35   0.817
marine GTDB         Sylph default       222  56   34   0.831
```

Same-sample tuned hybrid thresholds can improve individual datasets, but the
best thresholds differ enough that they should not be promoted to defaults.
Leave-one-dataset-out classifiers trained on only one sample also failed to
generalize.

Evidence:

- `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/NOTE.md`
- Claim `C005` in `research/claims.tsv`

## 2026-06-21: Six-Sample Mixed-Domain RF Call Model Improved F1

Added CAMI III Toy Human Gut samples 1 and 2 from the CAMI download site, then
evaluated CAMI II marine samples 0-2 plus Toy Human Gut samples 0-2 with
whole-sample holdout.

The fixed minco direct rule remains strong for marine but is too strict for Toy
Human Gut. A random-forest call model using joined unique/split readwise
features improved mixed-domain mean F1:

```text
method               samples  mean precision  mean recall  mean F1
loso_rf              6        0.855           0.684        0.745
sylph_default        6        0.810           0.682        0.725
minco_unique_direct  6        0.828           0.582        0.646
```

The top model features were split breadth, split depth CV, unique breadth, and
unique support. This supports the idea that split-mode recall is useful when
controlled by a trained call model.

Evidence:

- `research/experiments/2026-06-21_minco_multisample_call_calibration/NOTE.md`
- Claim `C006` in `research/claims.tsv`

## 2026-06-21: Third-Domain Plant Benchmark Supports No-Leak RF/HGB Ensemble

Added CAMI II plant-associated short-read samples 0-2 as a third domain,
scored in bacteria-only mode because the current comparison is GTDB/Sylph
bacterial references. The first mixed-domain RF result was rechecked and the
non-deployable dataset identity features were removed.

With no dataset labels, RF alone nearly tied Sylph. Adding deployable
strict/relaxed rule-evidence features and averaging RF+HGB probabilities gave
the best nine-sample leave-one-sample-out result:

```text
method               samples  mean precision  mean recall  mean F1  mean FP  mean FN
loso_rf_hgb_avg      9        0.733           0.700        0.698    30.7     33.0
loso_rf              9        0.776           0.642        0.677    23.0     38.4
sylph_default        9        0.719           0.665        0.676    32.0     36.4
minco_unique_direct  9        0.720           0.555        0.598    23.4     46.1
```

The mechanism is call-level calibration over minco evidence, not a new ANI
formula: breadth, depth CV, XnY support, min align fraction, and strict/relaxed
rule evidence recover additional true positives without the huge false-positive
cost of relaxed direct thresholds.

Evidence:

- `research/experiments/2026-06-21_minco_third_domain_plant_calibration/NOTE.md`
- Claim `C007` in `research/claims.tsv`

## 2026-06-21: External CAMI Strain Holdout Weakens Current Learned Filter

Added CAMI II strain-madness short-read samples 0-2 as a true external
holdout. The model was trained and threshold-tuned only on the existing
nine-sample marine/toy-gut/plant panel, then tested once on strain0-2.

The current RF/HGB averaged call model did not beat local Sylph on the new
strain-heavy domain:

```text
method               samples  mean precision  mean recall  mean F1  mean FP  mean FN
sylph_default        3        0.550           0.517        0.527    9.0      9.7
train9_rf_hgb_avg    3        0.543           0.417        0.466    7.3      11.7
minco_split_direct   3        0.601           0.383        0.463    5.3      12.3
train9_hgb           3        0.498           0.400        0.442    8.0      12.0
train9_rf            3        0.518           0.333        0.401    6.3      13.3
```

This is not mostly a sketch recovery failure: minco joined candidates contained
17 of 20 gold species for each strain sample, and the same three gold species
were absent from the current GTDB/taxmap mapping. An oracle threshold tuned on
the external holdout only raised RF/HGB average F1 from 0.466 to 0.478, so the
issue is feature ranking/filter generalization.

Conclusion: do not export the current RF/HGB ensemble as a production default.
Prefer a simpler deployable rule/model until external-domain validation shows a
repeatable gain. Keep strain-madness as a locked holdout while improving the
filter on other data.

Follow-up simple-rule diagnostics showed:

```text
method                         strain mean F1  mean FP+FN
u_or_s_direct_validtax          0.497        16.0
train9_rf_hgb_avg_validtax      0.488        17.3
u_or_s_direct                   0.477        17.3
train9_rf_hgb_avg               0.466        19.0
```

Generic-name suppression helps strain precision but slightly lowers the
nine-sample training F1, so it should be an optional reporting filter rather
than a default. Genus-capped relaxed rescue rules did not improve over
`u_or_s_direct_validtax`.

Evidence:

- `research/experiments/2026-06-21_minco_external_strain_holdout/NOTE.md`
- Claim `C008` in `research/claims.tsv`

## 2026-06-24: Integrated Full S2000 Dual Evidence and Aggregate EM Proxy

Implemented `minco ani --readwise-dual-evidence`, which appends marker-only
evidence columns while keeping the normal full-S2000 evidence columns. The
marker sizes matched the physical S2000 markerdb PSMP map for all comparable
emitted rows across the six-sample mouse+CAMI3 benchmark.

The integrated implementation confirmed the previous proxy result: it is useful
for storage/runtime and future experiments, but it did not beat the current
mixed MinCO rows or Sylph for F1/abundance. A follow-up aggregate shared-mass
assignment test also failed: unrestricted full add-back created huge FP counts,
and marker-present-only EM collapsed to marker-only presence.

```text
method                                      all_F1   all_L1
Sylph baseline                              0.7835   15.0384
MinCO current mixed F1-priority             0.7825   18.0146
integrated aggregate EM proxy               0.7752   18.0612
MinCO current mixed L1-priority             0.7710   17.8234
integrated marker-only                      0.7475   20.5137
```

Conclusion: the next EM attempt needs context-level or read-level ambiguity
groups. Aggregate per-reference full evidence totals are not enough.

Evidence:

- `research/experiments/2026-06-24_integrated_dual_evidence_benchmark/NOTE.md`
- Claim `C016` in `research/claims.tsv`

## 2026-06-24: Context-Level Ambiguity EM Prototype on Mouse0

Added an experimental MinCO readwise ambiguity edge dump controlled by
`MINCO_READWISE_EDGE_OUT` and related environment variables. The dump emits
read/context candidate-reference edges with `gid`, `diff`, `best_diff`,
`candidate_refs`, `selected_refs`, and `selected_by_mode`.

On full Toy Mouse sample0, the dump produced 4,952,131 edge rows across 671,619
ambiguous read/context groups. After mapping to GTDB and filtering to the
marker-supported species set, the EM scorer used 546,908 species-edge rows
across 434,448 groups.

The simple unfiltered result was negative. A stricter filtered follow-up found
a tiny abundance-L1 gain on sample0, but this did not transfer to samples 1-2:

```text
method                              F1        L1 pp
marker_l1_per_read_profile          0.923077  1.800456
raw best edge-EM blend              0.923077  1.801188
filtered selected_group_species2    0.923077  1.792004
edge-only selected EM               0.923077  42.135253
```

Mouse0-2 fixed-beta validation:

```text
method                              mean_F1   mean_L1
marker_l1_per_read_profile          0.934155  3.036063
filtered_group_species2_beta0.0001  0.934155  3.036827
filtered_group_species2_beta0.0045  0.934155  3.095701
```

Conclusion: real ambiguity groups do not automatically solve abundance.
Unfiltered ambiguous context mass is badly biased relative to true abundance.
The filtered `selected_group_species2` small-blend row is rejected as a fixed
default because the best mouse0-2 mean is beta `0`.

Evidence:

- `research/experiments/2026-06-24_context_level_ambiguity_em_prototype/NOTE.md`
- Claim `C017` in `research/claims.tsv`
