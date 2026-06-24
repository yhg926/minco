# Cross-Domain Spot Comparison

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO benchmarking

## Question

How do the currently relevant MinCO rows compare with Sylph on one sample each
from CAMI3 ToyGut, CAMI2 Marine, and CAMI2 Strain Madness, and does the
context-level edge-EM branch have a valid external score?

This note combines cached benchmark rows with new external edge-EM profiler and
NCBI-truth scorer runs.

## Inputs

The cached comparison uses existing result tables:

- CAMI3 ToyGut sample0:
  `research/experiments/2026-06-23_cami3_toygut_abundance_rescue/`
- CAMI2 Marine sample0:
  `research/experiments/2026-06-20_cami2_marine_s1000_choice/`
- CAMI2 Strain Madness sample0:
  `research/experiments/2026-06-21_minco_external_strain_holdout/`

The edge-EM comparison uses a new external scorer:

- Scorer:
  `research/experiments/2026-06-24_cross_domain_spot_comparison/score_external_edge_em.py`
- Edge-EM prototype source:
  `research/experiments/2026-06-24_context_level_ambiguity_em_prototype/search_mouse0_filtered_edge_em.py`
- S2000 marker reference:
  `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Edge/profile dumps:
  `/tmp/minco_context_em_external_20260624/`
- Scored outputs:
  `research/experiments/2026-06-24_cross_domain_spot_comparison/results/external_edge_em_summary.tsv`

## Edge-EM Scope Note

The previous caveat that edge-EM had no external score is now fixed. CAMI3
ToyGut sample0, CAMI2 Marine sample0, and CAMI2 Strain Madness sample0 were
scored against their NCBI species-truth namespaces.

The edge-EM rows are still diagnostic, not a promoted default:

- The scorer keeps the marker callset fixed, so TP/FP/FN/F1 do not change with
  beta; only abundance mass changes.
- The reported edge-EM row uses the best positive beta for each sample, so it is
  a per-sample diagnostic upper bound rather than a universal parameter.
- This test uses the `selected_group_species2` branch with `alpha=0` and
  `weight_mode=event`, inherited from the Toy Mouse prototype.

## Results

Summary table:

```text
dataset                 sample  winner by F1               best cached MinCO F1  Sylph F1  edge-EM best-positive F1
CAMI3 ToyGut            0       Sylph                       0.569892              0.606383  0.536313
CAMI2 Marine            0       MinCO S1000 ZIP-AAF         0.866019              0.831461  0.839187
CAMI2 Strain Madness    0       Sylph                       0.484848              0.578947  0.500000
```

Fair full-profile L1 comparison using one scorer:

```text
dataset               method                                F1        L1 pp     Pearson
CAMI3 ToyGut          Sylph                                 0.606383  27.2468   0.954895
CAMI3 ToyGut          current best MinCO ctxobj robust      0.569892  32.7134   0.948703
CAMI3 ToyGut          edge marker beta=0                    0.536313  32.9929   0.946652
CAMI3 ToyGut          adaptive edge-EM trigger              0.536313  32.9929   0.946652

CAMI2 Marine          current best MinCO S1000 ZIP-AAF      0.866019  33.4626   0.970537
CAMI2 Marine          Sylph, regenerated/profile rescored   0.844444  33.1027   0.950493
CAMI2 Marine          edge marker beta=0                    0.839187  15.3138   0.981656
CAMI2 Marine          adaptive edge-EM trigger              0.839187  15.0055   0.982462

CAMI2 Strain Madness  Sylph                                 0.578947  18.8281   0.990033
CAMI2 Strain Madness  current best MinCO RF/HGB average     0.484848  72.0718   0.779049
CAMI2 Strain Madness  edge marker beta=0                    0.500000  34.5766   0.937671
CAMI2 Strain Madness  adaptive edge-EM trigger              0.500000  34.5135   0.937632
```

The fair L1 scorer collapses predictions to NCBI species taxid, uses the same
GTDB r232 metadata accession map for regenerated marine/strain rows, normalizes
positive predicted abundance to sum 1, and computes L1 over the union of truth
and predicted species. Marine Sylph F1 differs from the earlier cached row
because this fair scorer regenerated the profile and remapped accessions through
GTDB r232 metadata instead of the lost temporary full-lineage taxmap.

Adaptive edge-EM trigger pilot:

```text
rule                                            beta
marker_raw_targets >= 150                       0.02
marker_raw_targets <= 25                        0.0015
otherwise                                       0
```

Observed adaptive effect versus the S2000 edge marker beta=0 baseline:

```text
dataset               beta    delta L1 pp  delta F1
CAMI3 ToyGut          0       0.0000       0.000000
CAMI2 Marine          0.02   -0.3083       0.000000
CAMI2 Strain Madness  0.0015 -0.0631       0.000000
```

This satisfies the narrow adaptive-EM test condition: it improves marine and
strain abundance versus the edge marker baseline and does not hurt CAMI3. It is
still a pilot rule, because it was designed on only these three samples and does
not improve the marker callset or F1.

Original cached/edge spot rows retained for provenance:

CAMI3 sample0, NCBI bacterial species truth:

```text
method                                TP  FP  FN  F1        L1 pp     Pearson
Sylph                                 57  13  61  0.606383  27.2468   0.954895
MinCO ctxobj robust depth             53  15  65  0.569892  32.7134   0.948703
MinCO S2000 edge-EM marker beta=0     48  13  70  0.536313  32.9929   0.946652
MinCO edge-EM best positive beta      48  13  70  0.536313  32.9950   0.946658
```

Marine sample0, NCBI species truth excluding unidentified/unidentified plasmid:

```text
method                                TP   FP  FN  F1        L1 pp     Pearson   TP MAE
MinCO S1000 unique ZIP                223  36  33  0.866019  NA        0.981267  0.003001
Sylph                                 222  56  34  0.831461  NA        0.984781  0.002551
MinCO S2000 edge-EM marker beta=0     227  58  29  0.839187  15.3138   0.981656  0.000418
MinCO edge-EM best positive beta      227  58  29  0.839187  15.0055   0.982462  0.000404
```

Strain-madness sample0, NCBI bacterial species truth:

```text
method                                TP  FP  FN  F1        L1 pp     Pearson
Sylph                                 11  7   9   0.578947  NA        NA
MinCO RF/HGB average                   8  5   12  0.484848  NA        NA
MinCO S2000 edge-EM marker beta=0     10  10  10  0.500000  34.5766   0.937671
MinCO edge-EM best positive beta      10  10  10  0.500000  34.5135   0.937632
```

Edge dump coverage:

```text
dataset          edge rows   cap hit  selected rows  selected groups  selected taxids
CAMI3 ToyGut     10,822,782  no       504,054        261,026          60
CAMI2 Marine      6,173,755  no        78,180         52,156          268
Strain Madness    4,564,256  no       202,597         80,278          20
```

## Interpretation

The ranking remains domain dependent:

- CAMI3 sample0: Sylph wins by both F1 and abundance L1. Edge-EM is below the
  cached MinCO ctxobj row. The adaptive trigger correctly leaves CAMI3 at
  beta=0, avoiding the positive-beta L1 regression.
- Marine sample0: cached MinCO S1000 unique ZIP-AAF still wins by F1. The
  adaptive edge-EM trigger gives the best fair full-profile L1 among rows tested,
  but does not overtake the current-best MinCO F1 row.
- Strain-madness sample0: Sylph still wins by F1. The S2000 edge marker callset
  is slightly better than the cached RF/HGB MinCO row by F1, and adaptive edge-EM
  gives a small L1 improvement over edge marker beta=0. Sylph is still much
  better by L1 on this sample.

The edge-EM branch now has valid external scores, but those scores do not
support promoting it as the default. It is useful as a diagnostic experiment for
shared-context mass assignment; it still needs a domain-independent trigger or
parameter policy before it can replace the current best MinCO strategy.

## Caveats

- This picks one sample per dataset, as requested.
- Truth namespaces differ: CAMI3 and strain use NCBI bacterial species taxids;
  marine excludes unidentified placeholder taxa.
- The fair L1 section reruns marine current-best MinCO and Sylph sample0, and
  remaps accessions through GTDB r232 metadata. It should be used for abundance
  comparison; the earlier cached marine abundance values were TP-only metrics.
- The edge-EM beta is selected per sample by best positive L1; this is not a
  fair fixed-parameter deployment setting.
- The adaptive edge-EM trigger is a three-sample pilot based only on
  `marker_raw_targets`; it should not be treated as validated across domains.
- Strain RF/HGB is a call model, not a native abundance estimator. Its fair L1
  row attaches abundance from the selected unique/split MinCO feature rows, so
  that abundance number is a diagnostic, not an optimized strain abundance
  model.
- The "best cached MinCO" row is chosen from rows already produced for that
  dataset, not from a new universal parameter search.

## Conclusion

The invalid placeholder caveat is fixed: edge-EM has now been externally scored
on CAMI3 ToyGut sample0, CAMI2 Marine sample0, and CAMI2 Strain Madness sample0.
The adaptive edge-EM trigger passes the narrow test of improving abundance
versus S2000 edge marker beta=0 without hurting CAMI3, but it is not strong
enough to promote edge-EM as the default. Current evidence says to keep
edge-EM as an experimental abundance module while preserving the stronger
non-edge MinCO detection strategy for F1.
