# MinCO Strategy Stage Summary

Date: 2026-06-24

This is the current quick-recall summary after the S1000/S2000 marker, coden15,
dual-evidence, abundance, ANI, and edge-EM experiments.

## Current Defaults

Use different MinCO strategies for different goals. There is not yet one
universal strategy that beats Sylph on F1, abundance L1, and ANI across all
tested domains.

```text
Goal                     Current first choice
Species F1, marine       S1000 GTDB unique ZIP-AAF
Toy Mouse abundance      ctx-marker robust effective depth + intra-genus rescue
Mixed mouse+CAMI3 F1     current mixed F1-priority ctxobj blend
Mixed mouse+CAMI3 L1     current mixed L1-priority ctxobj/cami blend
ANI reporting            coden15 Ref_zip_aaf_ani, not emitted naive ANI
Future refdb work        full S2000 dual evidence, not physical markerdb copies
Edge-EM                  diagnostic abundance module only
```

## Top Measured Rows

### GTDB Toy Mouse Abundance

Best MinCO abundance result remains the old ctx-marker robust-depth rescue:

```text
method                                      mean Pearson  mean L1 pp
MinCO ctx-marker robust depth + rescue      0.999945      1.7945
Sylph reported abundance                    0.999972      2.2870
```

Use this as the first baseline for new abundance work on Toy Mouse-like data.
Do not replace it with coden15 by default.

### CAMI II Marine Presence

The practical marine default remains S1000 unique ZIP-AAF:

```text
method                           sample0 F1  sample0 TP/FP/FN
MinCO S1000 unique ZIP-AAF        0.866019   223/36/33
Sylph                             0.831461   222/56/34
```

This result supports `S=1000` as the practical first reference size for marine
readwise profiling.

### Mixed Mouse+CAMI3 Abundance

On the six-sample mouse+CAMI3 mixed panel, Sylph still leads:

```text
method                            all F1    all L1 pp
Sylph                             0.783546  15.038388
MinCO mixed F1-priority           0.782475  18.014617
MinCO mixed L1-priority           0.770982  17.823364
```

This is the main reason no MinCO abundance strategy should be called universal
yet.

### Fair Spot Full-Profile L1

The latest fair full-profile L1 spot comparison used the same species-level
normalization/scorer for Sylph, current MinCO, and edge-EM where data were
available:

```text
dataset               method                               F1        L1 pp
CAMI3 ToyGut          Sylph                                0.606383  27.2468
CAMI3 ToyGut          current best MinCO ctxobj robust     0.569892  32.7134
CAMI3 ToyGut          adaptive edge-EM trigger             0.536313  32.9929

CAMI2 Marine          current best MinCO S1000 ZIP-AAF     0.866019  33.4626
CAMI2 Marine          Sylph                                0.844444  33.1027
CAMI2 Marine          adaptive edge-EM trigger             0.839187  15.0055

CAMI2 Strain Madness  Sylph                                0.578947  18.8281
CAMI2 Strain Madness  current best MinCO RF/HGB average    0.484848  72.0718
CAMI2 Strain Madness  adaptive edge-EM trigger             0.500000  34.5135
```

Interpretation: adaptive edge-EM is promising for abundance mass assignment on
marine/strain under the S2000 edge-marker branch, but it is not a default
species-detection strategy and does not beat Sylph on CAMI3 or strain.

## Adaptive Edge-EM Pilot

Only use this as a diagnostic abundance module after generating S2000 edge
marker evidence:

```text
if marker_raw_targets >= 150:
    beta = 0.02
elif marker_raw_targets <= 25:
    beta = 0.0015
else:
    beta = 0
```

Observed effect versus S2000 edge-marker beta 0:

```text
dataset               beta    delta L1 pp  delta F1
CAMI3 ToyGut          0       0.0000       0.000000
CAMI2 Marine          0.02   -0.3083       0.000000
CAMI2 Strain Madness  0.0015 -0.0631       0.000000
```

Status: passes the narrow "do not hurt CAMI3, help marine/strain abundance"
test, but is still only a three-sample pilot. It should not be promoted without
larger validation.

## Do Not Promote Yet

```text
Strategy                          Reason
coden15 abundance default          Better than Sylph in some Toy Mouse settings, but worse than old ctx-marker and not CAMI3-general.
edge-EM default                     Does not improve F1; positive beta hurts CAMI3 unless guarded.
aggregate per-reference EM          Not enough; needs true context/read ambiguity groups.
train9 RF/HGB production default    External strain holdout remains weaker than Sylph.
full S2000 dual evidence alone      Architecture is useful, but accuracy is not a standalone win.
emitted naive ANI                   Saturates at 1.0; use coden15 Ref_zip_aaf_ani for reporting.
```

## Next Decision Point

Validate the adaptive edge-EM trigger on a larger panel before changing any
default. The useful next benchmark is not another same-sample beta sweep; it is
multi-domain validation of the no-truth trigger policy with fixed thresholds.

Evidence:

- `research/CURRENT_BEST.md`
- `research/experiments/2026-06-24_cross_domain_spot_comparison/NOTE.md`
- `research/experiments/2026-06-24_context_level_ambiguity_em_prototype/NOTE.md`
- `research/experiments/2026-06-24_integrated_dual_evidence_benchmark/NOTE.md`
- `research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/NOTE.md`
