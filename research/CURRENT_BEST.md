# Current Best MinCO Baselines

Date pinned: 2026-06-24

## Stage Summary

Use `research/STAGE_SUMMARY_2026-06-24.md` as the quick-recall page for the
current top strategies.

Current practical choices:

```text
Species F1, marine       S1000 GTDB unique ZIP-AAF
Toy Mouse abundance      ctx-marker robust effective depth + intra-genus rescue
Mixed mouse+CAMI3 F1     current mixed F1-priority ctxobj blend
Mixed mouse+CAMI3 L1     current mixed L1-priority ctxobj/cami blend
ANI reporting            coden15 Ref_zip_aaf_ani
Future refdb work        full S2000 dual evidence
Edge-EM                  diagnostic abundance module only
```

Latest edge-EM status: the guarded trigger
`marker_raw_targets>=150 -> beta=0.02`, `marker_raw_targets<=25 -> beta=0.0015`,
otherwise `beta=0`, improved abundance versus the S2000 edge-marker baseline on
marine and strain while leaving CAMI3 unchanged. It still does not improve F1
and is not a default.

## Abundance First Try: GTDB Toy Mouse

Use this as the first baseline for abundance quantification work:

```text
ctx-marker + robust effective depth + intra-genus winner rescue
```

Evidence:

- Note: `research/experiments/2026-06-22_minco_sylph_abundance_model/NOTE.md`
- Mean metrics: `research/experiments/2026-06-22_minco_sylph_abundance_model/intragenus_winner_rescue_mean_metrics.tsv`
- Per-sample metrics: `research/experiments/2026-06-22_minco_sylph_abundance_model/intragenus_winner_rescue_sample_metrics.tsv`
- Claim: `C012` in `research/claims.tsv`

Mean GTDB Toy Mouse samples 0-2 abundance metrics:

| method | Pearson | Spearman | L1 percentage points |
|---|---:|---:|---:|
| MinCO ctx-marker robust depth + rescue | 0.999945 | 0.910342 | 1.7945 |
| Sylph reported abundance | 0.999972 | 0.987668 | 2.2870 |

Interpretation: this is the current best MinCO abundance-L1 result. Sylph still has slightly better mean Pearson and much better Spearman/presence recall, so do not describe it as a complete profiler win.

## Exact Logic

Run MinCO with the ctx-marker refdb and readwise split/naive/product filter profile output:

```text
minco_core/bin/minco ani -p16 -r <ctx_markerdb> --qraw <reads> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0
```

Use robust effective depth:

```text
if Reliable_Ref_hit_median_depth >= 20:
    Effective_abundance_depth = Reliable_Ref_hit_median_depth
else:
    Effective_abundance_depth = Reliable_Ref_mean_depth / Reliable_Ref_zip_af^1.05
```

Then apply conservative intra-genus winner rescue:

```text
Start from rows passing the current active gate.
For each rejected row:
  same genus must already have at least one active call;
  ANI_naive_calc >= 0.999;
  XnY_ctx >= 100;
  Effective_abundance_depth >= 5;
  Reliable_Ref_zip_af <= 0.25;
  Effective_abundance_depth >= 3 * max active Effective_abundance_depth in that genus.
Keep only the strongest rescued rejected row per genus.
Normalize abundance from Effective_abundance_depth over active plus rescued rows.
```

## Practical Reminder

For new abundance benchmarks, try this baseline before coden15, full/shared-context refdb, or new markerdb variants. Coden15 improved MinCO presence/F1 slightly in the later Toy Mouse test, but it did not reproduce this abundance advantage over Sylph.

Retest result: coden15 ctx-marker plus the same robust-depth rescue increased mean L1 from 1.7945 to 2.4022 percentage points on Toy Mouse samples 0-2. The old ctx-marker robust-rescue baseline remains the first try for abundance work.

Mechanism retune result: coden15 improves to L1 2.1857 with `Reliable_ztp_af >= 0.35` and VMR trigger lowered from `>50` to `>10`, but this is still worse than old ctx-marker L1 1.7945.

Best coden15 diagnostic so far uses a model-derived depth-adjusted AF floor: `Reliable_ztp_af >= 0.95^24 = 0.291989`, with `XnY_ctx >= 10`, delta trigger `VMR > 50`, delta cutoff `0.02`, robust-depth exponent `1.0`, and median cutoff `10`. It reaches mean L1 2.0091, better than Sylph's 2.2870 but still worse than old ctx-marker. Treat coden15 formula-AF as diagnostic until independently validated.

CAMI3 ToyGut transfer test: the Toy Mouse ranking did not reproduce. On CAMI3 bacterial NCBI species-taxid truth, Sylph remained best by mean L1 and Pearson (`27.7897`, `0.923617`), coden15 formula-AF was second (`33.7127`, `0.909277`), old ctx-marker robust-rescue was third (`39.2328`, `0.898686`), and coden15 original gate was fourth (`45.9801`, `0.896053`). Formula-AF transfers better than old ctx-marker on CAMI3, but does not beat Sylph there. Evidence: `research/experiments/2026-06-23_cami3_toygut_coden15_formula_af/NOTE.md`.

CAMI3 GTDB-mapped ANI check: GTDB truth transferred from CAMI3 NCBI species taxids is highly ambiguous, so use both unique and possible-GTDB views. In selected-call ANI scoring, Sylph had no high-ANI calls outside possible GTDB truth; coden15 formula-AF averaged 3.0 such calls/sample but had higher unique-GTDB recall than old ctx-marker and coden15 original. Evidence: `research/experiments/2026-06-23_cami3_toygut_gtdb_ani_accuracy/NOTE.md`.

CAMI3 source/ref ANI check: when read-based ANI is compared directly with skani source-genome-to-reference ANI, MinCO coden15 `Ref_zip_aaf_ani` slightly outperforms Sylph `Adjusted_ANI` on Pearson and MAE for the paired-source selected-call subset. MinCO emitted/naive ANI is saturated at `1.0` and should not be used as a continuous ANI accuracy signal here. Evidence: `research/experiments/2026-06-23_cami3_source_ref_ani_accuracy/NOTE.md`.

## 2026-06-24 Mixed Mouse+CAMI3 Abundance Check

On the six-sample mixed panel used for universal tuning (`mouse_gtdb` samples
0-2 plus `cami3_ncbi` samples 0-2), the best MinCO abundance result found so
far still does not beat Sylph.

Best raw MinCO L1 on this mixed panel:

```text
callset: marker_l1_ctxobj_cami_blend
abundance score: raw value_sum
all F1: 0.770982
all L1: 17.823364
mouse L1: 2.989748
CAMI3 L1: 32.656979
```

Best F1-priority MinCO raw rule on this mixed panel:

```text
callset: marker_l1_ctxobj_f1_blend
abundance score: raw value_sum
all F1: 0.782475
all L1: 18.014617
```

Sylph on the same mixed panel:

```text
all F1: 0.783546
all L1: 15.038388
mouse L1: 2.287032
CAMI3 L1: 27.789743
```

A 120-model leave-one-sample-out row-feature calibration did not improve over
raw `value_sum`; best held-out L1 was `18.302036`. Treat row-level calibration
as a negative result. The next abundance direction should be read-level or
context-level EM-style mass assignment, not more independent per-row regression.
Evidence: `research/experiments/2026-06-24_abundance_feature_calibration/NOTE.md`.

## 2026-06-24 Full S2000 Marked-Unique Benchmark

The full S2000 marked-unique design is now implemented in `minco ani` as
`--readwise-dual-evidence`. It remains useful as an implementation
architecture, but it is not yet the best accuracy strategy.

Evidence:

- Full-index marker-size equivalence:
  `research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/NOTE.md`
- Full F1/abundance/ANI comparison:
  `research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark/NOTE.md`
- Integrated dual-counter implementation benchmark:
  `research/experiments/2026-06-24_integrated_dual_evidence_benchmark/NOTE.md`

The integrated dual-counter output matched the physical markerdb size map for
all comparable emitted rows across the six-sample benchmark (`0` comparable
mismatches). The previous proxy result was therefore confirmed by the real
implementation.

Six-sample mouse+CAMI3 comparison:

```text
strategy                                           all_F1    all_L1
Sylph                                              0.783546  15.038388
MinCO current mixed F1-priority                    0.782475  18.014617
MinCO integrated aggregate EM proxy                0.775173  18.061240
MinCO current mixed L1-priority                    0.770982  17.823364
MinCO integrated full S2000 dual evidence XnY700   0.750078  20.478490
MinCO integrated marker-only                       0.747536  20.513677
```

Interpretation: the integrated full S2000 marked-unique path slightly improves old
ctx-marker inside MinCO (`+0.002543` F1 and `-0.035187` L1 pp), driven by Toy
Mouse, but it does not beat the current mixed MinCO rows or Sylph. The tested
aggregate EM proxy is also not sufficient: unrestricted full add-back creates
too many remote shared-context FPs, while marker-present-only redistribution
collapses to marker-only presence and worsens L1 relative to the old physical
ctx+obj markerdb blend.

ANI remains a separate reporting layer. On the CAMI3 source/ref ANI subset,
MinCO coden15 `Ref_zip_aaf_ani` slightly beat Sylph adjusted ANI by MAE
(`0.007627` versus `0.008038`). Keep the coden15 ANI reporter, not old emitted
naive ANI.

Practical status: use integrated dual evidence for future full/marker
optimization because it avoids the physical markerdb copy. Do not expect this
architecture alone to beat Sylph on F1 or abundance; the next work should
change the abundance mechanism. If EM is pursued, it should be context-level or
read-level EM over actual ambiguity groups, not aggregate per-reference full
totals.

## 2026-06-24 Context-Level Ambiguity EM Prototype

A first true read/context-level ambiguity edge dump was implemented and tested
on full Toy Mouse samples. The edge dump worked. Raw EM over ambiguous
full-S2000 context groups did not improve abundance. A stricter filtered
follow-up produced a tiny sample0-only L1 gain, but failed the mouse0-2
validation as a fixed global blend.

Mouse0 per-read profile result:

```text
method                              F1        L1
marker_l1_per_read_profile          0.923077  1.800456
raw best edge-EM blend              0.923077  1.801188
filtered selected_group_species2    0.923077  1.792004
edge-only selected EM               0.923077  42.135253
```

Mouse0-2 validation:

```text
method                              mean_F1   mean_L1
marker_l1_per_read_profile          0.934155  3.036063
filtered_group_species2_beta0.0001  0.934155  3.036827
filtered_group_species2_beta0.0045  0.934155  3.095701
```

Interpretation: unfiltered ambiguous context mass is not abundance-proportional
on these samples. The sample0 filtered gain does not transfer: the best
three-sample mean is beta `0`, meaning no edge-EM abundance blend. Keep edge
dumps for diagnostics only; the current default abundance path remains the
marker/robust-depth strategy, not context-level EM.

Evidence:

- `research/experiments/2026-06-24_context_level_ambiguity_em_prototype/NOTE.md`
