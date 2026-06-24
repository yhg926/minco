# Minco Readwise Correction Research

Date: 2026-06-20
Author/agent: Codex
Project: minco
Code checkpoint: `b55b4d3` plus this experiment script

## Question

Why does minco S10000 miss low-abundance CAMI marine sample0 species that
Sylph detects, and which correction direction is most promising?

Directions tested:

1. Coverage-corrected context AAF, including ZIP and overdispersed ZINB
   corrections.
2. Sequencing-error/distance inflation correction using Sylph adjusted ANI as
   pseudo-truth.
3. A learned feature model using raw ANI, AAF, MoE count features, support, and
   depth features.
4. Shared-context coverage reassignment, evaluated by code-path review because
   the current TSV output does not retain per-query-context hit lists.

## Dataset

- Query reads: `/tmp/cami_marine_sample0_reads.fq.gz`
- Gold profile: `/tmp/gs_marine_short.profile`
- Gold sample ID: `marmgCAMI2_short_read_sample_0`
- minco S10000 unfiltered detail table:
  `/tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_unfiltered.tsv`
- Sylph profile:
  `/tmp/sylph_marine_sample0/profile.tsv`
- Taxmap:
  `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`

The scoring path uses the same OPAL-like species filter as the earlier CAMI
notes: species-level rows for sample0, dropping `unidentified`,
`unidentified plasmid`, and `unidentified virus`. This gives 256 gold species.

## Methods

The analysis script is:

```text
research/experiments/2026-06-20_minco_readwise_correction_research/scripts/analyze_readwise_corrections.py
```

It joined minco rows, Sylph rows, taxmap species IDs, and CAMI gold species.
Large generated tables are under:

```text
/tmp/minco_readwise_correction_20260620
```

Context AAF ANI used the current minco context `k = Bitslen.ctx / 2 = 11`:

```text
AAF_ANI = 1 + log(AF) / 11
```

Coverage corrections tested:

- Raw AAF: `AF = Ref_breadth`.
- ZIP: solve `observed_breadth = pi * (1 - exp(-lambda))` and
  `Ref_mean_depth = pi * lambda` for latent `pi`.
- ZINB: same latent `pi`, but present contexts follow a negative-binomial
  distribution fitted by the observed mean, variance, and zero fraction.
- Oracle Sylph effective coverage: `pi = Ref_breadth / (1 - exp(-Sylph_Eff_cov))`.
  This is not deployable directly, but tests whether a correct effective
  coverage estimate would rescue the low-abundance taxa.

Learned correction tests:

- A row-level HGB classifier for whether a candidate row's species taxid is in
  the CAMI gold set. This is diagnostic only because it trains and tests on one
  sample.
- Ridge and HGB regressors from minco features to Sylph adjusted ANI on 288
  matched Sylph/minco accessions.
- A same-sample HGB Sylph-calibrated ANI column applied back to all minco rows.

## Results

Input sanity checks:

```text
minco rows:                         200,709
minco rows mapped to species:       199,923
gold species:                       256
Sylph mapped species predictions:   278
Sylph TP not in minco S10000 base:  25
```

Main species-level strategy results are in `summary.tsv` and
`/tmp/minco_readwise_correction_20260620/strategy_summary.tsv`.

Key rows:

```text
strategy                                  TP   FP    FN   F1      recovered 25 Sylph-only TP
Sylph default                            222   56    34   0.831   25
minco S10000 current opt                 205   100   51   0.731   0
raw ctx AAF best                          83   18   173   0.465   0
ZIP ctx AAF best                          99   18   157   0.531   0
ZINB ctx AAF best                        240  6178   16   0.072   24
oracle Sylph-effcov ctx AAF best         212   36    44   0.841   21
HGB Sylph-calibrated distance best       196   18    60   0.834   13
```

Interpretation by correction:

- Raw AAF is not enough. Low-abundance true species have small observed
  breadth, so raw AAF ANI is strongly underestimated.
- ZIP using minco `Ref_mean_depth` barely helps. The mean depth is inflated by
  uneven/shared hits, so the fitted latent AF remains close to raw breadth.
- ZINB sees the high variance/zero inflation and rescues 24 of 25 missed
  Sylph-only TPs, but it also rescues thousands of false positives. This is a
  useful diagnostic but not a usable filter alone.
- The oracle Sylph effective-coverage AAF is the strongest signal: F1 0.841,
  TP 212, FP 36, FN 44. This beats current minco S10000 F1 0.731 and is close
  to local Sylph F1 0.831. It implies that a correct effective coverage estimate
  can rescue much of the low-abundance failure mode.
- Learned distance correction is promising but incomplete. HGB prediction of
  Sylph ANI reduced CV MAE from raw minco-vs-Sylph 0.0406 to 0.00371 on 288
  matched rows, but as a caller it recovered only 13 of the 25 missed Sylph TPs.
  It improves precision but does not solve low-abundance recall by itself.

## Model Metrics

From `/tmp/minco_readwise_correction_20260620/model_summary.tsv`:

```text
analysis                          rows     metric
row gold classifier holdout       199923   ROC-AUC 0.9253, AP 0.3159
raw minco ANI vs Sylph ANI        288      MAE 0.04055, Pearson 0.427
Ridge feature calibration CV      288      MAE 0.00533, Pearson 0.613
HGB feature calibration CV        288      MAE 0.00371, Pearson 0.813
```

The classifier result should not be treated as a deployable accuracy estimate:
it is a one-sample diagnostic with row labels derived from the same gold set.

## Shared-Context Reassignment Review

Current code path:

- Query density units are processed in
  `minco_core/src/command_ani.c::ani_process_density_ctxobj_unit`.
- For each query context, the code walks every matching reference `ctxgid`.
- `ref_ctx_cov` is a packed `uint32_t`: 28 bits coverage plus 4 bits best
  object-difference code (`ANI_REFCOV_DIFF_SHIFT = 28`).
- `ani_ref_covdiff_add_hit()` increments coverage for every matching reference
  context and keeps the minimum object difference for that reference context.
- `ani_unique_best_features_from_ref_covdiff()` later converts that coverage
  array into one best-diff feature row per reference context.

This means late hits with fewer object differences are already kept as the best
diff for that same reference context. However, the current code does not decide
which reference should receive a shared query context before adding coverage.
Every matched reference gets coverage. After the run, the TSV only has aggregate
per-reference summaries, so true shared-context reassignment cannot be tested
post hoc.

Recommended implementation:

1. Inside each query-context run in `ani_process_density_ctxobj_unit`, first
   collect all matching reference `ctxgid` candidates and their minimum object
   difference.
2. Find the best object difference across candidates for that query context.
3. Add coverage only to candidates with the best difference, or split coverage
   among tied best candidates.
4. Prefer species-level or accession-cluster splitting when many references of
   the same species tie, otherwise strain-dense databases will overcount shared
   contexts.
5. Store coverage as fixed-point in the existing 28 coverage bits, for example
   Q8 or Q10, while keeping the 4-bit best-diff code. This allows fractional
   split coverage without doubling the `ref_ctx_cov` memory footprint.
6. Keep unweighted support counters separately from reassigned coverage so ANI
   and abundance can be calibrated independently.

This should be tested before installing a new default because it changes both
abundance and ANI features.

## Conclusion

The best answer from this experiment is that minco's low-abundance misses are
primarily an effective-coverage and shared-coverage assignment problem, not
only a raw distance-model problem.

The next implementation target should be shared-context best-diff reassignment
plus an effective-coverage estimator. The oracle Sylph-effcov AAF result shows
the target behavior is plausible. ZINB alone is too permissive, and the learned
distance model alone is too conservative.

## Caveats

- This is one CAMI sample, not a general benchmark.
- S10000 used GTDB-only minco references, while the taxmap includes virus rows
  and the earlier S1000 run used a plus-virus database.
- Sylph adjusted ANI is used as pseudo-truth for the learned distance model,
  not as independent ANIm ground truth.
- The oracle effective-coverage result uses Sylph `Eff_cov`, so it is an upper
  bound, not a deployable minco method.
- The Python run emitted a SciPy binary-size warning during sklearn import, but
  the analysis completed and reproduced expected counts.

## Next Experiment

Implement a guarded prototype of shared-context best-diff reassignment in the
readwise density path, emit both old and reassigned coverage/ANI columns, and
rerun CAMI marine sample0 S10000. Then train an effective-coverage model on
multiple CAMI samples, not only sample0.
