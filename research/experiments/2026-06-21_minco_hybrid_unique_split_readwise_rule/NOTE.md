# Minco hybrid unique split readwise rule

Date: 2026-06-21
Author/agent: Codex
Project: minco
Code commit: b55b4d3 plus dirty worktree changes in readwise split/ZIP experiment code
KSSD3A binary/tool version: minco 0.1

## Question

Can a hybrid readwise profiling rule use `best-diff-unique` features as the
precision anchor and `best-diff-split` features as a recall rescue, improving
species detection on both CAMI3 toy human-gut and CAMI2 marine sample0?

## Hypothesis

Split mode recovers true taxa missed by unique mode, but also adds shared-context
false positives. A conservative hybrid rule should keep unique mode precision
while rescuing only split-supported taxa with enough unique confirmation.

## Dataset

- Toy input reads:
  `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz`
- Toy gold profile:
  `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt`
- Marine input reads:
  `/tmp/cami_marine_sample0_reads.fq.gz`
- Marine gold profile:
  `/tmp/gs_marine_short.profile`
- Reference sketch:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
- Tax map:
  `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Samples: 2 benchmark samples, CAMI3 toy sample0 and CAMI2 marine sample0.
- Selection criteria: GTDB-only species-level comparison. Toy is restricted to
  bacterial gold taxa because the local Sylph database is GTDB-only.
- Storage location: repo note files plus large temporary outputs under
  `/tmp/minco_hybrid_20260621`.

## Methods

Key parameters:

```text
Candidate tables:
  minco ani --qraw READS --query-density ref --abundance-est depth
    --readwise-profile-only --readwise-assign best-diff-unique --readwise-ani zip-aaf
    -m0 -f0 -n0 -t0

  minco ani --qraw READS --query-density ref --abundance-est depth
    --readwise-profile-only --readwise-assign best-diff-split --readwise-ani zip-aaf
    -m0 -f0 -n0 -t0

Hybrid script:
  scripts/hybrid_unique_split_score.py

Taxid aggregation:
  one best row per species taxid per mode, sorted by ANI, Ref_breadth,
  XnY_ctx, then Ref_mean_depth.

Baseline rules:
  unique_current/split_current:
    XnY_ctx >= 10, ANI >= 0.94, Real_min_align_fraction >= 0.05

  unique_relaxed/split_relaxed:
    XnY_ctx >= 1, ANI >= 0.90, Ref_breadth >= 0.005

Hybrid rule family:
  primary unique rule OR split rescue rule.
  Joint-best parameters:
    u_xny=10, u_ani=0.90, u_breadth=0.005
    s_xny=25, s_ani=0.94, s_breadth=0.05
    u_confirm_xny=1

Classifier test:
  Logistic regression, random forest, and HistGradientBoosting on joined
  unique/split features, evaluated by leave-one-dataset-out.
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

Important comparison caveat: this experiment uses a conservative best-row
species-taxid aggregation so unique and split features can be joined. It is not
the same scoring path as the previous direct CLI marine benchmark. The earlier
direct CLI `S=1000 best-diff-unique + ZIP AAF` result still beat local Sylph on
CAMI marine sample0 (`F1 0.866` vs `0.831`) and across samples 0-2 (`mean F1
0.862` vs `0.832`). The table below evaluates whether a new hybrid unique/split
rule should replace that recipe.

Sanity check on the current unfiltered marine unique output with the original
row-level direct filter (`XnY_ctx>=10`, `ANI>=0.94`,
`Real_min_align_fraction>=0.05`) reproduced the previous result:

```text
current minco S1000 unique direct filter: TP 223 FP 36 FN 33 F1 0.866
Sylph default:                            TP 222 FP 56 FN 34 F1 0.831
```

```text
Toy GTDB bacteria:
  unique_relaxed       TP 57 FP 18 FN 51 F1 0.623
  split_current        TP 35 FP 10 FN 73 F1 0.458
  hybrid_joint_best    TP 47 FP 12 FN 61 F1 0.563
  Sylph default        TP 55 FP 13 FN 53 F1 0.625

Marine GTDB species:
  unique_current       TP 197 FP 34 FN 59 F1 0.809
  split_current        TP 213 FP 51 FN 43 F1 0.819
  hybrid_joint_best    TP 221 FP 64 FN 35 F1 0.817
  Sylph default        TP 222 FP 56 FN 34 F1 0.831

Posthoc per-dataset best hybrid rule:
  Toy GTDB bacteria:   TP 60 FP 20 FN 48 F1 0.638
  Marine GTDB species: TP 212 FP 40 FN 44 F1 0.835
  These are same-sample tuned results, so they are useful diagnostics but not
  valid default parameters.

Leave-one-dataset-out classifier:
  Train toy -> test marine:
    best test F1 among tested models: HGB 0.623, with 218 TP, 226 FP, 38 FN
  Train marine -> test toy:
    best test F1 among tested models: logistic 0.450, with 34 TP, 9 FP, 74 FN

Runtime for newly generated marine candidate tables:
  split unfiltered:    2:11 wall, 4.12 GB RSS
  unique unfiltered:   1:53 wall, 3.68 GB RSS
```

## Validation

- Checks run:
  - `minco ani` completed normally for marine unique and split candidate tables.
  - Hybrid script completed and wrote joined features, grid, best params, and summary.
  - Classifier script completed; both scripts pass `python3 -m py_compile`.
  - The scorer was corrected after detecting invalid max-per-column aggregation.
- Expected behavior:
  - Split mode should improve recall but risk more false positives.
  - Hybrid rule should reduce split FPs if unique confirmation is informative.
- Observed behavior:
  - Split improved strict/current recall on both datasets.
  - The tested hybrid family did not beat the best unique rule on toy and did
    not beat Sylph on either toy bacteria or marine sample0.
  - Same-sample tuning can beat the fixed baselines, but the winning thresholds
    differ enough that they should not be promoted to defaults.
- Failure modes checked:
  - Avoided mixing features from different rows of the same taxon by using one
    best row per taxid per mode.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, matrices, trees, figures, and sketches.

## Conclusion

The simple unique-primary/split-rescue threshold family is not an improvement
over the previous direct CLI marine recipe. It confirms that split mode has
useful recall information, but the current hybrid rules do not separate split
true positives from split false positives well enough. A quick classifier test
also failed to generalize across toy and marine when trained on only one
dataset. The next improvement should gather more samples and train a
domain-aware calibrated effective-coverage/call model, not install a model from
this two-sample experiment.

## Paper-Relevant Claim

No supported paper claim from this experiment. This weakens the hypothesis that
a small hand-written hybrid threshold rule can robustly combine unique and split
readwise ZIP features.

## Caveats

- Only CAMI3 toy sample0 and CAMI2 marine sample0 were used in this first
  hybrid grid.
- Sylph and minco use different GTDB releases.
- The scorer is species-taxid based and local, not official OPAL.
- Toy all-species scoring is not fair to Sylph because the local Sylph database
  is GTDB-only; therefore toy comparison uses bacterial taxa.
- Best-row taxid aggregation is conservative and does not exactly reproduce
  direct minco row-level filtering, but it avoids impossible feature mixtures.

## Next Experiment

Train a calibrated classifier using row-level unique and split features with
labels from multiple CAMI samples. Include domain indicators, abundance features,
ZIP features, naive ANI, and unique-vs-split deltas. Evaluate by holding out
entire samples, not by tuning and testing on the same sample.
