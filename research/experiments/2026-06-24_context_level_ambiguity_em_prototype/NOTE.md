# Context-Level Ambiguity EM Prototype

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: `b55b4d3-dirty`
- Working tree: dirty, with prior project changes and the new edge-dump patch
- Changed implementation file for this experiment: `minco_core/src/command_ani.c`
- New scorer: `score_mouse0_edge_em.py`
- Provenance snapshots:
  - `provenance/code_status.txt`
  - `provenance/code_diff.patch`
  - `provenance/code_diffstat.txt`

## Question

Can true read/context-level ambiguity groups improve abundance estimation over
the aggregate per-reference full-evidence add-back that failed in the integrated
dual-evidence benchmark?

## Hypothesis

The previous aggregate EM proxy failed because it only had per-reference totals.
If MinCO emits per-read/per-context candidate groups, each shared context can
contribute total mass `1` across its candidate refs instead of being counted
once per ref. This should reduce shared-context overcounting and might improve
abundance L1.

## Implementation

Added an experimental environment-controlled ambiguity edge dump to MinCO
readwise mode. It does not add a public CLI option.

Environment variables:

```text
MINCO_READWISE_EDGE_OUT
MINCO_READWISE_EDGE_MAX
MINCO_READWISE_EDGE_MAX_CANDIDATES
MINCO_READWISE_EDGE_AMBIGUOUS_ONLY
MINCO_READWISE_EDGE_SELECTED_ONLY
```

The emitted table has one row per candidate reference edge in a retained
read/context group:

```text
edge_id
read_id
unit_id
qctx
edge_rank
ref_begin
gid
diff
best_diff
candidate_refs
selected_refs
selected_by_mode
cov_inc
```

For this first prototype, only groups with `candidate_refs > 1` were emitted,
and groups with more than `64` candidates were skipped. The output was capped
at `5,000,000` edge rows.

## Dataset

One full CAMI II Toy Mouse Gut sample was tested:

```text
sample: mouse sample0
reads: /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
refdb: /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
truth: existing GTDB species truth from evaluate_abundance_estimators.load_truth_profile(0)
```

The MinCO profile was rerun with `--density-block-ctx 0` so each event retained
a read-level `read_id`. This differs from the earlier density-blocked
six-sample benchmark, so comparisons here are within this experiment unless
otherwise stated.

## Methods

Commands are fully recorded in `commands.sh`; parameters are recorded in
`parameters.tsv`.

Readwise profile and edge dump:

```text
minco ani -p4
  -r /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
  --qraw mouse_sample0_reads
  --query-density ref
  --abundance-est depth
  --readwise-profile-only
  --readwise-dual-evidence
  --readwise-assign best-diff-split
  --readwise-ani naive
  --readwise-ctx-filter product-topfrac-median
  --readwise-fake-threshold 0.25
  --density-block-ctx 0
  -m0 -f0 -n0 -t0
```

EM scorer:

1. Score `MARKER_L1` from the marker columns of the same per-read profile.
2. Map edge `gid` values to GTDB species using `minco.stat` sample names and
   the GTDB metadata loader.
3. Collapse candidate refs to species per `(unit_id, qctx, species)`.
4. Keep only species already present in the marker callset. This makes this
   first prototype an abundance-only test and avoids the FP explosion seen with
   unrestricted full evidence.
5. Run EM on retained ambiguity groups:

```text
weight(species | event) proportional to
  abundance(species) * exp(-alpha * (diff - best_diff))
```

6. Test two edge modes:

```text
selected: selected_by_mode == 1 only
allcand: all candidates in retained groups, with diff likelihood
```

7. Blend normalized marker abundance with normalized edge-EM abundance:

```text
prediction = (1 - beta) * marker_abundance + beta * edge_em_abundance
```

Tested values:

```text
alpha = 0, 0.5, 1, 2, 4
beta = 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.05, 0.1
```

## Validation

Smoke test on 1,000 reads:

```text
edge rows: 181
event groups: 28
selected edges: 71
runtime: 0:27.94 wall
max RSS: 4.91 GB
```

Full mouse0 edge dump:

```text
reads processed: 33,170,320
edge rows written: 4,952,131
edge cap: 5,000,000, not hit
event groups: 671,619
selected edges: 1,660,122
runtime: 1:53.81 wall
max RSS: 6.20 GB
edge TSV size: 356 MB
profile TSV size: 103 MB
```

After mapping and marker filtering:

```text
mapped edge rows: 4,949,196
marker-filtered edge rows: 546,908
marker-filtered groups: 434,448
species-level edge rows: 546,908
species-level groups: 434,448
marker edge species: 65
```

Scorer runtime:

```text
runtime: 1:01.09 wall
max RSS: 2.47 GB
stderr: empty
```

## Results: Simple Ambiguous EM

Initial unfiltered top rows by abundance L1:

```text
method                                  F1        L1 pp     Pearson
marker_l1_per_read_profile              0.923077  1.800456  0.999959
blend_selected_alpha0_beta0.0001        0.923077  1.801188  0.999959
blend_selected_alpha0.5_beta0.0001      0.923077  1.801188  0.999959
blend_selected_alpha1_beta0.0001        0.923077  1.801188  0.999959
blend_allcand_alpha0.5_beta0.0001       0.923077  1.801281  0.999959
```

The marker baseline from this per-read profile was best. Every tested edge-EM
blend worsened L1, even at `beta=0.0001`.

Edge-only abundance was very poor:

```text
method                         F1        L1 pp     Pearson
edge_only_selected_alpha0       0.923077  42.1353   0.9248
edge_only_allcand_alpha0.5      0.923077  44.3925   0.9203
edge_only_allcand_alpha0        0.923077  44.9090   0.9193
```

The presence F1 is unchanged in this prototype because the edge EM was
restricted to marker-supported species. The test is therefore a direct
abundance-signal test, not a rescue-call test.

## Results: Filtered Edge-EM Follow-up

The follow-up search reused the cached marker-filtered edge table and tested
118 filters over selected/all candidate edges. The main added diagnostics were
context event depth, `best_diff * qctx_depth`, and the number of
marker-supported species in each ambiguity group.

The best broad-grid row required selected edges from groups containing at least
two marker-supported species:

```text
filter: selected_group_species2
rows: 115,290
groups: 89,965
species: 58
alpha: 0
weight: event
beta: 0.003
F1: 0.923077
L1: 1.793005 pp
Pearson: 0.999961
```

A local beta refinement of the same filter found a slightly better blend:

```text
filter: selected_group_species2
alpha: 0
weight: event
beta: 0.0045
F1: 0.923077
L1: 1.792004 pp
Pearson: 0.999961
```

This is a small improvement over the same-profile marker baseline:

```text
marker baseline L1: 1.800456 pp
filtered edge-EM L1: 1.792004 pp
absolute gain: 0.008452 pp
relative L1 gain: 0.47%
```

Edge-table diagnostics:

```text
unique qctx count: 22,118
event depth median: 2
event depth q90/q99/max: 13 / 479 / 3030
product median/q90/q99/max: 0 / 1 / 5.83 / 2480
```

## Results: Toy Mouse 0-2 Validation

The proposed next step was run on Toy Mouse samples 1 and 2. Sample1 initially
hit the old `5,000,000` edge cap, so it was rerun with
`MINCO_READWISE_EDGE_MAX=20000000`. Sample2 also used the larger cap. Neither
larger-cap run was truncated.

MinCO edge/profile run details:

```text
sample  reads       edge rows  profile rows  wall time  max RSS
0       33170320    4952131    118722        1:53.81    6.20 GB
1       33270054    6453967    153588        1:57.39    6.53 GB
2       33165928    6776571    155309        1:55.51    6.50 GB
```

Marker-filtered edge diagnostics:

```text
sample  marker species  truth species  filtered edge rows  selected_group_species2 rows
0       65              65             546908              115290
1       75              80             556997              67141
2       68              73             716192              181500
```

Across Toy Mouse samples 0-2, the best mean L1 was the marker baseline
(`beta=0`). Every positive fixed beta worsened mean abundance L1:

```text
method                              mean F1   mean L1 pp  mean Pearson
marker_l1_per_read_profile          0.934155  3.036063    0.999953
filtered_group_species2_beta0.0001  0.934155  3.036827    0.999953
filtered_group_species2_beta0.003   0.934155  3.069947    0.999944
filtered_group_species2_beta0.0045  0.934155  3.095701    0.999938
filtered_group_species2_beta0.006   0.934155  3.126429    0.999931
```

Per-sample best rows:

```text
sample  marker L1 pp  best filtered beta  best filtered L1 pp  delta vs marker
0       1.800456      0.0045              1.792004             -0.008452
1       1.628925      none                1.628925              0.000000
2       5.678809      0.02                5.583216             -0.095593
```

The fixed beta that looked best on sample0 does not transfer. Sample1 is harmed
substantially by the same filtered ambiguous mass, while sample2 improves only
at a much larger beta. This means a global fixed small-blend rule is not a
default candidate.

## Interpretation

The true ambiguity-group dump worked, and it did prevent double-counting within
each retained group. However, raw ambiguous-context mass is not proportional to
true abundance on mouse sample0. It is enriched for shared, repetitive, or
related species contexts, so adding it shifts abundance away from the already
strong marker-depth estimate.

This falsifies the simplest version of the idea:

```text
marker abundance + small amount of context-level ambiguous EM mass
```

It does not falsify all possible EM strategies. It shows that EM over only the
ambiguous full-S2000 groups is not useful unless the ambiguous events are
filtered or modeled more carefully.

The filtered follow-up shows that a very small amount of carefully restricted
ambiguous mass can slightly improve mouse0 abundance, but the three-sample
validation rejects it as a global default. Context-level EM should remain a
diagnostic branch, but adoption requires an adaptive rule that can decide when
ambiguous mass is helpful.

## Conclusion

For Toy Mouse sample0, simple unfiltered context-level ambiguity EM did not
improve MinCO abundance. A filtered version gave a tiny abundance-L1 gain while
leaving presence F1 unchanged. However, the follow-up on Toy Mouse samples 0-2
showed that the best fixed setting is no edge-EM blend at all.

Same-profile marker-only baseline:

```text
F1=0.923077
L1=1.800456 percentage points
```

Best filtered edge-EM diagnostic row:

```text
F1=0.923077
L1=1.792004 percentage points
```

Decision: do not promote this EM prototype. Keep the edge-dump machinery for
diagnostics. The `selected_group_species2` small-blend rule is rejected as a
global default because it worsens mean Toy Mouse abundance L1.

## Caveats

- Filter search was broad only on sample0; samples 1-2 tested the best sample0
  filter and beta neighborhood, not all 118 filters.
- The profile used `--density-block-ctx 0`, so it is not numerically identical
  to previous density-blocked six-sample rows.
- Candidate groups with more than `64` refs were skipped.
- Edge EM was restricted to marker-supported species, so it cannot rescue
  missing true species in this prototype.
- The edge table contains only ambiguous groups (`candidate_refs > 1`), not all
  unique context events.

## Next Experiment

If continuing this direction:

```text
1. do not expand the fixed selected_group_species2 beta rule to CAMI3;
2. inspect why sample1 is harmed while sample2 benefits, using edge-mass
   species deltas and diagnostics;
3. if a reliable adaptive trigger appears, retest on mouse0-2 before CAMI3;
4. otherwise return to marker-only robust abundance for the default path.
```
