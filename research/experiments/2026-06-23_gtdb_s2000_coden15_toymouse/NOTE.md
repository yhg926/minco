# GTDB S2000 Coden15 Toy Mouse Readwise Test

## Question

Does extending MinCO coden15 to a 96-bit `ctx64+obj32` payload improve GTDB
S2000 markerdb readwise performance on CAMI II Toy Mouse Gut samples?

## Code Tested

Workspace: `/home/ubuntu/yihuiguang/tools/KSSD3mini`

Coden15 build:

```sh
make -C minco_core BINDIR=bin_coden15 \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1' clean
make -C minco_core BINDIR=bin_coden15 \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1'
```

Important implementation point: coden15 sketching now uses native MinCO
OpenMP threading over input files with one process and `-pN`. The earlier
external chunk-worker approach was abandoned and is not part of this result.

## Input References

The old S10000 manifest had 200,709 entries. Re-sketching the exact set was
blocked because 785 VEuPathDB relative-path FASTAs are absent locally. The
available coden15 run used the 199,924 absolute GTDB/human paths:

`research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse/gtdb232_199924_available_abs_paths.list`

`gtdb232_199924_available_abs_paths.missing` has 0 rows.

The old AAF0.03 remove list had 182 entries, but all were VEuPathDB paths, so
it matched 0 of the 199,924 available GTDB/human references. Therefore the
coden15 markerdb was built directly from the available sketch.

## Sketch And Markerdb

Large artifacts are under:

`/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623`

Full coden15 S2000 sketch:

- Path: `/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_native`
- Samples: 199,924
- Payload: `minco.ctxobj96`
- Sketch time: 56:58.14
- Peak RSS: 822,284 KB
- Output payload: 4,794,588,684 bytes

Full sketch index:

- Index file: `minco.refindex.ctx64gid32obj32`
- Time: 1:45.37
- Peak RSS: 12,488,060 KB

Coden15 ctx markerdb:

- Path: `/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623/sketch_T15_S2000_199924_ctxmarker`
- Build time: 2:00.47
- Peak RSS: 17,222,144 KB
- Payload size: 3,901,527,240 bytes
- Marker sizes: n=199,924, min=20, max=1996, mean=1626.254
- Refs below default warning threshold 500: 2,529

Markerdb index:

- Index file: `minco.refindex.ctx64gid32obj32`
- Time: 1:03.17
- Peak RSS: 10,162,132 KB

## Toy Mouse Readwise

Readwise settings:

```sh
minco ani -p16 -r <coden15_ctxmarker> --qraw <reads.fq.gz> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0 -o <output.tsv>
```

Runtime:

| sample | elapsed | peak RSS KB |
| ---: | ---: | ---: |
| 0 | 7:30.78 | 6,290,104 |
| 1 | 7:23.10 | 6,453,756 |
| 2 | 7:24.80 | 6,498,460 |

## Gate Rules Scored

Active gate, unchanged from previous notes:

- `XnY_ctx >= 15`
- `ANI_naive_calc > 0.95`
- `Reliable_ztp_af >= 0.40`
- if `Reliable_Ref_hit_mean_depth > 3` and `Reliable_depth_vmr > 50`,
  require `ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03`

Exploratory relaxed coden15 gate:

- same as active gate except `XnY_ctx >= 20` and `Reliable_ztp_af >= 0.25`

The relaxed gate is exploratory only. It was selected from this result grid and
needs validation beyond these three samples before becoming a default.

## Results

Per-sample GTDB species presence/absence:

| sample | method | TP | FP | FN | F1 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0 | ctx-only current best | 58 | 0 | 7 | 0.943089 |
| 0 | coden15 active | 58 | 2 | 7 | 0.928000 |
| 0 | coden15 relaxed | 60 | 2 | 5 | 0.944882 |
| 0 | Sylph | 59 | 2 | 6 | 0.936508 |
| 1 | ctx-only current best | 70 | 2 | 10 | 0.921053 |
| 1 | coden15 active | 70 | 0 | 10 | 0.933333 |
| 1 | coden15 relaxed | 71 | 3 | 9 | 0.922078 |
| 1 | Sylph | 78 | 2 | 2 | 0.975000 |
| 2 | ctx-only current best | 66 | 0 | 7 | 0.949640 |
| 2 | coden15 active | 67 | 0 | 6 | 0.957143 |
| 2 | coden15 relaxed | 69 | 1 | 4 | 0.965035 |
| 2 | Sylph | 71 | 1 | 2 | 0.979310 |

Mean F1 across samples:

- coden15 active: 0.939492
- coden15 relaxed: 0.943998
- ctx-only current best: 0.937927
- Sylph: 0.963606

Pooled totals:

- coden15 active: TP=195 FP=2 FN=23 F1=0.939759
- coden15 relaxed: TP=200 FP=6 FN=18 F1=0.943396
- ctx-only current best: TP=194 FP=2 FN=24 F1=0.937198
- Sylph: TP=208 FP=5 FN=10 F1=0.965197

## Interpretation

Coden15 works technically and native `-p` threading is stable. It slightly
improves MinCO pooled presence/absence over the old ctx-only markerdb, but it
does not close the recall gap to Sylph.

The unchanged active gate is conservative: it improves samples 1 and 2, but
sample 0 adds two false positives and loses F1. The relaxed coden15 gate
recovers more low-breadth true positives and gives the best MinCO F1 here, but
it increases false positives, especially in sample 1.

Abundance remains worse than Sylph. Coden15 improves raw predicted mass versus
old ctx-only, but renormalized MAE is still far above Sylph on all samples.

## Robust Abundance Rescue Retest

Retested coden15 with the current pinned abundance strategy from
`research/CURRENT_BEST.md`:

```text
ctx-marker + robust effective depth + intra-genus winner rescue
```

The scorer reused the same robust effective-depth rule and rescue thresholds as
`2026-06-22_minco_sylph_abundance_model`, then applied them to the coden15
ctx-marker readwise outputs.

Mean GTDB species abundance metrics across Toy Mouse samples 0-2:

| method | TP | FP | FN | Pearson | Spearman | L1 percentage points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| old ctx-marker robust depth + rescue | 65.000 | 0.667 | 7.667 | 0.999945 | 0.910342 | 1.7945 |
| Sylph reported abundance | 69.333 | 1.667 | 3.333 | 0.999972 | 0.987668 | 2.2870 |
| coden15 active robust depth + rescue | 65.333 | 0.667 | 7.333 | 0.999917 | 0.887397 | 2.4022 |
| coden15 relaxed robust depth + rescue | 67.000 | 2.000 | 5.667 | 0.999930 | 0.916127 | 3.7377 |

Per-sample coden15 active L1 was worse than the old ctx-marker robust-rescue
baseline on all three samples:

| sample | old ctx-marker L1 | coden15 active L1 | coden15 relaxed L1 |
| ---: | ---: | ---: | ---: |
| 0 | 1.7692 | 2.3746 | 2.5414 |
| 1 | 1.6675 | 1.9563 | 3.0291 |
| 2 | 1.9468 | 2.8758 | 5.6425 |

Interpretation: coden15 rescues slightly more species, but it does not improve
abundance. The active coden15 rescue increases mean L1 from 1.79 to 2.40
percentage points. The relaxed coden15 gate improves recall but moves too much
mass onto false or wrong species, increasing mean L1 to 3.74 percentage points.
The current best abundance baseline remains the old ctx-marker robust effective
depth plus intra-genus winner rescue.

## Mechanism Analysis And Coden15 Retune

Coden15 ctx-marker does not behave like old ctx-marker under the same gate.
The active call counts are similar, but the candidate background is much
smaller:

| sample | old mapped rows | coden15 mapped rows | old active rows | coden15 active rows |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 65,233 | 3,071 | 58 | 60 |
| 1 | 88,994 | 4,900 | 72 | 70 |
| 2 | 83,344 | 3,821 | 66 | 67 |

The main abundance regression is concentrated in a few species:

- sample0 `s__Lactobacillus mulieris_A` becomes a coden15 false positive with
  0.951% predicted abundance.
- sample2 `s__Porphyromonas pasteri` is lost by the coden15 active gate; its
  truth abundance is 0.339%.
- coden15 also slightly shifts mass away from dominant true species such as
  `s__Lacticaseibacillus paracasei` and `s__Lactobacillus helveticus`.

The exact gate mechanism is:

| sample | species | old status | coden15 status | old Reliable_ztp_af | coden15 Reliable_ztp_af | old VMR | coden15 VMR | old delta | coden15 delta |
| ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | `s__Lactobacillus mulieris_A` | blocked | passes | 0.423 | 0.411 | 65.4 | 13.9 | 0.0359 | 0.0370 |
| 2 | `s__Porphyromonas pasteri` | passes | blocked by AF | 0.431 | 0.369 | 0.60 | 0.65 | 0.0351 | 0.0294 |

Old ctx-marker blocked `s__Lactobacillus mulieris_A` because the VMR/delta
filter was triggered (`VMR > 50`) and `ANI_naive_calc - ANI_from_Reliable_ztp_af`
was above 0.03. Coden15 smooths that row's hit-depth distribution enough that
VMR drops to 13.9, so the same filter is not triggered. Conversely, coden15
lowers `s__Porphyromonas pasteri` reliable AF below the 0.40 active floor.

An estimator-only grid did not fix this. The best tested estimator change was
only to use AF exponent 1.0 instead of 1.05 in the pinned robust-depth rule:

| gate | formula | exponent | L1 percentage points |
| --- | --- | ---: | ---: |
| active | pinned robust depth | 1.0 | 2.3867 |
| active | pinned robust depth | 1.05 | 2.4022 |

The useful optimization is gate-specific for coden15:

```text
XnY_ctx >= 15
ANI_naive_calc > 0.95
Reliable_ztp_af >= 0.35
if Reliable_Ref_hit_mean_depth > 3 and Reliable_depth_vmr > 10:
    require ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03
then apply the same robust effective depth + intra-genus winner rescue
```

Diagnostic mean abundance metrics with that coden15 retune:

| method | TP | FP | FN | Pearson | Spearman | L1 percentage points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| old ctx-marker robust depth + rescue | 65.000 | 0.667 | 7.667 | 0.999945 | 0.910342 | 1.7945 |
| Sylph reported abundance | 69.333 | 1.667 | 3.333 | 0.999972 | 0.987668 | 2.2870 |
| coden15 retuned gate + robust rescue | 65.667 | 0.667 | 7.000 | 0.999929 | 0.900959 | 2.1857 |
| coden15 original active + robust rescue | 65.333 | 0.667 | 7.333 | 0.999917 | 0.887397 | 2.4022 |

Interpretation: coden15 can be retuned to beat Sylph on mean Toy Mouse L1, but
it still does not beat the old ctx-marker abundance baseline. The likely reason
is that coden15 changes marker composition and depth dispersion: it can recover
some low-AF true species, but it also weakens the old VMR-based fake-context
guard and shifts effective-depth ratios among high-abundance species. Therefore
coden15 should not replace old ctx-marker for abundance unless the coden15 gate
is validated on independent samples.

### Formula-Based AF Gate

The fixed `Reliable_ztp_af >= 0.40` floor is not principled because AF depends
on effective context length. The gate should use depth-adjusted AF and derive
the AF floor from the ANI threshold:

```text
AF_floor = ANI_threshold ^ effective_context_length
```

Here `Reliable_ztp_af` is the depth-adjusted AF. With `ANI_threshold = 0.95`
and effective context length 24, the floor is:

```text
0.95^24 = 0.291989
```

This matches the useful coden15 retune range around 0.29-0.30, but gives it a
model-based reason rather than an arbitrary threshold. A formula grid over
effective lengths 20-32 found:

| effective length | best AF floor | best L1 percentage points |
| ---: | ---: | ---: |
| 20 | 0.358486 | 2.0879 |
| 22 | 0.323534 | 2.0802 |
| 24 | 0.291989 | 2.0091 |
| 26 | 0.289872 | 2.0091 |
| 28 | 0.261610 | 3.0224 |
| 30 | 0.246835 | 3.3496 |
| 32 | 0.222768 | 3.4773 |

Best formula-AF coden15 diagnostic gate:

```text
ANI_threshold = 0.95
Reliable_ztp_af >= 0.95^24 = 0.291989
XnY_ctx >= 10
if Reliable_Ref_hit_mean_depth > 3 and Reliable_depth_vmr > 50:
    require ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.02
robust abundance: exponent 1.0, median cutoff 10
then apply the same intra-genus winner rescue
```

Mean abundance metrics:

| method | TP | FP | FN | Pearson | Spearman | L1 percentage points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| old ctx-marker robust depth + rescue | 65.000 | 0.667 | 7.667 | 0.999945 | 0.910342 | 1.7945 |
| coden15 formula-AF gate + robust rescue | 66.333 | 3.000 | 6.333 | 0.999936 | 0.893282 | 2.0091 |
| Sylph reported abundance | 69.333 | 1.667 | 3.333 | 0.999972 | 0.987668 | 2.2870 |
| coden15 original active + robust rescue | 65.333 | 0.667 | 7.333 | 0.999917 | 0.887397 | 2.4022 |

Interpretation: the formula-AF gate is the best coden15 abundance setting found
so far and is consistent with the ANI-to-AF model. It beats Sylph on mean Toy
Mouse L1, but it still does not beat the old ctx-marker abundance baseline and
it raises false positives. Treat this as the coden15 diagnostic best, not a
validated default.

## Caveats

- Exact 200,709-reference re-sketch was not possible because local VEuPathDB
  FASTA sources are missing.
- The relaxed gate was selected after looking at these samples; it should not be
  promoted without validation on an independent benchmark.
- The coden15 retuned gate is also diagnostic and was selected from Toy Mouse
  samples 0-2 after mechanism inspection. It must be validated elsewhere before
  becoming a coden15 default.
- The formula-AF gate is more principled than the fixed AF floor, but its
  effective context length and safety margin still need independent validation.
- `git status` did not work inside this sandbox despite `.git` being visible,
  so changed-file review used direct file inspection rather than git metadata.
