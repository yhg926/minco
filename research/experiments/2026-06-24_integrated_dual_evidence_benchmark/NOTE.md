# Integrated Full S2000 Dual-Evidence Benchmark

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory used: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: local dirty worktree
- Changed implementation files:
  - `minco_core/src/command_ani.c`
  - `minco_core/src/command_ani.h`
  - `minco_core/src/command_ani_wrapper.c`
- Build outputs:
  - coden11: `minco_core/bin/minco`
  - coden15 smoke-test binary: `minco_core/bin_coden15/minco`
- Provenance snapshots:
  - `provenance/code_status.txt`
  - `provenance/code_diffstat.txt`

## Question

After implementing integrated full S2000 dual counters in MinCO `ani`, does the
full-refdb output with marker-only evidence improve F1 or abundance over:

1. Sylph;
2. current mixed MinCO F1/L1 rows;
3. the previous full S2000 marked-unique proxy?

## Implementation

Added `--readwise-dual-evidence` to `minco ani`.

When enabled in readwise depth mode, MinCO keeps the normal full-refdb evidence
columns unchanged and appends marker-only columns computed from the same full
inverted index. A reference context is marker-only when its context occurs in
exactly one reference in the sorted full index. The same conflict-object rule
is preserved: `--ignoreconflict` still removes reference context groups with
multiple objects.

Appended columns include:

```text
Marker_XnY_ctx
Marker_Raw_XnY_ctx
Marker_N_diff_obj
Marker_N_diff_obj_section
Marker_N_mut2_ctx
Marker_Ref_ctx_total
ctx_marker_size
Marker_Ref_breadth
Marker_Ref_mean_depth
Marker_Ref_hit_mean_depth
Marker_Ref_depth_variance
Marker_Ref_depth_cv
Marker_Ref_zero_fraction
Marker_Relative_abundance_depth
Marker_Effective_abundance_depth
Marker_Ref_zip_af
Marker_Ref_zip_aaf_ani
Marker_Reliable_Ref_breadth
Marker_Reliable_Ref_mean_depth
Marker_Reliable_Ref_hit_ctx
Marker_Reliable_Ref_hit_mean_depth
Marker_Reliable_Ref_hit_median_depth
Marker_Reliable_Ref_hit_depth_variance
Marker_Reliable_Ref_zip_af
```

The marker object-difference columns are required so downstream gates can
recompute marker `ANI_naive_calc` exactly as for physical markerdb output.

## Dataset

F1 and abundance were scored on the same six-sample panel as the previous
mixed abundance benchmark:

```text
mouse_gtdb: CAMI II Toy Mouse Gut samples 0,1,2; GTDB species truth
cami3_ncbi: CAMI3 ToyGut samples 0,1,2; bacterial NCBI species-taxid truth
```

Reference:

```text
/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup
```

## Commands

Full commands are recorded in `commands.sh`.

Each sample used:

```text
minco_core/bin/minco ani -p4 \
  -r /tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup \
  --qraw <reads.fq.gz> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-dual-evidence \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0
```

## Scoring Rule

The scoring rule intentionally matches the previous hybrid proxy:

```text
marker branch:
  use Marker_* columns as if they came from a physical ctx-markerdb

full branch:
  use normal full-refdb columns

baseline:
  marker robust effective depth + intra-genus winner rescue

full fallback candidates:
  ctx_marker_size < 200
  full XnY_ctx >= 500, 700, or 1000
  abundance contribution = 0.4 * full robust effective depth
  do not add full fallback when marker branch already selected that species/taxid
```

## Results

Six-sample mouse+CAMI3 comparison:

```text
method                                                all_F1    all_L1_pp
Sylph baseline                                        0.783546  15.038388
MinCO current mixed F1-priority                       0.782475  18.014617
MinCO current mixed L1-priority                       0.770982  17.823364
integrated_hybrid_baseline_plus_full_xny700_scale04   0.750078  20.478490
integrated_hybrid_baseline_plus_full_xny1000_scale04  0.748810  20.498238
integrated_hybrid_baseline_plus_full_xny500_scale04   0.748574  20.451236
integrated_marker_only                                0.747536  20.513677
```

Best integrated full fallback:

```text
method: integrated_hybrid_baseline_plus_full_xny700_scale04
all F1: 0.750078
all L1: 20.478490 percentage points
mouse F1: 0.944844
mouse L1: 1.724130
CAMI3 F1: 0.555312
CAMI3 L1: 39.232850
```

This exactly reproduces the previous proxy row and does not improve over it.

## Validation

Builds:

```text
make -C minco_core
make -C minco_core OBJDIR=obj_coden15 BINDIR=bin_coden15 \
  CFLAGS='-std=gnu11 -Wno-format-overflow -Wno-unused-result -O3 -flto -fopenmp -DNUM_CODENS=15 -DMINCO_ENABLE_CTXOBJ96=1'
```

Smoke tests:

- `minco_core/bin/minco ani --help` shows `--readwise-dual-evidence`.
- `minco_core/bin_coden15/minco ani --help` shows `--readwise-dual-evidence`.
- A small coden15 readwise smoke test emitted the new marker columns including
  `Marker_N_diff_obj`, `Marker_N_diff_obj_section`, and `Marker_N_mut2_ctx`.

Runtime:

```text
mouse_s0: 1:50.95 wall, 6.20 GB max RSS
mouse_s1: 1:53.79 wall, 6.52 GB max RSS
mouse_s2: 1:51.47 wall, 6.49 GB max RSS
cami3_s0: 1:49.64 wall, 6.13 GB max RSS
cami3_s1: 1:49.61 wall, 6.19 GB max RSS
cami3_s2: 1:48.95 wall, 6.06 GB max RSS
```

Marker-size validation against the physical markerdb PSMP map:

```text
file               comparable_rows  missing_from_psmp  comparable_mismatches
cami3_s0_dual.tsv  102809           320                0
cami3_s1_dual.tsv  111041           361                0
cami3_s2_dual.tsv  101441           309                0
mouse_s0_dual.tsv  118403           284                0
mouse_s1_dual.tsv  153204           362                0
mouse_s2_dual.tsv  154923           337                0
```

Rows missing from the old PSMP map were not comparable; all comparable rows
matched exactly.

Scorer:

```text
score_integrated_dual_evidence.py elapsed 3:27.57, max RSS 4.53 GB, exit status 0
```

The only scorer stderr was a pre-existing SciPy binary compatibility warning.

## Add-back and EM Follow-up

After the first integrated benchmark showed the inherited hybrid fallback was
not enough, the scorer was extended to test the previously best raw
`value_sum` shape on the integrated output:

```text
marker_l1 full_f1/full_l1/full_cami direct value_sum
marker_l1 full_f1/full_l1/full_cami marker-present-only value_sum
EM-proxy redistribution grid:
  beta: 0.25,0.5,0.75,1.0
  seed_fraction: 0,0.001,0.01,0.05,0.1
  prior_power: 0.5,1.0,2.0
```

The direct unrestricted full-S2000 add-back failed badly because many remote
shared-context references pass the full evidence gate:

```text
method                                      all_F1    all_L1   mouse_F1  cami3_F1
integrated full_f1 direct value_sum          0.1371   88.9278   0.1567    0.1176
integrated full_l1 direct value_sum          0.1287   84.6045   0.1426    0.1148
integrated full_cami direct value_sum        0.0969   93.9318   0.1094    0.0843
```

The best EM-proxy and marker-present-only variants avoided the FP explosion,
but they collapsed to the marker-supported candidate set and did not beat the
current mixed MinCO or Sylph rows:

```text
method                                                     all_F1   all_L1
Sylph baseline                                             0.7835   15.0384
MinCO current mixed F1-priority                            0.7825   18.0146
integrated full_f1 marker-present EM b1 s0 p1              0.7752   18.0612
integrated full_l1 marker-present value_sum                0.7752   20.0815
MinCO current mixed L1-priority                            0.7710   17.8234
```

Interpretation: aggregate per-reference totals are not enough for useful shared
mass assignment. A real EM implementation needs context-level or read-level
ambiguity groups, so that each shared context/read can distribute its own mass
only among references that actually contain that context.

## Conclusion

The integrated dual-counter implementation works and removes the need for a
physical markerdb copy when testing marker/full logic from the full S2000
refdb.

Accuracy did not improve over the previous proxy, and the aggregate add-back/EM
prototypes did not recover the older physical ctx+obj markerdb blend. Therefore
the architecture is useful for storage/runtime and future optimization, but it
is not itself the missing abundance mechanism.

Current decision:

```text
F1 priority: Sylph still slightly leads the six-sample mixed panel.
MinCO F1-priority mixed row remains the best MinCO F1 candidate.
Abundance L1: Sylph still clearly leads the six-sample mixed panel.
Integrated full S2000 dual evidence only reproduces the prior proxy.
Aggregate full-S2000 shared-mass assignment is not sufficient.
```

## Caveats

- The benchmark panel has only six samples and mixes GTDB mouse truth with
  CAMI3 NCBI species-taxid truth.
- The integrated fallback thresholds are inherited from the previous proxy
  experiment, not retuned after implementation.
- The integrated code currently appends dual evidence but does not implement a
  production default gate inside MinCO.
- ANI accuracy was not rerun here; keep the independent coden15 ANI reporter
  conclusion from the prior ANI notes.

## Next Experiment

Use the integrated full/marker output as the default input for optimization,
but change the mechanism:

```text
1. implement context-level/read-level shared-context ambiguity groups;
2. run EM-like mass redistribution within each ambiguity group;
3. keep F1, abundance L1, and ANI reporting as separate objectives.
```
