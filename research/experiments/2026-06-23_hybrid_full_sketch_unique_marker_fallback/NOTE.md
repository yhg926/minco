# Hybrid Full Sketch Unique Marker Fallback

Date: 2026-06-23
Author/agent: ubuntu
Project: KSSD3mini

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Commit: `NA`
- Ref/describe: `NA`
- Working tree: `NA`
- Binary/script/tool: `minco-dev`
- Provenance files: `provenance/code_status.txt`

Git metadata was unavailable in this sandbox, so the experiment is tracked by
local script paths, command provenance, and result files.

## Question

Can MinCO keep the full S1000/S2000 sketch available, use ctx-marker logic as
the default, and add full-sketch evidence only for references with too few
unique markers, improving F1 and abundance without reintroducing shared-context
false positives?

## Hypothesis

The inverted reference index tells us which references are marker-poor. The safe
mode is not to replace markerdb scoring, but to keep the ctx-marker baseline and
only add full-reference evidence for marker-poor references when support is very
strong. Full fallback should add missing taxa; it should not overwrite marker
abundance for taxa already selected by marker evidence.

## Dataset

- CAMI II toy mouse samples 0-2 with GTDB species ground truth from previous
  experiment outputs.
- CAMI3 ToyGut samples 0-2 with NCBI species truth profiles.
- Reference sketches:
  - ctx-marker: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
  - full S2000 dedup: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Marker-size map: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv`

## Methods

Exploratory steps:

1. Sweeped naive marker-size fallback thresholds on mouse. Simple replacement
   was rejected because it created many false positives.
2. Added a full-fallback support threshold. `XnY_ctx >= 500` recovered mouse
   TPs but created CAMI3 Veillonella FPs.
3. Raised fallback support to `XnY_ctx >= 700`. This kept the mouse rescue and
   removed the CAMI3 extra FPs.
4. Scaled fallback abundance depth. A 0.4 multiplier gave the best mouse L1
   among the tested F1-equivalent rules.
5. Changed to baseline-plus logic: keep all ctx-marker selected taxa and add
   full fallback only for marker-poor taxa not already selected by marker.

Selected rule:

```text
Use coden11 S2000 ctx-marker robust depth + intra-genus rescue as baseline.
For full-reference fallback candidates:
  ctx_marker_size < 200
  XnY_ctx >= 700
  active gate passes under the existing product0-adjusted rule
  species/taxid is not already selected by ctx-marker evidence
  abundance contribution = 0.4 * robust_effective_depth
```

Commands are recorded in `commands.sh`; parameters are recorded in
`parameters.tsv`.

## Results

Mouse GTDB species truth, samples 0-2:

```text
old ctx-marker:         F1=0.939759  L1=1.794503
hybrid baseline-plus:   F1=0.944844  L1=1.724130
delta:                 +0.005085    -0.070374 pp
```

CAMI3 ToyGut NCBI species truth, samples 0-2:

```text
old ctx-marker:         F1=0.555312  L1=39.232850
hybrid XnY>=700:        F1=0.555312  L1=39.232850
hybrid XnY>=500:        F1=0.552304  L1=39.178343  rejected due extra FP
Sylph:                  F1=0.603485  L1=27.789743
```

CAMI3 full-reference readwise generation used 6.40 GB maximum RSS and
approximately 1:54-2:11 wall clock per sample with 8 threads.

## Validation

- `python3 -m py_compile` passed for the three scoring scripts.
- Mouse final scorer reran and reproduced the selected result table:
  `hybrid_mouse_baseline_plus_mean_metrics.tsv`.
- CAMI3 final scorer reran and reproduced the selected result table:
  `hybrid_cami3_summary.tsv`.
- CAMI3 `XnY>=500` failure mode was inspected: it added three sample0
  Veillonella false positives with support around 596-615. `XnY>=700` removed
  them.

## Important Artifacts

See `artifacts.md` for generated outputs, scripts, logs, and external
reference paths.

## Conclusion

The hybrid idea is worth keeping, but not as a replacement of markerdb scoring.
The current best candidate is a conservative baseline-plus fallback:

`hybrid_baseline_plus_full_xny700_scale04`.

It improves MinCO over the old ctx-marker baseline on the three mouse samples
for both F1 and abundance L1, and it is neutral on CAMI3 ToyGut. It does not
close the CAMI3 gap to Sylph.

## Paper-Relevant Claim

Full-sketch evidence can safely rescue marker-poor references when used as a
high-support additive fallback, but low-threshold full-context fallback
reintroduces shared-context false positives.

## Caveats

- Marine and plant full-reference hybrid-compatible panels were not regenerated
  in this experiment.
- The chosen `ctx_marker_size < 200`, `XnY_ctx >= 700`, and abundance scale 0.4
  are empirical from mouse plus CAMI3. They should be retested before becoming a
  hard default.
- ANI reporting was not retested here; keep the previous independent ANI
  strategy until a separate ANI-focused benchmark says otherwise.

## Next Experiment

Run the same baseline-plus fallback on the marine and plant calibration panels,
or implement it inside MinCO so the full-sketch inverted index and marker-size
logic are available without post-processing.
