# Fresh Raw-Run Profile Rescue Validation

## Question

Does the high-specificity uncalled-candidate profile rescue remain safe when
MinCO is rerun from raw reads with `--report-all`, rather than replayed from
cached profile tables?

## Scope

This experiment validates the cached best candidate:

- `profile_rescue_p002_x300_ani095_af06_b025_d1_top1`
- Current abundance guard: ANI floor `0.90` plus sparse high-depth cap.
- Fresh raw-run samples with local reads:
  - CAMI3 toy human gut samples 0, 1, and 2.
  - HMP air/skin sample 13.

Large fresh profiles are written under
`/tmp/minco_profile_rescue_fresh_raw_validation_20260701`.

## Decision Rule

The rescue is promotable only if fresh runs preserve the cached behavior:

- no added false species in the fresh non-CAMI panel checked here;
- L1 improves or stays neutral relative to the current abundance guard;
- F1 does not regress relative to the current called set.

If these conditions fail, keep the rule diagnostic and leave the current default
unchanged.

## Fresh Results

Fresh raw-read profiling completed for all four locally available samples:

| Panel | Samples | Current pooled F1 | Rescue pooled F1 | Current L1 pp | Rescue L1 pp | Added truth | Added FP |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CAMI3 toy human gut | 0,1,2 | 0.814570 | 0.826230 | 58.155186 | 39.708920 | 6 | 0 |
| HMP air/skin | 13 | 0.925000 | 0.925000 | 54.747984 | 54.747984 | 0 | 0 |

Overall across the two fresh panels, the high-specificity rescue reduced mean
official L1 from `56.451585` to `47.228452` pp and improved mean Pearson from
`0.877137` to `0.913351`. The rescue added six rows, all CAMI3 truth species,
with no fresh false additions.

The fresh profiles were small enough for repo-safe summaries only:

- CAMI3 sample0: 90.04 s, 3,189,548 KB peak RSS.
- CAMI3 sample1: 112.04 s, 3,218,352 KB peak RSS.
- CAMI3 sample2: 111.67 s, 3,150,288 KB peak RSS.
- HMP air/skin sample13: 88.13 s, 3,285,296 KB peak RSS.

Toy mouse samples 5-7 and HMP air/skin samples 6, 11, and 28 were not rerun
because their raw read files were no longer local after cleanup.

## Conclusion

The high-specificity rescue passed this fresh raw-run check on the available
mixed samples. It is a reasonable candidate for a controlled default trial, but
the evidence is still narrower than the cached cross-panel replay because the
toy mouse and additional HMP raw reads were not available locally.
