# Default Profile Rescue Productization Validation

## Question

After wiring the strict split-evidence rescue into `scripts/minco_profile`, does
the normal default entry point recover the same high-support candidates without
requiring users to pass the expert calibrated-wrapper switch?

## Methods

The validation runs the user-facing launcher `scripts/minco_profile` from saved
unique/split profile tables under
`/tmp/minco_profile_rescue_fresh_raw_validation_20260701/work`.

Compared methods:

- `default_candidate`: normal `scripts/minco_profile` candidate preset.
- `default_candidate_no_profile_rescue`: same command plus `--no-profile-rescue`.

Both are run with explicit `--taxmap`, `--model-cache`, `--scope bacteria`, and
`--report-all` because this validation starts from precomputed tables rather
than a packaged `--ref` sidecar.

Large profile outputs are written under
`/tmp/minco_profile_default_rescue_productization_validation_20260701`.

## Parameters

The default `candidate` preset injects the strict split-evidence rescue switch
`split-p002-x300-ani095-af06-b025-d1-top1` unless `--no-profile-rescue` is
passed. It also keeps the candidate abundance policy
`normalized-depth-alpha2` and the abundance-only sparse-depth guards:

- `abundance_ani_floor = 0.90`
- `abundance_sparse_depth_cap = poisson-breadth`
- `abundance_sparse_breadth_max = 0.15`
- `abundance_sparse_depth_ratio_min = 200`

The strict rescue keeps native raw abundance mass for newly rescued calls; the
normalized-depth rescue mass is only applied to older loose rescue rows.

## Results

The default entry point matched the strict-switch validation behavior.

Panel-level scores:

- CAMI3 source-readmap panel, samples 0-2:
  - default candidate: pooled F1 `0.826230`, TP `252`, FP `62`, FN `44`,
    official L1 `39.708920`, Pearson `0.916179`
  - opt-out: pooled F1 `0.814570`, TP `246`, FP `62`, FN `50`, official L1
    `58.155186`, Pearson `0.843750`
- HMP sample 13 panel:
  - default candidate and opt-out were unchanged: pooled F1 `0.925000`, TP
    `37`, FP `3`, FN `3`, official L1 `54.747984`, Pearson `0.910524`

Across the two panels, default candidate had mean official L1 `47.228452` and
mean Pearson `0.913351`; opt-out had mean official L1 `56.451585` and mean
Pearson `0.877137`.

The rescue added 6 rows across the CAMI3 panel. All 6 were truth-matching rows
and 0 were false additions. The added rows represented `50.120600` percentage
points of truth abundance mass in the benchmark mapping.

Table-mode runtime, starting from saved unique/split tables, was similar with
and without rescue:

- default candidate: `8.83-11.03 s`, peak RSS `1.41-1.55 GB`
- opt-out: `8.74-10.95 s`, peak RSS `1.41-1.55 GB`

Raw-input smoke test:

- command helper: `run_raw_default_smoke.py`
- temporary run root: `/tmp/minco_default_raw_smoke_20260701`
- output rows: `1`
- called rows: `1`
- candidate rescue switch: `split-p002-x300-ani095-af06-b025-d1-top1`
- candidate surface switch: `accession-ani90-xny100-br01-af70`
- candidate abundance policy: `normalized-depth-alpha2`
- sparse-depth abundance guard: `poisson-breadth`
- exact-density replay: used, source `density_cache`
- exit status: `0`
- elapsed time: `0:01.48`
- peak RSS: `226096 KB`

## Conclusion

The user-facing `scripts/minco_profile` default now exposes the validated strict
rescue without requiring users to know the calibrated-wrapper option. In this
validation it improved the CAMI3 source-readmap panel without adding false
rows, and it left the HMP sample unchanged.

## Caveats

This validation starts from saved unique/split tables, so it measures the
default profile decision layer and not raw input scanning. Candidate-surface
features were disabled in this table-mode run because the supplied taxmap is a
species-level table without the candidate-surface sidecar. Broad release claims
still rely on the larger cross-panel benchmark notes.
