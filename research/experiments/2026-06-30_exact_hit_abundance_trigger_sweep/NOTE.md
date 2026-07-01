# Exact-Hit Abundance Trigger Sweep

Date: 2026-06-30

## Question

Can the abundance-only exact-split sidecar be enabled by default only on samples
where the block-mode profile indicates enough unresolved abundance mass, while
avoiding samples where exact-hit abundance hurts?

## Code Provenance

- Repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git commit: `87f426e7d4739e7d1664ae8c520722c264050e17`
- Branch: `main`
- Git status: dirty. Relevant uncommitted changes include
  `scripts/minco_profile_calibrated.py`,
  `tests/test_minco_profile_calibrated_auto_exact.py`, and this experiment
  folder.
- Scripts:
  - `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/sweep_exact_hit_abundance_trigger.py`
  - `research/experiments/2026-06-30_exact_hit_abundance_trigger_sweep/score_trigger_validation_reruns.py`

## Method

The first pass was an offline cached-profile replay across CAMI2 toy mouse gut,
HMP airskin, HMP gastrooral, and CAMI3 toy human gut profiles. It held calls
fixed and changed only abundance mass when an exact-split abundance sidecar
would be requested.

Trigger rule tested:

- Sidecar requested when
  `auto_exact_split_block_p_extra_mass_ratio >= threshold`
  and `auto_exact_split_guard_passed`
  and `auto_exact_split_low_extra_gate_passed`.
- Offline thresholds: disabled/current, `0.50`, `0.40`, `0.30`, `0.25`,
  `0.20`, `0.15`, `0.10`, `0.05`, `0.00`.
- Exact-hit abundance uses `Ref_hit_mean_depth` from the exact split table
  matched by `s_best_accession`; missing exact hits fall back to the normal
  profile abundance raw value.
- Calls are fixed, so F1 changes are not expected.

The cached sweep suggested `>=0.10`, but four triggered samples needed raw-read
validation because local cached exact tables or accession columns were
incomplete. Those samples were restored or regenerated under
`/tmp/minco_exact_trigger_validation_20260630` and rerun with
`--exact-split-abundance-trigger 0.10`.

## Results

Raw rerun deltas for the four threshold-0.10 triggered samples:

| panel | sample | baseline L1 pp | trigger 0.10 L1 pp | delta L1 pp | baseline Pearson | trigger Pearson | delta Pearson |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CAMI3 toy human gut | 0 | 62.229113 | 58.074209 | -4.154904 | 0.916810 | 0.939886 | 0.023076 |
| CAMI3 toy human gut | 1 | 67.193605 | 66.923142 | -0.270463 | 0.721659 | 0.704242 | -0.017417 |
| CAMI3 toy human gut | 2 | 43.818985 | 107.850907 | 64.031922 | 0.954493 | 0.692052 | -0.262442 |
| HMP airskin | 13 | 54.747984 | 34.686507 | -20.061477 | 0.910524 | 0.956996 | 0.046471 |

Validated cross-panel abundance summary:

| method | mean official L1 pp | delta vs current L1 pp | max panel L1 regression pp | improved panels | worsened panels | mean Pearson | delta Pearson | decision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| validated current, no abundance exact sidecar | 35.760147 | 0.000000 | 0.000000 | 0 | 0 | 0.928874 | 0.000000 | baseline |
| exact-hit abundance trigger `>=0.10` | 32.993713 | -2.766434 | 19.868851 | 2 | 1 | 0.930787 | 0.001912 | reject: raw CAMI3 regression |
| exact-hit abundance trigger `>=0.30` | 28.235474 | -7.524673 | 0.000000 | 1 | 0 | 0.951701 | 0.022827 | promote default candidate |
| Sylph external baseline | 12.446060 | NA | NA | NA | NA | 0.984206 | NA | external abundance baseline |

Panel-level summary for the promoted `>=0.30` trigger:

| panel | current L1 pp | trigger 0.30 L1 pp | Sylph L1 pp | current Pearson | trigger 0.30 Pearson | Sylph Pearson |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CAMI2 toy mouse gut | 13.767837 | 13.767837 | 6.731523 | 0.987977 | 0.987977 | 0.994040 |
| CAMI3 toy human gut | 57.747234 | 57.747234 | 32.632841 | 0.864321 | 0.864321 | 0.945708 |
| HMP airskin | 22.538526 | 22.538526 | 5.356745 | 0.980924 | 0.980924 | 0.997950 |
| HMP gastrooral | 48.986991 | 18.888299 | 5.063129 | 0.882277 | 0.973584 | 0.999125 |

Runtime and memory for the four raw validation reruns with trigger `0.10`:

| sample | wall time | max RSS KB |
| --- | ---: | ---: |
| CAMI3 sample 0 | 3:22.72 | 3,627,124 |
| CAMI3 sample 1 | 3:15.63 | 3,658,032 |
| CAMI3 sample 2 | 3:05.47 | 3,585,840 |
| HMP airskin sample 13 | 3:06.76 | 3,691,952 |

Restored raw-read inputs occupy about 18 GB under
`/tmp/minco_exact_trigger_validation_20260630/reads`; rerun outputs occupy
about 522 MB under `/tmp/minco_exact_trigger_validation_20260630/run`.

## Interpretation

The raw reruns reject `0.10`: CAMI3 sample 2 has a large abundance regression
when abundance exact-hit mode is triggered there. The stricter `0.30` trigger
does not fire on the CAMI3 triggered samples, preserves current CAMI3 behavior,
and still captures the large HMP gastrooral sample 6 improvement seen in the
cached validation.

## Conclusion

Promote `--exact-split-abundance-trigger 0.30` as the calibrated profile
default. Keep `--exact-split-abundance-trigger inf` as the explicit opt-out.

This improves MinCO's validated abundance relative to the previous current
default in this panel, but Sylph remains the stronger external abundance
baseline on these datasets.

## Caveats

- The trigger changes abundance only; it does not change species calls/F1.
- The validated panel is still limited to the available CAMI2 toy mouse gut,
  CAMI3 toy human gut, HMP airskin, and HMP gastrooral profiles.
- `0.30` was selected to avoid the observed raw CAMI3 regression; more holdout
  panels should continue to monitor high-extra-mass samples.
- The `/tmp` raw reads and rerun outputs are temporary and may be removed.
