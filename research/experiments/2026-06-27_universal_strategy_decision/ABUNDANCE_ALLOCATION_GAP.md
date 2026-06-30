# Abundance Allocation Gap Manifest

Date: 2026-06-29

This tracked note is generated from cached decision TSVs. It records the
remaining abundance gap for the selected MinCO default without rerunning
raw profiling jobs.

## Source Tables

- `results/abundance_error_decomposition_delta.tsv`
- `results/cross_panel_abundance_variant_overall.tsv`
- `results/candidate_default_vs_sylph.tsv`
- `results/candidate_default_vs_current.tsv`
- `results/candidate_preset_genus_xny_blend_audit.tsv`
- `results/candidate_preset_genus_xny_blend_panel_delta.tsv`
- `results/candidate_preset_genus_xny_blend_guard_audit.tsv`
- `results/candidate_preset_genus_xny_alpha_sweep_audit.tsv`
- `results/candidate_preset_genus_xny_alpha_sweep_overall.tsv`

## Panel Gap

| Panel | Main gap | Selected default - Sylph L1 pp | Legacy matched-gap pp | Legacy missing-gap pp | Legacy extra-gap pp |
|---|---:|---:|---:|---:|---:|
| cami2_toy_mouse_gut | matched_abs_error | 6.798669 | 5.481572 | 1.554743 | -0.147326 |
| hmp_airskin_gtdb_source_abundance | matched_abs_error | 17.056351 | 11.310500 | 0.949829 | 4.921452 |
| hmp_gastrooral_gtdb_source_abundance | matched_abs_error | 41.294813 | 38.274541 | 4.474012 | 1.175309 |
| cami3_toy_human_gut_gtdb_source_readmap | missing_truth_mass | 25.114394 | 5.585748 | 19.641662 | 4.871778 |

## Cached Allocation Candidate

| Method | Mean L1 delta vs current pp | Max worse pp | Improved panels | Status |
|---|---:|---:|---:|---|
| blend_current_genus_realloc_xny_a0.25 | -0.399949 | -0.151769 | 4 | cached_signal_not_promoted_without_raw_validation |

## Selected-Default Blend Audit

| Method | Mean L1 delta vs selected default pp | Max worse pp | Sample direction | Decision |
|---|---:|---:|---|---|
| default_candidate_preset_genus_xny_blend_a0.25 | -0.780908 | 0.993420 | improved=26;worsened=6 | do_not_promote_abundance_blend |

| Panel | Mean L1 delta pp | Max worse pp | Improved samples | Worsened samples |
|---|---:|---:|---:|---:|
| cami2_toy_mouse_gut | -0.156283 | 0.019185 | 2 | 1 |
| cami3_toy_human_gut_gtdb_source_readmap | -0.599122 | -0.178298 | 3 | 0 |
| hmp_airskin_gtdb_source_abundance | -0.932102 | 0.993420 | 20 | 4 |
| hmp_gastrooral_gtdb_source_abundance | -0.176193 | 0.261544 | 1 | 1 |

## Adaptive Blend Guard Audit

| Tested rules | Sample-safe in-panel | LOPO result | Decision |
|---:|---:|---|---|
| 3722 | 394 | mean_delta=-0.505883;holdouts_with_mean_regression=0;holdouts_with_sample_regression=2;max_sample_worse=0.220479 | sample_safe_in_panel_but_fails_lopo |

## Conservative Alpha Sweep

| Variants | Sample-safe variants | Best ranked variant | Decision |
|---:|---:|---|---|
| 6 | 0 | genus_xny_a0.15;mean_panel_delta=-0.309376;worsened_samples=4;max_worse=0.580806 | do_not_promote_alpha_sweep |

| Method | Mean panel L1 delta pp | Worsened samples | Max sample worse pp |
|---|---:|---:|---:|
| genus_xny_a0.15 | -0.309376 | 4 | 0.580806 |
| genus_xny_a0.10 | -0.222767 | 4 | 0.374499 |
| genus_xny_a0.05 | -0.122719 | 4 | 0.168191 |
| genus_xny_a0.02 | -0.050842 | 4 | 0.044407 |
| genus_xny_a0.25 | -0.465925 | 6 | 0.993420 |
| genus_xny_a0.25_cap0.005 | -0.413174 | 6 | 0.949474 |

## Decision

- Selected default status: `selected_default_is_best_minco_candidate_not_abundance_release_claim`.
- The candidate preset improves MinCO versus its previous default on the
  matched panels, but it still loses abundance L1 to Sylph on all four
  cached comparison panels.
- The fixed-call allocator sweep has a small safe cached signal, but the
  selected-default posthoc audit has individual-sample regressions, so
  it is not promoted.
- Output-derived blend guards can be sample-safe in-panel, but the
  leave-one-panel-out audit still has held-out sample regressions, so
  no adaptive blend guard is promoted.
- Lower-alpha and capped genus-XnY sweeps improve all panel means, but
  no nonzero tested variant is sample-safe, so this family remains
  off by default.
- Next abundance work should target matched-call mass allocation first.
  The selected-call oracle feasibility audit in
  `ABUNDANCE_NEXT_TARGET.md` narrows call recovery to one cached
  insufficient-call sample under the selected candidate preset.
