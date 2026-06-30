# Abundance Strategy Decision Matrix

Date: 2026-06-29

This tracked note is generated from cached summary TSVs. It records which
abundance-related strategies are part of the selected MinCO default and
which ones remain diagnostic only. It does not rerun raw profiling jobs.

## Source Tables

- `results/adaptive_abundance_switch_audit.tsv`
- `results/candidate_preset_genus_xny_alpha_sweep_overall.tsv`
- `results/candidate_preset_genus_xny_alpha_sweep_audit.tsv`
- `results/candidate_preset_genus_xny_blend_audit.tsv`
- `results/candidate_preset_genus_xny_blend_guard_audit.tsv`
- `results/cami3_source_profile_binomial_extension_audit.tsv`
- `results/cami3_source_readmap_recovery_inputs_audit.tsv`
- `results/candidate_callset_oracle_feasibility_audit.tsv`
- `results/candidate_callset_oracle_feasibility_summary.tsv`
- `results/candidate_default_decision_audit.tsv`
- `results/cross_panel_candidate_abundance_wrapper_audit.tsv`
- `results/edge_em_cross_domain_policy_summary.tsv`
- `results/feature_allocator_external_exactsplit_audit.tsv`
- `results/feature_allocator_combined_guard_audit.tsv`
- `results/feature_allocator_refined_independent_holdout_inventory_audit.tsv`
- `results/feature_allocator_refined_marine_diagnostic_audit.tsv`
- `results/feature_allocator_refined_guard_parity_audit.tsv`
- `results/feature_allocator_refined_holdout_recovery_audit.tsv`
- `results/feature_allocator_refined_score_replay_audit.tsv`
- `results/feature_allocator_release_candidate_audit.tsv`
- `results/feature_allocator_wrapper_parity_audit.tsv`
- `results/abundance_variant_safety_audit.tsv`
- `results/abundance_oracle_bounds_audit.tsv`
- `results/selected_call_feature_allocator_audit.tsv`
- `results/selected_call_feature_allocator_guard_audit.tsv`
- `results/selected_call_mass_transform_audit.tsv`
- `results/selected_call_mass_transform_guard_audit.tsv`
- `results/supervised_abundance_external_exactsplit_audit.tsv`
- `results/supervised_abundance_calibrator_audit.tsv`

## Decision Matrix

| Strategy family | Scope | Best signal | Failure mode | Decision |
|---|---|---|---|---|
| selected candidate preset | default wrapper | F1=0.007199;L1=-1.994230;Pearson=0.013705 | F1_wins=1/4;L1_wins=0/4;Pearson_wins=0/4 | promote_candidate_for_minco_default_study_not_release_claim |
| candidate row abundance | candidate-only rows | mean_L1=-0.722425;max_worse_L1=0.000000;mean_Pearson=0.004448;worsened_samples=0 | counts=0;F1=1.11e-16;L1=4.2;Pearson=0.0302 | keep_as_candidate_preset_component_not_general_allocator |
| truth-aware allocation upper bound | current calls only | recoverable_L1_pp=22.82540327656117 | oracle_vs_sylph_panels=2/4 | allocation_model_can_close_some_not_all_sylph_gap |
| fixed-call formula sweep | 32 cached scored profiles | genus_realloc_xny;mean_delta=-1.865556;worsened_samples=9;max_worse=4.177992 | sample_safe=0 | reject_fixed_call_abundance_replacement |
| adaptive output switch | output-derived thresholds | genus_realloc_xny;breadth_ratio_median>=1.11111;mean_delta=-1.195622;worsened_panels=0;max_worse=0.000000;switched=25 | mean_delta=-0.737243;worsened=1;max_worse=1.218914 | reject_as_overfit_by_leave_one_panel_out |
| supervised row calibrator | leave-one-panel-out | lopo_hgb_log_l2_1_blend0.75;mean_delta=-0.090593;worsened_panels=1;max_worse=0.062803 | candidate_requires_independent_holdout | candidate_requires_independent_holdout |
| supervised external stress test | cached external exact-split profiles | rf_log_leaf3_blend0.75;mean_delta_l1=-0.000068;worsened_datasets=1;max_dataset_worse=0.001128 | external_safe_blends=0 | reject_supervised_abundance_external_default |
| edge-level mass redistribution | three spot panels | vs_current_L1_delta=-18.5786236509 | vs_current_F1_delta=-0.0150869449894; vs_sylph_F1_delta=-0.0514250837048 | not promotable: F1 is lower on CAMI3 and marine |
| selected-default genus-XnY blend | posthoc selected-default profiles | mean_L1_delta_pp=-0.780907827013158 | sample_direction=improved=26;worsened=6; max_worse_pp=0.9934199183847312 | do_not_promote_abundance_blend |
| selected-default blend guard | output-derived guard rules | two_feature_and;s_xny_mean<=783.804 AND blend_perturb_max_frac>=0.00184827;mean_delta=-0.798743;switched=22 | mean_delta=-0.505883;holdouts_with_mean_regression=0;holdouts_with_sample_regression=2;max_sample_worse=0.220479 | sample_safe_in_panel_but_fails_lopo |
| selected-default alpha/cap sweep | conservative nonzero blends | genus_xny_a0.15;mean_panel_delta=-0.309376;worsened_samples=4;max_worse=0.580806 | sample_safe=0 | do_not_promote_alpha_sweep |
| selected call-set oracle feasibility | 32 selected-candidate cached samples | allocation_only_can_close_samples=31/32 | call_recovery_required_samples=1 | prioritize_allocator_plus_call_recovery |
| selected-call simple mass transforms | base called rows; candidate-row mass preserved | base_genus_xny_a0.10;mean_sample_delta=-0.382481;worsened_samples=4;max_worse=0.161705 | sample_safe=0 | do_not_promote_selected_mass_transform |
| selected-call guarded mass transforms | output-derived guards over selected-call profiles | base_global_xny_a0.10;single;n_calls>=71.4;mean_delta=-0.491859;switched=13 | mean_delta=-0.129226;holdouts_with_mean_regression=1;holdouts_with_sample_regression=2;max_sample_worse=0.825999 | sample_safe_in_panel_but_fails_lopo |
| selected-call feature allocators | fixed calls; output depth/quality target features | genus_raw_hit_breadth_a0.02;mean_sample_delta=-0.318875;worsened_samples=1;max_worse=0.040714 | sample_safe=0 | do_not_promote_feature_allocator |
| selected-call guarded feature allocators | output-derived guards over feature allocators | genus_raw_hit_breadth_a0.02;single;base_mass_multi_genus_frac>=0.375015;mean_delta=-0.312324;switched=28 | mean_delta=-0.146115;holdouts_with_mean_regression=0;holdouts_with_sample_regression=0;max_sample_worse=0.000000 | candidate_guard_needs_wrapper_and_independent_holdout |
| implemented guarded feature allocator | wrapper function parity on cached selected profiles | mean_L1_delta_pp=-0.318866 | worsened_samples=0; needs independent holdout | candidate_needs_independent_holdout |
| implemented guarded feature allocator external | plant/strain cached exact-split diagnostic profiles | mean_delta_l1=-0.004134;improved_datasets=2;worsened_datasets=0;improved_samples=4;worsened_samples=2;max_sample_worse=0.000249;switched_samples=6;adjusted_rows=130 | diagnostic_nonrelease; counts=0;F1=5.55e-17;L1=2.22e-16;Pearson=2.22e-16 | keep_feature_allocator_experimental_off_by_default |
| guarded feature allocator release candidate | wrapper parity plus external exact-split | samples=38;switched=36;improved=34;worsened=2;mean_L1_delta_pp=-0.269171;max_worse_L1_delta_pp=0.000249;mean_F1_delta=0 | plant_holdout3:delta_L1=0.000248910;base_multi=0.794009;rows=40;plant_holdout4:delta_L1=0.000060987;base_multi=0.762779;rows=33 | keep_feature_allocator_experimental_off_by_default |
| refined guarded feature allocator | combined output-only guard, wrapper parity, score replay | samples=38;switched=33;improved=33;worsened=0;mean_L1_delta_pp=-0.262955;selected_mean_L1_delta_pp=-0.311474;external_mean_L1_delta_pp=-0.004186;max_worse_L1_delta_pp=0.000000 | guard_match=max_abs_mean_delta=3.10619365229e-07; independent_holdout_missing | refined_guard_score_replay_pass_candidate_needs_independent_holdout |
| refined guarded feature allocator marine diagnostic | independent cached marine exact-split profiles; nonrelease GTDB transfer truth | samples=2;switched=2;improved=2;worsened=0;mean_L1_delta=-0.000521634;max_worse_L1=-0.000503001;mean_F1_delta=0.000000000 | samples=marine0,marine2;min_truth_mass_mapped_pct_bacteria_archaea=74.120670 | diagnostic_supports_refined_allocator_holdout_candidate |
| refined guarded feature allocator holdout inventory | release-grade manifest overlap and omitted HMP input check | release_panels=4;release_samples=32;selected_cache_overlap=32;independent_release_ready_samples=0 | samples=2,8,12,26,27;truth=5/5;minco_profile=0/5;sylph_profile=0/5;ready_pairs=0/5 | release_grade_independent_holdout_missing |
| refined guarded feature allocator recovery plan | ranked independent holdout recovery routes | post_recovery_scored=True;candidate_vs_sylph=minco_F1_wins=0/3;minco_L1_wins=0/3;refined_effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000 | extension_result_does_not_support_default_promotion | do_not_promote_refined_allocator_without_independent_release_holdout |
| CAMI3 source-profile extension diagnostic | samples3-5 local taxonomic-profile source rows | samples=3,4,5;min_mapped_profile_scope_pct=62.728832;ready_samples= | minco_F1_wins=0/3;minco_L1_wins=0/3 | keep_refined_allocator_opt_in |
| CAMI3 source-readmap recovery inputs | samples3-5 exact per-read truth recovery | candidate_profiles=3/3;raw_tables=3/3;sylph_profiles=3/3 |  | post_recovery_scoring_ready |

## Current Default

- Keep `scripts/minco_profile_default.py` default preset as `candidate`.
- Candidate preset components remain: `universal-auto-exact`,
  `emitted-ani90-xny100-br01-af70` rescue,
  `accession-ani90-xny100-br01-af70` surface, and
  `normalized-depth-alpha2` candidate abundance.
- Do not promote any additional abundance allocator yet. Every tested
  nonzero allocator family either has sample regressions, weak external
  validation, or a call-level F1 cost.
- This is a best-current-MinCO default, not a claim that MinCO broadly
  beats Sylph on abundance.

## Next Abundance Work

- Prioritize matched-call mass allocation, because that is the largest
  known gap on most cached panels and the selected-call oracle can
  close the selected-candidate L1 gap on 31/32 cached samples.
- Pair allocation work with targeted call recovery for the remaining
  sample where the selected call set is still insufficient.
- Avoid promoting panel-mean-only improvements unless they are also
  sample-safe on independent holdouts.
- The selected-call mass-transform sweep confirms that small
  unguarded base-row transforms are not sample-safe, even when panel
  means improve.
- The selected-call guard audit found in-panel sample-safe rules, but
  leave-one-panel-out still regressed held-out samples, so guarded
  transforms also remain diagnostic.
- The selected-call feature allocator sweep found useful
  panel-level signals from existing depth features, but no strict
  sample-safe variant, so it is not a default change.
- The corrected feature-allocator guard audit excludes cached
  F1/L1/Pearson score columns and is the first allocator branch to
  pass leave-one-panel-out without held-out sample regressions.
- The implemented guarded feature allocator switch reproduces the
  sample-safe cached signal when accession taxmap labels are used
  for genus grouping. The stricter release-candidate audit improves
  34/38 combined samples with unchanged F1, but two external
  diagnostic samples still have small L1 regressions.
- The refined `guarded-genus-hit-breadth-a002-xny230` switch adds an
  output-only split-support guard that removes those cached external
  regressions while preserving 97.7% of the selected-panel gain. Its
  wrapper guard parity and fixed-call score replay are validated, but
  it still needs release-grade independent holdout evidence before
  default promotion. A marine exact-split diagnostic replay is
  supportive, but its GTDB transfer truth is not release-grade.
- The independent-holdout inventory confirms the current cache has
  no unused release-grade profile set for this allocator: all
  release-grade manifest samples are already in the guard-selection
  cache, and omitted HMP IDs lack cached MinCO/Sylph profile pairs.
- The recovery plan ranks CAMI3 source-readmap samples3-5 as the
  best local route: selected-default replay profiles, raw tables,
  and baseline profiles are present for 3/3 samples; the remaining
  blocker is source-readmap truth cache.
- A source-profile fallback using local CAMI taxonomic profiles was
  tested for samples3-5 and rejected as a release substitute because
  in-scope GTDB truth mapping remains below threshold.
- The exact source-readmap recovery audit shows remote archive URLs
  are present in the local manifest and rescoring artifacts are ready,
  but the extracted source-readmap truth files for samples3-5 are
  missing. Anonymous reads are optional for this rescore because the
  MinCO and Sylph profiling artifacts already exist.
- A post-recovery CAMI3 extension scorer is now wired into the
  evidence path. It is currently blocked only by `reads_mapping.tsv.gz`
  files and will score selected-default MinCO, the refined allocator,
  and Sylph once those files are restored.

Machine-readable matrix: `results/abundance_strategy_decision_matrix.tsv`
