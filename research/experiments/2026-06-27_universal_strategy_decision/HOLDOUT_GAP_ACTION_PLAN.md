# Holdout Gap Action Plan

Date: 2026-06-29

This generated note turns the holdout component of the release-gate
expected gaps into concrete next actions. It uses cached manifests and
summaries only.

## Summary

- Holdout status: `release=4;diagnostic=0;all_release=false`
- Panel counts: `release=4;diagnostic=0;nonrelease=4`
- Release gate: `pass=15;fail=0;expected_gap=2`
- Abundance gate: `selected_default_L1_losses=4/4;posthoc_blend=do_not_promote_abundance_blend;oracle=prioritize_allocator_plus_call_recovery;feature_allocator=do_not_promote_until_external_sample_safe;refined_score_replay=refined_guard_score_replay_pass_candidate_needs_independent_holdout;marine_diagnostic=diagnostic_supports_refined_allocator_holdout_candidate;holdout_inventory=release_grade_independent_holdout_missing;recovery_plan=do_not_promote_refined_allocator_without_independent_release_holdout;source_profile_extension=keep_refined_allocator_opt_in;source_readmap_inputs=post_recovery_scoring_ready` (abundance_expected_gap_not_release_claim)

## Priority Actions

| Rank | Panel | Grade | Samples | Action |
|---:|---|---|---:|---|
| 50 | mixed_readiness_26 | nonrelease | 26 | Do not use for release claims. Rebuild this dataset with clean GTDB species truth and same-release MinCO/Sylph profiles, or keep it diagnostic only. |
| 60 | cami3_toy_human_gut_0_5 | nonrelease | 6 | Do not repeat this completed CAMI3 extension unless inputs or strategy change. Keep the samples0-2 accepted truth policy as release evidence, treat samples3-5 as a negative holdout for allocator promotion, and use another independent route for the next release-gap closure. Post-recovery extension is scored for samples3-5: minco_F1_wins=0/3;minco_L1_wins=0/3; allocator effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000. |
| 60 | cami3_toy_human_gut_0_5_gtdb_taxid_transfer | nonrelease | 6 | Do not repeat this completed CAMI3 extension unless inputs or strategy change. Keep the samples0-2 accepted truth policy as release evidence, treat samples3-5 as a negative holdout for allocator promotion, and use another independent route for the next release-gap closure. Post-recovery extension is scored for samples3-5: minco_F1_wins=0/3;minco_L1_wins=0/3; allocator effect=improved=0;worsened=0;mean_L1_delta_pp=0.000000000. |
| 65 | cami2_marine_0_3_5_gtdb_taxid_transfer | nonrelease | 4 | Marine truth mapping and selected-default profile scoring are complete for the available same-namespace samples, and the result is negative for default promotion. Do not repeat this marine route unless the strategy or input profiles change. Exact-binomial fallback audit rescues truth mass but remains below release threshold: all-sample min mapped B/A %=91.005806, scored-sample min=91.005806; decision=truth_transfer_still_partial. Candidate-rule audit also fails: rescued_rows=51;rescued_mass_pct_all=3.146300;min_mapped_pct_bacteria_archaea=91.292431;ready_samples=2/10; decision=truth_ambiguity_not_resolved_by_current_candidate_rules. Source-readmap feasibility audit sampled 6000000 rows and 314110 unresolved rows, but found no direct source-sequence-to-GTDB assembly matches; decision=contig_or_otu_to_assembly_mapping_missing. Setup metadata plus local assembly-summary strict source rule clears the truth threshold, but same-namespace selected-default profiles are still needed: min_mapped_pct_bacteria_archaea=95.029253;ready_samples=10/10; decision=truth_mapping_threshold_cleared_but_profiles_need_same_namespace_rescore. Selected-default same-namespace profiles are now scored and do not support promotion: minco_mean_F1=0.781358;sylph_mean_F1=0.838716;minco_mean_L1_union_pp=41.299694;sylph_mean_L1_union_pp=40.997305;minco_mean_Pearson_union=0.886706;sylph_mean_Pearson_union=0.923667; decision=route_profile_evidence_ready_for_review. |
| 75 | cami2_toy_mouse_gut_5_7 | release | 3 | Keep as clean GTDB species release evidence; useful for confirming F1 regressions when candidate call/rescue rules change. |
| 80 | hmp_airskin_source_abundance_24sample | release | 24 | Release panel is refreshed to the 24-sample same-release r232 source-abundance run. Use it as a strong counterexample/accuracy stress panel, not as a Sylph-beating claim. |
| 85 | hmp_gastrooral_source_abundance_0_6 | release | 2 | Keep as a small but high-quality same-release source-abundance counterexample; use it to test abundance improvements before promotion. |
| 90 | cami3_toy_human_gut_0_2_gtdb_source_readmap | release | 3 | Keep as release-grade evidence and rerun only after default-strategy changes. |

## Decision

- The current default remains the candidate preset.
- The abundance expected gap is tracked separately in
  `ABUNDANCE_RELEASE_BLOCKER.md`; this plan only ranks
  holdout/truth/profile coverage actions.
- The next evidence target is not another local threshold sweep; it is
  cleaner release-grade GTDB holdout coverage in panel families that
  still lack high-coverage GTDB truth or selected-default profiles.
- The 24-sample HMP airskin release panel is now reflected in the
  manifest, but it is one dataset family and remains a counterexample
  where the external baseline is stronger for F1 and abundance.

## Outputs

- `results/holdout_gap_action_plan.tsv`
- `results/holdout_gap_action_plan_audit.tsv`
