# Achievable Release Goal

Date: 2026-06-30

This note records the narrowed goal that can be achieved from the current evidence: make MinCO profiling release-ready as a simple, documented default with clear claim boundaries. It does not claim universal superiority over external profilers.

## Goal

MinCO is release-ready when a user can run one default command, validate a packaged reference with preflight, and rely on documented defaults that are supported by cached cross-panel evidence without manually choosing research modes.

## Decision

- Status: `pass`
- Decision: `achievable_release_goal_supported`
- Boundary: the original broad universal/Sylph-beating research goal remains incomplete; this release goal is narrower and user-facing.

## Checks

| Check | Status | Observed | Decision |
| --- | --- | --- | --- |
| goal_scope | `pass` | default-ready goal; broad Sylph-beating claim stays rejected | achievable_goal_defined |
| simple_entrypoint | `pass` | scripts/minco_profile -r ref.minco --reads reads.fq.gz -p16 -o calibrated.profile.tsv | user_can_run_scripts_minco_profile |
| preflight | `pass` | scripts/minco_profile --check-ref -r ref.minco --reads reads.fq.gz | preflight_documented_and_manifested |
| packaged_defaults | `pass` | optional minco_profile_defaults.tsv sidecar supports whitelisted domain defaults such as scope=bacteria | sidecar_defaults_supported |
| selected_default | `pass` | preset=candidate;strategy=universal-auto-exact | candidate_preset_is_release_default |
| current_minco_support | `pass` | F1_regressions=0;L1_regressions=0;Pearson_regressions=0 | selected_default_supported_vs_previous_minco |
| external_claim_boundary | `pass` | F1_wins=1/4;L1_wins=0/4;Pearson_wins=0/4 | claim_boundary_is_explicit |
| expected_gaps_accounted | `pass` | fail=0;expected_gap=2 | expected_gaps_are_scope_boundaries |
| original_goal_boundary | `pass` | goal_not_complete_expected_gaps_remain | do_not_mark_broad_goal_complete |
| summary | `pass` | pass=9;fail=0 | achievable_release_goal_supported |

Machine-readable output: `results/achievable_release_goal.tsv`.
