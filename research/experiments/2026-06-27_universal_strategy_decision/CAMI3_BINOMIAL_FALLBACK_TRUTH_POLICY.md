# CAMI3 Binomial Fallback Truth Policy

Date: 2026-06-29

This cached audit decides whether the exact GTDB-binomial fallback is
acceptable as CAMI3 source-readmap GTDB-species truth. It does not rerun
profiling and does not change MinCO calls.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `fallback_rows_reviewed` | 51 | nonzero_fallback_set |
| `fallback_added_rows` | 1732020 | truth_rows_added_by_policy |
| `deterministic_unique_binomial` | 51/51 | pass |
| `source_mapping_priority` | priority_conflicts=0;same_species_remainders=2 | pass |
| `profile_scope_coverage` | ready_samples=0,1,2;min_adjusted_pct=95.941 | pass |
| `selected_default_scored` | 3 | pass |
| `truth_policy_decision` | accept_exact_binomial_fallback_for_release_truth | release_truth_policy_pass |

## Detail Summary

| Fallback rows | Added rows | Unique-binomial rows | Prior-mapped violations |
|---:|---:|---:|---:|
| 51 | 1732020 | 51 | 0 |

## Decision

- Accept the exact GTDB-binomial fallback as clean CAMI3
  source-readmap GTDB-species truth if all audit checks pass.
- This policy fills rows left unmapped by the stronger WGS-prefix
  and unique-taxid transfer routes. If a source already mapped by
  WGS/unique-taxid, the fallback is allowed only when it assigns the
  unmapped remainder to that same GTDB species.
- The policy does not imply MinCO abundance beats Sylph; it only decides
  whether the CAMI3 panel is clean enough to use as release evidence.

## Outputs

- `results/cami3_binomial_fallback_truth_policy_detail.tsv`
- `results/cami3_binomial_fallback_truth_policy_audit.tsv`
