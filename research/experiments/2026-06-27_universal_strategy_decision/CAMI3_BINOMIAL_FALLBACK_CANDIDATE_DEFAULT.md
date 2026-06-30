# CAMI3 Binomial Fallback Candidate Default

Date: 2026-06-29

This cached-only audit scores the selected `candidate` default MinCO
profiles against the CAMI3 exact-binomial fallback truth table. It
uses cached profile outputs and does not change MinCO calls.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `candidate_profiles_scored` | 3 | selected_default_profile_set |
| `minco_vs_sylph_mean_F1` | 0.862183 vs 0.762024 | minco_higher_F1 |
| `minco_vs_sylph_mean_L1` | 55.919344 vs 31.103234 | sylph_lower_L1 |
| `minco_vs_sylph_mean_Pearson` | 0.866593 vs 0.947203 | sylph_higher_Pearson |
| `candidate_delta_vs_prior_minco_rescore` | mean_F1_delta=0.009125;mean_L1_delta=-5.964824;mean_Pearson_delta=0.040308 | candidate_default_rescore_replaces_prior_minco_rescore |
| `promotion_decision` | release_evidence_candidate_if_truth_policy_accepts_binomial_fallback | truth_policy_is_remaining_gate |

## Summary Scores

| Method | Mean F1 | Pooled F1 | Mean L1 pp | Mean Pearson |
|---|---:|---:|---:|---:|
| minco_candidate_default_binomial_fallback | 0.862183 | 0.860759 | 55.919344 | 0.866593 |
| sylph_binomial_fallback | 0.762024 | 0.762617 | 31.103234 | 0.947203 |

## Delta Versus Prior MinCO Rescore

| Sample | Method | Delta F1 | Delta L1 pp | Delta Pearson |
|---:|---|---:|---:|---:|
| 0 | minco_candidate_default_binomial_fallback | 0.011338 | -1.387836 | 0.006011 |
| 0 | sylph_binomial_fallback | 0.000000 | 0.000000 | -0.000000 |
| 1 | minco_candidate_default_binomial_fallback | 0.010703 | -8.471789 | 0.071252 |
| 1 | sylph_binomial_fallback | 0.000000 | 0.000000 | 0.000000 |
| 2 | minco_candidate_default_binomial_fallback | 0.005333 | -8.034847 | 0.043659 |
| 2 | sylph_binomial_fallback | 0.000000 | 0.000000 | -0.000000 |

## Decision

- This is the CAMI3 binomial-fallback score that matches the selected
  no-expertise MinCO default.
- It can replace the older MinCO rescore only if the exact-binomial
  fallback truth policy is accepted.
- It supports MinCO on F1 for this panel, but abundance remains weaker
  than Sylph.

## Outputs

- `results/cami3_binomial_fallback_candidate_default_scores.tsv`
- `results/cami3_binomial_fallback_candidate_default_summary.tsv`
- `results/cami3_binomial_fallback_candidate_default_delta.tsv`
- `results/cami3_binomial_fallback_candidate_default_audit.tsv`
