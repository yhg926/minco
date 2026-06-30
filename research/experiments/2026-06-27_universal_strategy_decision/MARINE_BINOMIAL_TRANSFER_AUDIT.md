# Marine Binomial Transfer Audit

Date: 2026-06-29

This cached audit tests whether exact GTDB-binomial fallback can upgrade
CAMI II marine GTDB truth transfer. It uses gold-profile species labels
and cached GTDB metadata only; no profilers are rerun.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `samples_audited` | 10 | all_gold_samples |
| `profile_scored_samples` | 0,3,4,5 | cached_profile_subset |
| `fallback_mapped_mass_pct_all_sum` | 59.034500 | binomial_fallback_rescues_truth_mass |
| `all_samples_min_mapped_pct_bacteria_archaea` | 91.005806 | below_release_threshold |
| `scored_samples_min_mapped_pct_bacteria_archaea` | 91.005806 | below_release_threshold |
| `promotion_decision` | do_not_promote | truth_transfer_still_partial |

## Quality

| Sample | Mapped B/A % | Unique-taxid mass | Binomial fallback mass | Remaining ambiguous mass | Ready | Scored |
|---:|---:|---:|---:|---:|---|---|
| 0 | 91.006 | 20.3749 | 4.2955 | 2.4123 | false | true |
| 1 | 92.494 | 23.6389 | 7.5347 | 2.4129 | false | false |
| 2 | 93.515 | 22.1520 | 5.7964 | 1.8834 | false | false |
| 3 | 92.573 | 23.4667 | 5.8918 | 2.2382 | false | true |
| 4 | 94.088 | 22.0485 | 7.0171 | 1.7400 | false | true |
| 5 | 94.096 | 28.7660 | 6.3444 | 2.0668 | false | true |
| 6 | 94.525 | 27.6543 | 6.3855 | 1.7794 | false | false |
| 7 | 93.169 | 23.8250 | 5.5261 | 2.0559 | false | false |
| 8 | 93.734 | 24.4522 | 4.6188 | 1.8088 | false | false |
| 9 | 92.245 | 24.7979 | 5.6242 | 2.3634 | false | false |

## Decision

- Exact-binomial fallback is useful for marine truth transfer.
- Do not promote the marine panel unless both truth coverage is high
  and selected-default MinCO/Sylph profiles are available in the same
  benchmark namespace.

## Outputs

- `results/marine_binomial_transfer_truth.tsv`
- `results/marine_binomial_transfer_quality.tsv`
- `results/marine_binomial_transfer_detail.tsv`
- `results/marine_binomial_transfer_audit.tsv`
