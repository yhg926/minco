# Abundance Next-Target Feasibility

Date: 2026-06-29

This note is generated from cached selected-candidate profile TSVs and
cached truth/comparator scores. It keeps the selected candidate call set
fixed, then computes a truth-aware abundance oracle over detected true
species. The oracle is not an implementable strategy; it is a feasibility
bound for deciding the next algorithmic target.

## Summary

| Panel | Samples | Selected-Sylph L1 pp | Oracle-Sylph L1 pp | Allocation headroom pp | Detected truth % | Already not worse | Allocation-only | Call recovery first | Residual work |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cami2_toy_mouse_gut | 3 | 6.798669 | -3.896734 | 10.695402 | 98.582605 | 0 | 3 | 0 | 0 |
| cami3_toy_human_gut_gtdb_source_readmap | 3 | 25.114394 | -16.281277 | 41.395671 | 91.824218 | 0 | 3 | 0 | 0 |
| hmp_airskin_gtdb_source_abundance | 24 | 17.064849 | -3.543929 | 20.608778 | 99.093592 | 0 | 23 | 1 | 0 |
| hmp_gastrooral_gtdb_source_abundance | 2 | 41.294813 | -5.058279 | 46.353092 | 99.997575 | 0 | 2 | 0 | 0 |
| all | 32 | 18.371412 | -4.865778 | 23.237190 | 98.420682 | 0 | 31 | 1 | 0 |

## Decision

- Baseline validation max delta: `4.54747350886e-13`.
- Allocation-only can close the selected-candidate L1 gap on `31` of `32` cached samples.
- Call recovery is required first on `1` cached samples.
- Residual work remains after allocation on `0` cached samples.
- Therefore the next default-strategy work should not be another
  panel-mean-only abundance rescaling. It should combine a sample-safe
  matched-call allocator with high-confidence call recovery for samples
  whose detected truth mass is not enough.

## Outputs

- `results/candidate_callset_oracle_feasibility.tsv`
- `results/candidate_callset_oracle_feasibility_summary.tsv`
- `results/candidate_callset_oracle_feasibility_audit.tsv`
