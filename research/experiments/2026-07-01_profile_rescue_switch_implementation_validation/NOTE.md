# Profile Rescue Switch Implementation Validation

## Question

Does the implemented `--candidate-rescue-switch
split-p002-x300-ani095-af06-b025-d1-top1` reproduce the fresh offline
high-specificity rescue result when run through the real calibrated profile
wrapper?

## Methods

The validation reuses saved unique/split tables from the fresh raw-run
validation in `/tmp/minco_profile_rescue_fresh_raw_validation_20260701/work`.
This avoids another raw-read scan while still exercising
`scripts/minco_profile_calibrated.py` output generation.

Both outputs use:

- `--strategy universal-auto-exact`
- `--scope bacteria`
- `--abundance-ani-floor 0.90`
- `--abundance-sparse-depth-cap poisson-breadth`
- `--abundance-sparse-breadth-max 0.15`
- `--abundance-sparse-depth-ratio-min 200`
- `--report-all`

Compared methods:

- `current_calibrated_abundance`: `--candidate-rescue-switch off`
- `profile_rescue_p002_x300_ani095_af06_b025_d1_top1`:
  `--candidate-rescue-switch split-p002-x300-ani095-af06-b025-d1-top1`

## Status

Run `run_switch_validation.py` from the repository root to regenerate the
summary tables. Large profile outputs are written under
`/tmp/minco_profile_rescue_switch_implementation_validation_20260701`.

## Results

The implemented switch reproduced the fresh offline validation exactly on the
available four-sample panel:

| Panel | Samples | Current pooled F1 | Switch pooled F1 | Current L1 pp | Switch L1 pp | Added truth | Added FP |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CAMI3 toy human gut | 0,1,2 | 0.814570 | 0.826230 | 58.155186 | 39.708920 | 6 | 0 |
| HMP air/skin | 13 | 0.925000 | 0.925000 | 54.747984 | 54.747984 | 0 | 0 |

Overall mean official L1 improved from `56.451585` to `47.228452` pp and mean
official Pearson improved from `0.877137` to `0.913351`. The maximum panel L1
worsening was `0.0` pp.

The switch added six called rows, all CAMI3 truth species:

- sample0: `s__Blautia_A faecis`, `s__Faecalibacterium taiwanense`,
  `s__Bacteroides sp947646015`
- sample1: `s__Faecalibacterium taiwanense`, `s__Blautia_A faecis`,
  `s__Bacteroides sp947646015`

Table-based wrapper runtime was `8.53-11.00` seconds per profile, with peak RSS
`1.40-1.57` GB. This does not include the original raw-read scan time, which was
measured in the fresh raw-run validation.

## Conclusion

The implementation is correct for the validated rule: real wrapper output
matches the offline rescue result, rescued rows retain native profile raw
abundance mass, and the default remains unchanged unless the user explicitly
passes `--candidate-rescue-switch split-p002-x300-ani095-af06-b025-d1-top1`.
