# CAMI3 Source-Readmap Recovery Inputs

Date: 2026-06-29

This audit records the exact local and remote inputs needed to recover
per-read source mapping for CAMI3 samples3-5. It does not download or
extract data.

## Audit

| Metric | Value | Decision |
|---|---:|---|
| `remote_manifest_urls` | 3/3 | download_urls_available |
| `local_archives` | 0/3 | archives_missing_locally |
| `local_readmaps` | 3/3 | readmaps_present |
| `local_anonymous_reads` | 0/3 | not_required_for_truth_rescore |
| `rescoring_artifacts` | candidate_profiles=3/3;raw_tables=3/3;sylph_profiles=3/3 | rescoring_artifacts_ready |
| `compressed_archive_size_estimate` | mean_existing_archive_bytes=4864188990;samples3_5_estimated_bytes=14592566970 | size_estimate_from_samples0_2 |
| `readmap_only_recovery_plan` | commands=3/3 | stream_extract_readmaps_without_storing_archives |
| `recovery_blockers` |  | no_recovery_blockers |
| `promotion_decision` | source_readmap_truth_recovery_inputs_ready | post_recovery_scoring_ready |

## Detail

| Sample | URL present | Archive | Readmap | Reads optional | Candidate profile | Raw tables | Baseline profile |
|---:|---|---|---|---|---|---|---|
| 3 | True | False | True | True | True | True | True |
| 4 | True | False | True | True | True | True | True |
| 5 | True | False | True | True | True | True | True |

## Decision

- Samples3-5 now have extracted `reads_mapping.tsv.gz` files,
  selected-default profiles, raw tables, and baseline profiles ready
  for rescoring.
- Anonymous reads are not needed for this truth-rescore path because
  profiling artifacts already exist.
- The generated command file streams each archive URL through `tar` and
  extracts only `reads_mapping.tsv.gz`, avoiding storage of the full
  sample read archives.
- The source-profile fallback has already been tested and is not a release
  substitute.

## Outputs

- `results/cami3_source_readmap_recovery_inputs.tsv`
- `results/cami3_source_readmap_recovery_inputs_audit.tsv`
- `results/cami3_source_readmap_readmap_only_recovery_commands.sh`

Promotion decision is `post_recovery_scoring_ready`.
