# CAMI3 Source-Readmap Extension Cache Audit

Date: 2026-06-29

This audit checks whether the accepted source-readmap release panel can
be extended from local cached files alone. It only checks file presence;
it does not rerun profilers or infer release truth from weaker tables.

This is now a historical pre-recovery cache snapshot. The current route state
is recorded in `CAMI3_SOURCE_READMAP_EXTENSION_AFTER_RECOVERY.md`: readmaps are
present for samples3-5, the extension has been scored, and the result is
negative for refined-allocator promotion.

## Result

- Accepted subset ready: `3/3`.
- Extension subset release-ready: `0/3`.
- Extension blockers: `missing_source_readmap_truth`.

Older profile outputs, raw best-reference tables, baseline profiles,
and selected-default replay outputs exist for the later samples. The
historical cache-only blocker was source-readmap truth cache in the same
namespace (`requires_source_readmap_truth_cache`; selected replay status:
`selected_default_cache_present_for_extension`). That blocker was later removed
by readmap-only stream extraction.

## Detail

| Sample | Scope | Truth Cache | Selected Default Profile | Older Profile | Raw Tables | Baseline Profile | Ready |
|---:|---|---|---|---|---|---|---|
| 0 | accepted_release_subset | True | True | True | True | True | True |
| 1 | accepted_release_subset | True | True | True | True | True | True |
| 2 | accepted_release_subset | True | True | True | True | True | True |
| 3 | extension_candidate_subset | False | True | True | True | True | False |
| 4 | extension_candidate_subset | False | True | True | True | True | False |
| 5 | extension_candidate_subset | False | True | True | True | True | False |

## Outputs

- `results/cami3_source_readmap_extension_cache_detail.tsv`
- `results/cami3_source_readmap_extension_cache_audit.tsv`
