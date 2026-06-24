# Readwise Tracking Benchmark on Toy Mouse sample0

Date: 2026-06-24

## Question

Does the new `minco ani --readwise-track` sidecar output work on a full Toy Mouse gut sample, and does it make readwise profiling significantly slower?

## Inputs

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Code commit: `efa6fafdbd4df4521431555c4d5a41db438eff95`
- Binary: `/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/minco`
- Sample: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- Refdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Generated GTDB taxmap: `/tmp/minco_readwise_tracking_toymouse0_20260624/gtdb_r232_accession.taxmap.tsv`
- Output directory: `/tmp/minco_readwise_tracking_toymouse0_20260624`

## Methods

Four full-sample runs were compared:

1. Normal readwise profile without tracking, default `density_block_ctx=100`.
2. Exact per-read readwise profile without tracking, `--density-block-ctx 0`.
3. Tracking with taxonomy disabled.
4. Tracking with `--readwise-taxonomy gtdb` using a generated GTDB r232 accession taxmap.

The exact no-tracking output is the fair baseline for tracking overhead, because `--readwise-track` forces `--density-block-ctx 0` to preserve exact read IDs and context offsets.

## Results

See `summary.tsv` for the timing table.

The main profile output from both tracking runs was byte-identical to the exact no-tracking profile:

```text
baseline_exact.tsv == track_none.profile.tsv
baseline_exact.tsv == track_gtdb.profile.tsv
sha256 ab20ff23cb7c166a90e937936b182b0b4a3ad6b8d006199905e69a7e6fa8a615
```

The GTDB tracking summary reported:

```text
total_reads                         33170320
reads_with_density_ctx              26731332
reads_with_ref_hit                  1163752
reads_with_multi_ref_hit            7603
tracked_read_pct                    3.508413546
total_density_ctx                   54088053
total_matched_ctx                   1188657
total_selected_ctx                  1188657
total_selected_ref_events           1188657
sampled_ctx_ref_hit_pct             2.197633182
estimated_ref_absent_ctx_pct        97.80236682
```

The original run was made before the tracking summary was renamed from a
read-level no-hit diagnostic to a context-level absent-context metric. The
correct context-level calculation for this markerdb run is:

```text
estimated_ref_absent_ctx_pct =
  100 * (1 - 1188657 / 54088053) = 97.80236682
```

Because the tested refdb is a ctx-markerdb, this means absent from retained
marker contexts. It should not be interpreted as absent from the original whole
GTDB genomes.

Tracking table distribution:

```text
target_ref_count=1                  1156149 reads
target_ref_count=2                  7455 reads
target_ref_count=3                  148 reads
target_ref_count>=4                 0 reads
offset_lists_truncated              0 reads
```

GTDB LCA rank distribution:

```text
species                             1155872
superkingdom                        4539
genus                               1344
phylum                              822
root                                480
NA                                  292
family                              170
class                               154
order                               79
```

## Conclusion

The feature worked on the full Toy Mouse sample and did not alter the main profile output. Runtime overhead is meaningful but not catastrophic: tracking without taxonomy was 1.49x the exact per-read baseline, and tracking with GTDB LCA was 1.54x. Peak RSS was unchanged for no-taxonomy tracking and increased from 4.11 GB to 4.74 GB with the 377 MB GTDB taxmap loaded.

The observed slowdown is mostly from per-read offset extraction plus writing 1,163,752 sidecar rows. GTDB LCA adds only about 5.3 seconds and 0.63 GB RSS in this run.

## Caveats

- This is one full Toy Mouse sample only.
- The generated GTDB taxmap uses assembly accessions as stable identifiers because GTDB does not provide numeric species taxids in the same CAMI taxmap schema.
- The absent-context percentage in this run is markerdb absence, not whole-genome reference absence. Whole-genome absence requires a full context refdb or an exact whole-reference context membership index.
- Runtime order may benefit later runs from filesystem cache; the most relevant comparison is exact no-tracking versus tracking in the same run series.
