# virtual_s2000_marker_index_benchmark

Date: 2026-06-24
Author/agent: Codex
Project: KSSD3mini / MinCO metagenomic profiling

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Git metadata directory used: `/home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git`
- Commit: `b55b4d3f4c986de29098d2f1a092ee28251008aa`
- Ref/describe: `main`, `b55b4d3-dirty`
- Working tree: dirty; status and diffstat are recorded in `provenance/code_status.txt`.
- Main script: `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/benchmark_virtual_markers.py`

## Question

Have we already tested a full S2000 refdb where unique-marker contexts are marked virtually from the inverted index, instead of physically building a separate ctx-markerdb? If not, can the full S2000 inverted index reproduce the physical ctx-markerdb marker-size table exactly?

## Prior Status

The exact virtual marked-unique full-refdb mode had not been tested end-to-end inside MinCO. The closest previous experiment was `2026-06-23_hybrid_full_sketch_unique_marker_fallback`, which combined physical ctx-markerdb rows with full-refdb fallback rows. That is not the same as computing marker/shared status directly inside a full-refdb scan.

This experiment tests the critical prerequisite: whether the full S2000 sorted inverted index can reproduce the physical ctx-markerdb marker sizes exactly.

## Dataset

- Full S2000 dedup reference:
  `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Physical S2000 ctx-markerdb:
  `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Physical ctx-marker size table:
  `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv`
- Number of refs: `200,527`
- Full S2000 entries: `401,054,000`
- Full sorted inverted index:
  `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup/minco.refindex.ctxgid64obj32`
- Availability status: available in `/tmp` on 2026-06-24.
- Retention risk: temporary; preserve elsewhere if needed.

## Methods

The script memory-maps the full S2000 sorted inverted index. Each record is interpreted as:

```text
ctxgid: uint64, sorted by (ctx << 20) | gid
obj:    uint32
```

For each contiguous context run:

```text
ctx = ctxgid >> 20
gid = ctxgid & ((1 << 20) - 1)

if first_gid == last_gid:
    context is unique to exactly one reference
    marker_size[gid] += run_length
else:
    context is shared
```

Because the index is sorted by ctx then gid, `first_gid == last_gid` is sufficient to test whether the run contains only one distinct reference gid. The input sketch is non-conflict, so unique context runs have length 1 in this S2000 build.

The resulting virtual marker sizes were compared against:

1. Physical ctx-markerdb offsets from `sketch_T_S2000_aaf003_dedup_ctxmarker/minco.ctxobj64.offsets`.
2. Existing physical size table `markerdb_ctx_s2000_dedup.psmp.tsv`.

## Results

The virtual full-index marker computation exactly reproduced the physical ctx-markerdb marker sizes:

```text
refs                                  200,527
full entries                          401,054,000
virtual marker entries                239,132,175
physical marker entries               239,132,175
virtual vs physical mismatched refs   0
virtual vs psmp mismatched refs       0
```

Index scan summary:

```text
total ctx runs       281,442,561
unique ctx runs      239,132,175
shared ctx runs       42,310,386
unique ctx fraction   0.849666
max ctx run length    64,417
```

Marker-size distribution:

```text
min       11
p01       218
p05       501
p25       918
median   1253
p75      1506
p95      1715
p99      1797
max      1934
mean     1192.52
```

Low-marker references:

```text
marker_size < 500:  9,936 refs
marker_size < 200:  1,878 refs
marker_size = 0:        0 refs
```

Prototype runtime:

```text
full index records      401,054,000
chunk size               10,000,000 records
internal scan time        15.79 s
command wall time         33.80 s
peak RSS                   8.67 GB
```

The Python prototype is not memory-optimized. A C implementation in the existing readwise index scan should need only the current ctx run, per-reference counters, and optionally a small `ref_marker_size[gid]` array.

Storage implication:

```text
full S2000 refdb directory       8,146,213,104 bytes
physical ctx-markerdb directory  4,886,921,796 bytes
```

Using virtual markers avoids storing the extra physical ctx-markerdb copy for this build. A `uint32_t ref_marker_size[n_refs]` sidecar would be about `802,108` bytes for 200,527 refs. A per-entry marker flag would be about `401 MB`, but this benchmark shows it is not needed for correctness.

## Validation

- `python3 -m py_compile` passed for `benchmark_virtual_markers.py`.
- The first prototype scan showed mismatches that changed with chunk size; this identified a chunk-boundary carry bug.
- After fixing carry handling, the virtual marker sizes exactly matched both the physical ctx-markerdb offsets and `markerdb_ctx_s2000_dedup.psmp.tsv`.
- Final mismatch table has only a header row.

## Relationship To Previous Accuracy Results

This benchmark does not rerun readwise profiling. It proves that the full S2000 inverted index can reproduce the physical ctx-markerdb marker-size semantics exactly.

Therefore, once implemented inside MinCO, marker-only metrics from virtual marked-unique full S2000 should match the physical ctx-markerdb branch when using the same context/object and conflict semantics.

Existing physical ctx-marker/full-fallback results remain the accuracy reference:

```text
Mouse GTDB 0-2, hybrid baseline-plus:
  F1 = 0.944844
  L1 = 1.724130

CAMI3 ToyGut 0-2, hybrid XnY>=700:
  F1 = 0.555312
  L1 = 39.232850
```

## Conclusion

The virtual full S2000 marked-unique design is feasible and should replace a separate physical ctx-markerdb copy in the next implementation.

The recommended architecture is:

```text
Use full S2000 refdb as the only reference payload.
Use the sorted inverted index to classify each ctx hit as unique or shared.
Maintain marker and full counters separately during readwise analysis.
Keep a tiny ref_marker_size[gid] array for marker-breadth denominators and warnings.
Do not store a physical ctx-markerdb unless needed for backward-compatible benchmarking.
```

## Paper-Relevant Claim

The full S2000 inverted index is sufficient to reproduce physical ctx-markerdb unique-context marker sizes exactly, avoiding a duplicate 4.6 GB markerdb copy for this benchmark build.

## Caveats

- This is an index-equivalence and storage/runtime benchmark, not an integrated MinCO readwise benchmark.
- End-to-end F1/L1 must still be tested after implementing marker/full counters inside `ani`.
- The prototype scanned the whole index to compute all marker sizes. In readwise profiling, marker/shared classification should be done only for query-hit context runs.
- The Python prototype has high peak RSS because it uses large vectorized chunks; C can be much smaller.
- The reference data are under `/tmp`, so paths are temporary.

## Next Experiment

Implement virtual marker counters in MinCO `ani` over the full S2000 refdb:

```text
full_XnY, full_depth, full_breadth
marker_XnY, marker_depth, marker_breadth
marker_size
marker_ztp_af, marker_effective_depth
```

Then rerun the mouse and CAMI3 benchmark against physical ctx-markerdb outputs to verify identical marker-only calls and measure the cost of marker/full dual counters in the real readwise path.
