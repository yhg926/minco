# CAMI II Toy Mouse Density Block Argument Check

## Question

Does changing MinCO readwise `--density-block-ctx` from the default `100` to
`0` or `1` change the current best ctx-only S2000 markerdb result on CAMI II Toy
Mouse GTDB-species scoring?

## Setup

Reference:

- `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`

Profiler settings:

- `--query-density ref`
- `--abundance-est depth`
- `--readwise-profile-only`
- `--readwise-assign best-diff-split`
- `--readwise-ani naive`
- `--readwise-ctx-filter product-topfrac-median`
- `--readwise-fake-threshold 0.25`
- active product0 gate from the Toy Mouse GTDB benchmark

`--density-block-ctx 100` is the existing default and batches consecutive reads
until at least 100 retained density contexts. `--density-block-ctx 0` forces
exact per-read lookup. In the current implementation, block mode is active only
when `N > 1`, so `N=1` follows the same path as `N=0`.

## Results

`N=0` and `N=1` on sample0 produced byte-identical outputs:

```text
bf5f2d44d57520ecdefa52d883ea3437c3855738fca1fd05cbfb86a4cadb95b3  sample0_ctxmarker_block0.tsv
bf5f2d44d57520ecdefa52d883ea3437c3855738fca1fd05cbfb86a4cadb95b3  sample0_ctxmarker_block1.tsv
```

GTDB-species presence/absence:

| sample | variant | TP | FP | FN | F1 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 0 | N100 default | 58 | 0 | 7 | 0.943089 |
| 0 | N0 per-read | 57 | 0 | 8 | 0.934426 |
| 1 | N100 default | 70 | 2 | 10 | 0.921053 |
| 1 | N0 per-read | 69 | 2 | 11 | 0.913907 |
| 2 | N100 default | 66 | 0 | 7 | 0.949640 |
| 2 | N0 per-read | 66 | 0 | 7 | 0.949640 |

Pooled totals:

| variant | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| N100 default | 194 | 2 | 24 | 0.937198 |
| N0 per-read | 192 | 2 | 26 | 0.932039 |

Mean abundance metrics over samples 0-2:

| variant | renorm | Pearson | Spearman | MAE pct points | L1 pct points |
| --- | --- | ---: | ---: | ---: | ---: |
| N100 default | no | 0.988697 | 0.883994 | 0.486323 | 35.540141 |
| N0 per-read | no | 0.988685 | 0.885508 | 0.485777 | 35.500466 |
| N100 default | yes | 0.988697 | 0.883994 | 0.230133 | 17.008147 |
| N0 per-read | yes | 0.988685 | 0.885508 | 0.230285 | 17.019358 |

Changed true-positive calls:

- sample0 loses `s__Bacteroides fragilis`: `Reliable_ztp_af` changes from
  `0.444906` at `N=100` to `0.379301` at `N=0`, falling below the active `0.40`
  floor.
- sample1 loses `s__Porphyromonas gingivalis`: `Reliable_ztp_af` changes from
  `1.0` at `N=100` to `0.141666` at `N=0`, also below the floor.
- sample2 call set is unchanged.

Runtime and memory were similar. `N=0` did not reduce memory materially; peak RSS
remained about 4.3-4.5 GB for ctx-only markerdb.

## Conclusion

`N=0` and `N=1` are equivalent for the current implementation. The exact
per-read mode is not better than the default `N=100` pseudo-read block mode on
this benchmark. It slightly improves raw abundance MAE but loses two true
species across samples 0-2, lowering pooled F1 from `0.937198` to `0.932039`.

Keep the default `--density-block-ctx 100` for the current active Toy Mouse
strategy.
