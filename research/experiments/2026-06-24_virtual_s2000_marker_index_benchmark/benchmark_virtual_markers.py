#!/usr/bin/env python3
"""Benchmark virtual unique-marker classification from a full S2000 ref index.

The physical context markerdb keeps full ctxobj entries whose context occurs in
exactly one reference. This script computes the same marker size per reference
directly from the full sorted inverted index, without building a markerdb copy.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


GID_BITS = 20
GID_MASK = (1 << GID_BITS) - 1
CTXGIDOBJ_DTYPE = np.dtype([("ctxgid", "<u8"), ("obj", "<u4")])


def load_offsets(path: Path) -> np.ndarray:
    arr = np.fromfile(path, dtype="<u8")
    if arr.ndim != 1 or len(arr) < 2:
        raise ValueError(f"bad offsets file: {path}")
    return arr


def load_psmp(path: Path) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t", header=None, names=["psmp_size", "sample"])
    rows["psmp_size"] = pd.to_numeric(rows["psmp_size"], errors="raise").astype(np.uint64)
    return rows


def finalize_runs(
    marker_entries: np.ndarray,
    unique_run_lengths: list[int],
    first_gids: np.ndarray,
    last_gids: np.ndarray,
    run_lengths: np.ndarray,
) -> None:
    if len(first_gids) == 0:
        return
    unique = first_gids == last_gids
    if not bool(unique.any()):
        return
    gids = first_gids[unique].astype(np.int64, copy=False)
    lens = run_lengths[unique].astype(np.uint64, copy=False)
    counts = np.bincount(gids, weights=lens, minlength=len(marker_entries))
    marker_entries += counts.astype(np.uint64)
    unique_run_lengths.extend(lens.astype(np.int64).tolist())


def compute_virtual_marker_sizes(index_path: Path, n_refs: int, chunk_records: int) -> tuple[np.ndarray, dict[str, object]]:
    file_size = index_path.stat().st_size
    if file_size % CTXGIDOBJ_DTYPE.itemsize != 0:
        raise ValueError(f"{index_path} size {file_size} is not a multiple of {CTXGIDOBJ_DTYPE.itemsize}")
    n_records = file_size // CTXGIDOBJ_DTYPE.itemsize
    index = np.memmap(index_path, dtype=CTXGIDOBJ_DTYPE, mode="r", shape=(n_records,))

    marker_entries = np.zeros(n_refs, dtype=np.uint64)
    unique_run_lengths: list[int] = []
    total_ctx_runs = 0
    shared_ctx_runs = 0
    max_ctx_run_len = 0
    max_shared_gid_span = 0

    carry_ctx: int | None = None
    carry_first_gid = 0
    carry_last_gid = 0
    carry_len = 0

    start_time = time.time()
    for start in range(0, n_records, chunk_records):
        stop = min(start + chunk_records, n_records)
        ctxgid = index[start:stop]["ctxgid"]
        ctx = ctxgid >> GID_BITS
        gid = ctxgid & GID_MASK
        if len(ctx) == 0:
            continue

        run_starts = np.empty(len(ctx), dtype=np.int64)
        run_starts[0] = 0
        n_starts = 1
        changed = np.flatnonzero(ctx[1:] != ctx[:-1]) + 1
        if len(changed):
            run_starts[1 : len(changed) + 1] = changed
            n_starts += len(changed)
        run_starts = run_starts[:n_starts]
        run_ends = np.empty_like(run_starts)
        run_ends[:-1] = run_starts[1:]
        run_ends[-1] = len(ctx)

        run_ctx = ctx[run_starts].astype(np.uint64, copy=False)
        run_first_gid = gid[run_starts].astype(np.uint64, copy=True)
        run_last_gid = gid[run_ends - 1].astype(np.uint64, copy=True)
        run_len = (run_ends - run_starts).astype(np.uint64, copy=True)

        if carry_ctx is not None:
            if int(run_ctx[0]) == carry_ctx:
                run_first_gid[0] = carry_first_gid
                run_len[0] += np.uint64(carry_len)
            else:
                finalize_runs(
                    marker_entries,
                    unique_run_lengths,
                    np.array([carry_first_gid], dtype=np.uint64),
                    np.array([carry_last_gid], dtype=np.uint64),
                    np.array([carry_len], dtype=np.uint64),
                )
                total_ctx_runs += 1
                if carry_first_gid != carry_last_gid:
                    shared_ctx_runs += 1
                    max_shared_gid_span = max(max_shared_gid_span, carry_last_gid - carry_first_gid + 1)
                max_ctx_run_len = max(max_ctx_run_len, carry_len)

        # Hold the final run as carry because it may continue in the next chunk.
        complete_stop = len(run_starts) - 1
        if complete_stop > 0:
            fgid = run_first_gid[:complete_stop]
            lgid = run_last_gid[:complete_stop]
            lens = run_len[:complete_stop]
            finalize_runs(marker_entries, unique_run_lengths, fgid, lgid, lens)
            total_ctx_runs += len(fgid)
            shared = fgid != lgid
            shared_ctx_runs += int(shared.sum())
            if len(lens):
                max_ctx_run_len = max(max_ctx_run_len, int(lens.max()))
            if bool(shared.any()):
                spans = (lgid[shared] - fgid[shared] + 1).astype(np.uint64)
                max_shared_gid_span = max(max_shared_gid_span, int(spans.max()))

        last_i = len(run_starts) - 1
        carry_ctx = int(run_ctx[last_i])
        carry_first_gid = int(run_first_gid[last_i])
        carry_last_gid = int(run_last_gid[last_i])
        carry_len = int(run_len[last_i])

    if carry_ctx is not None:
        finalize_runs(
            marker_entries,
            unique_run_lengths,
            np.array([carry_first_gid], dtype=np.uint64),
            np.array([carry_last_gid], dtype=np.uint64),
            np.array([carry_len], dtype=np.uint64),
        )
        total_ctx_runs += 1
        if carry_first_gid != carry_last_gid:
            shared_ctx_runs += 1
            max_shared_gid_span = max(max_shared_gid_span, carry_last_gid - carry_first_gid + 1)
        max_ctx_run_len = max(max_ctx_run_len, carry_len)

    elapsed = time.time() - start_time
    unique_ctx_runs = len(unique_run_lengths)
    stats = {
        "index_path": str(index_path),
        "index_bytes": file_size,
        "index_records": int(n_records),
        "chunk_records": int(chunk_records),
        "elapsed_seconds": elapsed,
        "records_per_second": float(n_records / elapsed) if elapsed > 0 else None,
        "n_refs": int(n_refs),
        "total_ctx_runs": int(total_ctx_runs),
        "unique_ctx_runs": int(unique_ctx_runs),
        "shared_ctx_runs": int(shared_ctx_runs),
        "unique_ctx_fraction": float(unique_ctx_runs / total_ctx_runs) if total_ctx_runs else None,
        "max_ctx_run_len": int(max_ctx_run_len),
        "max_shared_gid_span": int(max_shared_gid_span),
        "unique_run_len_min": int(np.min(unique_run_lengths)) if unique_run_lengths else 0,
        "unique_run_len_max": int(np.max(unique_run_lengths)) if unique_run_lengths else 0,
    }
    return marker_entries, stats


def summarize_sizes(full_size: np.ndarray, virtual_size: np.ndarray, physical_size: np.ndarray) -> dict[str, object]:
    def quantiles(arr: np.ndarray) -> dict[str, float]:
        return {
            "min": float(np.min(arr)),
            "p01": float(np.quantile(arr, 0.01)),
            "p05": float(np.quantile(arr, 0.05)),
            "p25": float(np.quantile(arr, 0.25)),
            "median": float(np.median(arr)),
            "p75": float(np.quantile(arr, 0.75)),
            "p95": float(np.quantile(arr, 0.95)),
            "p99": float(np.quantile(arr, 0.99)),
            "max": float(np.max(arr)),
            "mean": float(np.mean(arr)),
        }

    ratio = np.divide(
        virtual_size.astype(float),
        full_size.astype(float),
        out=np.zeros(len(full_size), dtype=float),
        where=full_size > 0,
    )
    return {
        "n_refs": int(len(full_size)),
        "full_entries_total": int(full_size.sum()),
        "virtual_marker_entries_total": int(virtual_size.sum()),
        "physical_marker_entries_total": int(physical_size.sum()),
        "virtual_vs_physical_mismatched_refs": int(np.count_nonzero(virtual_size != physical_size)),
        "virtual_vs_physical_max_abs_diff": int(np.max(np.abs(virtual_size.astype(np.int64) - physical_size.astype(np.int64)))),
        "refs_marker_lt_500": int(np.count_nonzero(virtual_size < 500)),
        "refs_marker_lt_200": int(np.count_nonzero(virtual_size < 200)),
        "refs_marker_eq_0": int(np.count_nonzero(virtual_size == 0)),
        "marker_size": quantiles(virtual_size),
        "full_size": quantiles(full_size),
        "marker_to_full_ratio": quantiles(ratio),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-dir", type=Path, default=Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup"))
    parser.add_argument("--ctxmarker-dir", type=Path, default=Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker"))
    parser.add_argument("--psmp", type=Path, default=Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results"))
    parser.add_argument("--chunk-records", type=int, default=20_000_000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    full_offsets = load_offsets(args.full_dir / "minco.ctxobj64.offsets")
    physical_offsets = load_offsets(args.ctxmarker_dir / "minco.ctxobj64.offsets")
    if len(full_offsets) != len(physical_offsets):
        raise ValueError("full and ctxmarker offsets have different ref counts")
    n_refs = len(full_offsets) - 1
    full_size = np.diff(full_offsets).astype(np.uint64)
    physical_size = np.diff(physical_offsets).astype(np.uint64)

    virtual_size, index_stats = compute_virtual_marker_sizes(
        args.full_dir / "minco.refindex.ctxgid64obj32",
        n_refs,
        args.chunk_records,
    )
    psmp = load_psmp(args.psmp)
    if len(psmp) != n_refs:
        raise ValueError(f"psmp rows {len(psmp)} != refs {n_refs}")
    psmp_size = psmp["psmp_size"].to_numpy(dtype=np.uint64)

    detail = pd.DataFrame(
        {
            "sample_id": np.arange(n_refs, dtype=np.uint32),
            "sample": psmp["sample"],
            "full_size": full_size,
            "virtual_marker_size": virtual_size,
            "physical_marker_size": physical_size,
            "psmp_marker_size": psmp_size,
            "virtual_minus_physical": virtual_size.astype(np.int64) - physical_size.astype(np.int64),
            "virtual_minus_psmp": virtual_size.astype(np.int64) - psmp_size.astype(np.int64),
            "marker_to_full_ratio": np.divide(
                virtual_size.astype(float),
                full_size.astype(float),
                out=np.zeros(n_refs, dtype=float),
                where=full_size > 0,
            ),
        }
    )
    detail.to_csv(args.outdir / "virtual_marker_sizes.tsv", sep="\t", index=False)

    summary = summarize_sizes(full_size, virtual_size, physical_size)
    summary.update(
        {
            "psmp_mismatched_refs": int(np.count_nonzero(virtual_size != psmp_size)),
            "psmp_max_abs_diff": int(np.max(np.abs(virtual_size.astype(np.int64) - psmp_size.astype(np.int64)))),
        }
    )
    with open(args.outdir / "virtual_marker_summary.json", "w") as fp:
        json.dump({"index_stats": index_stats, "size_summary": summary}, fp, indent=2, sort_keys=True)

    flat_rows = []
    for key, value in index_stats.items():
        flat_rows.append({"section": "index", "metric": key, "value": value})
    for key, value in summary.items():
        if isinstance(value, dict):
            for subkey, subvalue in value.items():
                flat_rows.append({"section": key, "metric": subkey, "value": subvalue})
        else:
            flat_rows.append({"section": "size", "metric": key, "value": value})
    pd.DataFrame(flat_rows).to_csv(args.outdir / "virtual_marker_summary.tsv", sep="\t", index=False)

    below = detail.loc[detail["virtual_marker_size"] < 500].sort_values(["virtual_marker_size", "sample"])
    below.head(200).to_csv(args.outdir / "low_marker_refs_head200.tsv", sep="\t", index=False)
    mismatches = detail.loc[(detail["virtual_marker_size"] != detail["physical_marker_size"]) | (detail["virtual_marker_size"] != detail["psmp_marker_size"])]
    mismatches.to_csv(args.outdir / "virtual_marker_mismatches.tsv", sep="\t", index=False)

    print(pd.read_csv(args.outdir / "virtual_marker_summary.tsv", sep="\t").to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
