#!/usr/bin/env python3
"""Create a smaller bottom-k minco sketch from an existing larger sketch."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import struct
import sys
import time

import numpy as np


PATHLEN = 256
MINCO_SEED = 0x9E3779B97F4A7C15
MASK64 = (1 << 64) - 1

STAT_PREFIX = struct.Struct("<IBBxxiiiiii")
STAT_EXT = struct.Struct("<IHHIIIII4xQIIIIQQQQIIIIIIII")
CTXMETA = struct.Struct("<BBHIQQQQQQQQ")

MINCO_STAT_EXT_MAGIC = 0x4D434F53
MINCO_STAT_EXT_VERSION = 1
MINCO_STAT_SELECTION_BOTTOMK = 1
MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK = 1
MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE = 1
MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE = 2
MINCO_STAT_DENSITY_SAMPLE_ID_NONE = 0xFFFFFFFF
MINCO_STAT_FLAG_SPARSE_CTX_HASH = 0x08

SKETCH_FILE = "minco.ctxobj64"
OFFSETS_FILE = "minco.ctxobj64.offsets"
STAT_FILE = "minco.stat"
CTXMETA_FILE = "minco.ctxmeta"
OPTIONAL_COPY_FILES = (
    "minco.anno",
    "minco.infilemeta",
    "minco.qc",
)


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def mix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return (x ^ (x >> 31)) & MASK64


def fold_id64(x: int) -> int:
    value = (x ^ (x >> 32)) & 0xFFFFFFFF
    return value or 1


def sketch_id(feature_id: int, target_size: int, flags: int) -> int:
    x = 0x6D696E636F2D5331  # "minco-S1"
    payload_flags = flags & MINCO_STAT_FLAG_SPARSE_CTX_HASH
    x = mix64(x ^ feature_id)
    x = mix64(x ^ (target_size << 1))
    x = mix64(x ^ (MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK << 41))
    x = mix64(x ^ (payload_flags << 49))
    return fold_id64(x)


def estimate_unique(observed: int, threshold: int, hash_bits: int) -> int:
    if observed <= 0 or hash_bits <= 0:
        return max(observed, 0)
    denom = threshold + 1
    if denom <= 0:
        return 0
    return (observed * (1 << hash_bits) + denom // 2) // denom


def read_stat(path: Path) -> tuple[bytearray, tuple[int, ...], tuple[int, ...], int, int]:
    data = bytearray(path.read_bytes())
    if len(data) < STAT_PREFIX.size:
        raise RuntimeError(f"{path}: stat file is too small")
    prefix = STAT_PREFIX.unpack_from(data, 0)
    infile_num = prefix[-1]
    if infile_num < 0:
        raise RuntimeError(f"{path}: negative sample count {infile_num}")
    base = STAT_PREFIX.size + infile_num * PATHLEN
    if len(data) < base + STAT_EXT.size:
        raise RuntimeError(f"{path}: missing minco stat extension")
    ext = STAT_EXT.unpack_from(data, base)
    if ext[0] != MINCO_STAT_EXT_MAGIC or ext[1] != MINCO_STAT_EXT_VERSION:
        raise RuntimeError(f"{path}: unsupported minco stat extension")
    if ext[2] < STAT_EXT.size:
        raise RuntimeError(f"{path}: short minco stat extension size {ext[2]}")
    return data, prefix, ext, infile_num, base


def write_stat(
    out_path: Path,
    stat_data: bytearray,
    ext: tuple[int, ...],
    base: int,
    target_size: int,
    density_summary: dict[str, int],
) -> int:
    feature_id = ext[4]
    flags = ext[11]
    hash_bits = ext[7]
    new_sketch_id = sketch_id(feature_id, target_size, flags)

    prefix = list(STAT_PREFIX.unpack_from(stat_data, 0))
    prefix[0] = new_sketch_id
    STAT_PREFIX.pack_into(stat_data, 0, *prefix)
    new_ext = (
        MINCO_STAT_EXT_MAGIC,
        MINCO_STAT_EXT_VERSION,
        STAT_EXT.size,
        target_size,
        feature_id,
        new_sketch_id,
        ext[6],
        hash_bits,
        ext[8],
        ext[9],
        MINCO_STAT_SELECTION_BOTTOMK,
        flags,
        0,
        MASK64,
        density_summary["min_threshold"],
        density_summary["max_threshold"],
        density_summary["universal_threshold"],
        density_summary["min_sample_id"],
        density_summary["max_sample_id"],
        density_summary["universal_sample_id"],
        density_summary["valid_count"],
        hash_bits,
        density_summary["universal_policy"],
        density_summary["flags"],
        0,
    )
    STAT_EXT.pack_into(stat_data, base, *new_ext)
    out_path.write_bytes(stat_data)
    return new_sketch_id


def prepare_output(dst: Path, force: bool) -> None:
    if dst.exists():
        if not force:
            raise RuntimeError(f"{dst} already exists; use --force to replace it")
        shutil.rmtree(dst)
    dst.mkdir(parents=True)


def read_offsets(path: Path, infile_num: int) -> np.ndarray:
    offsets = np.fromfile(path, dtype="<u8")
    if offsets.size != infile_num + 1:
        raise RuntimeError(
            f"{path}: found {offsets.size} offsets, expected {infile_num + 1}"
        )
    return offsets


def load_ctxmeta(path: Path, infile_num: int) -> list[tuple[int, ...]] | None:
    if not path.exists():
        return None
    data = path.read_bytes()
    expected = infile_num * CTXMETA.size
    if len(data) != expected:
        raise RuntimeError(f"{path}: {len(data)} bytes, expected {expected}")
    return [CTXMETA.unpack_from(data, i * CTXMETA.size) for i in range(infile_num)]


def downsample(src: Path, dst: Path, target_size: int, force: bool) -> None:
    stat_data, _prefix, ext, infile_num, base = read_stat(src / STAT_FILE)
    hash_bits = ext[7]
    n_obj_bits = 64 - hash_bits
    if n_obj_bits < 0 or n_obj_bits > 64:
        raise RuntimeError(f"{src}: invalid hash_bits={hash_bits}")

    old_target = ext[3]
    if old_target and target_size > old_target:
        raise RuntimeError(
            f"target size {target_size} is larger than source target {old_target}"
        )

    offsets = read_offsets(src / OFFSETS_FILE, infile_num)
    entries_path = src / SKETCH_FILE
    if entries_path.stat().st_size != int(offsets[-1]) * 8:
        raise RuntimeError(f"{entries_path}: size does not match offsets")
    old_ctxmeta = load_ctxmeta(src / CTXMETA_FILE, infile_num)

    prepare_output(dst, force)
    for name in OPTIONAL_COPY_FILES:
        source = src / name
        if source.exists():
            shutil.copy2(source, dst / name)

    new_offsets = np.zeros(infile_num + 1, dtype="<u8")
    ctxmeta_records: list[tuple[int, ...]] = []
    valid_count = 0
    min_threshold = MASK64
    max_threshold = 0
    min_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE
    max_sample_id = MINCO_STAT_DENSITY_SAMPLE_ID_NONE

    log(
        f"downsample {src} -> {dst}; samples={infile_num}; "
        f"S={target_size}; hash_bits={hash_bits}"
    )
    with entries_path.open("rb") as inp, (dst / SKETCH_FILE).open("wb") as out:
        for i in range(infile_num):
            start = int(offsets[i])
            end = int(offsets[i + 1])
            count = end - start
            keep = min(target_size, count)
            new_offsets[i + 1] = int(new_offsets[i]) + keep

            if keep > 0:
                inp.seek(start * 8)
                data = inp.read(keep * 8)
                if len(data) != keep * 8:
                    raise RuntimeError(f"{entries_path}: short read at sample {i}")
                out.write(data)
                values = np.frombuffer(data, dtype="<u8")
                ctx = values >> n_obj_bits
                threshold = int(ctx[-1])
                unique_ctx = int(np.count_nonzero(ctx[1:] != ctx[:-1]) + 1) if keep > 1 else 1
                estimate = estimate_unique(unique_ctx, threshold, hash_bits)
                valid = 1
                valid_count += 1
                if threshold < min_threshold:
                    min_threshold = threshold
                    min_sample_id = i
                if threshold >= max_threshold:
                    max_threshold = threshold
                    max_sample_id = i
            else:
                threshold = 0
                unique_ctx = 0
                estimate = 0
                valid = 0

            mode = old_ctxmeta[i][0] if old_ctxmeta else 3
            ctxmeta_records.append(
                (
                    mode,
                    valid,
                    0,
                    hash_bits,
                    threshold,
                    keep,
                    unique_ctx,
                    estimate,
                    unique_ctx,
                    estimate,
                    unique_ctx,
                    estimate,
                )
            )
            if (i + 1) % 5000 == 0 or i + 1 == infile_num:
                log(
                    f"processed {i + 1}/{infile_num} samples; "
                    f"entries={int(new_offsets[i + 1])}"
                )

    new_offsets.tofile(dst / OFFSETS_FILE)
    with (dst / CTXMETA_FILE).open("wb") as handle:
        for rec in ctxmeta_records:
            handle.write(CTXMETA.pack(*rec))

    if valid_count == 0:
        density_summary = {
            "min_threshold": MASK64,
            "max_threshold": MASK64,
            "universal_threshold": MASK64,
            "min_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "max_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "universal_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "valid_count": 0,
            "universal_policy": 0,
            "flags": 0,
        }
    else:
        density_summary = {
            "min_threshold": min_threshold,
            "max_threshold": max_threshold,
            "universal_threshold": max_threshold,
            "min_sample_id": min_sample_id,
            "max_sample_id": max_sample_id,
            "universal_sample_id": max_sample_id,
            "valid_count": valid_count,
            "universal_policy": MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE
            if valid_count == 1
            else MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE,
            "flags": 0,
        }
    new_id = write_stat(
        dst / STAT_FILE,
        stat_data,
        ext,
        base,
        target_size,
        density_summary,
    )
    log(
        f"done: samples={infile_num}; entries={int(new_offsets[-1])}; "
        f"sketch_id={new_id}; universal_threshold={density_summary['universal_threshold']}"
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("src", type=Path, help="source minco sketch directory")
    parser.add_argument("dst", type=Path, help="output minco sketch directory")
    parser.add_argument("-S", "--sketch-size", type=int, required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if args.sketch_size < 1:
        raise RuntimeError("--sketch-size must be positive")
    downsample(args.src, args.dst, args.sketch_size, args.force)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
