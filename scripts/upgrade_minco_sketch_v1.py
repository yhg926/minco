#!/usr/bin/env python3
"""Upgrade development minco sketches to the v1 metadata layout.

This is an in-place metadata upgrader for sketches produced during minco
development. It canonicalizes legacy minco filenames, converts
minco.ctxmeta.tsv into binary minco.ctxmeta records, and appends/refreshes the
v1 extension in minco.stat with the current density-summary semantics.
"""

from __future__ import annotations

import argparse
import csv
import os
import struct
import sys
from collections import Counter
from pathlib import Path


PATHLEN = 256
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1
MINCO_SEED = 0x9E3779B97F4A7C15

STAT_PREFIX = struct.Struct("<IBB2xiiiiii")
STAT_EXT = struct.Struct("<IHHIIIII4xQIIIIQQQQIIIIIIII")
CTXMETA_RECORD = struct.Struct("<BBHIQQQQQQQQ")

MINCO_STAT_EXT_MAGIC = 0x4D434F53
MINCO_STAT_EXT_VERSION = 1
MINCO_STAT_HASH_FUNCTION_SPLITMIX64 = 1
MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK = 1
MINCO_STAT_SELECTION_BOTTOMK = 1
MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE = 1
MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE = 2
MINCO_STAT_DENSITY_SAMPLE_ID_NONE = UINT32_MAX
MINCO_STAT_DENSITY_FLAG_MIXED_HASH_BITS = 0x01
MINCO_STAT_FLAG_STREAM_BOTTOMK = 0x02
MINCO_STAT_FLAG_MINCO_HASH_BOTTOMK = 0x01
MINCO_STAT_FLAG_SPARSE_CTX_HASH = 0x08

LEGACY_NAMES = {
    "lcofiles.stat": "minco.stat",
    "lcofiles.qc": "minco.qc",
    "lcofiles.anno": "minco.anno",
    "lcofiles.infilemeta": "minco.infilemeta",
    "comblco": "minco.ctxobj64",
    "comblco.index": "minco.ctxobj64.offsets",
    "comblco.a": "minco.ctxobj64.abund",
    "comblco.position": "minco.ctxobj64.position",
    "sortedcomb_ctxgid64obj32": "minco.refindex.ctxgid64obj32",
}

CTXMETA_MODES = {
    "none": 0,
    "preconflict": 1,
    "postconflict": 2,
    "both": 3,
}


def u64(x: int) -> int:
    return x & UINT64_MAX


def mix64(x: int) -> int:
    x = u64(x + 0x9E3779B97F4A7C15)
    x = u64((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9)
    x = u64((x ^ (x >> 27)) * 0x94D049BB133111EB)
    return u64(x ^ (x >> 31))


def fold_id64(x: int) -> int:
    ident = ((x & UINT32_MAX) ^ (x >> 32)) & UINT32_MAX
    return ident or 1


def hash_bits_from_stat(stat: dict[str, int | bool]) -> int:
    ctx_bits = 4 * int(stat["coden_len"]) if int(stat["coden_len"]) > 0 else 4 * int(stat["hclen"])
    obj_bits = 2 * int(stat["klen"]) - ctx_bits
    if obj_bits < 0 or obj_bits > 64:
        raise ValueError(f"invalid context/object dimensions: {stat}")
    return 64 - obj_bits


def feature_id_from_stat(stat: dict[str, int | bool]) -> int:
    iolen = 0 if int(stat["coden_len"]) > 0 else int(stat["klen"]) - 2 * (int(stat["hclen"]) + int(stat["holen"]))
    x = 0x6D696E636F2D4631
    x = mix64(x ^ (int(stat["coden_len"]) & UINT32_MAX))
    x = mix64(x ^ ((int(stat["klen"]) & UINT32_MAX) << 8))
    x = mix64(x ^ ((int(stat["hclen"]) & UINT32_MAX) << 16))
    x = mix64(x ^ ((int(stat["holen"]) & UINT32_MAX) << 24))
    x = mix64(x ^ ((iolen & UINT32_MAX) << 32))
    x = mix64(x ^ (MINCO_STAT_HASH_FUNCTION_SPLITMIX64 << 48))
    x = mix64(x ^ MINCO_SEED)
    x = mix64(x ^ hash_bits_from_stat(stat))
    return fold_id64(x)


def sketch_id_from_stat(stat: dict[str, int | bool], target_sketch_size: int, flags: int) -> int:
    if target_sketch_size <= 0:
        raise ValueError(f"invalid target sketch size: {target_sketch_size}")
    payload_flags = flags & MINCO_STAT_FLAG_SPARSE_CTX_HASH
    x = 0x6D696E636F2D5331
    x = mix64(x ^ feature_id_from_stat(stat))
    x = mix64(x ^ (target_sketch_size << 1))
    x = mix64(x ^ (MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK << 41))
    x = mix64(x ^ (payload_flags << 49))
    return fold_id64(x)


def parse_stat(path: Path) -> tuple[dict[str, int | bool], bytes, dict[str, int] | None]:
    data = path.read_bytes()
    if len(data) < STAT_PREFIX.size:
        raise ValueError(f"{path}: stat too small")
    fields = STAT_PREFIX.unpack_from(data, 0)
    stat = {
        "hash_id": fields[0],
        "koc": bool(fields[1]),
        "conflict": bool(fields[2]),
        "coden_len": fields[3],
        "klen": fields[4],
        "hclen": fields[5],
        "holen": fields[6],
        "compat_filter_shift": fields[7],
        "infile_num": fields[8],
    }
    infile_num = int(stat["infile_num"])
    if infile_num < 0:
        raise ValueError(f"{path}: negative sample count")
    base = STAT_PREFIX.size + infile_num * PATHLEN
    if len(data) < base:
        raise ValueError(f"{path}: truncated stat name table")
    names = data[STAT_PREFIX.size:base]
    ext = None
    if len(data) >= base + STAT_EXT.size:
        values = STAT_EXT.unpack_from(data, base)
        if values[0] == MINCO_STAT_EXT_MAGIC and values[1] == MINCO_STAT_EXT_VERSION and values[2] >= STAT_EXT.size:
            ext = {
                "target_sketch_size": values[3],
                "feature_id": values[4],
                "sketch_id": values[5],
                "hash_function": values[6],
                "hash_bits": values[7],
                "hash_seed": values[8],
                "sketch_model": values[9],
                "selection_mode": values[10],
                "flags": values[11],
                "density_threshold": values[13],
                "density_min_threshold": values[14],
                "density_max_threshold": values[15],
                "density_universal_threshold": values[16],
                "density_min_sample_id": values[17],
                "density_max_sample_id": values[18],
                "density_universal_sample_id": values[19],
                "density_valid_sample_count": values[20],
                "density_hash_bits": values[21],
                "density_universal_policy": values[22],
                "density_flags": values[23],
            }
    return stat, names, ext


def canonicalize_names(sketch_dir: Path, dry_run: bool) -> list[str]:
    actions: list[str] = []
    for legacy, current in LEGACY_NAMES.items():
        src = sketch_dir / legacy
        dst = sketch_dir / current
        if not src.exists():
            continue
        if dst.exists():
            if src.stat().st_size == dst.stat().st_size:
                actions.append(f"keep {dst.name}; legacy {src.name} still present")
                continue
            raise FileExistsError(f"{dst} already exists and differs from {src}")
        actions.append(f"rename {src.name} -> {dst.name}")
        if not dry_run:
            src.rename(dst)
    return actions


def parse_ctxmeta_tsv(path: Path, expected_samples: int) -> tuple[bytes, dict[str, int], int]:
    rows: list[tuple[int, tuple[int, ...]]] = []
    target_counts: Counter[int] = Counter()
    valid_count = 0
    min_info: tuple[float, int, int, int] | None = None
    max_info: tuple[float, int, int, int] | None = None
    hash_bits_seen: set[int] = set()

    with path.open(newline="") as fp:
        reader = csv.DictReader(fp, delimiter="\t")
        required = [
            "sample_id",
            "mode",
            "valid",
            "hash_bits",
            "threshold",
            "sketch_entries",
            "selected_observed_ctx",
            "selected_estimated_unique_ctx",
            "preconflict_observed_ctx",
            "preconflict_estimated_unique_ctx",
            "postconflict_observed_ctx",
            "postconflict_estimated_unique_ctx",
        ]
        missing = [name for name in required if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(missing)}")

        records = [None] * expected_samples
        for line_no, row in enumerate(reader, start=2):
            sample_id = int(row["sample_id"])
            if sample_id < 0 or sample_id >= expected_samples:
                raise ValueError(f"{path}:{line_no}: sample_id out of range: {sample_id}")
            if records[sample_id] is not None:
                raise ValueError(f"{path}:{line_no}: duplicate sample_id {sample_id}")
            mode = CTXMETA_MODES.get(row["mode"])
            if mode is None:
                raise ValueError(f"{path}:{line_no}: unknown mode {row['mode']!r}")
            valid = int(row["valid"])
            hash_bits = int(row["hash_bits"])
            threshold = int(row["threshold"])
            sketch_entries = int(row["sketch_entries"])
            selected_observed = int(row["selected_observed_ctx"])
            selected_estimated = int(row["selected_estimated_unique_ctx"])
            pre_observed = int(row["preconflict_observed_ctx"])
            pre_estimated = int(row["preconflict_estimated_unique_ctx"])
            post_observed = int(row["postconflict_observed_ctx"])
            post_estimated = int(row["postconflict_estimated_unique_ctx"])
            packed_fields = (
                mode,
                1 if valid else 0,
                0,
                hash_bits,
                threshold,
                sketch_entries,
                selected_observed,
                selected_estimated,
                pre_observed,
                pre_estimated,
                post_observed,
                post_estimated,
            )
            records[sample_id] = packed_fields
            if sketch_entries > 0:
                target_counts[sketch_entries] += 1
            if valid:
                valid_count += 1
                hash_bits_seen.add(hash_bits)
                density = (threshold + 1) / float(1 << hash_bits)
                info = (density, sample_id, hash_bits, threshold)
                if min_info is None or density < min_info[0]:
                    min_info = info
                if max_info is None or density > max_info[0]:
                    max_info = info

    if any(record is None for record in records):
        missing_count = sum(record is None for record in records)
        raise ValueError(f"{path}: missing {missing_count} sample records")
    binary = b"".join(CTXMETA_RECORD.pack(*record) for record in records if record is not None)
    if len(binary) != expected_samples * CTXMETA_RECORD.size:
        raise AssertionError("internal ctxmeta packing size mismatch")

    target = target_counts.most_common(1)[0][0] if target_counts else 10000
    if valid_count == 0 or min_info is None or max_info is None:
        density = {
            "min_threshold": UINT64_MAX,
            "max_threshold": UINT64_MAX,
            "universal_threshold": UINT64_MAX,
            "min_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "max_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "universal_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "valid_sample_count": 0,
            "hash_bits": 0,
            "universal_policy": 0,
            "flags": 0,
        }
    else:
        density = {
            "min_threshold": min_info[3],
            "max_threshold": max_info[3],
            "universal_threshold": max_info[3],
            "min_sample_id": min_info[1],
            "max_sample_id": max_info[1],
            "universal_sample_id": max_info[1],
            "valid_sample_count": valid_count,
            "hash_bits": max_info[2],
            "universal_policy": MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE
            if expected_samples <= 1
            else MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE,
            "flags": MINCO_STAT_DENSITY_FLAG_MIXED_HASH_BITS if len(hash_bits_seen) > 1 else 0,
        }
    return binary, density, target


def summarize_ctxmeta_records(records: list[tuple[int, ...]], expected_samples: int) -> tuple[dict[str, int], int]:
    target_counts: Counter[int] = Counter()
    valid_count = 0
    min_info: tuple[float, int, int, int] | None = None
    max_info: tuple[float, int, int, int] | None = None
    hash_bits_seen: set[int] = set()

    for sample_id, record in enumerate(records):
        valid = int(record[1])
        hash_bits = int(record[3])
        threshold = int(record[4])
        sketch_entries = int(record[5])
        if sketch_entries > 0:
            target_counts[sketch_entries] += 1
        if valid:
            valid_count += 1
            hash_bits_seen.add(hash_bits)
            density_value = (threshold + 1) / float(1 << hash_bits)
            info = (density_value, sample_id, hash_bits, threshold)
            if min_info is None or density_value < min_info[0]:
                min_info = info
            if max_info is None or density_value > max_info[0]:
                max_info = info

    target = target_counts.most_common(1)[0][0] if target_counts else 10000
    if valid_count == 0 or min_info is None or max_info is None:
        density = {
            "min_threshold": UINT64_MAX,
            "max_threshold": UINT64_MAX,
            "universal_threshold": UINT64_MAX,
            "min_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "max_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "universal_sample_id": MINCO_STAT_DENSITY_SAMPLE_ID_NONE,
            "valid_sample_count": 0,
            "hash_bits": 0,
            "universal_policy": 0,
            "flags": 0,
        }
    else:
        density = {
            "min_threshold": min_info[3],
            "max_threshold": max_info[3],
            "universal_threshold": max_info[3],
            "min_sample_id": min_info[1],
            "max_sample_id": max_info[1],
            "universal_sample_id": max_info[1],
            "valid_sample_count": valid_count,
            "hash_bits": max_info[2],
            "universal_policy": MINCO_STAT_DENSITY_POLICY_SINGLE_SAMPLE
            if expected_samples <= 1
            else MINCO_STAT_DENSITY_POLICY_LARGEST_SAMPLE,
            "flags": MINCO_STAT_DENSITY_FLAG_MIXED_HASH_BITS if len(hash_bits_seen) > 1 else 0,
        }
    return density, target


def parse_ctxmeta_bin(path: Path, expected_samples: int) -> tuple[dict[str, int], int]:
    data = path.read_bytes()
    expected = expected_samples * CTXMETA_RECORD.size
    if len(data) != expected:
        raise ValueError(f"{path}: {len(data)} bytes, expected {expected}")
    records = list(CTXMETA_RECORD.iter_unpack(data))
    return summarize_ctxmeta_records(records, expected_samples)


def write_atomic(path: Path, data: bytes, dry_run: bool) -> None:
    if dry_run:
        return
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def write_stat(
    path: Path,
    stat: dict[str, int | bool],
    names: bytes,
    old_ext: dict[str, int] | None,
    density: dict[str, int],
    inferred_target: int,
    dry_run: bool,
) -> tuple[int, int]:
    target = old_ext["target_sketch_size"] if old_ext and old_ext.get("target_sketch_size", 0) else inferred_target
    selection_mode = old_ext["selection_mode"] if old_ext and old_ext.get("selection_mode", 0) else MINCO_STAT_SELECTION_BOTTOMK
    flags = old_ext["flags"] if old_ext is not None else (MINCO_STAT_FLAG_MINCO_HASH_BOTTOMK | MINCO_STAT_FLAG_STREAM_BOTTOMK)
    density_threshold = old_ext["density_threshold"] if old_ext is not None else UINT64_MAX
    feature_id = feature_id_from_stat(stat)
    sketch_id = sketch_id_from_stat(stat, target, flags)
    hash_bits = hash_bits_from_stat(stat)

    prefix = STAT_PREFIX.pack(
        sketch_id,
        1 if stat["koc"] else 0,
        1 if stat["conflict"] else 0,
        int(stat["coden_len"]),
        int(stat["klen"]),
        int(stat["hclen"]),
        int(stat["holen"]),
        0,
        int(stat["infile_num"]),
    )
    ext = STAT_EXT.pack(
        MINCO_STAT_EXT_MAGIC,
        MINCO_STAT_EXT_VERSION,
        STAT_EXT.size,
        target,
        feature_id,
        sketch_id,
        MINCO_STAT_HASH_FUNCTION_SPLITMIX64,
        hash_bits,
        MINCO_SEED,
        MINCO_STAT_SKETCH_MODEL_CTX_BOTTOMK,
        selection_mode,
        flags,
        0,
        density_threshold,
        density["min_threshold"],
        density["max_threshold"],
        density["universal_threshold"],
        density["min_sample_id"],
        density["max_sample_id"],
        density["universal_sample_id"],
        density["valid_sample_count"],
        density["hash_bits"] or hash_bits,
        density["universal_policy"],
        density["flags"],
        0,
    )
    write_atomic(path, prefix + names + ext, dry_run)
    return target, sketch_id


def remove_legacy_text(sketch_dir: Path, dry_run: bool) -> list[str]:
    actions = []
    for name in ("minco.ctxmeta.tsv", "minco.ctxsetmeta.tsv"):
        path = sketch_dir / name
        if path.exists():
            actions.append(f"remove {name}")
            if not dry_run:
                path.unlink()
    return actions


def upgrade_one(sketch_dir: Path, dry_run: bool, remove_text: bool) -> None:
    if not sketch_dir.is_dir():
        raise NotADirectoryError(sketch_dir)
    print(f"[{sketch_dir}]")
    for action in canonicalize_names(sketch_dir, dry_run):
        print(f"  {action}")

    stat_path = sketch_dir / "minco.stat"
    stat_read_path = stat_path
    if not stat_read_path.exists() and (sketch_dir / "lcofiles.stat").exists():
        stat_read_path = sketch_dir / "lcofiles.stat"
    ctxmeta_tsv = sketch_dir / "minco.ctxmeta.tsv"
    ctxmeta_bin = sketch_dir / "minco.ctxmeta"
    if not stat_read_path.exists():
        raise FileNotFoundError(f"{stat_path} not found after canonicalization")
    if not ctxmeta_tsv.exists() and not ctxmeta_bin.exists():
        raise FileNotFoundError(f"{sketch_dir}: neither minco.ctxmeta.tsv nor minco.ctxmeta exists")

    stat, names, old_ext = parse_stat(stat_read_path)
    sample_count = int(stat["infile_num"])
    if ctxmeta_tsv.exists():
        binary, density, inferred_target = parse_ctxmeta_tsv(ctxmeta_tsv, sample_count)
        write_atomic(ctxmeta_bin, binary, dry_run)
        print(f"  write minco.ctxmeta ({sample_count} records, {len(binary)} bytes)")
    else:
        density, inferred_target = parse_ctxmeta_bin(ctxmeta_bin, sample_count)
        print(f"  read existing minco.ctxmeta ({sample_count} records)")

    target, sketch_id = write_stat(stat_path, stat, names, old_ext, density, inferred_target, dry_run)
    print(
        "  refresh minco.stat "
        f"(samples={sample_count}, S={target}, sketch_id={sketch_id}, "
        f"valid_density={density['valid_sample_count']}, "
        f"universal_sample={density['universal_sample_id']}, "
        f"universal_threshold={density['universal_threshold']})"
    )
    if remove_text:
        for action in remove_legacy_text(sketch_dir, dry_run):
            print(f"  {action}")


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sketch_dir", nargs="+", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--remove-legacy-text", action="store_true")
    args = parser.parse_args(argv)
    for sketch_dir in args.sketch_dir:
        upgrade_one(sketch_dir, args.dry_run, args.remove_legacy_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
