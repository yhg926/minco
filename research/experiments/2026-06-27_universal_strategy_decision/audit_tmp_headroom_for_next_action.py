#!/usr/bin/env python3
"""Audit /tmp disk headroom before running the next release-grade sample."""

from __future__ import annotations

import csv
import os
import re
import shutil
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
NEXT_ACTIONS = RESULTS / "next_release_grade_actions.tsv"
HEADROOM_OUT = RESULTS / "tmp_headroom_audit.tsv"
CANDIDATES_OUT = RESULTS / "tmp_cleanup_candidates.tsv"

TMP = Path("/tmp")
ALTERNATE_WORK_ROOT = Path("/mnt/new3T/minco_release_holdouts_20260628")
UNSEEN_ROOT = TMP / "cami2_hmp_unseen_transfer_20260626"
AIRSKIN_ROOT = TMP / "cami2_hmp_airskin_20260625"
SAFETY_BUFFER_BYTES = 2 * 1024**3
MIN_RECOMMENDED_FREE_BYTES = 12 * 1024**3


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def file_size(path: Path) -> int:
    try:
        return path.stat().st_size if path.is_file() else 0
    except OSError:
        return 0


def tree_size(path: Path) -> int:
    if path.is_file():
        return file_size(path)
    total = 0
    if not path.exists():
        return 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            total += file_size(Path(root) / name)
    return total


def gib(n_bytes: int) -> float:
    return n_bytes / 1024**3


def disk_usage_for_target(path: Path) -> shutil._ntuple_diskusage:
    probe = path
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe)


def next_sample_id() -> int:
    for row in read_tsv(NEXT_ACTIONS):
        if row["action_id"].startswith("restore_profile_hmp_airskin_sample"):
            match = re.search(r"sample([0-9]+)", row["action_id"])
            if match:
                return int(match.group(1))
    return -1


def existing_sample_fastqs() -> list[Path]:
    paths: list[Path] = []
    for root in (UNSEEN_ROOT, AIRSKIN_ROOT):
        paths.extend((root / "reads").glob("airskinurogenital_sample*.nonzero.fastq.gz"))
    return [path for path in paths if path.exists()]


def existing_sample_sketches() -> list[Path]:
    paths: list[Path] = []
    for root in (UNSEEN_ROOT, AIRSKIN_ROOT):
        paths.extend((root / "run_sylph_r232").glob("sylph_sample*/*.sylsp"))
    return [path for path in paths if path.exists()]


def existing_minco_dirs() -> list[Path]:
    return [path for path in TMP.glob("minco_current_code_hmp_airskin*_20260627") if path.is_dir()]


def cleanup_rows() -> list[dict[str, object]]:
    candidates: list[tuple[Path, str, str]] = []
    completed_samples = [0, 5, 22, 28]
    for sample in completed_samples:
        candidates.extend(
            [
                (
                    UNSEEN_ROOT / f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
                    "completed_hmp_fastq",
                    "completed sample; rerunnable from recorded command and archive if needed",
                ),
                (
                    AIRSKIN_ROOT / f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
                    "completed_hmp_fastq",
                    "completed sample; rerunnable from recorded command and archive if needed",
                ),
                (
                    UNSEEN_ROOT
                    / f"run_sylph_r232/sylph_sample{sample}/airskinurogenital_sample{sample}.nonzero.fastq.gz.sylsp",
                    "completed_hmp_sylph_sketch",
                    "profile output is recorded; sketch can be regenerated from FASTQ if needed",
                ),
                (
                    AIRSKIN_ROOT
                    / f"run_sylph_r232/sylph_sample{sample}/airskinurogenital_sample{sample}.nonzero.fastq.gz.sylsp",
                    "completed_hmp_sylph_sketch",
                    "profile output is recorded; sketch can be regenerated from FASTQ if needed",
                ),
            ]
        )

    candidates.extend(
        [
            (
                TMP / "cami_marine_sample1_reads.fq.gz",
                "old_large_raw_read_cache",
                "large temporary read cache outside current HMP next-action path",
            ),
            (
                TMP / "cami_marine_sample2_reads.fq.gz",
                "old_large_raw_read_cache",
                "large temporary read cache outside current HMP next-action path",
            ),
            (
                TMP / "minco_context_em_external_20260624/cami_marine_sample0_reads.fq.gz",
                "old_large_raw_read_cache",
                "large temporary read cache outside current HMP next-action path",
            ),
            (
                TMP / "minco_context_em_external_20260624/cami3_s0_edges.tsv",
                "old_large_diagnostic_table",
                "large temporary diagnostic table outside current HMP next-action path",
            ),
            (
                TMP / "minco_universal_controlled_rescue_20260625/full_panel_hybrid_max_support.tsv",
                "old_large_diagnostic_table",
                "large temporary diagnostic table outside current HMP next-action path",
            ),
        ]
    )

    rows: list[dict[str, object]] = []
    for path, category, rationale in candidates:
        size = tree_size(path)
        if size <= 0:
            continue
        rows.append(
            {
                "path": str(path),
                "bytes": size,
                "GiB": f"{gib(size):.3f}",
                "category": category,
                "delete_status": "manual_review_required",
                "rationale": rationale,
            }
        )
    rows.sort(key=lambda row: int(row["bytes"]), reverse=True)
    return rows


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    usage = shutil.disk_usage(TMP)
    alternate_usage = disk_usage_for_target(ALTERNATE_WORK_ROOT)
    sample = next_sample_id()

    fastq_sizes = [file_size(path) for path in existing_sample_fastqs()]
    sketch_sizes = [file_size(path) for path in existing_sample_sketches()]
    minco_sizes = [tree_size(path) for path in existing_minco_dirs()]
    expected_fastq = max(fastq_sizes, default=4 * 1024**3)
    expected_sketch = max(sketch_sizes, default=128 * 1024**2)
    expected_minco = max(minco_sizes, default=512 * 1024**2)
    estimated_required = expected_fastq + expected_sketch + expected_minco + SAFETY_BUFFER_BYTES
    required_floor = max(estimated_required, MIN_RECOMMENDED_FREE_BYTES)
    headroom_pass = usage.free >= required_floor
    alternate_headroom_pass = alternate_usage.free >= required_floor

    rows: list[dict[str, object]] = [
        {"metric": "next_sample", "value": sample, "unit": "sample_id"},
        {"metric": "tmp_total_bytes", "value": usage.total, "unit": "bytes"},
        {"metric": "tmp_used_bytes", "value": usage.used, "unit": "bytes"},
        {"metric": "tmp_free_bytes", "value": usage.free, "unit": "bytes"},
        {"metric": "tmp_free_GiB", "value": f"{gib(usage.free):.3f}", "unit": "GiB"},
        {"metric": "observed_max_hmp_fastq_bytes", "value": expected_fastq, "unit": "bytes"},
        {"metric": "observed_max_hmp_sylph_sketch_bytes", "value": expected_sketch, "unit": "bytes"},
        {"metric": "observed_max_hmp_minco_work_bytes", "value": expected_minco, "unit": "bytes"},
        {"metric": "safety_buffer_bytes", "value": SAFETY_BUFFER_BYTES, "unit": "bytes"},
        {"metric": "estimated_required_free_bytes", "value": estimated_required, "unit": "bytes"},
        {
            "metric": "estimated_required_free_GiB",
            "value": f"{gib(estimated_required):.3f}",
            "unit": "GiB",
        },
        {
            "metric": "min_recommended_free_bytes",
            "value": MIN_RECOMMENDED_FREE_BYTES,
            "unit": "bytes",
        },
        {
            "metric": "headroom_pass",
            "value": headroom_pass,
            "unit": "boolean",
        },
        {"metric": "alternate_work_root", "value": str(ALTERNATE_WORK_ROOT), "unit": "path"},
        {
            "metric": "alternate_free_bytes",
            "value": alternate_usage.free,
            "unit": "bytes",
        },
        {
            "metric": "alternate_free_GiB",
            "value": f"{gib(alternate_usage.free):.3f}",
            "unit": "GiB",
        },
        {
            "metric": "alternate_headroom_pass",
            "value": alternate_headroom_pass,
            "unit": "boolean",
        },
        {
            "metric": "recommended_action",
            "value": (
                "run_next_sample"
                if headroom_pass
                else "run_with_alternate_work_root"
                if alternate_headroom_pass
                else "free_tmp_space_before_next_sample"
            ),
            "unit": "decision",
        },
    ]
    write_tsv(HEADROOM_OUT, rows, ["metric", "value", "unit"])
    cleanup = cleanup_rows()
    write_tsv(CANDIDATES_OUT, cleanup, ["path", "bytes", "GiB", "category", "delete_status", "rationale"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
