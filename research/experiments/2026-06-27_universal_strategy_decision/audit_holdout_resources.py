#!/usr/bin/env python3
"""Audit local resources needed to promote diagnostic panels to release grade."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

SYLPH = Path("/home/ubuntu/yihuiguang/bin/sylph")
SYLPH_DB_ROOT = Path("/mnt/new3T/sylph_db")
SYLPH_R226 = SYLPH_DB_ROOT / "gtdb-r226-c200-dbv1.syldb"
SYLPH_R232_CANDIDATES = [
    SYLPH_DB_ROOT / "gtdb-r232-c200-dbv1.syldb",
    SYLPH_DB_ROOT / "gtdb-r232-c200.syldb",
    Path("/mnt/new3T/gtdbr220/gtdb232/gtdb-r232-c200-dbv1.syldb"),
    Path("/mnt/new3T/gtdbr220/sylph/gtdb-r232-c200-dbv1.syldb"),
]
SYLPH_R232_CHUNK_LIST = SYLPH_DB_ROOT / "gtdb-r232-c200-dbv1_chunks/gtdb-r232-c200-dbv1.chunk_syldb.list"

HMP_TRUTH = Path("/tmp/cami2_hmp_airskin_20260625/truth")
HMP_RUN = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run")
HMP_R232_RUN = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232")
HMP_SAMPLE_PATHS = {
    "sample6_minco_default": HMP_RUN / "minco_sample6_universal_strategy.tsv",
    "sample6_minco_unique": HMP_RUN / "minco_sample6_unique_zip_unfiltered.tsv",
    "sample6_minco_split": HMP_RUN / "minco_sample6_split_zip_unfiltered.tsv",
    "sample6_sylph_r226_profile": HMP_RUN / "sylph_sample6/profile.tsv",
    "sample11_minco_default": HMP_RUN / "minco_sample11_universal_strategy.tsv",
    "sample11_minco_unique": HMP_RUN / "minco_sample11_unique_zip_unfiltered.tsv",
    "sample11_minco_split": HMP_RUN / "minco_sample11_split_zip_unfiltered.tsv",
    "sample11_sylph_r226_profile": HMP_RUN / "sylph_sample11/profile.tsv",
    "truth_genome_to_id": HMP_TRUTH / "genome_to_id.tsv",
    "truth_abundance6": HMP_TRUTH / "abundance6.tsv",
    "truth_abundance11": HMP_TRUTH / "abundance11.tsv",
}
HMP_R232_PROFILE_PATHS = {
    "sample6_sylph_r232_profile": HMP_R232_RUN / "sylph_sample6/profile.tsv",
    "sample11_sylph_r232_profile": HMP_R232_RUN / "sylph_sample11/profile.tsv",
}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sylph_version() -> str:
    if not SYLPH.exists():
        return ""
    proc = subprocess.run([str(SYLPH), "--version"], check=False, text=True, capture_output=True)
    return (proc.stdout or proc.stderr).strip()


def r232_candidates() -> list[Path]:
    candidates = list(SYLPH_R232_CANDIDATES) + [SYLPH_R232_CHUNK_LIST]
    if SYLPH_DB_ROOT.is_dir():
        candidates.extend(sorted(SYLPH_DB_ROOT.glob("*r232*.syldb")))
        candidates.extend(sorted(SYLPH_DB_ROOT.glob("*232*.syldb")))
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in candidates:
        if path not in seen:
            seen.add(path)
            unique.append(path)
    return unique


def main() -> int:
    resources: list[dict[str, object]] = []

    def add(group: str, name: str, path: Path, required_for: str, note: str = "") -> None:
        resources.append(
            {
                "group": group,
                "name": name,
                "path": str(path),
                "exists": str(path.exists()).lower(),
                "required_for": required_for,
                "note": note,
            }
        )

    add("tool", "sylph_binary", SYLPH, "rerun_sylph_same_release")
    add("database", "sylph_gtdb_r226_syldb", SYLPH_R226, "cached_diagnostic_baseline")
    for path in r232_candidates():
        add("database", "sylph_gtdb_r232_candidate", path, "release_grade_hmp_airskin")
    for name, path in HMP_SAMPLE_PATHS.items():
        add("hmp_airskin", name, path, "release_grade_hmp_airskin")
    for name, path in HMP_R232_PROFILE_PATHS.items():
        add("hmp_airskin_r232", name, path, "release_grade_hmp_airskin")

    r232_exists = any(path.exists() for path in r232_candidates())
    r232_chunk_list_exists = SYLPH_R232_CHUNK_LIST.exists()
    cached_hmp_complete = all(path.exists() for path in HMP_SAMPLE_PATHS.values())
    r232_profiles_complete = all(path.exists() for path in HMP_R232_PROFILE_PATHS.values())
    same_release_ready = r232_exists and cached_hmp_complete and r232_profiles_complete
    version = sylph_version()

    summary = [
        {"metric": "sylph_binary_exists", "value": str(SYLPH.exists()).lower()},
        {"metric": "sylph_version", "value": version},
        {"metric": "sylph_gtdb_r226_syldb_exists", "value": str(SYLPH_R226.exists()).lower()},
        {"metric": "sylph_gtdb_r232_syldb_exists", "value": str(r232_exists).lower()},
        {"metric": "sylph_gtdb_r232_chunk_list_exists", "value": str(r232_chunk_list_exists).lower()},
        {"metric": "hmp_airskin_cached_inputs_complete", "value": str(cached_hmp_complete).lower()},
        {"metric": "hmp_airskin_r232_profiles_complete", "value": str(r232_profiles_complete).lower()},
        {"metric": "hmp_airskin_same_release_sylph_ready", "value": str(same_release_ready).lower()},
        {
            "metric": "main_blocker",
            "value": ""
            if same_release_ready
            else "GTDB r232 Sylph DB/profile is missing; cached HMP Sylph outputs use GTDB r226",
        },
    ]

    write_tsv(
        RESULTS / "holdout_resource_audit.tsv",
        resources,
        ["group", "name", "path", "exists", "required_for", "note"],
    )
    write_tsv(RESULTS / "holdout_resource_audit_summary.tsv", summary, ["metric", "value"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
