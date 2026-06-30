#!/usr/bin/env python3
"""Audit cheap prefix-read preflight for exact-sidecar decisions."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

EXACT_TRIGGER = 0.10
MIN_BASE_N = 10
MIN_RAW_UNIQUE_N = 10

PREFLIGHTS = [
    {
        "dataset": "hmp_gastrooral",
        "sample": "0",
        "reads": "50000",
        "profile": "/tmp/minco_exact_preflight_20260627/sample0_preflight_universal.tsv",
        "time_log": "/tmp/minco_exact_preflight_20260627/sample0_preflight_universal.time.log",
        "full_profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv",
        "full_policy": "exact_candidate",
    },
    {
        "dataset": "hmp_gastrooral",
        "sample": "6",
        "reads": "50000",
        "profile": "/tmp/minco_exact_preflight_20260627/sample6_preflight_universal.tsv",
        "time_log": "/tmp/minco_exact_preflight_20260627/sample6_preflight_universal.time.log",
        "full_profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv",
        "full_policy": "avoid_exact_high_extra",
    },
    {
        "dataset": "hmp_gastrooral",
        "sample": "0",
        "reads": "200000",
        "profile": "/tmp/minco_exact_preflight_20260627/sample0_preflight200k_universal.tsv",
        "time_log": "/tmp/minco_exact_preflight_20260627/sample0_preflight200k_universal.time.log",
        "full_profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv",
        "full_policy": "exact_candidate",
    },
    {
        "dataset": "hmp_gastrooral",
        "sample": "6",
        "reads": "200000",
        "profile": "/tmp/minco_exact_preflight_20260627/sample6_preflight200k_universal.tsv",
        "time_log": "/tmp/minco_exact_preflight_20260627/sample6_preflight200k_universal.time.log",
        "full_profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv",
        "full_policy": "avoid_exact_high_extra",
    },
]


def read_first(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", nrows=1)
    if df.empty:
        return {}
    return df.iloc[0].to_dict()


def val(row: dict[str, object], key: str, default: object = "") -> object:
    out = row.get(key, default)
    if pd.isna(out):
        return default
    return out


def fnum(row: dict[str, object], key: str, default: float = math.nan) -> float:
    try:
        return float(val(row, key, default))
    except (TypeError, ValueError):
        return default


def parse_seconds(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "", ""
    wall = ""
    rss = ""
    for line in path.read_text(errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            wall = line.rsplit(": ", 1)[1].strip()
        elif "Maximum resident set size" in line:
            match = re.search(r"(\d+)", line)
            if match:
                rss = f"{int(match.group(1)) / 1024.0 / 1024.0:.4f}"
    if not wall:
        return "", rss
    parts = wall.split(":")
    if len(parts) == 2:
        seconds = float(parts[0]) * 60.0 + float(parts[1])
    elif len(parts) == 3:
        seconds = float(parts[0]) * 3600.0 + float(parts[1]) * 60.0 + float(parts[2])
    else:
        seconds = math.nan
    return f"{seconds:.2f}", rss


def raw_policy(row: dict[str, object]) -> str:
    p_extra = fnum(row, "probability_extra_mass_ratio")
    guard = fnum(row, "joined_base_median_uaf") < 0.35 or fnum(row, "raw_unique_to_base_ratio") >= 0.80
    low_added = int(fnum(row, "low_extra_split_rescue_added_n", 0.0))
    if p_extra > EXACT_TRIGGER:
        return "avoid_exact_high_extra"
    if not guard:
        return "avoid_exact_guard"
    if low_added > 0:
        return "avoid_exact_low_extra_rescue"
    return "exact_candidate"


def support_aware_policy(row: dict[str, object]) -> str:
    base_n = int(fnum(row, "joined_base_n", 0.0))
    raw_n = int(fnum(row, "raw_unique95_n", 0.0))
    if base_n < MIN_BASE_N or raw_n < MIN_RAW_UNIQUE_N:
        return "insufficient_support"
    return raw_policy(row)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    rows: list[dict[str, object]] = []
    for spec in PREFLIGHTS:
        profile = Path(spec["profile"])
        first = read_first(profile)
        seconds, rss = parse_seconds(Path(spec["time_log"]))
        raw = raw_policy(first) if first else "missing_profile"
        support = support_aware_policy(first) if first else "missing_profile"
        full_first = read_first(Path(spec["full_profile"]))
        full_block_p_extra = fnum(full_first, "auto_exact_split_block_p_extra_mass_ratio", fnum(full_first, "probability_extra_mass_ratio"))
        row = {
            "dataset": spec["dataset"],
            "sample": spec["sample"],
            "reads": spec["reads"],
            "profile_exists": str(bool(first)).lower(),
            "seconds": seconds,
            "peak_rss_gib": rss,
            "joined_base_n": val(first, "joined_base_n"),
            "probability_gate_n": val(first, "probability_gate_n"),
            "probability_extra_n": val(first, "probability_extra_n"),
            "probability_extra_mass_ratio": val(first, "probability_extra_mass_ratio"),
            "joined_base_median_uaf": val(first, "joined_base_median_uaf"),
            "raw_unique95_n": val(first, "raw_unique95_n"),
            "raw_unique_to_base_ratio": val(first, "raw_unique_to_base_ratio"),
            "low_extra_split_rescue_added_n": val(first, "low_extra_split_rescue_added_n"),
            "raw_preflight_policy": raw,
            "support_aware_preflight_policy": support,
            "full_policy": spec["full_policy"],
            "raw_policy_matches_full": str(raw == spec["full_policy"]).lower(),
            "support_aware_policy_matches_full": str(support == spec["full_policy"]).lower(),
            "full_joined_base_n": val(full_first, "joined_base_n"),
            "full_probability_extra_mass_ratio": full_block_p_extra,
            "full_raw_unique_to_base_ratio": val(full_first, "raw_unique_to_base_ratio"),
            "profile": str(profile),
            "time_log": spec["time_log"],
        }
        rows.append(row)

    raw_decisive = [row for row in rows if row["raw_preflight_policy"] != "insufficient_support"]
    support_decisive = [row for row in rows if row["support_aware_preflight_policy"] != "insufficient_support"]
    support_insufficient = [row for row in rows if row["support_aware_preflight_policy"] == "insufficient_support"]
    raw_mismatches = [row for row in rows if row["raw_policy_matches_full"] == "false"]
    support_mismatches = [
        row for row in support_decisive
        if row["support_aware_policy_matches_full"] == "false"
    ]
    audit_rows = [
        {
            "metric": "evaluated_preflights",
            "value": len(rows),
            "evidence": "exact_preflight_policy_instances.tsv",
            "decision": "diagnostic_small_hmp_prefix_test",
        },
        {
            "metric": "raw_policy_mismatches",
            "value": len(raw_mismatches),
            "evidence": "raw prefix policy without support floor",
            "decision": "raw_prefix_policy_unstable",
        },
        {
            "metric": "support_aware_decisive_preflights",
            "value": len(support_decisive),
            "evidence": f"joined_base_n>={MIN_BASE_N};raw_unique95_n>={MIN_RAW_UNIQUE_N}",
            "decision": "support_floor_removes_unstable_decisions",
        },
        {
            "metric": "support_aware_insufficient_preflights",
            "value": len(support_insufficient),
            "evidence": "exact_preflight_policy_instances.tsv",
            "decision": "first50k_and_first200k_too_sparse_here",
        },
        {
            "metric": "support_aware_mismatches",
            "value": len(support_mismatches),
            "evidence": "only counts decisive support-aware rows",
            "decision": "no_decisive_rows_to_validate",
        },
        {
            "metric": "promotion_decision",
            "value": "reject_prefix_preflight_default",
            "evidence": "first50k/first200k HMP gastrooral sample0/sample6",
            "decision": "cheap_prefix_preflight_not_stable_or_decisive",
        },
    ]

    fields = [
        "dataset",
        "sample",
        "reads",
        "profile_exists",
        "seconds",
        "peak_rss_gib",
        "joined_base_n",
        "probability_gate_n",
        "probability_extra_n",
        "probability_extra_mass_ratio",
        "joined_base_median_uaf",
        "raw_unique95_n",
        "raw_unique_to_base_ratio",
        "low_extra_split_rescue_added_n",
        "raw_preflight_policy",
        "support_aware_preflight_policy",
        "full_policy",
        "raw_policy_matches_full",
        "support_aware_policy_matches_full",
        "full_joined_base_n",
        "full_probability_extra_mass_ratio",
        "full_raw_unique_to_base_ratio",
        "profile",
        "time_log",
    ]
    write_tsv(RESULTS / "exact_preflight_policy_instances.tsv", rows, fields)
    write_tsv(RESULTS / "exact_preflight_policy_audit.tsv", audit_rows, ["metric", "value", "evidence", "decision"])
    print(pd.DataFrame(rows)[fields].to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
