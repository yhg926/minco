#!/usr/bin/env python3
"""Summarize raw-read current-default HMP gastrooral reruns."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
RUN = Path("/tmp/minco_current_code_hmp_gastrooral_20260627")

RAW_SUMMARY = RESULTS / "hmp_gastrooral_raw_default_r232_source_abundance_summary.tsv"
TABLE_SUMMARY = RESULTS / "hmp_gastrooral_r232_source_abundance_summary.tsv"

SAMPLES = {
    "0": {
        "time_log": RUN / "sample0_current_default.time.log",
        "profile": RUN / "minco_sample0_current_default.tsv",
    },
    "6": {
        "time_log": RUN / "sample6_current_default.time.log",
        "profile": RUN / "minco_sample6_current_default.tsv",
    },
}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def parse_elapsed_seconds(text: str) -> float:
    parts = text.strip().split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return float(minutes) * 60.0 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)
    return float("nan")


def parse_time_log(path: Path) -> dict[str, object]:
    out: dict[str, object] = {"time_log": str(path)}
    if not path.exists():
        out.update({"wall_time": "", "seconds": "", "peak_rss_kb": "", "peak_rss_gib": ""})
        return out
    for line in path.read_text(errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            value = line.rsplit(": ", 1)[1].strip()
            out["wall_time"] = value
            out["seconds"] = f"{parse_elapsed_seconds(value):.2f}"
        elif "Maximum resident set size" in line:
            match = re.search(r"(\d+)", line)
            if match:
                kb = int(match.group(1))
                out["peak_rss_kb"] = kb
                out["peak_rss_gib"] = f"{kb / 1024.0 / 1024.0:.4f}"
        elif "User time (seconds)" in line:
            out["user_time_s"] = line.split(":", 1)[1].strip()
        elif "System time (seconds)" in line:
            out["system_time_s"] = line.split(":", 1)[1].strip()
    return out


def first_profile_fields(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"profile": str(path), "exists": "false"}
    profile = pd.read_csv(path, sep="\t", nrows=1)
    out: dict[str, object] = {"profile": str(path), "exists": "true"}
    for col in [
        "profile_strategy",
        "adaptive_mode",
        "auto_exact_split_requested",
        "auto_exact_split_used",
        "auto_exact_split_source",
        "auto_exact_split_unavailable_reason",
        "low_extra_split_rescue_added_n",
        "probability_extra_mass_ratio",
    ]:
        value = profile[col].iloc[0] if col in profile.columns and len(profile) else ""
        out[col] = "" if pd.isna(value) else value
    return out


def method_row(path: Path, method: str) -> dict[str, str]:
    if not path.exists():
        return {}
    rows = list(csv.DictReader(path.open(newline=""), delimiter="\t"))
    for row in rows:
        if row.get("method") == method:
            return row
    return {}


def main() -> int:
    runtime_rows = []
    for sample, paths in SAMPLES.items():
        row = {"sample": sample}
        row.update(parse_time_log(paths["time_log"]))
        row.update(first_profile_fields(paths["profile"]))
        runtime_rows.append(row)

    raw_minco = method_row(RAW_SUMMARY, "minco_current_raw_default_gastrooral_source_abundance")
    raw_sylph = method_row(RAW_SUMMARY, "sylph_gtdb_r232_gastrooral_source_abundance")
    table_minco = method_row(TABLE_SUMMARY, "minco_current_universal_gastrooral_source_abundance")
    audit_rows = [
        {
            "metric": "raw_default_minco_mean_F1",
            "value": raw_minco.get("mean_F1", ""),
            "evidence": str(RAW_SUMMARY),
        },
        {
            "metric": "raw_default_minco_mean_L1_union_pp",
            "value": raw_minco.get("mean_L1_union_pp", ""),
            "evidence": str(RAW_SUMMARY),
        },
        {
            "metric": "raw_default_sylph_mean_F1",
            "value": raw_sylph.get("mean_F1", ""),
            "evidence": str(RAW_SUMMARY),
        },
        {
            "metric": "raw_default_sylph_mean_L1_union_pp",
            "value": raw_sylph.get("mean_L1_union_pp", ""),
            "evidence": str(RAW_SUMMARY),
        },
    ]
    if raw_minco and table_minco:
        raw_l1 = float(raw_minco["mean_L1_union_pp"])
        table_l1 = float(table_minco["mean_L1_union_pp"])
        audit_rows.append(
            {
                "metric": "raw_default_L1_delta_vs_table_mode_pp",
                "value": f"{raw_l1 - table_l1:.6f}",
                "evidence": f"{RAW_SUMMARY};{TABLE_SUMMARY}",
            }
        )

    runtime_fields = [
        "sample",
        "wall_time",
        "seconds",
        "peak_rss_kb",
        "peak_rss_gib",
        "user_time_s",
        "system_time_s",
        "profile",
        "exists",
        "profile_strategy",
        "adaptive_mode",
        "auto_exact_split_requested",
        "auto_exact_split_used",
        "auto_exact_split_source",
        "auto_exact_split_unavailable_reason",
        "low_extra_split_rescue_added_n",
        "probability_extra_mass_ratio",
        "time_log",
    ]
    write_tsv(RESULTS / "hmp_gastrooral_raw_default_runtime.tsv", runtime_rows, runtime_fields)
    write_tsv(RESULTS / "hmp_gastrooral_raw_default_audit.tsv", audit_rows, ["metric", "value", "evidence"])
    print(pd.DataFrame(runtime_rows).to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
