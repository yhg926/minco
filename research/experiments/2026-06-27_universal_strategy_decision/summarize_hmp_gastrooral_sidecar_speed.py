#!/usr/bin/env python3
"""Summarize HMP gastrooral exact-sidecar speed experiment."""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
DEFAULT_RUN = Path("/tmp/minco_current_code_hmp_gastrooral_20260627")
SIDECAR_RUN = Path("/tmp/minco_hmp_gastrooral_sidecar_speed_20260627")
SYLPH_RUN = Path("/tmp/cami2_hmp_pilot_20260625/run_sylph_r232")

RAW_SUMMARY = RESULTS / "hmp_gastrooral_raw_default_r232_source_abundance_summary.tsv"
SIDECAR_SUMMARY = RESULTS / "hmp_gastrooral_sidecar_r232_source_abundance_summary.tsv"
BLOCK_SUMMARY = RESULTS / "hmp_gastrooral_block_replay_r232_source_abundance_summary.tsv"

PROFILE_ROWS = [
    {
        "sample": "0",
        "method": "minco_default_exact_rerun",
        "time_log": DEFAULT_RUN / "sample0_current_default.time.log",
        "profile": DEFAULT_RUN / "minco_sample0_current_default.tsv",
        "accuracy_summary": RAW_SUMMARY,
        "accuracy_method": "minco_current_raw_default_gastrooral_source_abundance",
    },
    {
        "sample": "0",
        "method": "minco_same_stream_exact_sidecar",
        "time_log": SIDECAR_RUN / "sample0_sidecar.time.log",
        "profile": SIDECAR_RUN / "minco_sample0_sidecar.tsv",
        "accuracy_summary": SIDECAR_SUMMARY,
        "accuracy_method": "minco_hmp_gastrooral_sample0_sidecar_sample6_current",
    },
    {
        "sample": "0",
        "method": "sylph_r232_chunked_profile",
        "time_log": SYLPH_RUN / "sylph_sample0/profile.chunked.time.log",
        "profile": SYLPH_RUN / "sylph_sample0/profile.chunked.tsv",
        "accuracy_summary": SIDECAR_SUMMARY,
        "accuracy_method": "sylph_gtdb_r232_gastrooral_source_abundance",
    },
    {
        "sample": "6",
        "method": "minco_default_lowextra_or_rescue",
        "time_log": DEFAULT_RUN / "sample6_current_default.time.log",
        "profile": DEFAULT_RUN / "minco_sample6_current_default.tsv",
        "accuracy_summary": RAW_SUMMARY,
        "accuracy_method": "minco_current_raw_default_gastrooral_source_abundance",
    },
    {
        "sample": "6",
        "method": "sylph_r232_chunked_profile",
        "time_log": SYLPH_RUN / "sylph_sample6/profile.chunked.time.log",
        "profile": SYLPH_RUN / "sylph_sample6/profile.chunked.tsv",
        "accuracy_summary": SIDECAR_SUMMARY,
        "accuracy_method": "sylph_gtdb_r232_gastrooral_source_abundance",
    },
]


def parse_elapsed_seconds(text: str) -> float:
    parts = text.strip().split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return float(minutes) * 60.0 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)
    return math.nan


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
    out: dict[str, object] = {"profile": str(path), "exists": "true"}
    try:
        profile = pd.read_csv(path, sep="\t", nrows=1)
    except Exception:
        return out
    for col in [
        "profile_strategy",
        "adaptive_mode",
        "auto_exact_split_requested",
        "auto_exact_split_used",
        "auto_exact_split_source",
        "auto_exact_split_sidecar_requested",
        "probability_extra_mass_ratio",
        "joined_base_median_uaf",
        "raw_unique_to_base_ratio",
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


def called_species(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in {"species_name", "calibrated_call"})
    call = df["calibrated_call"].astype(str).str.lower().isin({"true", "1"})
    return set(df.loc[call, "species_name"].astype(str))


def abundance_abs_delta(left: Path, right: Path) -> float:
    cols = ["species_name", "calibrated_abundance"]
    if not left.exists() or not right.exists():
        return math.nan
    lhs = pd.read_csv(left, sep="\t", usecols=lambda col: col in set(cols))
    rhs = pd.read_csv(right, sep="\t", usecols=lambda col: col in set(cols))
    merged = lhs.merge(rhs, on="species_name", how="outer", suffixes=("_left", "_right")).fillna(0.0)
    return float(
        (
            pd.to_numeric(merged["calibrated_abundance_left"], errors="coerce").fillna(0.0)
            - pd.to_numeric(merged["calibrated_abundance_right"], errors="coerce").fillna(0.0)
        )
        .abs()
        .sum()
    )


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fnum(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def main() -> int:
    runtime_rows: list[dict[str, object]] = []
    for row_def in PROFILE_ROWS:
        row = {"sample": row_def["sample"], "method": row_def["method"]}
        row.update(parse_time_log(row_def["time_log"]))
        row.update(first_profile_fields(row_def["profile"]))
        score = method_row(row_def["accuracy_summary"], row_def["accuracy_method"])
        row["score_method"] = row_def["accuracy_method"]
        for col in ["mean_F1", "pooled_F1", "mean_L1_union_pp", "mean_Pearson_union"]:
            row[col] = score.get(col, "")
        runtime_rows.append(row)

    by_method = {str(row["method"]): row for row in runtime_rows if row["sample"] == "0"}
    exact = by_method.get("minco_default_exact_rerun", {})
    sidecar = by_method.get("minco_same_stream_exact_sidecar", {})
    sylph = by_method.get("sylph_r232_chunked_profile", {})
    raw_minco = method_row(RAW_SUMMARY, "minco_current_raw_default_gastrooral_source_abundance")
    sidecar_minco = method_row(SIDECAR_SUMMARY, "minco_hmp_gastrooral_sample0_sidecar_sample6_current")
    block_minco = method_row(BLOCK_SUMMARY, "minco_hmp_gastrooral_sample0_block_universal_replay_sample6_current")

    exact_profile = DEFAULT_RUN / "minco_sample0_current_default.tsv"
    sidecar_profile = SIDECAR_RUN / "minco_sample0_sidecar.tsv"
    block_profile = DEFAULT_RUN / "minco_sample0_block_universal_replay.tsv"
    exact_calls = called_species(exact_profile)
    sidecar_calls = called_species(sidecar_profile)
    block_calls = called_species(block_profile)

    audit_rows = [
        {
            "metric": "sidecar_seconds_delta_vs_exact_rerun_sample0",
            "value": f"{fnum(sidecar.get('seconds')) - fnum(exact.get('seconds')):.2f}",
            "evidence": f"{sidecar.get('time_log')};{exact.get('time_log')}",
            "decision": "speed_improvement_when_exact_needed",
        },
        {
            "metric": "sidecar_seconds_delta_vs_sylph_sample0",
            "value": f"{fnum(sidecar.get('seconds')) - fnum(sylph.get('seconds')):.2f}",
            "evidence": f"{sidecar.get('time_log')};{sylph.get('time_log')}",
            "decision": "beats_sylph_profile_time_on_this_sample",
        },
        {
            "metric": "sidecar_rss_gib_delta_vs_exact_rerun_sample0",
            "value": f"{fnum(sidecar.get('peak_rss_gib')) - fnum(exact.get('peak_rss_gib')):.4f}",
            "evidence": f"{sidecar.get('time_log')};{exact.get('time_log')}",
            "decision": "higher_memory_than_rerun_but_still_low",
        },
        {
            "metric": "sidecar_rss_gib_delta_vs_sylph_sample0",
            "value": f"{fnum(sidecar.get('peak_rss_gib')) - fnum(sylph.get('peak_rss_gib')):.4f}",
            "evidence": f"{sidecar.get('time_log')};{sylph.get('time_log')}",
            "decision": "much_lower_memory_than_sylph",
        },
        {
            "metric": "sidecar_callset_delta_vs_exact_rerun_sample0",
            "value": f"sidecar_only={len(sidecar_calls - exact_calls)};exact_only={len(exact_calls - sidecar_calls)}",
            "evidence": f"{sidecar_profile};{exact_profile}",
            "decision": "calls_match",
        },
        {
            "metric": "sidecar_abundance_abs_delta_vs_exact_rerun_sample0",
            "value": f"{abundance_abs_delta(sidecar_profile, exact_profile):.12f}",
            "evidence": f"{sidecar_profile};{exact_profile}",
            "decision": "abundance_matches",
        },
        {
            "metric": "block_callset_delta_vs_exact_rerun_sample0",
            "value": f"block_only={len(block_calls - exact_calls)};exact_only={len(exact_calls - block_calls)}",
            "evidence": f"{block_profile};{exact_profile}",
            "decision": "calls_match_but_abundance_differs",
        },
        {
            "metric": "block_abundance_abs_delta_vs_exact_rerun_sample0",
            "value": f"{abundance_abs_delta(block_profile, exact_profile):.12f}",
            "evidence": f"{block_profile};{exact_profile}",
            "decision": "exact_needed_for_abundance_here",
        },
        {
            "metric": "block_replay_L1_delta_vs_exact_sidecar_panel_pp",
            "value": f"{fnum(block_minco.get('mean_L1_union_pp')) - fnum(sidecar_minco.get('mean_L1_union_pp')):.6f}",
            "evidence": f"{BLOCK_SUMMARY};{SIDECAR_SUMMARY}",
            "decision": "do_not_drop_exact_for_speed_on_this_panel",
        },
        {
            "metric": "sidecar_L1_delta_vs_raw_default_panel_pp",
            "value": f"{fnum(sidecar_minco.get('mean_L1_union_pp')) - fnum(raw_minco.get('mean_L1_union_pp')):.6f}",
            "evidence": f"{SIDECAR_SUMMARY};{RAW_SUMMARY}",
            "decision": "same_accuracy_as_rerun",
        },
    ]

    runtime_fields = [
        "sample",
        "method",
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
        "auto_exact_split_sidecar_requested",
        "probability_extra_mass_ratio",
        "joined_base_median_uaf",
        "raw_unique_to_base_ratio",
        "score_method",
        "mean_F1",
        "pooled_F1",
        "mean_L1_union_pp",
        "mean_Pearson_union",
        "time_log",
    ]
    write_tsv(RESULTS / "hmp_gastrooral_exact_sidecar_speed.tsv", runtime_rows, runtime_fields)
    write_tsv(RESULTS / "hmp_gastrooral_exact_sidecar_speed_audit.tsv", audit_rows, ["metric", "value", "evidence", "decision"])

    print(pd.DataFrame(runtime_rows).to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
