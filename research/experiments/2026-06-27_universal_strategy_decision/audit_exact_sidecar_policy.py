#!/usr/bin/env python3
"""Audit whether exact sidecar generation can be promoted by policy.

The default wrapper can only know the full block-mode exact trigger after the
initial split/unique passes. Eager exact-sidecar generation is faster when
exact is truly needed, but it adds work when exact is skipped. This audit keeps
the current policy boundary explicit and reproducible from available profile
outputs.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
RESULTS = EXP / "results"

EXACT_TRIGGER = 0.10

PROFILE_INSTANCES = [
    {
        "dataset": "hmp_gastrooral",
        "sample": "0",
        "profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv",
        "block_profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_block_universal_replay.tsv",
        "class": "exact_needed_with_sidecar_timing",
        "evidence": "hmp_gastrooral_exact_sidecar_speed_audit.tsv",
    },
    {
        "dataset": "hmp_gastrooral",
        "sample": "6",
        "profile": "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv",
        "class": "exact_not_requested_high_extra",
        "evidence": "hmp_gastrooral_exact_sidecar_speed.tsv",
    },
    {
        "dataset": "cami2_toy_mouse",
        "sample": "6",
        "profile": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_lowextra_skip_p16.tsv",
        "sidecar_profile": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16.tsv",
        "class": "exact_skipped_low_extra_rescue",
        "evidence": "exact_split_lowextra_skip_sample6_score.tsv",
    },
    {
        "dataset": "hmp_airskin",
        "sample": "28",
        "profile": "/tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv",
        "class": "exact_not_requested_guard",
        "evidence": "hmp_airskin28_r232_source_abundance_summary.tsv",
    },
    {
        "dataset": "hmp_airskin",
        "sample": "6",
        "profile": "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv",
        "class": "exact_supplied_table_no_raw_sidecar_timing",
        "evidence": "hmp_current_refresh_r232_source_abundance_summary.tsv",
    },
    {
        "dataset": "hmp_airskin",
        "sample": "11",
        "profile": "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv",
        "class": "exact_supplied_table_no_raw_sidecar_timing",
        "evidence": "hmp_current_refresh_r232_source_abundance_summary.tsv",
    },
]


def read_first_row(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", nrows=1)
    if df.empty:
        return {}
    return df.iloc[0].to_dict()


def value(row: dict[str, object], key: str, default: object = "") -> object:
    val = row.get(key, default)
    if pd.isna(val):
        return default
    return val


def fnum(row: dict[str, object], key: str, default: float = math.nan) -> float:
    try:
        return float(value(row, key, default))
    except (TypeError, ValueError):
        return default


def bval(row: dict[str, object], key: str) -> bool | None:
    raw = str(value(row, key, "")).strip().lower()
    if raw in {"true", "1", "yes"}:
        return True
    if raw in {"false", "0", "no"}:
        return False
    return None


def called_species(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in {"species_name", "calibrated_call"})
    if "calibrated_call" not in df.columns:
        return set()
    calls = df["calibrated_call"].astype(str).str.lower().isin({"true", "1"})
    return set(df.loc[calls, "species_name"].astype(str))


def abundance_abs_delta(left: Path, right: Path) -> float:
    cols = {"species_name", "calibrated_abundance"}
    if not left.exists() or not right.exists():
        return math.nan
    lhs = pd.read_csv(left, sep="\t", usecols=lambda col: col in cols)
    rhs = pd.read_csv(right, sep="\t", usecols=lambda col: col in cols)
    merged = lhs.merge(rhs, on="species_name", how="outer", suffixes=("_left", "_right")).fillna(0.0)
    return float(
        (
            pd.to_numeric(merged["calibrated_abundance_left"], errors="coerce").fillna(0.0)
            - pd.to_numeric(merged["calibrated_abundance_right"], errors="coerce").fillna(0.0)
        )
        .abs()
        .sum()
    )


def read_key_value(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {row["metric"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_method_rows(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    key = "method" if rows and "method" in rows[0] else "mode"
    return {row[key]: row for row in rows}


def parse_elapsed_seconds(text: str) -> float:
    parts = str(text).strip().split(":")
    if len(parts) == 2:
        minutes, seconds = parts
        return float(minutes) * 60.0 + float(seconds)
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)
    return math.nan


def parse_time_log(path: Path) -> dict[str, object]:
    out: dict[str, object] = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            wall = line.rsplit(": ", 1)[1].strip()
            out["wall_time"] = wall
            out["seconds"] = f"{parse_elapsed_seconds(wall):.2f}"
        elif "Maximum resident set size" in line:
            match = re.search(r"(\d+)", line)
            if match:
                kb = int(match.group(1))
                out["peak_rss_kb"] = str(kb)
                out["peak_rss_gib"] = f"{kb / 1024.0 / 1024.0:.4f}"
    return out


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def classify_policy(row: dict[str, object]) -> str:
    if row["profile_exists"] != "true":
        return "missing_profile"
    p_extra = float(row["block_p_extra_mass_ratio"])
    guard = row["block_guard_passed"] == "true"
    low_extra_added = int(float(row["low_extra_split_rescue_added_n"]))
    low_extra_pass = low_extra_added == 0
    if p_extra > EXACT_TRIGGER:
        return "avoid_exact_high_extra"
    if not guard:
        return "avoid_exact_guard"
    if not low_extra_pass:
        return "avoid_exact_low_extra_rescue"
    return "exact_candidate"


def main() -> int:
    hmp_sidecar_audit = read_key_value(RESULTS / "hmp_gastrooral_exact_sidecar_speed_audit.tsv")
    hmp_sidecar_speed = read_method_rows(RESULTS / "hmp_gastrooral_exact_sidecar_speed.tsv")
    toy_runtime = read_method_rows(RESULTS / "exact_split_lowextra_skip_runtime.tsv")
    toy_score = read_method_rows(RESULTS / "exact_split_lowextra_skip_sample6_score.tsv")

    rows: list[dict[str, object]] = []
    for spec in PROFILE_INSTANCES:
        profile = Path(spec["profile"])
        first = read_first_row(profile)
        block_first = read_first_row(Path(spec.get("block_profile", spec["profile"])))
        block_p_extra = fnum(first, "auto_exact_split_block_p_extra_mass_ratio", fnum(block_first, "probability_extra_mass_ratio"))
        if math.isnan(block_p_extra):
            block_p_extra = fnum(first, "probability_extra_mass_ratio")
        out: dict[str, object] = {
            "dataset": spec["dataset"],
            "sample": spec["sample"],
            "class": spec["class"],
            "profile": str(profile),
            "profile_exists": str(bool(first)).lower(),
            "profile_strategy": value(first, "profile_strategy"),
            "adaptive_mode": value(block_first or first, "adaptive_mode"),
            "block_p_extra_mass_ratio": f"{block_p_extra:.12g}" if not math.isnan(block_p_extra) else "",
            "probability_count_ratio": value(block_first or first, "probability_count_ratio"),
            "joined_base_n": value(block_first or first, "joined_base_n"),
            "probability_gate_n": value(block_first or first, "probability_gate_n"),
            "probability_extra_n": value(block_first or first, "probability_extra_n"),
            "joined_base_median_uaf": value(block_first or first, "joined_base_median_uaf"),
            "raw_unique95_n": value(block_first or first, "raw_unique95_n"),
            "raw_unique_to_base_ratio": value(block_first or first, "raw_unique_to_base_ratio"),
            "low_extra_split_rescue_added_n": value(block_first or first, "low_extra_split_rescue_added_n", 0),
            "block_guard_passed": str(
                bool(
                    fnum(block_first or first, "joined_base_median_uaf") < 0.35
                    or fnum(block_first or first, "raw_unique_to_base_ratio") >= 0.80
                )
            ).lower(),
            "auto_exact_split_requested": bval(first, "auto_exact_split_requested"),
            "auto_exact_split_used": bval(first, "auto_exact_split_used"),
            "auto_exact_split_source": value(first, "auto_exact_split_source"),
            "auto_exact_split_unavailable_reason": value(first, "auto_exact_split_unavailable_reason"),
            "evidence": spec["evidence"],
        }
        out["policy_from_block_features"] = classify_policy(out)

        if spec["dataset"] == "hmp_gastrooral" and spec["sample"] == "0":
            out["sidecar_seconds_delta_vs_exact_rerun"] = hmp_sidecar_audit.get(
                "sidecar_seconds_delta_vs_exact_rerun_sample0", {}
            ).get("value", "")
            out["sidecar_seconds_delta_vs_sylph"] = hmp_sidecar_audit.get(
                "sidecar_seconds_delta_vs_sylph_sample0", {}
            ).get("value", "")
            out["block_L1_delta_vs_exact_sidecar_pp"] = hmp_sidecar_audit.get(
                "block_replay_L1_delta_vs_exact_sidecar_panel_pp", {}
            ).get("value", "")
            out["exact_sidecar_effect"] = "sidecar_matches_exact_and_block_abundance_worse"
        elif spec["dataset"] == "cami2_toy_mouse":
            low = toy_runtime.get("wrapper_lowextra_skip_default", {})
            side = toy_runtime.get("wrapper_same_stream_exact_sidecar", {})
            low_score = toy_score.get("minco_lowextra_skip_autoexact_sample6", {})
            side_score = toy_score.get("minco_sidecar_autoexact_sample6", {})
            if low and side:
                out["sidecar_seconds_delta_vs_exact_skip"] = f"{float(side['seconds']) - float(low['seconds']):.2f}"
                out["sidecar_rss_gib_delta_vs_exact_skip"] = f"{float(side['peak_rss_gib']) - float(low['peak_rss_gib']):.4f}"
            if low_score and side_score:
                out["sidecar_L1_delta_vs_exact_skip_pp"] = (
                    f"{float(side_score['l1_pct_points']) - float(low_score['l1_pct_points']):.6f}"
                )
            out["exact_sidecar_effect"] = "sidecar_same_calls_slower_and_slightly_worse_L1"
        elif spec["dataset"] == "hmp_gastrooral" and spec["sample"] == "6":
            minco = hmp_sidecar_speed.get("minco_default_lowextra_or_rescue", {})
            sylph = hmp_sidecar_speed.get("sylph_r232_chunked_profile", {})
            if minco and sylph:
                out["minco_seconds_delta_vs_sylph"] = f"{float(minco['seconds']) - float(sylph['seconds']):.2f}"
                out["minco_rss_gib_delta_vs_sylph"] = f"{float(minco['peak_rss_gib']) - float(sylph['peak_rss_gib']):.4f}"
            out["exact_sidecar_effect"] = "exact_not_requested_no_sidecar_timing"
        else:
            out["exact_sidecar_effect"] = "no_raw_sidecar_pair_available"
        rows.append(out)

    exact_candidates = [row for row in rows if row["policy_from_block_features"] == "exact_candidate"]
    timed_exact_candidates = [
        row for row in exact_candidates
        if row.get("sidecar_seconds_delta_vs_exact_rerun", "") not in {"", "nan"}
    ]
    avoid_rows = [row for row in rows if row["policy_from_block_features"].startswith("avoid_exact")]
    contradicted_avoids = [
        row for row in avoid_rows
        if row.get("exact_sidecar_effect") == "sidecar_matches_exact_and_block_abundance_worse"
    ]
    toy_overhead = next(
        (
            row for row in rows
            if row["dataset"] == "cami2_toy_mouse" and row["sample"] == "6"
        ),
        {},
    )
    audit_rows = [
        {
            "metric": "evaluated_instances",
            "value": len(rows),
            "evidence": "exact_sidecar_policy_instances.tsv",
            "decision": "diagnostic_dataset_small",
        },
        {
            "metric": "block_policy_exact_candidates",
            "value": len(exact_candidates),
            "evidence": "policy_from_block_features==exact_candidate",
            "decision": "exact_candidates_require_sidecar_or_rerun",
        },
        {
            "metric": "timed_exact_candidates_with_sidecar_win",
            "value": sum(float(row.get("sidecar_seconds_delta_vs_exact_rerun", "nan")) < 0 for row in timed_exact_candidates),
            "evidence": "HMP gastrooral sample0 sidecar timing",
            "decision": "sidecar_helps_when_exact_needed_in_available_timed_case",
        },
        {
            "metric": "avoid_exact_rows",
            "value": len(avoid_rows),
            "evidence": "policy_from_block_features starts avoid_exact",
            "decision": "avoid_sidecar_when_exact_not_needed",
        },
        {
            "metric": "avoid_exact_contradictions",
            "value": len(contradicted_avoids),
            "evidence": "exact_sidecar_policy_instances.tsv",
            "decision": "no_observed_contradiction_but_small_n",
        },
        {
            "metric": "toy_mouse_lowextra_eager_sidecar_penalty",
            "value": (
                f"seconds_delta={toy_overhead.get('sidecar_seconds_delta_vs_exact_skip', '')};"
                f"rss_gib_delta={toy_overhead.get('sidecar_rss_gib_delta_vs_exact_skip', '')};"
                f"L1_delta={toy_overhead.get('sidecar_L1_delta_vs_exact_skip_pp', '')}"
            ),
            "evidence": "exact_split_lowextra_skip_runtime.tsv;exact_split_lowextra_skip_sample6_score.tsv",
            "decision": "unconditional_sidecar_rejected",
        },
        {
            "metric": "promotion_decision",
            "value": "do_not_promote_auto_sidecar_default_without_preflight",
            "evidence": "exact_sidecar_policy_instances.tsv;hmp_gastrooral_exact_sidecar_speed_audit.tsv",
            "decision": "need validated preflight or cheaper conditional sidecar before default",
        },
    ]

    fields = [
        "dataset",
        "sample",
        "class",
        "profile_exists",
        "profile_strategy",
        "adaptive_mode",
        "block_p_extra_mass_ratio",
        "probability_count_ratio",
        "joined_base_n",
        "probability_gate_n",
        "probability_extra_n",
        "joined_base_median_uaf",
        "raw_unique95_n",
        "raw_unique_to_base_ratio",
        "low_extra_split_rescue_added_n",
        "block_guard_passed",
        "auto_exact_split_requested",
        "auto_exact_split_used",
        "auto_exact_split_source",
        "auto_exact_split_unavailable_reason",
        "policy_from_block_features",
        "exact_sidecar_effect",
        "sidecar_seconds_delta_vs_exact_rerun",
        "sidecar_seconds_delta_vs_sylph",
        "sidecar_seconds_delta_vs_exact_skip",
        "sidecar_rss_gib_delta_vs_exact_skip",
        "sidecar_L1_delta_vs_exact_skip_pp",
        "block_L1_delta_vs_exact_sidecar_pp",
        "minco_seconds_delta_vs_sylph",
        "minco_rss_gib_delta_vs_sylph",
        "evidence",
        "profile",
    ]
    write_tsv(RESULTS / "exact_sidecar_policy_instances.tsv", rows, fields)
    write_tsv(RESULTS / "exact_sidecar_policy_audit.tsv", audit_rows, ["metric", "value", "evidence", "decision"])
    print(pd.DataFrame(rows)[fields].to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
