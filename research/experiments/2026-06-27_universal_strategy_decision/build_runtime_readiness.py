#!/usr/bin/env python3
"""Build a runtime readiness gate for the current MinCO default.

This audit is deliberately separate from the accuracy-first release gate.  It
answers a narrower question: whether the no-expertise current default can be
claimed faster than Sylph on the timed comparable cases we have already
recorded.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

HMP_SPEED = RESULTS / "hmp_gastrooral_exact_sidecar_speed.tsv"
TOYMOUSE_RUNTIME = RESULTS / "exact_split_lowextra_skip_runtime.tsv"
TOYMOUSE_SCORE = RESULTS / "exact_split_lowextra_skip_sample6_score.tsv"
SIDECAR_POLICY = RESULTS / "exact_sidecar_policy_audit.tsv"
BLOCK_EXACT_SEMANTICS = RESULTS / "block_exact_split_semantics_audit.tsv"
CANDIDATE_RESTRICTED_EXACT = RESULTS / "candidate_restricted_exact_feasibility_audit.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def find_row(rows: list[dict[str, str]], predicate, label: str) -> dict[str, str]:
    matches = [row for row in rows if predicate(row)]
    if len(matches) != 1:
        raise SystemExit(f"expected exactly one row for {label}, found {len(matches)}")
    return matches[0]


def fnum(row: dict[str, str], key: str, default: float = math.nan) -> float:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return default
    return float(value)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def case_row(
    *,
    case_id: str,
    dataset: str,
    sample: str,
    minco_method: str,
    sylph_method: str,
    minco_seconds: float,
    sylph_seconds: float,
    minco_rss_gib: float,
    sylph_rss_gib: float,
    default_status: str,
    exact_behavior: str,
    accuracy_note: str,
    decision: str,
    source: str,
) -> dict[str, object]:
    seconds_delta = minco_seconds - sylph_seconds
    rss_delta = minco_rss_gib - sylph_rss_gib
    return {
        "case_id": case_id,
        "dataset": dataset,
        "sample": sample,
        "minco_method": minco_method,
        "sylph_method": sylph_method,
        "default_status": default_status,
        "exact_behavior": exact_behavior,
        "minco_seconds": f"{minco_seconds:.2f}",
        "sylph_seconds": f"{sylph_seconds:.2f}",
        "seconds_delta_minco_minus_sylph": f"{seconds_delta:.2f}",
        "minco_peak_rss_gib": f"{minco_rss_gib:.4f}",
        "sylph_peak_rss_gib": f"{sylph_rss_gib:.4f}",
        "rss_delta_minco_minus_sylph_gib": f"{rss_delta:.4f}",
        "faster_than_sylph": bool_text(seconds_delta < 0.0),
        "lower_memory_than_sylph": bool_text(rss_delta < 0.0),
        "accuracy_note": accuracy_note,
        "decision": decision,
        "source": source,
    }


def hmp_cases() -> list[dict[str, object]]:
    rows = read_tsv(HMP_SPEED)
    s0_sylph = find_row(
        rows,
        lambda row: row["sample"] == "0" and row["method"] == "sylph_r232_chunked_profile",
        "HMP gastrooral sample0 Sylph",
    )
    s6_sylph = find_row(
        rows,
        lambda row: row["sample"] == "6" and row["method"] == "sylph_r232_chunked_profile",
        "HMP gastrooral sample6 Sylph",
    )
    s0_default = find_row(
        rows,
        lambda row: row["sample"] == "0" and row["method"] == "minco_default_exact_rerun",
        "HMP gastrooral sample0 MinCO default",
    )
    s6_default = find_row(
        rows,
        lambda row: row["sample"] == "6" and row["method"] == "minco_default_lowextra_or_rescue",
        "HMP gastrooral sample6 MinCO default",
    )
    s0_sidecar = find_row(
        rows,
        lambda row: row["sample"] == "0" and row["method"] == "minco_same_stream_exact_sidecar",
        "HMP gastrooral sample0 MinCO sidecar",
    )
    panel_note = (
        "same two-sample panel metrics: MinCO F1 "
        f"{s0_default['mean_F1']} L1 {s0_default['mean_L1_union_pp']} pp; "
        f"Sylph F1 {s0_sylph['mean_F1']} L1 {s0_sylph['mean_L1_union_pp']} pp"
    )
    return [
        case_row(
            case_id="hmp_gastrooral_sample0_current_default_exact_rerun",
            dataset="hmp_gastrooral_r232",
            sample="0",
            minco_method="MinCO current default exact rerun",
            sylph_method="Sylph r232 chunked profile",
            minco_seconds=fnum(s0_default, "seconds"),
            sylph_seconds=fnum(s0_sylph, "seconds"),
            minco_rss_gib=fnum(s0_default, "peak_rss_gib"),
            sylph_rss_gib=fnum(s0_sylph, "peak_rss_gib"),
            default_status="current_default",
            exact_behavior="exact_rerun_used",
            accuracy_note=panel_note,
            decision="default_slower_than_sylph_here",
            source=str(HMP_SPEED.relative_to(NOTE_DIR)),
        ),
        case_row(
            case_id="hmp_gastrooral_sample6_current_default_no_exact",
            dataset="hmp_gastrooral_r232",
            sample="6",
            minco_method="MinCO current default no exact",
            sylph_method="Sylph r232 chunked profile",
            minco_seconds=fnum(s6_default, "seconds"),
            sylph_seconds=fnum(s6_sylph, "seconds"),
            minco_rss_gib=fnum(s6_default, "peak_rss_gib"),
            sylph_rss_gib=fnum(s6_sylph, "peak_rss_gib"),
            default_status="current_default",
            exact_behavior="exact_not_requested",
            accuracy_note=panel_note,
            decision="default_slightly_slower_than_sylph_here",
            source=str(HMP_SPEED.relative_to(NOTE_DIR)),
        ),
        case_row(
            case_id="hmp_gastrooral_sample0_optin_exact_sidecar",
            dataset="hmp_gastrooral_r232",
            sample="0",
            minco_method="MinCO opt-in same-stream exact sidecar",
            sylph_method="Sylph r232 chunked profile",
            minco_seconds=fnum(s0_sidecar, "seconds"),
            sylph_seconds=fnum(s0_sylph, "seconds"),
            minco_rss_gib=fnum(s0_sidecar, "peak_rss_gib"),
            sylph_rss_gib=fnum(s0_sylph, "peak_rss_gib"),
            default_status="opt_in_speed_substrate",
            exact_behavior="sidecar_used_and_matches_exact",
            accuracy_note=panel_note + "; sidecar call/abundance delta versus exact rerun is zero",
            decision="sidecar_faster_than_sylph_but_not_default",
            source=str(HMP_SPEED.relative_to(NOTE_DIR)),
        ),
    ]


def toymouse_cases() -> list[dict[str, object]]:
    runtime = read_tsv(TOYMOUSE_RUNTIME)
    score = read_tsv(TOYMOUSE_SCORE)
    default_runtime = find_row(
        runtime,
        lambda row: row["mode"] == "wrapper_lowextra_skip_default",
        "Toy Mouse sample6 MinCO low-extra default",
    )
    sidecar_runtime = find_row(
        runtime,
        lambda row: row["mode"] == "wrapper_same_stream_exact_sidecar",
        "Toy Mouse sample6 MinCO sidecar",
    )
    sylph_runtime = find_row(
        runtime,
        lambda row: row["mode"] == "Sylph sketch+profile",
        "Toy Mouse sample6 Sylph",
    )
    default_score = find_row(
        score,
        lambda row: row["method"] == "minco_lowextra_skip_autoexact_sample6",
        "Toy Mouse sample6 MinCO low-extra score",
    )
    sidecar_score = find_row(
        score,
        lambda row: row["method"] == "minco_sidecar_autoexact_sample6",
        "Toy Mouse sample6 MinCO sidecar score",
    )
    sylph_score = find_row(
        score,
        lambda row: row["method"] == "sylph_gtdb_profile_sample6",
        "Toy Mouse sample6 Sylph score",
    )
    default_note = (
        f"MinCO F1 {default_score['F1']} L1 {default_score['l1_pct_points']} pp; "
        f"Sylph F1 {sylph_score['F1']} L1 {sylph_score['l1_pct_points']} pp"
    )
    sidecar_note = (
        f"MinCO sidecar F1 {sidecar_score['F1']} L1 {sidecar_score['l1_pct_points']} pp; "
        f"Sylph F1 {sylph_score['F1']} L1 {sylph_score['l1_pct_points']} pp"
    )
    return [
        case_row(
            case_id="toymouse_sample6_current_default_lowextra_skip",
            dataset="cami2_toy_mouse",
            sample="6",
            minco_method="MinCO current default low-extra exact skip",
            sylph_method="Sylph GTDB profile",
            minco_seconds=fnum(default_runtime, "seconds"),
            sylph_seconds=fnum(sylph_runtime, "seconds"),
            minco_rss_gib=fnum(default_runtime, "peak_rss_gib"),
            sylph_rss_gib=fnum(sylph_runtime, "peak_rss_gib"),
            default_status="current_default",
            exact_behavior="exact_skipped_by_low_extra_rule",
            accuracy_note=default_note,
            decision="default_faster_but_accuracy_worse",
            source=str(TOYMOUSE_RUNTIME.relative_to(NOTE_DIR)),
        ),
        case_row(
            case_id="toymouse_sample6_optin_exact_sidecar",
            dataset="cami2_toy_mouse",
            sample="6",
            minco_method="MinCO opt-in same-stream exact sidecar",
            sylph_method="Sylph GTDB profile",
            minco_seconds=fnum(sidecar_runtime, "seconds"),
            sylph_seconds=fnum(sylph_runtime, "seconds"),
            minco_rss_gib=fnum(sidecar_runtime, "peak_rss_gib"),
            sylph_rss_gib=fnum(sylph_runtime, "peak_rss_gib"),
            default_status="opt_in_speed_substrate",
            exact_behavior="sidecar_used_but_exact_not_needed_by_current_default",
            accuracy_note=sidecar_note,
            decision="sidecar_slower_than_sylph_when_exact_not_needed",
            source=str(TOYMOUSE_RUNTIME.relative_to(NOTE_DIR)),
        ),
    ]


def policy_value(metric: str) -> str:
    if not SIDECAR_POLICY.exists():
        return ""
    rows = {row["metric"]: row for row in read_tsv(SIDECAR_POLICY)}
    return rows.get(metric, {}).get("value", "")


def block_exact_value(metric: str) -> str:
    if not BLOCK_EXACT_SEMANTICS.exists():
        return ""
    rows = {row["metric"]: row for row in read_tsv(BLOCK_EXACT_SEMANTICS)}
    return rows.get(metric, {}).get("value", "")


def candidate_exact_value(metric: str) -> str:
    if not CANDIDATE_RESTRICTED_EXACT.exists():
        return ""
    rows = {row["metric"]: row for row in read_tsv(CANDIDATE_RESTRICTED_EXACT)}
    return rows.get(metric, {}).get("value", "")


def main() -> int:
    cases = hmp_cases() + toymouse_cases()
    default_cases = [row for row in cases if row["default_status"] == "current_default"]
    opt_in_cases = [row for row in cases if row["default_status"] == "opt_in_speed_substrate"]

    def count_true(rows: list[dict[str, object]], field: str) -> int:
        return sum(1 for row in rows if row.get(field) == "true")

    default_faster = count_true(default_cases, "faster_than_sylph")
    default_lower_mem = count_true(default_cases, "lower_memory_than_sylph")
    opt_in_faster = count_true(opt_in_cases, "faster_than_sylph")
    opt_in_lower_mem = count_true(opt_in_cases, "lower_memory_than_sylph")
    sidecar_policy = policy_value("promotion_decision")
    block_exact_pairs = block_exact_value("not_interchangeable_pairs")
    candidate_exact_decision = candidate_exact_value("candidate_restricted_exact_decision")
    candidate_exact_keep = candidate_exact_value("min_exact_called_taxid_keep_fraction")

    audit_rows = [
        {
            "metric": "timed_default_cases",
            "value": len(default_cases),
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "small_current_timed_set",
        },
        {
            "metric": "timed_default_faster_than_sylph_cases",
            "value": f"{default_faster}/{len(default_cases)}",
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "default_speed_claim_not_supported",
        },
        {
            "metric": "timed_default_lower_memory_than_sylph_cases",
            "value": f"{default_lower_mem}/{len(default_cases)}",
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "memory_advantage_consistent_in_timed_cases",
        },
        {
            "metric": "timed_opt_in_sidecar_faster_than_sylph_cases",
            "value": f"{opt_in_faster}/{len(opt_in_cases)}",
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "sidecar_is_conditional_speed_substrate_not_default",
        },
        {
            "metric": "timed_opt_in_sidecar_lower_memory_than_sylph_cases",
            "value": f"{opt_in_lower_mem}/{len(opt_in_cases)}",
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "sidecar_memory_still_below_sylph_in_timed_cases",
        },
        {
            "metric": "sidecar_default_promotion_decision",
            "value": sidecar_policy or "missing_policy_audit",
            "evidence": "exact_sidecar_policy_audit.tsv",
            "decision": "do_not_promote_unconditional_sidecar_default",
        },
        {
            "metric": "block_exact_reuse_decision",
            "value": (
                f"not_interchangeable_pairs={block_exact_pairs}"
                if block_exact_pairs
                else "missing_block_exact_semantics_audit"
            ),
            "evidence": "block_exact_split_semantics_audit.tsv",
            "decision": "do_not_synthesize_block_split_from_exact_per_read",
        },
        {
            "metric": "candidate_restricted_exact_decision",
            "value": candidate_exact_decision or "missing_candidate_restricted_exact_audit",
            "evidence": "candidate_restricted_exact_feasibility_audit.tsv",
            "decision": "do_not_promote_candidate_only_exact_rerun",
        },
        {
            "metric": "candidate_restricted_exact_min_keep_fraction",
            "value": candidate_exact_keep or "missing_candidate_restricted_exact_audit",
            "evidence": "candidate_restricted_exact_feasibility_audit.tsv",
            "decision": "posthoc_output_reduction_not_scan_cost_reduction",
        },
        {
            "metric": "default_faster_than_sylph_claim_supported",
            "value": bool_text(default_faster == len(default_cases) and len(default_cases) > 0),
            "evidence": "runtime_strategy_cases.tsv",
            "decision": "known_gap_until_conditional_exact_sidecar_or_other_speed_fix",
        },
        {
            "metric": "runtime_next_target",
            "value": "conditional_exact_sidecar_or_lower_cost_dual_channel_processing",
            "evidence": (
                "HMP sample0 exact sidecar win; Toy Mouse sample6 eager sidecar "
                "penalty; block/exact split rows are not interchangeable; "
                "candidate-only exact rerun is not proven safe from aggregate outputs"
            ),
            "decision": "preserve current default accuracy while removing exact-rerun latency",
        },
    ]

    fields = [
        "case_id",
        "dataset",
        "sample",
        "minco_method",
        "sylph_method",
        "default_status",
        "exact_behavior",
        "minco_seconds",
        "sylph_seconds",
        "seconds_delta_minco_minus_sylph",
        "minco_peak_rss_gib",
        "sylph_peak_rss_gib",
        "rss_delta_minco_minus_sylph_gib",
        "faster_than_sylph",
        "lower_memory_than_sylph",
        "accuracy_note",
        "decision",
        "source",
    ]
    write_tsv(RESULTS / "runtime_strategy_cases.tsv", cases, fields)
    write_tsv(RESULTS / "runtime_readiness.tsv", audit_rows, ["metric", "value", "evidence", "decision"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
