#!/usr/bin/env python3
"""Audit whether exact per-read split output can replace block split output.

The exact sidecar speed path tempts a second shortcut: compute exact per-read
best-diff-split evidence and derive the normal block-mode split evidence from
it.  This script checks that assumption against cached real tables.  A mismatch
does not mean either table is wrong; it means density-block grouping is a real
part of the current default evidence, so exact-per-read rows are not a safe
drop-in replacement for block-mode rows.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

PAIRS = [
    {
        "dataset": "hmp_gastrooral_r232",
        "sample": "0",
        "block": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work/minco.best_diff_split.unfiltered.tsv",
        "exact": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work/minco.best_diff_split.exact.unfiltered.tsv",
    },
    {
        "dataset": "cami2_toy_mouse",
        "sample": "6",
        "block": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work/minco.best_diff_split.unfiltered.tsv",
        "exact": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work/minco.best_diff_split.exact.unfiltered.tsv",
    },
    {
        "dataset": "toy_mouse_first50k",
        "sample": "0_first50k",
        "block": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sidecar_work/minco.best_diff_split.unfiltered.tsv",
        "exact": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sidecar_work/minco.best_diff_split.exact.unfiltered.tsv",
    },
]

COMPARE_COLS = [
    "ANI",
    "XnY_ctx",
    "Raw_XnY_ctx",
    "N_diff_obj",
    "N_diff_obj_section",
    "Real_min_align_fraction",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
    "Normalized_abundance_depth",
]


def read_table(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ref = row.get("Ref", "")
            if ref:
                rows[ref] = row
    return rows


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return 0.0
    try:
        out = float(value)
    except ValueError:
        return 0.0
    if not math.isfinite(out):
        return 0.0
    return out


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def summarize_pair(spec: dict[str, str]) -> tuple[dict[str, object], list[dict[str, object]]]:
    block_path = Path(spec["block"])
    exact_path = Path(spec["exact"])
    block = read_table(block_path)
    exact = read_table(exact_path)
    shared = sorted(set(block) & set(exact))
    block_only = sorted(set(block) - set(exact))
    exact_only = sorted(set(exact) - set(block))
    if not block or not exact:
        return (
            {
                "dataset": spec["dataset"],
                "sample": spec["sample"],
                "block_exists": str(bool(block)).lower(),
                "exact_exists": str(bool(exact)).lower(),
                "block_rows": len(block),
                "exact_rows": len(exact),
                "shared_refs": len(shared),
                "block_only_refs": len(block_only),
                "exact_only_refs": len(exact_only),
                "decision": "missing_input",
                "block_path": str(block_path),
                "exact_path": str(exact_path),
            },
            [],
        )

    mismatches: dict[str, int] = {col: 0 for col in COMPARE_COLS}
    max_abs: dict[str, float] = {col: 0.0 for col in COMPARE_COLS}
    norm_abs_delta = 0.0
    mean_depth_abs_delta = 0.0
    top_rows: list[dict[str, object]] = []
    for ref in shared:
        b = block[ref]
        e = exact[ref]
        ref_norm_delta = abs(
            number(b, "Normalized_abundance_depth")
            - number(e, "Normalized_abundance_depth")
        )
        ref_mean_delta = abs(number(b, "Ref_mean_depth") - number(e, "Ref_mean_depth"))
        norm_abs_delta += ref_norm_delta
        mean_depth_abs_delta += ref_mean_delta
        for col in COMPARE_COLS:
            delta = abs(number(b, col) - number(e, col))
            if delta > 1e-12:
                mismatches[col] += 1
                if delta > max_abs[col]:
                    max_abs[col] = delta
        top_rows.append(
            {
                "dataset": spec["dataset"],
                "sample": spec["sample"],
                "Ref": ref,
                "block_XnY_ctx": number(b, "XnY_ctx"),
                "exact_XnY_ctx": number(e, "XnY_ctx"),
                "block_Ref_breadth": number(b, "Ref_breadth"),
                "exact_Ref_breadth": number(e, "Ref_breadth"),
                "block_Ref_mean_depth": number(b, "Ref_mean_depth"),
                "exact_Ref_mean_depth": number(e, "Ref_mean_depth"),
                "block_Normalized_abundance_depth": number(b, "Normalized_abundance_depth"),
                "exact_Normalized_abundance_depth": number(e, "Normalized_abundance_depth"),
                "normalized_abundance_abs_delta": ref_norm_delta,
                "ref_mean_depth_abs_delta": ref_mean_delta,
            }
        )

    top_rows.sort(
        key=lambda row: (
            float(row["normalized_abundance_abs_delta"]),
            float(row["ref_mean_depth_abs_delta"]),
        ),
        reverse=True,
    )
    top_rows = top_rows[:20]

    compared = max(len(shared), 1)
    summary: dict[str, object] = {
        "dataset": spec["dataset"],
        "sample": spec["sample"],
        "block_exists": "true",
        "exact_exists": "true",
        "block_rows": len(block),
        "exact_rows": len(exact),
        "shared_refs": len(shared),
        "block_only_refs": len(block_only),
        "exact_only_refs": len(exact_only),
        "norm_abundance_abs_delta_shared": f"{norm_abs_delta:.12g}",
        "ref_mean_depth_abs_delta_shared": f"{mean_depth_abs_delta:.12g}",
        "any_compared_col_mismatch": str(any(v > 0 for v in mismatches.values())).lower(),
        "decision": (
            "not_interchangeable"
            if block_only or exact_only or any(v > 0 for v in mismatches.values())
            else "interchangeable_on_this_pair"
        ),
        "block_path": str(block_path),
        "exact_path": str(exact_path),
    }
    for col in COMPARE_COLS:
        summary[f"{col}_mismatches"] = mismatches[col]
        summary[f"{col}_mismatch_fraction"] = f"{mismatches[col] / compared:.8f}"
        summary[f"{col}_max_abs_delta"] = f"{max_abs[col]:.12g}"
    return summary, top_rows


def main() -> int:
    summaries: list[dict[str, object]] = []
    top_rows: list[dict[str, object]] = []
    for spec in PAIRS:
        summary, top = summarize_pair(spec)
        summaries.append(summary)
        top_rows.extend(top)

    evaluated = [row for row in summaries if row["decision"] != "missing_input"]
    rejected = [row for row in evaluated if row["decision"] == "not_interchangeable"]
    audit_rows = [
        {
            "metric": "evaluated_pairs",
            "value": len(evaluated),
            "evidence": "block_exact_split_semantics.tsv",
            "decision": "cached_real_tables",
        },
        {
            "metric": "not_interchangeable_pairs",
            "value": f"{len(rejected)}/{len(evaluated)}",
            "evidence": "block_exact_split_semantics.tsv",
            "decision": "do_not_synthesize_block_split_from_exact_per_read",
        },
        {
            "metric": "max_shared_normalized_abundance_abs_delta",
            "value": (
                f"{max(float(row.get('norm_abundance_abs_delta_shared', 0.0)) for row in evaluated):.12g}"
                if evaluated
                else "NA"
            ),
            "evidence": "block_exact_split_semantics.tsv",
            "decision": "grouping_changes_abundance_evidence",
        },
        {
            "metric": "speed_optimization_boundary",
            "value": "exact sidecar can replace exact rerun, but exact per-read rows cannot replace block split rows",
            "evidence": "block_exact_split_semantics.tsv; command_ani.c density-unit grouping",
            "decision": "next_speed_fix_must_preserve_both_evidence_channels_or_predict_exact_need",
        },
    ]

    summary_fields = [
        "dataset",
        "sample",
        "block_exists",
        "exact_exists",
        "block_rows",
        "exact_rows",
        "shared_refs",
        "block_only_refs",
        "exact_only_refs",
        "norm_abundance_abs_delta_shared",
        "ref_mean_depth_abs_delta_shared",
        "any_compared_col_mismatch",
        "decision",
        *[
            field
            for col in COMPARE_COLS
            for field in (
                f"{col}_mismatches",
                f"{col}_mismatch_fraction",
                f"{col}_max_abs_delta",
            )
        ],
        "block_path",
        "exact_path",
    ]
    top_fields = [
        "dataset",
        "sample",
        "Ref",
        "block_XnY_ctx",
        "exact_XnY_ctx",
        "block_Ref_breadth",
        "exact_Ref_breadth",
        "block_Ref_mean_depth",
        "exact_Ref_mean_depth",
        "block_Normalized_abundance_depth",
        "exact_Normalized_abundance_depth",
        "normalized_abundance_abs_delta",
        "ref_mean_depth_abs_delta",
    ]
    write_tsv(RESULTS / "block_exact_split_semantics.tsv", summaries, summary_fields)
    write_tsv(RESULTS / "block_exact_split_top_deltas.tsv", top_rows, top_fields)
    write_tsv(
        RESULTS / "block_exact_split_semantics_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
