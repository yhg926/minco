#!/usr/bin/env python3
"""Import cross-domain edge-EM evidence into the universal-strategy note."""

from __future__ import annotations

import csv
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"
SOURCE_EXP = REPO_ROOT / "research/experiments/2026-06-24_cross_domain_spot_comparison"
SOURCE_RESULTS = SOURCE_EXP / "results"
FAIR_METHODS = SOURCE_RESULTS / "fair_full_l1_methods.tsv"
ADAPTIVE_DELTA = SOURCE_RESULTS / "adaptive_edge_em_delta.tsv"


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


def fnum(value: object) -> float:
    value = str(value)
    if value in {"", "NA", "nan", "None"}:
        return float("nan")
    return float(value)


def mean(values: list[float]) -> float:
    vals = [value for value in values if value == value]
    return sum(vals) / len(vals) if vals else float("nan")


def fmt(value: float) -> str:
    return "" if value != value else f"{value:.12g}"


def summarize_methods() -> list[dict[str, object]]:
    wanted = {
        "Sylph",
        "MinCO current best",
        "MinCO S2000 edge marker beta=0",
        "MinCO adaptive edge-EM trigger",
    }
    rows = []
    for row in read_tsv(FAIR_METHODS):
        method = row["method"]
        if not any(method.startswith(prefix) for prefix in wanted):
            continue
        rows.append(
            {
                "dataset": row["dataset"],
                "sample": row["sample"],
                "truth_scope": row["truth_scope"],
                "method": method,
                "method_role": row["method_role"],
                "F1": row["F1"],
                "L1_pp": row["l1_pct_points"],
                "Pearson": row["pearson"],
                "TP": row["TP"],
                "FP": row["FP"],
                "FN": row["FN"],
                "adaptive_beta": row.get("adaptive_beta", ""),
                "adaptive_reason": row.get("adaptive_reason", ""),
                "source": str(FAIR_METHODS.relative_to(REPO_ROOT)),
            }
        )
    return rows


def summarize_policy() -> list[dict[str, object]]:
    delta = read_tsv(ADAPTIVE_DELTA)
    comparisons = [
        (
            "adaptive_edge_em_vs_edge_marker",
            "adaptive_minus_edge_marker_l1",
            "adaptive_minus_edge_marker_F1",
            "narrow positive: improves edge-marker abundance without F1 cost",
        ),
        (
            "adaptive_edge_em_vs_current_best_minco",
            "adaptive_minus_current_best_l1",
            "adaptive_minus_current_best_F1",
            "not promotable: F1 is lower on CAMI3 and marine",
        ),
        (
            "adaptive_edge_em_vs_sylph",
            "adaptive_minus_sylph_l1",
            "adaptive_minus_sylph_F1",
            "not promotable: F1 is lower on all three spot samples",
        ),
    ]
    rows: list[dict[str, object]] = []
    for label, l1_field, f1_field, decision in comparisons:
        l1 = [fnum(row[l1_field]) for row in delta]
        f1 = [fnum(row[f1_field]) for row in delta]
        rows.append(
            {
                "comparison": label,
                "datasets": ",".join(row["dataset"] for row in delta),
                "sample_count": len(delta),
                "mean_delta_L1_pp": fmt(mean(l1)),
                "max_delta_L1_pp": fmt(max(l1)),
                "min_delta_L1_pp": fmt(min(l1)),
                "improved_L1_count": sum(1 for value in l1 if value < 0.0),
                "unchanged_L1_count": sum(1 for value in l1 if value == 0.0),
                "worsened_L1_count": sum(1 for value in l1 if value > 0.0),
                "mean_delta_F1": fmt(mean(f1)),
                "max_delta_F1": fmt(max(f1)),
                "min_delta_F1": fmt(min(f1)),
                "improved_F1_count": sum(1 for value in f1 if value > 0.0),
                "unchanged_F1_count": sum(1 for value in f1 if value == 0.0),
                "worsened_F1_count": sum(1 for value in f1 if value < 0.0),
                "decision": decision,
                "source": str(ADAPTIVE_DELTA.relative_to(REPO_ROOT)),
            }
        )
    return rows


def main() -> int:
    methods = summarize_methods()
    write_tsv(
        RESULTS / "edge_em_cross_domain_methods.tsv",
        methods,
        [
            "dataset",
            "sample",
            "truth_scope",
            "method",
            "method_role",
            "F1",
            "L1_pp",
            "Pearson",
            "TP",
            "FP",
            "FN",
            "adaptive_beta",
            "adaptive_reason",
            "source",
        ],
    )

    policy = summarize_policy()
    write_tsv(
        RESULTS / "edge_em_cross_domain_policy_summary.tsv",
        policy,
        [
            "comparison",
            "datasets",
            "sample_count",
            "mean_delta_L1_pp",
            "max_delta_L1_pp",
            "min_delta_L1_pp",
            "improved_L1_count",
            "unchanged_L1_count",
            "worsened_L1_count",
            "mean_delta_F1",
            "max_delta_F1",
            "min_delta_F1",
            "improved_F1_count",
            "unchanged_F1_count",
            "worsened_F1_count",
            "decision",
            "source",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
