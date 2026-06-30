#!/usr/bin/env python3
"""Summarize whether fixed-call abundance variants are safe defaults."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
SCORES = RESULTS / "cross_panel_abundance_variant_scores.tsv"
PANEL_SUMMARY = RESULTS / "cross_panel_abundance_variant_overall.tsv"
OUT = RESULTS / "abundance_variant_safety_audit.tsv"

CURRENT = "current_calibrated_abundance"
TRIVIAL_EQUIVALENTS = {"current_calibrated_raw", "current_raw_power_p1"}


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = ["metric", "value", "evidence", "decision"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    scores = pd.read_csv(SCORES, sep="\t")
    panel = pd.read_csv(PANEL_SUMMARY, sep="\t")

    current = scores.loc[
        scores["method"].eq(CURRENT),
        ["panel", "sample", "L1_union_pp", "L1_truth_only_pp"],
    ].rename(columns={"L1_union_pp": "current_L1_union_pp", "L1_truth_only_pp": "current_L1_truth_only_pp"})
    work = scores.merge(current, on=["panel", "sample"], how="left")
    work["official_L1_pp"] = work.apply(
        lambda row: row["L1_truth_only_pp"]
        if row["panel"] == "cami2_toy_mouse_gut"
        else row["L1_union_pp"],
        axis=1,
    )
    work["current_official_L1_pp"] = work.apply(
        lambda row: row["current_L1_truth_only_pp"]
        if row["panel"] == "cami2_toy_mouse_gut"
        else row["current_L1_union_pp"],
        axis=1,
    )
    work["delta_current_L1_pp"] = work["official_L1_pp"] - work["current_official_L1_pp"]
    noncurrent = work.loc[~work["method"].eq(CURRENT)].copy()
    nontrivial = noncurrent.loc[~noncurrent["method"].isin(TRIVIAL_EQUIVALENTS)].copy()

    sample_agg = (
        nontrivial.groupby("method", as_index=False)
        .agg(
            mean_sample_delta_L1_pp=("delta_current_L1_pp", "mean"),
            max_sample_worse_L1_pp=("delta_current_L1_pp", "max"),
            improved_samples=("delta_current_L1_pp", lambda s: int((s < -1e-9).sum())),
            worsened_samples=("delta_current_L1_pp", lambda s: int((s > 1e-9).sum())),
            sample_count=("delta_current_L1_pp", "size"),
        )
        .sort_values(["mean_sample_delta_L1_pp", "max_sample_worse_L1_pp"], kind="mergesort")
    )
    best_sample = sample_agg.iloc[0]
    strict_sample = sample_agg.loc[
        sample_agg["worsened_samples"].eq(0) & sample_agg["improved_samples"].gt(0)
    ]

    panel_noncurrent = panel.loc[
        ~panel["method"].isin({CURRENT, *TRIVIAL_EQUIVALENTS})
    ].copy()
    best_panel = panel_noncurrent.sort_values(
        ["mean_official_L1_pp", "max_worse_current_L1_pp"],
        kind="mergesort",
    ).iloc[0]
    strict_panel = panel_noncurrent.loc[
        panel_noncurrent["worsened_panel_count"].eq(0)
        & panel_noncurrent["improved_panel_count"].gt(0)
    ]

    rows = [
        {
            "metric": "evaluated_samples",
            "value": str(int(work[["panel", "sample"]].drop_duplicates().shape[0])),
            "evidence": str(SCORES.relative_to(NOTE_DIR)),
            "decision": "fixed_calls_only",
        },
        {
            "metric": "nontrivial_variants",
            "value": str(int(sample_agg.shape[0])),
            "evidence": str(SCORES.relative_to(NOTE_DIR)),
            "decision": "excludes_current_raw_equivalents",
        },
        {
            "metric": "strict_sample_safe_variants",
            "value": str(int(strict_sample.shape[0])),
            "evidence": str(SCORES.relative_to(NOTE_DIR)),
            "decision": "must_be_zero_to_reject_simple_default_replacement",
        },
        {
            "metric": "strict_panel_safe_variants",
            "value": str(int(strict_panel.shape[0])),
            "evidence": str(PANEL_SUMMARY.relative_to(NOTE_DIR)),
            "decision": "must_be_zero_to_reject_simple_default_replacement",
        },
        {
            "metric": "best_sample_mean_variant",
            "value": (
                f"{best_sample['method']};mean_delta={best_sample['mean_sample_delta_L1_pp']:.6f};"
                f"worsened_samples={int(best_sample['worsened_samples'])};"
                f"max_worse={best_sample['max_sample_worse_L1_pp']:.6f}"
            ),
            "evidence": str(SCORES.relative_to(NOTE_DIR)),
            "decision": "tradeoff_not_default",
        },
        {
            "metric": "best_panel_mean_variant",
            "value": (
                f"{best_panel['method']};mean_delta={float(best_panel['mean_delta_current_L1_pp']):.6f};"
                f"worsened_panels={int(best_panel['worsened_panel_count'])};"
                f"max_worse={float(best_panel['max_worse_current_L1_pp']):.6f}"
            ),
            "evidence": str(PANEL_SUMMARY.relative_to(NOTE_DIR)),
            "decision": "tradeoff_not_default",
        },
        {
            "metric": "promotion_decision",
            "value": "reject_fixed_call_abundance_replacement",
            "evidence": f"{SCORES.relative_to(NOTE_DIR)};{PANEL_SUMMARY.relative_to(NOTE_DIR)}",
            "decision": "current_calibrated_abundance_remains_default",
        },
    ]
    write_tsv(OUT, rows)
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
