#!/usr/bin/env python3
"""Summarize cached exact-split diagnostic score files.

The large exact-split profile tables remain in /tmp. This script keeps only
small score summaries in the note folder and marks them diagnostic because the
panels use mixed local scorers, cached /tmp outputs, and mixed reference-release
baselines.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
SOURCE_DIR = Path("/tmp/minco_exactsplit_universal_20260626")

CURRENT_METHOD = "minco_universal_exactsplit"
SYLPH_METHOD = "sylph"


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def load_scores() -> pd.DataFrame:
    files = sorted(SOURCE_DIR.glob("*_exactsplit_score.tsv"))
    if not files:
        raise SystemExit(f"no *_exactsplit_score.tsv files under {SOURCE_DIR}")
    frames = []
    for path in files:
        rows = pd.read_csv(path, sep="\t")
        rows["source_file"] = str(path)
        frames.append(rows)
    scores = pd.concat(frames, ignore_index=True, sort=False)
    scores = scores.drop_duplicates(
        subset=[
            "method",
            "sample_key",
            "dataset",
            "gold_taxa",
            "pred_taxa",
            "TP",
            "FP",
            "FN",
            "precision",
            "recall",
            "F1",
            "FP_plus_FN",
            "L1",
            "Pearson",
        ],
        keep="first",
    )
    scores = numeric(
        scores,
        [
            "gold_taxa",
            "pred_taxa",
            "TP",
            "FP",
            "FN",
            "precision",
            "recall",
            "F1",
            "FP_plus_FN",
            "L1",
            "Pearson",
            "tp_abundance_n",
        ],
    )
    return scores.sort_values(["dataset", "sample_key", "method"], kind="mergesort")


def summarize_methods(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, method), sub in scores.groupby(["dataset", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        pooled_precision = tp / (tp + fp) if tp + fp else 0.0
        pooled_recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = (
            2.0 * pooled_precision * pooled_recall / (pooled_precision + pooled_recall)
            if pooled_precision + pooled_recall
            else 0.0
        )
        rows.append(
            {
                "dataset": dataset,
                "method": method,
                "samples": ",".join(sorted(sub["sample_key"].astype(str).unique())),
                "sample_count": int(sub["sample_key"].nunique()),
                "mean_precision": float(sub["precision"].mean()),
                "mean_recall": float(sub["recall"].mean()),
                "mean_F1": float(sub["F1"].mean()),
                "mean_L1": float(sub["L1"].mean()),
                "mean_Pearson": float(sub["Pearson"].mean()),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": pooled_f1,
            }
        )
    return pd.DataFrame(rows).sort_values(["dataset", "method"], kind="mergesort")


def pairwise_current_vs_sylph(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    current = scores.loc[scores["method"].eq(CURRENT_METHOD)].copy()
    sylph = scores.loc[scores["method"].eq(SYLPH_METHOD)].copy()
    merged = current.merge(
        sylph,
        on=["dataset", "sample_key"],
        suffixes=("_minco", "_sylph"),
        how="inner",
    )
    if merged.empty:
        return merged, pd.DataFrame()
    merged["delta_F1_minco_minus_sylph"] = merged["F1_minco"] - merged["F1_sylph"]
    merged["delta_L1_minco_minus_sylph"] = merged["L1_minco"] - merged["L1_sylph"]
    merged["delta_Pearson_minco_minus_sylph"] = merged["Pearson_minco"] - merged["Pearson_sylph"]
    rows = []
    for dataset, sub in merged.groupby("dataset", sort=True):
        rows.append(
            {
                "dataset": dataset,
                "samples": ",".join(sorted(sub["sample_key"].astype(str).unique())),
                "comparable_samples": int(sub["sample_key"].nunique()),
                "mean_delta_F1_minco_minus_sylph": float(sub["delta_F1_minco_minus_sylph"].mean()),
                "mean_delta_L1_minco_minus_sylph": float(sub["delta_L1_minco_minus_sylph"].mean()),
                "mean_delta_Pearson_minco_minus_sylph": float(sub["delta_Pearson_minco_minus_sylph"].mean()),
                "F1_wins": int((sub["delta_F1_minco_minus_sylph"] > 0.0).sum()),
                "F1_losses": int((sub["delta_F1_minco_minus_sylph"] < 0.0).sum()),
                "L1_wins": int((sub["delta_L1_minco_minus_sylph"] < 0.0).sum()),
                "L1_losses": int((sub["delta_L1_minco_minus_sylph"] > 0.0).sum()),
                "Pearson_wins": int((sub["delta_Pearson_minco_minus_sylph"] > 0.0).sum()),
                "Pearson_losses": int((sub["delta_Pearson_minco_minus_sylph"] < 0.0).sum()),
            }
        )
    return merged, pd.DataFrame(rows).sort_values("dataset", kind="mergesort")


def audit(pair_summary: pd.DataFrame) -> list[dict[str, object]]:
    if pair_summary.empty:
        return [
            {
                "metric": "cached_exactsplit_diagnostic",
                "value": "no_comparable_sylph_samples",
                "evidence": str(SOURCE_DIR),
                "decision": "not_promotable",
            }
        ]
    f1_win_datasets = int((pair_summary["F1_wins"] > pair_summary["F1_losses"]).sum())
    f1_loss_datasets = int((pair_summary["F1_losses"] > pair_summary["F1_wins"]).sum())
    l1_win_datasets = int((pair_summary["L1_wins"] > pair_summary["L1_losses"]).sum())
    l1_loss_datasets = int((pair_summary["L1_losses"] > pair_summary["L1_wins"]).sum())
    return [
        {
            "metric": "cached_exactsplit_sources",
            "value": int(len(list(SOURCE_DIR.glob("*_exactsplit_score.tsv")))),
            "evidence": str(SOURCE_DIR),
            "decision": "diagnostic_nonrelease_cached_outputs",
        },
        {
            "metric": "comparable_datasets",
            "value": int(pair_summary["dataset"].nunique()),
            "evidence": "cached_exactsplit_pairwise_minco_vs_sylph.tsv",
            "decision": "diagnostic_minco_vs_sylph_subset",
        },
        {
            "metric": "current_exactsplit_vs_sylph",
            "value": (
                f"F1_win_datasets={f1_win_datasets};F1_loss_datasets={f1_loss_datasets};"
                f"L1_win_datasets={l1_win_datasets};L1_loss_datasets={l1_loss_datasets}"
            ),
            "evidence": "cached_exactsplit_pairwise_summary.tsv",
            "decision": "mixed_diagnostic_support_not_release_grade",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_release_grade",
            "evidence": "cached_exactsplit_pairwise_summary.tsv",
            "decision": "do_not_promote_or_reject_default_from_cached_local_scorers",
        },
    ]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    scores = load_scores()
    method_summary = summarize_methods(scores)
    pairwise, pair_summary = pairwise_current_vs_sylph(scores)
    audit_rows = audit(pair_summary)

    scores.to_csv(RESULTS / "cached_exactsplit_diagnostic_scores.tsv", sep="\t", index=False)
    method_summary.to_csv(RESULTS / "cached_exactsplit_diagnostic_method_summary.tsv", sep="\t", index=False)
    pairwise.to_csv(RESULTS / "cached_exactsplit_pairwise_minco_vs_sylph.tsv", sep="\t", index=False)
    pair_summary.to_csv(RESULTS / "cached_exactsplit_pairwise_summary.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "cached_exactsplit_diagnostic_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )

    print(pd.DataFrame(audit_rows).to_string(index=False))
    print("\nPAIRWISE SUMMARY")
    print(pair_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
