#!/usr/bin/env python3
"""Summarize cached current exact-split plant holdout scores.

The source score files live under /tmp because their raw exact-split tables are
large. This script copies only compact score summaries into the note folder so
the decision table can consume them reproducibly while keeping large artifacts
out of git.
"""

from __future__ import annotations

import csv
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
SOURCE_DIR = Path("/tmp/minco_exactsplit_universal_20260626")
SAMPLES = ["plant_holdout3", "plant_holdout4", "plant_holdout5"]
METHODS = ["minco_universal_exactsplit", "sylph"]


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


def sample_file(sample: str) -> Path:
    return SOURCE_DIR / f"{sample}_exactsplit_score.tsv"


def mean(rows: list[dict[str, str]], field: str) -> float:
    values = [float(row[field]) for row in rows if row.get(field, "") not in {"", "NA"}]
    if not values:
        return 0.0
    return sum(values) / len(values)


def total(rows: list[dict[str, str]], field: str) -> int:
    return sum(int(float(row[field])) for row in rows if row.get(field, "") not in {"", "NA"})


def unique_method_row(rows: list[dict[str, str]], method: str, path: Path) -> dict[str, str]:
    matches = [row for row in rows if row.get("method") == method]
    if len(matches) == 1:
        return matches[0]
    unique_rows = {
        tuple((key, row.get(key, "")) for key in sorted(row))
        for row in matches
    }
    if len(unique_rows) == 1 and matches:
        return matches[0]
    raise SystemExit(f"expected one {method} row in {path}, found {len(matches)}")


def main() -> int:
    score_rows: list[dict[str, object]] = []
    missing = [str(sample_file(sample)) for sample in SAMPLES if not sample_file(sample).exists()]
    if missing:
        raise SystemExit("missing cached plant exact-split score files: " + ", ".join(missing))

    for sample in SAMPLES:
        path = sample_file(sample)
        rows = read_tsv(path)
        for method in METHODS:
            row = dict(unique_method_row(rows, method, path))
            row["source_file"] = str(path)
            score_rows.append(row)

    score_fields = [
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
        "tp_abundance_n",
        "raw_unique_fallback",
        "source_file",
    ]
    write_tsv(RESULTS / "plant_holdout_exactsplit_current_scores.tsv", score_rows, score_fields)

    summary_rows: list[dict[str, object]] = []
    for method in METHODS:
        method_rows = [row for row in score_rows if row["method"] == method]
        tp = total(method_rows, "TP")
        fp = total(method_rows, "FP")
        fn = total(method_rows, "FN")
        summary_rows.append(
            {
                "method": method,
                "samples": "3,4,5",
                "mean_precision": mean(method_rows, "precision"),
                "mean_recall": mean(method_rows, "recall"),
                "mean_F1": mean(method_rows, "F1"),
                "mean_L1": mean(method_rows, "L1"),
                "mean_Pearson": mean(method_rows, "Pearson"),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "source_files": ",".join(str(sample_file(sample)) for sample in SAMPLES),
                "caveat": (
                    "nonrelease: cached current exact-split MinCO outputs, "
                    "plant bacteria-scope scorer, Sylph r226/default comparison"
                ),
            }
        )

    write_tsv(
        RESULTS / "plant_holdout_exactsplit_current_summary.tsv",
        summary_rows,
        [
            "method",
            "samples",
            "mean_precision",
            "mean_recall",
            "mean_F1",
            "mean_L1",
            "mean_Pearson",
            "pooled_TP",
            "pooled_FP",
            "pooled_FN",
            "source_files",
            "caveat",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
