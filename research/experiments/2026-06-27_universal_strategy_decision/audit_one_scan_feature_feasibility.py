#!/usr/bin/env python3
"""Audit whether split-pass outputs can replace the unique-pass table.

The current calibrated default runs separate best-diff-unique and
best-diff-split MinCO passes. A tempting speed shortcut is to derive the unique
channel from columns already present in the split table. This script tests that
shortcut on available real cached tables and records whether it is safe.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
OUT = EXP / "results/one_scan_feature_feasibility.tsv"

PAIRS = [
    (
        "hmp_airskin_sample11_block",
        Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample11_unique_zip_unfiltered.tsv"),
        Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample11_split_zip_unfiltered.tsv"),
    ),
    (
        "hmp_airskin_sample6_block",
        Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample6_unique_zip_unfiltered.tsv"),
        Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample6_split_zip_unfiltered.tsv"),
    ),
    (
        "toymouse_sample6_current",
        Path("/tmp/cami2_toymouse_current_default_20260626/run/sample6/minco.best_diff_unique.unfiltered.tsv"),
        Path("/tmp/cami2_toymouse_current_default_20260626/run/sample6/minco.best_diff_split.unfiltered.tsv"),
    ),
]

COMPARE_COLS = [
    ("XnY_ctx", "Raw_XnY_ctx"),
    ("XnY_ctx", "XnY_ctx"),
    ("ANI", "ANI"),
    ("Real_min_align_fraction", "Real_min_align_fraction"),
    ("Ref_breadth", "Ref_breadth"),
    ("Ref_mean_depth", "Ref_mean_depth"),
    ("Ref_zip_af", "Ref_zip_af"),
]


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def summarize_pair(label: str, unique_path: Path, split_path: Path) -> list[dict[str, object]]:
    if not unique_path.exists() or not split_path.exists():
        return [
            {
                "sample": label,
                "unique_path": str(unique_path),
                "split_path": str(split_path),
                "metric_pair": "all",
                "available": "no",
                "unique_rows": "",
                "split_rows": "",
                "shared_refs": "",
                "exact_matches": "",
                "mismatches": "",
                "mismatch_fraction": "",
                "max_abs_diff": "",
                "pearson": "",
                "decision": "input table missing; not scored",
            }
        ]

    unique = load(unique_path)
    split = load(split_path)
    if "Ref" not in unique.columns or "Ref" not in split.columns:
        raise SystemExit(f"missing Ref column in {label}")
    merged = unique.merge(split, on="Ref", how="inner", suffixes=("_u", "_s"))
    rows: list[dict[str, object]] = []
    for u_col, s_col in COMPARE_COLS:
        left_name = f"{u_col}_u"
        right_name = f"{s_col}_s" if s_col != "Raw_XnY_ctx" else "Raw_XnY_ctx_s"
        if left_name not in merged.columns or right_name not in merged.columns:
            continue
        left = numeric(merged, left_name)
        right = numeric(merged, right_name)
        diff = (left - right).abs()
        mismatches = int((diff > 1e-9).sum())
        shared = int(len(merged))
        pearson = float(left.corr(right)) if shared > 1 else 0.0
        max_abs = float(diff.max()) if shared else 0.0
        safe = mismatches == 0
        rows.append(
            {
                "sample": label,
                "unique_path": str(unique_path),
                "split_path": str(split_path),
                "metric_pair": f"unique.{u_col} vs split.{s_col}",
                "available": "yes",
                "unique_rows": int(len(unique)),
                "split_rows": int(len(split)),
                "shared_refs": shared,
                "exact_matches": int(shared - mismatches),
                "mismatches": mismatches,
                "mismatch_fraction": f"{(mismatches / shared) if shared else 0.0:.8f}",
                "max_abs_diff": f"{max_abs:.8g}",
                "pearson": f"{pearson:.8f}" if np.isfinite(pearson) else "NA",
                "decision": (
                    "safe_remap" if safe else "not_safe; requires real dual accumulator"
                ),
            }
        )
    return rows


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for label, unique_path, split_path in PAIRS:
        rows.extend(summarize_pair(label, unique_path, split_path))
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
