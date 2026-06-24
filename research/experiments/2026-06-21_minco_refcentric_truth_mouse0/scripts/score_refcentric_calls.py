#!/usr/bin/env python3
"""Score minco outputs against closest-reference truth accessions."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np
import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def accession(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def score(pred: Iterable[str], truth: set[str]) -> Dict[str, object]:
    pred_set = {str(x) for x in pred if str(x)}
    tp = len(pred_set & truth)
    fp = len(pred_set - truth)
    fn = len(truth - pred_set)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "truth_refs": len(truth),
        "pred_refs": len(pred_set),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "FP_plus_FN": fp + fn,
    }


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def load_minco(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    if "Ref" not in df.columns:
        raise SystemExit(f"{path} does not look like minco detail output: missing Ref")
    df["ref_accession"] = df["Ref"].map(accession)
    missing = df["ref_accession"].eq("")
    if missing.any() and "Ref_annotation" in df.columns:
        df.loc[missing, "ref_accession"] = df.loc[missing, "Ref_annotation"].map(accession)
    for col in [
        "ANI",
        "XnY_ctx",
        "Real_min_align_fraction",
        "Ref_breadth",
        "Ref_mean_depth",
        "Fake_ctx_fraction",
        "Fake_ctx_prob_mean",
    ]:
        df[col] = numeric(df, col)
    return df.loc[df["ref_accession"].astype(bool)].copy()


def direct_mask(df: pd.DataFrame) -> pd.Series:
    return (
        (df["XnY_ctx"] >= 10.0)
        & (df["ANI"] >= 0.94)
        & (df["Real_min_align_fraction"] >= 0.05)
    )


def relaxed_mask(df: pd.DataFrame) -> pd.Series:
    return (
        (df["XnY_ctx"] >= 1.0)
        & (df["ANI"] >= 0.90)
        & (df["Ref_breadth"] >= 0.005)
    )


def best_sweep(df: pd.DataFrame, truth: set[str]) -> Dict[str, object]:
    best: Dict[str, object] | None = None
    best_key: tuple[float, int, int, int] | None = None
    best_params = ""
    for ani in np.arange(0.90, 0.991, 0.005):
        for xny in [1, 3, 5, 10, 20, 50, 100, 150, 200, 300, 500, 1000]:
            for af in [0.0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]:
                mask = (
                    (df["XnY_ctx"] >= xny)
                    & (df["ANI"] >= ani)
                    & (df["Real_min_align_fraction"] >= af)
                )
                row = score(df.loc[mask, "ref_accession"], truth)
                key = (float(row["F1"]), -int(row["FP"]), int(row["TP"]), -int(row["FN"]))
                if best_key is None or key > best_key:
                    best_key = key
                    best = row
                    best_params = f"ani>={ani:.3f};xny>={xny};af>={af:.2f}"
    assert best is not None
    best["params"] = best_params
    return best


def parse_prediction_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("prediction must be LABEL=PATH")
    label, path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("empty prediction label")
    return label, Path(path)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-refs", type=Path, required=True)
    ap.add_argument("--prediction", action="append", type=parse_prediction_arg, required=True)
    ap.add_argument("--out-summary", type=Path, required=True)
    args = ap.parse_args(argv)

    truth_df = pd.read_csv(args.truth_refs, sep="\t", keep_default_na=False)
    if "truth_ref_accession" not in truth_df.columns:
        raise SystemExit(f"{args.truth_refs} missing truth_ref_accession")
    truth = set(truth_df["truth_ref_accession"].astype(str))

    rows = []
    for label, path in args.prediction:
        df = load_minco(path)
        for mode, mask in [
            ("direct_ani94_xny10_af05", direct_mask(df)),
            ("relaxed_ani90_xny1_breadth005", relaxed_mask(df)),
        ]:
            row = {"prediction": label, "mode": mode, "params": ""}
            row.update(score(df.loc[mask, "ref_accession"], truth))
            rows.append(row)
        row = {"prediction": label, "mode": "best_threshold_sweep"}
        row.update(best_sweep(df, truth))
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(
        ["mode", "F1", "precision", "recall"], ascending=[True, False, False, False]
    )
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_summary, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
