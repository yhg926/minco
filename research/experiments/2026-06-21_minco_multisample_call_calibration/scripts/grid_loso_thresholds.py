#!/usr/bin/env python3
"""Leave-one-sample-out grid search for interpretable minco direct thresholds."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
HELPER_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(HELPER_DIR))
THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

from analyze_readwise_corrections import load_minco, parse_gold_profile, parse_species_taxmap, score_taxids  # noqa: E402
from calibrate_multisample_calls import filter_scope_taxids, gold_scope, taxid_scope_map  # noqa: E402


def _num(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def score_pred(pred: Iterable[str], gold: Set[str]) -> Dict[str, object]:
    pred_set = {str(t) for t in pred if str(t)}
    out: Dict[str, object] = {"pred_taxa": len(pred_set)}
    out.update(score_taxids(pred_set, gold))
    out["FP_plus_FN"] = out["FP"] + out["FN"]
    return out


def params_grid() -> List[Dict[str, float]]:
    params: List[Dict[str, float]] = []
    for xny in [1, 5, 10, 25, 50, 100]:
        for ani in [0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98]:
            for minaf in [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20]:
                params.append({"xny": float(xny), "ani": float(ani), "minaf": float(minaf)})
    return params


def mask_for(rows: pd.DataFrame, p: Mapping[str, float]) -> pd.Series:
    return (_num(rows, "XnY_ctx") >= p["xny"]) & (_num(rows, "ANI") >= p["ani"]) & (_num(rows, "Real_min_align_fraction") >= p["minaf"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--taxmap", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    taxmap = parse_species_taxmap(args.taxmap)
    scope_by_taxid = taxid_scope_map(taxmap)
    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)

    samples: Dict[str, Dict[str, object]] = {}
    for _, rec in manifest.iterrows():
        key = str(rec["sample_key"])
        scope = str(rec["scope"])
        gold = gold_scope(parse_gold_profile(Path(rec["gold_profile"]), str(rec["gold_sample_id"])), scope)
        rows = load_minco(Path(rec["unique_unfiltered"]), taxmap, 11.0)
        rows = rows.loc[rows["taxid"].astype(bool)].copy()
        samples[key] = {
            "dataset": str(rec["dataset"]),
            "scope": scope,
            "gold": set(gold["taxid"].astype(str)),
            "rows": rows,
        }

    params = params_grid()
    grid_records: List[Dict[str, object]] = []
    loso_records: List[Dict[str, object]] = []

    def evaluate_sample(sample_key: str, p: Mapping[str, float]) -> Dict[str, object]:
        sample = samples[sample_key]
        rows = sample["rows"]
        pred = filter_scope_taxids(rows.loc[mask_for(rows, p), "taxid"].astype(str), str(sample["scope"]), scope_by_taxid)
        out = {
            "sample_key": sample_key,
            "dataset": sample["dataset"],
            "scope": sample["scope"],
            "gold_taxa": len(sample["gold"]),
            **p,
        }
        out.update(score_pred(pred, sample["gold"]))
        return out

    for i, p in enumerate(params):
        for sample_key in sorted(samples):
            row = {"param_id": i, **evaluate_sample(sample_key, p)}
            grid_records.append(row)
    grid = pd.DataFrame(grid_records)
    grid.to_csv(args.outdir / "all_threshold_grid.tsv", sep="\t", index=False)

    for heldout in sorted(samples):
        train = grid.loc[grid["sample_key"] != heldout].copy()
        mean_train = (
            train.groupby("param_id", as_index=False)
            .agg(mean_F1=("F1", "mean"), mean_FP_plus_FN=("FP_plus_FN", "mean"), min_F1=("F1", "min"))
            .sort_values(["mean_F1", "mean_FP_plus_FN", "min_F1"], ascending=[False, True, False])
        )
        best_id = int(mean_train.iloc[0]["param_id"])
        best = params[best_id]
        row = evaluate_sample(heldout, best)
        row.update({"method": "loso_threshold", "param_id": best_id, "train_samples": ",".join(sorted(k for k in samples if k != heldout))})
        loso_records.append(row)

        fixed = evaluate_sample(heldout, {"xny": 10.0, "ani": 0.94, "minaf": 0.05})
        fixed.update({"method": "fixed_default", "param_id": -1, "train_samples": ""})
        loso_records.append(fixed)

    loso = pd.DataFrame(loso_records)
    loso.to_csv(args.outdir / "loso_threshold_summary.tsv", sep="\t", index=False)
    mean = (
        loso.groupby("method", as_index=False)
        .agg(
            samples=("sample_key", "count"),
            mean_precision=("precision", "mean"),
            mean_recall=("recall", "mean"),
            mean_F1=("F1", "mean"),
            mean_FP=("FP", "mean"),
            mean_FN=("FN", "mean"),
        )
        .sort_values("mean_F1", ascending=False)
    )
    mean.to_csv(args.outdir / "loso_threshold_mean.tsv", sep="\t", index=False)

    global_best = (
        grid.groupby("param_id", as_index=False)
        .agg(mean_F1=("F1", "mean"), mean_FP_plus_FN=("FP_plus_FN", "mean"), min_F1=("F1", "min"))
        .sort_values(["mean_F1", "mean_FP_plus_FN", "min_F1"], ascending=[False, True, False])
        .head(20)
    )
    global_best.to_csv(args.outdir / "global_best_thresholds.tsv", sep="\t", index=False)
    print(mean.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
