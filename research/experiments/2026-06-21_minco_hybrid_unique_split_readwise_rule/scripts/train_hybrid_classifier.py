#!/usr/bin/env python3
"""Leave-one-dataset-out classifier test for joined unique/split minco features."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Set

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


FEATURES = [
    "u_ANI",
    "u_Real_min_align_fraction",
    "u_XnY_ctx",
    "u_Ref_breadth",
    "u_Ref_mean_depth",
    "u_Ref_hit_mean_depth",
    "u_Ref_depth_cv",
    "u_Relative_abundance_depth",
    "u_Ref_zip_af",
    "u_Ref_zip_aaf_ani",
    "s_ANI",
    "s_Real_min_align_fraction",
    "s_XnY_ctx",
    "s_Ref_breadth",
    "s_Ref_mean_depth",
    "s_Ref_hit_mean_depth",
    "s_Ref_depth_cv",
    "s_Relative_abundance_depth",
    "s_Ref_zip_af",
    "s_Ref_zip_aaf_ani",
    "split_unique_xny_ratio",
    "split_unique_breadth_ratio",
    "split_unique_ani_delta",
]


def load_gold_taxids(path: Path) -> Set[str]:
    df = pd.read_csv(path, sep="\t")
    return set(df["taxid"].astype(str))


def score(pred: Iterable[str], truth: Set[str]) -> Dict[str, float]:
    pred_set = set(str(x) for x in pred if str(x))
    tp = len(pred_set & truth)
    fp = len(pred_set - truth)
    fn = len(truth - pred_set)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "pred_taxa": len(pred_set),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "FP_plus_FN": fp + fn,
    }


def best_threshold(prob: np.ndarray, taxids: pd.Series, truth: Set[str]) -> float:
    best_t = 0.5
    best_key = (-1.0, -1e18)
    for t in np.linspace(0.05, 0.95, 91):
        pred = taxids[prob >= t]
        s = score(pred, truth)
        key = (s["F1"], -s["FP_plus_FN"])
        if key > best_key:
            best_key = key
            best_t = float(t)
    return best_t


def prep(path: Path, truth: Set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    for col in FEATURES:
        if col not in df.columns:
            df[col] = 0.0
        df[col] = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    df["label"] = df["taxid"].astype(str).isin(truth).astype(int)
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--toy-features", type=Path, required=True)
    ap.add_argument("--toy-gold-taxids", type=Path, required=True)
    ap.add_argument("--marine-features", type=Path, required=True)
    ap.add_argument("--marine-gold-taxids", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    toy_truth = load_gold_taxids(args.toy_gold_taxids)
    marine_truth = load_gold_taxids(args.marine_gold_taxids)
    data = {
        "toy_gtdb_bacteria": prep(args.toy_features, toy_truth),
        "marine_gtdb_species": prep(args.marine_features, marine_truth),
    }
    truth = {"toy_gtdb_bacteria": toy_truth, "marine_gtdb_species": marine_truth}

    models = {
        "logreg": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=1000, class_weight="balanced", random_state=1),
        ),
        "rf": RandomForestClassifier(
            n_estimators=300,
            min_samples_leaf=3,
            class_weight="balanced_subsample",
            random_state=1,
            n_jobs=8,
        ),
        "hgb": HistGradientBoostingClassifier(
            max_iter=200,
            learning_rate=0.05,
            l2_regularization=0.1,
            random_state=1,
        ),
    }

    records: List[Dict[str, object]] = []
    for train_name, test_name in [
        ("toy_gtdb_bacteria", "marine_gtdb_species"),
        ("marine_gtdb_species", "toy_gtdb_bacteria"),
    ]:
        train = data[train_name]
        test = data[test_name]
        X_train = train[FEATURES].to_numpy(dtype=float)
        y_train = train["label"].to_numpy(dtype=int)
        X_test = test[FEATURES].to_numpy(dtype=float)
        for model_name, model in models.items():
            model.fit(X_train, y_train)
            train_prob = model.predict_proba(X_train)[:, 1]
            test_prob = model.predict_proba(X_test)[:, 1]
            threshold = best_threshold(train_prob, train["taxid"].astype(str), truth[train_name])
            train_score = score(train.loc[train_prob >= threshold, "taxid"], truth[train_name])
            test_score = score(test.loc[test_prob >= threshold, "taxid"], truth[test_name])
            records.append(
                {
                    "train_dataset": train_name,
                    "test_dataset": test_name,
                    "model": model_name,
                    "threshold": threshold,
                    "split": "train",
                    **train_score,
                }
            )
            records.append(
                {
                    "train_dataset": train_name,
                    "test_dataset": test_name,
                    "model": model_name,
                    "threshold": threshold,
                    "split": "test",
                    **test_score,
                }
            )

    out = pd.DataFrame(records)
    out.to_csv(args.outdir / "classifier_leave_one_dataset_out.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
