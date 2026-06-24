#!/usr/bin/env python3
"""Retune saved RF/HGB probabilities with one global threshold."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


OUT = Path("research/experiments/2026-06-23_universal_strategy_search/results")
TRAIN_PRED = Path("/mnt/new3T/minco_cami2_plant_20260621/calibration_marine_toy_plant_20260621/model_noleak_rules_ensemble/loso_predictions.tsv")
STRAIN_PRED = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/model_predictions.tsv")
TRAIN_SUMMARY = Path("/mnt/new3T/minco_cami2_plant_20260621/calibration_marine_toy_plant_20260621/model_noleak_rules_ensemble/summary.tsv")
STRAIN_SUMMARY = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/summary.tsv")
TRAIN_FEATURES = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/train.joined_features.tsv")
STRAIN_FEATURES = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/test.joined_features.tsv")


def load_gold_counts() -> dict[str, int]:
    out: dict[str, int] = {}
    for path in [TRAIN_SUMMARY, STRAIN_SUMMARY]:
        df = pd.read_csv(path, sep="\t")
        for row in df[["sample_key", "gold_taxa"]].drop_duplicates().itertuples(index=False):
            out[str(row.sample_key)] = int(row.gold_taxa)
    return out


def score(df: pd.DataFrame, pred: np.ndarray, gold: dict[str, int]) -> pd.DataFrame:
    rows = []
    sample_values = df["sample_key"].astype(str).to_numpy()
    label = df["label"].to_numpy(dtype=int)
    for sample_key, gold_count in sorted(gold.items()):
        mask = sample_values == sample_key
        if not mask.any():
            continue
        tp = int(np.sum(pred[mask] & (label[mask] == 1)))
        fp = int(np.sum(pred[mask] & (label[mask] == 0)))
        fn = max(int(gold_count) - tp, 0)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "sample_key": sample_key,
                "dataset": str(df.loc[mask, "dataset"].iloc[0]),
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": precision,
                "recall": recall,
                "F1": f1,
            }
        )
    return pd.DataFrame(rows)


def summarize(sample: pd.DataFrame) -> dict[str, float]:
    by_dataset = sample.groupby("dataset")["F1"].mean()
    out = {
        "mean_F1": float(sample["F1"].mean()),
        "mean_precision": float(sample["precision"].mean()),
        "mean_recall": float(sample["recall"].mean()),
        "mean_FP": float(sample["FP"].mean()),
        "mean_FN": float(sample["FN"].mean()),
        "mean_FP_plus_FN": float((sample["FP"] + sample["FN"]).mean()),
        "min_dataset_F1": float(by_dataset.min()),
    }
    for dataset, f1 in by_dataset.items():
        out[f"{dataset}_F1"] = float(f1)
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    train = pd.read_csv(TRAIN_PRED, sep="\t")
    train = train.loc[train["model"] == "rf_hgb_avg", ["sample_key", "taxid", "probability", "label"]].copy()
    strain = pd.read_csv(STRAIN_PRED, sep="\t")
    strain = strain.loc[strain["method"] == "train9_rf_hgb_avg", ["sample_key", "taxid", "probability", "label"]].copy()
    probs = pd.concat([train, strain], ignore_index=True, sort=False)
    probs["sample_key"] = probs["sample_key"].astype(str)
    probs["taxid"] = probs["taxid"].astype(str)

    features = pd.concat([pd.read_csv(TRAIN_FEATURES, sep="\t"), pd.read_csv(STRAIN_FEATURES, sep="\t")], ignore_index=True, sort=False)
    features["sample_key"] = features["sample_key"].astype(str)
    features["taxid"] = features["taxid"].astype(str)
    features["direct"] = (
        pd.to_numeric(features.get("u_direct_call", 0), errors="coerce").fillna(0).to_numpy() > 0
    ) | (
        pd.to_numeric(features.get("s_direct_call", 0), errors="coerce").fillna(0).to_numpy() > 0
    )
    df = probs.merge(features[["sample_key", "taxid", "dataset", "scope", "direct"]], on=["sample_key", "taxid"], how="left")
    df["direct"] = df["direct"].fillna(False).astype(bool)

    rows = []
    sample_rows = []
    gold = load_gold_counts()
    thresholds = np.unique(np.r_[np.linspace(0.01, 0.1, 10), np.linspace(0.11, 0.95, 85), 0.975, 0.99])
    for mode in ["prob", "direct_or_prob"]:
        for threshold in thresholds:
            pred = df["probability"].to_numpy(dtype=float) >= threshold
            if mode == "direct_or_prob":
                pred = pred | df["direct"].to_numpy(dtype=bool)
            sample = score(df, pred, gold)
            rec = summarize(sample)
            rec["mode"] = mode
            rec["threshold"] = float(threshold)
            rows.append(rec)
            if mode == "prob" and abs(float(threshold) - 0.35) < 1e-12:
                sample["mode"] = mode
                sample["threshold"] = float(threshold)
                sample_rows.append(sample)

    summary = pd.DataFrame(rows).sort_values(["mean_F1", "mean_FP_plus_FN", "min_dataset_F1"], ascending=[False, True, False])
    summary.to_csv(OUT / "previous_rf_hgb_threshold_sweep.tsv", sep="\t", index=False)
    if sample_rows:
        pd.concat(sample_rows, ignore_index=True).to_csv(OUT / "previous_rf_hgb_t035_sample_metrics.tsv", sep="\t", index=False)
    print(summary.head(30).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
