#!/usr/bin/env python3
"""Train minco call models on one panel and test an external holdout panel."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
CALIB_DIR = REPO_ROOT / "research/experiments/2026-06-21_minco_multisample_call_calibration/scripts"
sys.path.insert(0, str(CALIB_DIR))

from calibrate_multisample_calls import (  # noqa: E402
    MODEL_COLS,
    direct_mask,
    filter_scope_taxids,
    gold_scope,
    joined_features,
    load_numeric_minco,
    model_suite,
    predict_probability,
    relaxed_mask,
    score_prediction,
    score_sylph,
    taxid_scope_map,
    tune_threshold,
)

HELPER_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(HELPER_DIR))

from analyze_readwise_corrections import parse_gold_profile, parse_species_taxmap  # noqa: E402


def _metric_row(
    *,
    sample_key: str,
    dataset: str,
    scope: str,
    method: str,
    pred_taxids: Iterable[str],
    gold_taxids: Set[str],
    split: str,
    threshold: float | None = None,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "split": split,
        "sample_key": sample_key,
        "dataset": dataset,
        "scope": scope,
        "method": method,
        "gold_taxa": len(gold_taxids),
    }
    if threshold is not None:
        row["threshold"] = threshold
    row.update(score_prediction(pred_taxids, gold_taxids))
    row["FN"] = len(gold_taxids) - int(row["TP"])
    precision = row["TP"] / (row["TP"] + row["FP"]) if row["TP"] + row["FP"] else 0.0
    recall = row["TP"] / (row["TP"] + row["FN"]) if row["TP"] + row["FN"] else 0.0
    row["precision"] = precision
    row["recall"] = recall
    row["F1"] = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    row["FP_plus_FN"] = row["FP"] + row["FN"]
    return row


def _read_manifest(path: Path) -> pd.DataFrame:
    required = {
        "sample_key",
        "dataset",
        "scope",
        "gold_profile",
        "gold_sample_id",
        "unique_unfiltered",
        "split_unfiltered",
        "sylph_profile",
    }
    manifest = pd.read_csv(path, sep="\t", keep_default_na=False)
    missing = sorted(required.difference(manifest.columns))
    if missing:
        raise SystemExit(f"{path} is missing required columns: {', '.join(missing)}")
    return manifest


def load_panel(
    manifest: pd.DataFrame,
    *,
    split: str,
    taxmap: Mapping[str, Mapping[str, str]],
    scope_by_taxid: Mapping[str, str],
    feature_dir: Path,
) -> tuple[pd.DataFrame, Dict[str, Set[str]], Dict[str, int], List[Dict[str, object]]]:
    frames: List[pd.DataFrame] = []
    gold_by_sample: Dict[str, Set[str]] = {}
    gold_counts: Dict[str, int] = {}
    baseline_rows: List[Dict[str, object]] = []

    for _, rec in manifest.iterrows():
        sample_key = str(rec["sample_key"])
        dataset = str(rec["dataset"])
        scope = str(rec["scope"])
        gold = gold_scope(parse_gold_profile(Path(rec["gold_profile"]), str(rec["gold_sample_id"])), scope)
        gold_taxids = set(gold["taxid"].astype(str))
        gold_by_sample[sample_key] = gold_taxids
        gold_counts[sample_key] = len(gold_taxids)

        unique_rows = load_numeric_minco(Path(rec["unique_unfiltered"]), taxmap)
        split_rows = load_numeric_minco(Path(rec["split_unfiltered"]), taxmap)

        unique_direct = filter_scope_taxids(unique_rows.loc[direct_mask(unique_rows), "taxid"].astype(str), scope, scope_by_taxid)
        split_direct = filter_scope_taxids(split_rows.loc[direct_mask(split_rows), "taxid"].astype(str), scope, scope_by_taxid)
        unique_relaxed = filter_scope_taxids(unique_rows.loc[relaxed_mask(unique_rows), "taxid"].astype(str), scope, scope_by_taxid)

        for method, pred_taxids in [
            ("minco_unique_direct", unique_direct),
            ("minco_split_direct", split_direct),
            ("minco_unique_relaxed", unique_relaxed),
        ]:
            baseline_rows.append(
                _metric_row(
                    sample_key=sample_key,
                    dataset=dataset,
                    scope=scope,
                    method=method,
                    pred_taxids=pred_taxids,
                    gold_taxids=gold_taxids,
                    split=split,
                )
            )

        sylph = score_sylph(Path(rec["sylph_profile"]), taxmap, scope, scope_by_taxid, gold_taxids)
        sylph["FN"] = len(gold_taxids) - int(sylph["TP"])
        precision = sylph["TP"] / (sylph["TP"] + sylph["FP"]) if sylph["TP"] + sylph["FP"] else 0.0
        recall = sylph["TP"] / (sylph["TP"] + sylph["FN"]) if sylph["TP"] + sylph["FN"] else 0.0
        sylph["precision"] = precision
        sylph["recall"] = recall
        sylph["F1"] = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        sylph["FP_plus_FN"] = sylph["FP"] + sylph["FN"]
        baseline_rows.append(
            {
                "split": split,
                "sample_key": sample_key,
                "dataset": dataset,
                "scope": scope,
                "method": "sylph_default",
                "gold_taxa": len(gold_taxids),
                **sylph,
            }
        )

        features = joined_features(unique_rows, split_rows)
        if scope != "all":
            features = features.loc[features["taxid"].astype(str).map(lambda t: scope_by_taxid.get(t, "other") == scope)].copy()
        features["split"] = split
        features["sample_key"] = sample_key
        features["dataset"] = dataset
        features["scope"] = scope
        features["label"] = features["taxid"].astype(str).isin(gold_taxids).astype(int)
        frames.append(features)
        features.to_csv(feature_dir / f"{split}_{sample_key}.joined_features.tsv", sep="\t", index=False)

    all_features = pd.concat(frames, ignore_index=True, sort=False).fillna(0.0)
    return all_features, gold_by_sample, gold_counts, baseline_rows


def evaluate_model(
    *,
    method: str,
    threshold: float,
    prob: np.ndarray,
    test_features: pd.DataFrame,
    gold_by_sample: Mapping[str, Set[str]],
) -> tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    metric_rows: List[Dict[str, object]] = []
    prediction_rows: List[Dict[str, object]] = []
    taxids = test_features["taxid"].astype(str).to_numpy()
    samples = test_features["sample_key"].astype(str).to_numpy()

    for sample_key in sorted(set(samples)):
        mask = samples == sample_key
        sample_taxids = taxids[mask]
        sample_prob = prob[mask]
        pred_taxids = set(sample_taxids[sample_prob >= threshold])
        first = test_features.loc[mask].iloc[0]
        metric_rows.append(
            _metric_row(
                sample_key=sample_key,
                dataset=str(first["dataset"]),
                scope=str(first["scope"]),
                method=method,
                pred_taxids=pred_taxids,
                gold_taxids=gold_by_sample[sample_key],
                split="test",
                threshold=threshold,
            )
        )
        labels = test_features.loc[mask, "label"].astype(int).to_numpy()
        for taxid, p, label in zip(sample_taxids, sample_prob, labels):
            prediction_rows.append(
                {
                    "sample_key": sample_key,
                    "method": method,
                    "threshold": threshold,
                    "taxid": taxid,
                    "probability": float(p),
                    "predicted": int(p >= threshold),
                    "label": int(label),
                }
            )
    return metric_rows, prediction_rows


def summarize(summary: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for (split, method), group in summary.groupby(["split", "method"], dropna=False):
        rows.append(
            {
                "split": split,
                "method": method,
                "samples": int(len(group)),
                "mean_precision": float(pd.to_numeric(group["precision"], errors="coerce").mean()),
                "mean_recall": float(pd.to_numeric(group["recall"], errors="coerce").mean()),
                "mean_F1": float(pd.to_numeric(group["F1"], errors="coerce").mean()),
                "mean_FP": float(pd.to_numeric(group["FP"], errors="coerce").mean()),
                "mean_FN": float(pd.to_numeric(group["FN"], errors="coerce").mean()),
                "mean_FP_plus_FN": float(pd.to_numeric(group["FP_plus_FN"], errors="coerce").mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["split", "mean_F1"], ascending=[True, False])


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-manifest", type=Path, required=True)
    ap.add_argument("--test-manifest", type=Path, required=True)
    ap.add_argument("--taxmap", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    feature_dir = args.outdir / "features"
    feature_dir.mkdir(exist_ok=True)

    taxmap = parse_species_taxmap(args.taxmap)
    scope_by_taxid = taxid_scope_map(taxmap)
    train_manifest = _read_manifest(args.train_manifest)
    test_manifest = _read_manifest(args.test_manifest)

    train_features, _, train_gold_counts, train_baselines = load_panel(
        train_manifest,
        split="train",
        taxmap=taxmap,
        scope_by_taxid=scope_by_taxid,
        feature_dir=feature_dir,
    )
    test_features, test_gold_by_sample, _, test_baselines = load_panel(
        test_manifest,
        split="test",
        taxmap=taxmap,
        scope_by_taxid=scope_by_taxid,
        feature_dir=feature_dir,
    )

    train_features.to_csv(args.outdir / "train.joined_features.tsv", sep="\t", index=False)
    test_features.to_csv(args.outdir / "test.joined_features.tsv", sep="\t", index=False)

    x_train = train_features[MODEL_COLS].to_numpy(dtype=float)
    y_train = train_features["label"].to_numpy(dtype=int)
    sample_train = train_features["sample_key"].astype(str).to_numpy()
    x_test = test_features[MODEL_COLS].to_numpy(dtype=float)

    summary_rows: List[Dict[str, object]] = []
    summary_rows.extend(train_baselines)
    summary_rows.extend(test_baselines)
    prediction_rows: List[Dict[str, object]] = []
    train_prob_by_model: Dict[str, np.ndarray] = {}
    test_prob_by_model: Dict[str, np.ndarray] = {}

    for model_name, model in model_suite().items():
        model.fit(x_train, y_train)
        train_prob = predict_probability(model, x_train)
        test_prob = predict_probability(model, x_test)
        train_prob_by_model[model_name] = train_prob
        test_prob_by_model[model_name] = test_prob
        threshold = tune_threshold(train_prob, y_train, sample_train, train_gold_counts)
        rows, preds = evaluate_model(
            method=f"train9_{model_name}",
            threshold=threshold,
            prob=test_prob,
            test_features=test_features,
            gold_by_sample=test_gold_by_sample,
        )
        summary_rows.extend(rows)
        prediction_rows.extend(preds)

    if "rf" in train_prob_by_model and "hgb" in train_prob_by_model:
        threshold = tune_threshold(
            0.5 * (train_prob_by_model["rf"] + train_prob_by_model["hgb"]),
            y_train,
            sample_train,
            train_gold_counts,
        )
        rows, preds = evaluate_model(
            method="train9_rf_hgb_avg",
            threshold=threshold,
            prob=0.5 * (test_prob_by_model["rf"] + test_prob_by_model["hgb"]),
            test_features=test_features,
            gold_by_sample=test_gold_by_sample,
        )
        summary_rows.extend(rows)
        prediction_rows.extend(preds)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    summarize(summary).to_csv(args.outdir / "mean_summary.tsv", sep="\t", index=False)
    pd.DataFrame(prediction_rows).to_csv(args.outdir / "model_predictions.tsv", sep="\t", index=False)

    rf_model = model_suite().get("rf")
    if rf_model is not None:
        rf_model.fit(x_train, y_train)
        final = rf_model.steps[-1][1] if hasattr(rf_model, "steps") else rf_model
        if hasattr(final, "feature_importances_"):
            pd.DataFrame({"feature": MODEL_COLS, "importance": final.feature_importances_}).sort_values(
                "importance", ascending=False
            ).to_csv(args.outdir / "rf_feature_importance.tsv", sep="\t", index=False)

    print(summarize(summary).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
