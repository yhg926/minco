#!/usr/bin/env python3
"""Cross-domain MinCO model scoring on cached joined-feature tables.

This script uses the same row-level features and threshold objective as the
2026-06-21 calibration work, but evaluates two stricter questions:

1. leave-one-sample-out over all currently cached domains;
2. leave-one-dataset-out, where an entire domain is unseen during training.

Dataset identity is never used as a model feature.
"""

from __future__ import annotations

import sys
import os
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP = Path(__file__).resolve().parent
OUT = EXP / "results"
PANEL = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2")
TRAIN_FEATURES = PANEL / "train.joined_features.tsv"
TEST_FEATURES = PANEL / "test.joined_features.tsv"
TRAIN9_SUMMARY = Path(
    "/mnt/new3T/minco_cami2_plant_20260621/"
    "calibration_marine_toy_plant_20260621/model_noleak_rules_ensemble/summary.tsv"
)
STRAIN_SUMMARY = PANEL / "summary.tsv"

CALIBRATION_SCRIPTS = ROOT / "research/experiments/2026-06-21_minco_multisample_call_calibration/scripts"
sys.path.insert(0, str(CALIBRATION_SCRIPTS))

from calibrate_multisample_calls import MODEL_COLS, model_suite, predict_probability  # noqa: E402


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def load_features() -> pd.DataFrame:
    frames = []
    for path in [TRAIN_FEATURES, TEST_FEATURES]:
        df = pd.read_csv(path, sep="\t")
        for col in df.columns:
            if col not in {"taxid", "sample_key", "dataset", "scope", "split"}:
                df[col] = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True, sort=False).fillna(0.0)
    out["taxid"] = out["taxid"].astype(str)
    out["sample_key"] = out["sample_key"].astype(str)
    out["dataset"] = out["dataset"].astype(str)
    for col in MODEL_COLS:
        if col not in out.columns:
            out[col] = 0.0
    return out


def load_gold_counts() -> dict[str, int]:
    rows = []
    for path in [TRAIN9_SUMMARY, STRAIN_SUMMARY]:
        df = pd.read_csv(path, sep="\t")
        rows.append(df[["sample_key", "gold_taxa"]])
    merged = pd.concat(rows, ignore_index=True)
    merged["sample_key"] = merged["sample_key"].astype(str)
    merged["gold_taxa"] = pd.to_numeric(merged["gold_taxa"], errors="coerce").fillna(0).astype(int)
    return dict(merged.drop_duplicates("sample_key")[["sample_key", "gold_taxa"]].itertuples(index=False))


def load_baseline_summary() -> pd.DataFrame:
    train9 = pd.read_csv(TRAIN9_SUMMARY, sep="\t")
    strain = pd.read_csv(STRAIN_SUMMARY, sep="\t")
    rows = []
    train_map = {
        "sylph_default": "sylph_default",
        "minco_unique_direct": "minco_unique_direct",
        "minco_split_direct": "minco_split_direct",
        "loso_rf_hgb_avg": "previous_rf_hgb_avg",
        "loso_rf": "previous_rf",
        "loso_hgb": "previous_hgb",
    }
    strain_map = {
        "sylph_default": "sylph_default",
        "minco_unique_direct": "minco_unique_direct",
        "minco_split_direct": "minco_split_direct",
        "train9_rf_hgb_avg": "previous_rf_hgb_avg",
        "train9_rf": "previous_rf",
        "train9_hgb": "previous_hgb",
    }
    for raw, canonical in train_map.items():
        sub = train9.loc[train9["method"] == raw].copy()
        sub["method"] = canonical
        rows.append(sub)
    for raw, canonical in strain_map.items():
        sub = strain.loc[(strain["split"] == "test") & (strain["method"] == raw)].copy()
        sub["method"] = canonical
        rows.append(sub)
    sample = pd.concat(rows, ignore_index=True, sort=False)
    sample["source"] = "baseline"
    return sample


def direct_mask(df: pd.DataFrame) -> np.ndarray:
    return (
        (numeric(df, "u_direct_call").to_numpy() > 0)
        | (numeric(df, "s_direct_call").to_numpy() > 0)
    )


def feature_matrix(df: pd.DataFrame, mode: str) -> tuple[np.ndarray, list[str]]:
    if mode == "base":
        cols = list(MODEL_COLS)
        return df[cols].to_numpy(dtype=float), cols

    if mode != "enriched":
        raise ValueError(f"unsupported feature mode: {mode}")

    out = df[list(MODEL_COLS)].copy()
    pairs = [
        ("ANI", "u_ANI_max", "s_ANI_max"),
        ("XnY", "u_XnY_ctx_max", "s_XnY_ctx_max"),
        ("minAF", "u_Real_min_align_fraction_max", "s_Real_min_align_fraction_max"),
        ("breadth", "u_Ref_breadth_max", "s_Ref_breadth_max"),
        ("mean_depth", "u_Ref_mean_depth_max", "s_Ref_mean_depth_max"),
        ("hit_depth", "u_Ref_hit_mean_depth_max", "s_Ref_hit_mean_depth_max"),
        ("norm_depth", "u_Normalized_abundance_depth_max", "s_Normalized_abundance_depth_max"),
        ("zip_af", "u_Ref_zip_af_max", "s_Ref_zip_af_max"),
        ("zip_ani", "u_Ref_zip_aaf_ani_max", "s_Ref_zip_aaf_ani_max"),
    ]
    for name, u_col, s_col in pairs:
        u = numeric(df, u_col)
        s = numeric(df, s_col)
        max_v = np.maximum(u.to_numpy(), s.to_numpy())
        min_v = np.minimum(u.to_numpy(), s.to_numpy())
        out[f"max_{name}"] = max_v
        out[f"min_nonzero_{name}"] = np.where(min_v > 0, min_v, max_v)
        out[f"delta_s_minus_u_{name}"] = s.to_numpy() - u.to_numpy()
        if name in {"XnY", "mean_depth", "hit_depth", "norm_depth"}:
            out[f"log1p_max_{name}"] = np.log1p(np.maximum(max_v, 0.0))

    u_cv = numeric(df, "u_Ref_depth_cv_min").to_numpy()
    s_cv = numeric(df, "s_Ref_depth_cv_min").to_numpy()
    u_cv = np.where(u_cv > 0, u_cv, np.nan)
    s_cv = np.where(s_cv > 0, s_cv, np.nan)
    min_cv = np.nanmin(np.vstack([u_cv, s_cv]), axis=0)
    min_cv = np.where(np.isfinite(min_cv), min_cv, 0.0)
    out["min_nonzero_cv"] = min_cv
    out["log1p_min_nonzero_cv"] = np.log1p(np.maximum(min_cv, 0.0))
    out["either_direct"] = direct_mask(df).astype(float)
    out["both_direct"] = (
        (numeric(df, "u_direct_call").to_numpy() > 0) & (numeric(df, "s_direct_call").to_numpy() > 0)
    ).astype(float)
    out["either_relaxed_not_direct"] = (
        (
            (numeric(df, "u_relaxed_call").to_numpy() > 0)
            | (numeric(df, "s_relaxed_call").to_numpy() > 0)
        )
        & (~direct_mask(df))
    ).astype(float)
    out = out.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return out.to_numpy(dtype=float), list(out.columns)


def f1_from_counts(tp: int, fp: int, gold: int) -> tuple[int, float, float, float, int]:
    fn = max(gold - tp, 0)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return fn, precision, recall, f1, fp + fn


def tune_threshold(
    prob: np.ndarray,
    y: np.ndarray,
    sample: np.ndarray,
    gold_counts: Mapping[str, int],
    *,
    base_mask: np.ndarray | None = None,
) -> float:
    if base_mask is None:
        base_mask = np.zeros(len(prob), dtype=bool)
    best_t = 0.5
    best_key = (-1.0, -1e18, -1e18)
    thresholds = np.unique(np.r_[np.linspace(0.01, 0.10, 10), np.linspace(0.11, 0.95, 85), 0.975, 0.99])
    for threshold in thresholds:
        f1_values = []
        fpfn_values = []
        pred = base_mask | (prob >= threshold)
        for sample_key in sorted(set(sample)):
            mask = sample == sample_key
            tp = int(np.sum(pred[mask] & (y[mask] == 1)))
            fp = int(np.sum(pred[mask] & (y[mask] == 0)))
            fn, _, _, f1, fpfn = f1_from_counts(tp, fp, int(gold_counts[sample_key]))
            _ = fn
            f1_values.append(f1)
            fpfn_values.append(fpfn)
        key = (float(np.mean(f1_values)), -float(np.mean(fpfn_values)), -abs(float(threshold) - 0.5))
        if key > best_key:
            best_key = key
            best_t = float(threshold)
    return best_t


def score_sample(
    *,
    method: str,
    sample_key: str,
    dataset: str,
    scope: str,
    pred: np.ndarray,
    label: np.ndarray,
    gold_counts: Mapping[str, int],
    validation: str,
    threshold: float,
    train_group: str,
) -> dict[str, object]:
    tp = int(np.sum(pred & (label == 1)))
    fp = int(np.sum(pred & (label == 0)))
    fn, precision, recall, f1, fpfn = f1_from_counts(tp, fp, int(gold_counts[sample_key]))
    return {
        "method": method,
        "validation": validation,
        "sample_key": sample_key,
        "dataset": dataset,
        "scope": scope,
        "gold_taxa": int(gold_counts[sample_key]),
        "pred_taxa": int(np.sum(pred)),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "FP_plus_FN": fpfn,
        "threshold": threshold,
        "train_group": train_group,
    }


def summarize(sample: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, group in sample.groupby("method"):
        by_dataset = group.groupby("dataset", as_index=False).agg(dataset_F1=("F1", "mean"))
        rec = {
            "method": method,
            "samples": int(len(group)),
            "mean_precision": float(pd.to_numeric(group["precision"], errors="coerce").mean()),
            "mean_recall": float(pd.to_numeric(group["recall"], errors="coerce").mean()),
            "mean_F1": float(pd.to_numeric(group["F1"], errors="coerce").mean()),
            "mean_FP": float(pd.to_numeric(group["FP"], errors="coerce").mean()),
            "mean_FN": float(pd.to_numeric(group["FN"], errors="coerce").mean()),
            "mean_FP_plus_FN": float(pd.to_numeric(group["FP_plus_FN"], errors="coerce").mean()),
            "min_dataset_F1": float(by_dataset["dataset_F1"].min()),
        }
        for row in by_dataset.itertuples(index=False):
            rec[f"{row.dataset}_F1"] = float(row.dataset_F1)
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["mean_F1", "mean_FP_plus_FN"], ascending=[False, True])


def evaluate_folded(df: pd.DataFrame, *, feature_mode: str, validation: str) -> pd.DataFrame:
    x, _ = feature_matrix(df, feature_mode)
    y = numeric(df, "label").to_numpy(dtype=int)
    sample = df["sample_key"].astype(str).to_numpy()
    dataset = df["dataset"].astype(str).to_numpy()
    scope = df["scope"].astype(str).to_numpy()
    gold_counts = load_gold_counts()
    base = direct_mask(df)
    requested_models = {
        name.strip()
        for name in os.environ.get("MINCO_MODEL_NAMES", "logreg,hgb").split(",")
        if name.strip()
    }
    models = {name: model for name, model in model_suite().items() if name in requested_models}
    if not models:
        raise SystemExit(f"no requested models available: {sorted(requested_models)}")

    if validation == "loso":
        groups = [(heldout, sample != heldout, sample == heldout) for heldout in sorted(set(sample))]
    elif validation == "lodo":
        groups = [(heldout, dataset != heldout, dataset == heldout) for heldout in sorted(set(dataset))]
    else:
        raise ValueError(f"unsupported validation: {validation}")

    records: list[dict[str, object]] = []
    for heldout, train_idx, test_idx in groups:
        print(f"[{validation}/{feature_mode}] held out {heldout}; models={','.join(models)}", file=sys.stderr, flush=True)
        train_prob_by_model: dict[str, np.ndarray] = {}
        test_prob_by_model: dict[str, np.ndarray] = {}
        for model_name, model in models.items():
            model.fit(x[train_idx], y[train_idx])
            train_prob_by_model[model_name] = predict_probability(model, x[train_idx])
            test_prob_by_model[model_name] = predict_probability(model, x[test_idx])

        model_probs: dict[str, tuple[np.ndarray, np.ndarray]] = dict(train_prob_by_model)
        model_probs = {name: (train_prob_by_model[name], test_prob_by_model[name]) for name in train_prob_by_model}
        if "rf" in train_prob_by_model and "hgb" in train_prob_by_model:
            model_probs["rf_hgb_avg"] = (
                0.5 * (train_prob_by_model["rf"] + train_prob_by_model["hgb"]),
                0.5 * (test_prob_by_model["rf"] + test_prob_by_model["hgb"]),
            )

        for model_name, (train_prob, test_prob) in model_probs.items():
            for hybrid, base_train, base_test in [
                ("model", None, None),
                ("direct_or_model", base[train_idx], base[test_idx]),
            ]:
                threshold = tune_threshold(
                    train_prob,
                    y[train_idx],
                    sample[train_idx],
                    gold_counts,
                    base_mask=base_train,
                )
                pred_test = test_prob >= threshold
                if base_test is not None:
                    pred_test = base_test | pred_test
                method = f"{validation}_{feature_mode}_{hybrid}_{model_name}"
                for sample_key in sorted(set(sample[test_idx])):
                    smask = test_idx & (sample == sample_key)
                    records.append(
                        score_sample(
                            method=method,
                            sample_key=sample_key,
                            dataset=str(dataset[smask][0]),
                            scope=str(scope[smask][0]),
                            pred=pred_test[sample[test_idx] == sample_key],
                            label=y[smask],
                            gold_counts=gold_counts,
                            validation=validation,
                            threshold=threshold,
                            train_group=str(heldout),
                        )
                    )
    return pd.DataFrame(records)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_features()
    sample_rows = []
    validations = [v.strip() for v in os.environ.get("MINCO_VALIDATIONS", "loso,lodo").split(",") if v.strip()]
    feature_modes = [v.strip() for v in os.environ.get("MINCO_FEATURE_MODES", "base,enriched").split(",") if v.strip()]
    for validation in validations:
        for feature_mode in feature_modes:
            sample_rows.append(evaluate_folded(df, feature_mode=feature_mode, validation=validation))

    model_samples = pd.concat(sample_rows, ignore_index=True, sort=False)
    model_samples.to_csv(OUT / "model_cv_sample_metrics.tsv", sep="\t", index=False)
    model_summary = summarize(model_samples)
    model_summary["source"] = "model_cv"
    model_summary.to_csv(OUT / "model_cv_summary.tsv", sep="\t", index=False)

    baselines = load_baseline_summary()
    baseline_summary = summarize(baselines)
    baseline_summary["source"] = "baseline"
    combined = pd.concat([model_summary, baseline_summary], ignore_index=True, sort=False)
    combined = combined.sort_values(["mean_F1", "mean_FP_plus_FN", "min_dataset_F1"], ascending=[False, True, False])
    combined.to_csv(OUT / "model_cv_combined_summary.tsv", sep="\t", index=False)
    print(combined.head(60).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
