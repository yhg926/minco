#!/usr/bin/env python3
"""Fit drop-in-style MoE linear coefficients with current denominator params."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler


P_GLOBAL = np.array(
    [
        1.0178887,
        1.0135146,
        1.0302181,
        1.1205122,
        1.0136223,
        0.9803396,
        1.0021964,
        0.9839395,
        0.9988935,
        0.9690589,
        0.9888171,
        0.9817371,
    ],
    dtype=float,
)

P_N95 = np.array(
    [
        1.0145656,
        1.2708196,
        1.2169829,
        0.3904416,
        0.8147527,
        1.3389150,
        2.0581364,
        -1.4624736,
        1.6877678,
        2.4390067,
        0.7878482,
        0.6697437,
    ],
    dtype=float,
)

FEATURE_NAMES = (
    "intercept",
    "XnY_ctx",
    "N_diff_obj_section",
    "N_mut2_ctx",
    "N_diff_obj",
    "denom1",
    "denom2",
    "denom3",
    "XnY_ctx:denom1",
    "N_diff_obj_section:denom1",
    "N_mut2_ctx:denom1",
    "XnY_ctx:denom2",
    "N_diff_obj_section:denom2",
    "N_mut2_ctx:denom2",
    "XnY_ctx:denom3",
    "N_diff_obj_section:denom3",
    "N_mut2_ctx:denom3",
)


def as_float(row: Dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return float("nan")


def design_from_counts(xny: np.ndarray, nsection: np.ndarray, nmut: np.ndarray, ndiff: np.ndarray, p: np.ndarray) -> np.ndarray:
    eps = 1e-8
    denom1 = 1.0 / (p[0] * xny + p[1] * nsection + p[2] * nmut + p[3] * ndiff + eps)
    denom2 = 1.0 / (p[4] * xny + p[5] * nsection + p[6] * nmut + p[7] * ndiff + eps)
    denom3 = 1.0 / (p[8] * xny + p[9] * nsection + p[10] * nmut + p[11] * ndiff + eps)
    return np.column_stack(
        [
            np.ones_like(xny),
            xny,
            nsection,
            nmut,
            ndiff,
            denom1,
            denom2,
            denom3,
            xny * denom1,
            nsection * denom1,
            nmut * denom1,
            xny * denom2,
            nsection * denom2,
            nmut * denom2,
            xny * denom3,
            nsection * denom3,
            nmut * denom3,
        ]
    )


def read_rows(path: Path) -> Tuple[List[Dict[str, str]], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows: List[Dict[str, str]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "ok":
                continue
            rows.append(row)
    y_ani = np.array([as_float(row, "anim_truth") for row in rows], dtype=float)
    y_dist = 1.0 - y_ani
    xny = np.array([as_float(row, "minco_XnY_ctx") for row in rows], dtype=float)
    nsection = np.array([as_float(row, "minco_N_diff_obj_section") for row in rows], dtype=float)
    nmut = np.array([as_float(row, "minco_N_mut2_ctx") for row in rows], dtype=float)
    ndiff = np.array([as_float(row, "minco_N_diff_obj") for row in rows], dtype=float)
    bins = np.array([row["ani_bin"] for row in rows])
    return rows, y_ani, y_dist, xny, nsection, nmut, ndiff, bins


def pearson_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2:
        return float("nan")
    y0 = y_true - y_true.mean()
    y1 = y_pred - y_pred.mean()
    den = math.sqrt(float((y0 * y0).sum() * (y1 * y1).sum()))
    return float((y0 * y1).sum() / den) if den > 0 else float("nan")


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    y_pred = np.clip(y_pred, 0.0, 1.0)
    err = y_pred - y_true
    truth_pos = y_true >= 0.95
    pred_pos = y_pred >= 0.95
    return {
        "n": float(len(y_true)),
        "bias": float(err.mean()),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "max_abs": float(np.max(np.abs(err))),
        "pearson_r": pearson_r(y_true, y_pred),
        "tp95": float(np.sum(truth_pos & pred_pos)),
        "fp95": float(np.sum(~truth_pos & pred_pos)),
        "fn95": float(np.sum(truth_pos & ~pred_pos)),
        "tn95": float(np.sum(~truth_pos & ~pred_pos)),
    }


def fit_scaled_ridge(X: np.ndarray, y: np.ndarray) -> Tuple[RidgeCV, StandardScaler, np.ndarray]:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X[:, 1:])
    model = RidgeCV(alphas=(1e-8, 1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0), fit_intercept=True)
    model.fit(X_scaled, y)
    coef_raw_no_intercept = model.coef_ / scaler.scale_
    intercept_raw = model.intercept_ - float(np.sum(model.coef_ * scaler.mean_ / scaler.scale_))
    coef_raw = np.concatenate([[intercept_raw], coef_raw_no_intercept])
    return model, scaler, coef_raw


def predict_scaled(model: RidgeCV, scaler: StandardScaler, X: np.ndarray) -> np.ndarray:
    return model.predict(scaler.transform(X[:, 1:]))


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260618)
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows, y_ani, y_dist, xny, nsection, nmut, ndiff, bins = read_rows(args.features)
    idx = np.arange(len(rows))
    train_idx, test_idx = train_test_split(idx, test_size=0.30, random_state=args.seed, stratify=bins)

    X_global = design_from_counts(xny, nsection, nmut, ndiff, P_GLOBAL)
    X_n95 = design_from_counts(xny, nsection, nmut, ndiff, P_N95)
    global_model, global_scaler, global_coef_raw = fit_scaled_ridge(X_global[train_idx], y_dist[train_idx])
    n95_train_mask = np.isin(idx, train_idx) & (y_ani >= 0.95)
    n95_model, n95_scaler, n95_coef_raw = fit_scaled_ridge(X_n95[n95_train_mask], y_dist[n95_train_mask])

    global_dist = np.clip(predict_scaled(global_model, global_scaler, X_global), 0.0, 1.0)
    two_stage_dist = global_dist.copy()
    n95_gate = global_dist < 0.05
    if np.any(n95_gate):
        two_stage_dist[n95_gate] = np.clip(predict_scaled(n95_model, n95_scaler, X_n95[n95_gate]), 0.0, 1.0)
    preds = {
        "moe_fixedp_global": 1.0 - global_dist,
        "moe_fixedp_twostage": 1.0 - two_stage_dist,
    }

    metrics_path = args.outdir / "fixedp_moe_metrics.tsv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("model", "split", "n", "bias", "mae", "rmse", "max_abs", "pearson_r", "tp95", "fp95", "fn95", "tn95"),
            delimiter="\t",
        )
        writer.writeheader()
        masks = {
            "train": np.isin(idx, train_idx),
            "test": np.isin(idx, test_idx),
            "all": np.ones_like(idx, dtype=bool),
        }
        for model_name, pred in preds.items():
            for split, mask in masks.items():
                out = {"model": model_name, "split": split}
                out.update(metrics(y_ani[mask], pred[mask]))
                writer.writerow(out)

    coef_path = args.outdir / "fixedp_moe_coefficients.tsv"
    with coef_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("model", "feature", "coefficient"))
        for name, coef in zip(FEATURE_NAMES, global_coef_raw):
            writer.writerow(("global", name, coef))
        writer.writerow(("global", "alpha", global_model.alpha_))
        for name, coef in zip(FEATURE_NAMES, n95_coef_raw):
            writer.writerow(("n95", name, coef))
        writer.writerow(("n95", "alpha", n95_model.alpha_))
    print(f"metrics -> {metrics_path}")
    print(f"coefficients -> {coef_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
