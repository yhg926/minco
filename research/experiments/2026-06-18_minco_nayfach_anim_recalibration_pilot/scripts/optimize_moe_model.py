#!/usr/bin/env python3
"""Optimize minco MoE denominator parameters and linear coefficients."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.optimize import differential_evolution, minimize
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import train_test_split


OLD_GLOBAL_P = np.array(
    [1.0178887, 1.0135146, 1.0302181, 1.1205122,
     1.0136223, 0.9803396, 1.0021964, 0.9839395,
     0.9988935, 0.9690589, 0.9888171, 0.9817371],
    dtype=float,
)
OLD_N95_P = np.array(
    [1.0145656, 1.2708196, 1.2169829, 0.3904416,
     0.8147527, 1.3389150, 2.0581364, -1.4624736,
     1.6877678, 2.4390067, 0.7878482, 0.6697437],
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


@dataclass
class MoeFit:
    p: np.ndarray
    coef: np.ndarray
    alpha: float


def as_float(row: Dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError):
        return float("nan")


def rel_to_p(z: np.ndarray) -> np.ndarray:
    return np.array([1.0, z[0], z[1], z[2],
                     1.0, z[3], z[4], z[5],
                     1.0, z[6], z[7], z[8]], dtype=float)


def p_to_rel(p: np.ndarray) -> np.ndarray:
    return np.array([p[1] / p[0], p[2] / p[0], p[3] / p[0],
                     p[5] / p[4], p[6] / p[4], p[7] / p[4],
                     p[9] / p[8], p[10] / p[8], p[11] / p[8]], dtype=float)


def read_arrays(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows: List[Dict[str, str]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") == "ok":
                rows.append(row)
    y_ani = np.array([as_float(row, "anim_truth") for row in rows], dtype=float)
    xny = np.array([as_float(row, "minco_XnY_ctx") for row in rows], dtype=float)
    nsection = np.array([as_float(row, "minco_N_diff_obj_section") for row in rows], dtype=float)
    nmut = np.array([as_float(row, "minco_N_mut2_ctx") for row in rows], dtype=float)
    ndiff = np.array([as_float(row, "minco_N_diff_obj") for row in rows], dtype=float)
    bins = np.array([row["ani_bin"] for row in rows])
    finite = np.isfinite(y_ani) & np.isfinite(xny) & np.isfinite(nsection) & np.isfinite(nmut) & np.isfinite(ndiff)
    return y_ani[finite], xny[finite], nsection[finite], nmut[finite], ndiff[finite], bins[finite]


def design_matrix(xny: np.ndarray, nsection: np.ndarray, nmut: np.ndarray, ndiff: np.ndarray, p: np.ndarray) -> np.ndarray | None:
    eps = 1e-8
    den1 = p[0] * xny + p[1] * nsection + p[2] * nmut + p[3] * ndiff
    den2 = p[4] * xny + p[5] * nsection + p[6] * nmut + p[7] * ndiff
    den3 = p[8] * xny + p[9] * nsection + p[10] * nmut + p[11] * ndiff
    if not (np.all(np.isfinite(den1)) and np.all(np.isfinite(den2)) and np.all(np.isfinite(den3))):
        return None
    if min(np.min(np.abs(den1)), np.min(np.abs(den2)), np.min(np.abs(den3))) < 1e-5:
        return None
    d1 = 1.0 / (den1 + eps)
    d2 = 1.0 / (den2 + eps)
    d3 = 1.0 / (den3 + eps)
    return np.column_stack(
        [
            np.ones_like(xny),
            xny,
            nsection,
            nmut,
            ndiff,
            d1,
            d2,
            d3,
            xny * d1,
            nsection * d1,
            nmut * d1,
            xny * d2,
            nsection * d2,
            nmut * d2,
            xny * d3,
            nsection * d3,
            nmut * d3,
        ]
    )


def fit_scaled_linear(X: np.ndarray, y: np.ndarray, alpha: float) -> Tuple[np.ndarray, np.ndarray]:
    mean = X[:, 1:].mean(axis=0)
    scale = X[:, 1:].std(axis=0)
    scale[scale == 0.0] = 1.0
    Xs = np.column_stack([np.ones(X.shape[0]), (X[:, 1:] - mean) / scale])
    reg = np.eye(Xs.shape[1]) * alpha
    reg[0, 0] = 0.0
    beta_scaled = np.linalg.solve(Xs.T @ Xs + reg, Xs.T @ y)
    coef = np.empty_like(beta_scaled)
    coef[1:] = beta_scaled[1:] / scale
    coef[0] = beta_scaled[0] - np.sum(beta_scaled[1:] * mean / scale)
    pred = X @ coef
    return coef, pred


def fit_best_alpha(X: np.ndarray, y: np.ndarray, alphas: Sequence[float]) -> Tuple[np.ndarray, float, np.ndarray]:
    best = None
    for alpha in alphas:
        coef, pred = fit_scaled_linear(X, y, alpha)
        mse = float(np.mean((pred - y) ** 2))
        if best is None or mse < best[0]:
            best = (mse, coef, alpha, pred)
    assert best is not None
    return best[1], best[2], best[3]


def objective_factory(xny: np.ndarray, nsection: np.ndarray, nmut: np.ndarray, ndiff: np.ndarray, y_dist: np.ndarray, alpha: float):
    def objective(z: np.ndarray) -> float:
        p = rel_to_p(z)
        X = design_matrix(xny, nsection, nmut, ndiff, p)
        if X is None:
            return 1e6
        try:
            _, pred = fit_scaled_linear(X, y_dist, alpha)
        except np.linalg.LinAlgError:
            return 1e6
        residual = pred - y_dist
        mse = float(np.mean(residual * residual))
        out_penalty = float(np.mean(np.maximum(-pred, 0.0) ** 2 + np.maximum(pred - 1.0, 0.0) ** 2))
        return mse + 0.05 * out_penalty
    return objective


def optimize_p(label: str,
               xny: np.ndarray,
               nsection: np.ndarray,
               nmut: np.ndarray,
               ndiff: np.ndarray,
               y_dist: np.ndarray,
               start_p: np.ndarray,
               seed: int,
               maxiter: int,
               popsize: int,
               polish_iter: int) -> np.ndarray:
    bounds = [(-3.0, 4.0)] * 9
    obj = objective_factory(xny, nsection, nmut, ndiff, y_dist, alpha=1e-6)
    start = p_to_rel(start_p)
    start_score = obj(start)
    print(f"{label}: start objective {start_score:.8g}", flush=True)
    result = differential_evolution(
        obj,
        bounds=bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        polish=False,
        updating="immediate",
        workers=1,
        x0=start,
        tol=1e-5,
    )
    print(f"{label}: DE objective {result.fun:.8g}", flush=True)
    local = minimize(
        obj,
        result.x,
        method="Powell",
        bounds=bounds,
        options={"maxiter": polish_iter, "xtol": 1e-5, "ftol": 1e-7, "disp": False},
    )
    best = local.x if local.fun < result.fun else result.x
    print(f"{label}: final objective {min(local.fun, result.fun):.8g}", flush=True)
    return rel_to_p(best)


def predict_dist(fit: MoeFit,
                 xny: np.ndarray,
                 nsection: np.ndarray,
                 nmut: np.ndarray,
                 ndiff: np.ndarray) -> np.ndarray:
    X = design_matrix(xny, nsection, nmut, ndiff, fit.p)
    if X is None:
        raise RuntimeError("invalid design matrix")
    return X @ fit.coef


def pearson_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
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


def write_c_arrays(path: Path, global_fit: MoeFit, n95_fit: MoeFit) -> None:
    def fmt_arr(name: str, values: np.ndarray, comments: Sequence[str] | None = None) -> List[str]:
        lines = [f"static const double {name}[{len(values)}] =", "    {"]
        for i, value in enumerate(values):
            suffix = "," if i + 1 < len(values) else ""
            comment = f" // {comments[i]}" if comments else ""
            lines.append(f"        {float(value):.15g}{suffix}{comment}")
        lines.append("    };")
        return lines

    lines: List[str] = []
    lines.extend(fmt_arr("MINCO_OPT_DENOM_PARAMS4X3", global_fit.p, [
        "denom1", "denom1", "denom1", "denom1",
        "denom2", "denom2", "denom2", "denom2",
        "denom3", "denom3", "denom3", "denom3",
    ]))
    lines.append("")
    lines.extend(fmt_arr("MINCO_N95OPT_DENOM_PARAMS4X3", n95_fit.p, [
        "denom1", "denom1", "denom1", "denom1",
        "denom2", "denom2", "denom2", "denom2",
        "denom3", "denom3", "denom3", "denom3",
    ]))
    lines.append("")
    lines.extend(fmt_arr("MINCO_LINEAR_COEFFS_3WAY", global_fit.coef, FEATURE_NAMES))
    lines.append("")
    lines.extend(fmt_arr("MINCO_N95LINEAR_COEFFS_3WAY", n95_fit.coef, FEATURE_NAMES))
    lines.append("")
    path.write_text("\n".join(lines))


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--maxiter", type=int, default=18)
    parser.add_argument("--popsize", type=int, default=7)
    parser.add_argument("--polish-iter", type=int, default=80)
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    y_ani, xny, nsection, nmut, ndiff, bins = read_arrays(args.features)
    y_dist = 1.0 - y_ani
    idx = np.arange(len(y_ani))
    train_idx, test_idx = train_test_split(idx, test_size=0.30, random_state=args.seed, stratify=bins)

    p_global = optimize_p(
        "global",
        xny[train_idx],
        nsection[train_idx],
        nmut[train_idx],
        ndiff[train_idx],
        y_dist[train_idx],
        OLD_GLOBAL_P,
        args.seed,
        args.maxiter,
        args.popsize,
        args.polish_iter,
    )
    n95_train = train_idx[y_ani[train_idx] >= 0.95]
    p_n95 = optimize_p(
        "n95",
        xny[n95_train],
        nsection[n95_train],
        nmut[n95_train],
        ndiff[n95_train],
        y_dist[n95_train],
        OLD_N95_P,
        args.seed + 95,
        args.maxiter,
        args.popsize,
        args.polish_iter,
    )

    alphas = (1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 1.0)
    X_global_train = design_matrix(xny[train_idx], nsection[train_idx], nmut[train_idx], ndiff[train_idx], p_global)
    assert X_global_train is not None
    coef_global, alpha_global, _ = fit_best_alpha(X_global_train, y_dist[train_idx], alphas)
    X_n95_train = design_matrix(xny[n95_train], nsection[n95_train], nmut[n95_train], ndiff[n95_train], p_n95)
    assert X_n95_train is not None
    coef_n95, alpha_n95, _ = fit_best_alpha(X_n95_train, y_dist[n95_train], alphas)

    global_fit = MoeFit(p=p_global, coef=coef_global, alpha=alpha_global)
    n95_fit = MoeFit(p=p_n95, coef=coef_n95, alpha=alpha_n95)
    dist_global = np.clip(predict_dist(global_fit, xny, nsection, nmut, ndiff), 0.0, 1.0)
    dist_twostage = dist_global.copy()
    gate = dist_global < 0.05
    if np.any(gate):
        dist_twostage[gate] = np.clip(predict_dist(n95_fit, xny[gate], nsection[gate], nmut[gate], ndiff[gate]), 0.0, 1.0)

    preds = {
        "moe_opt_global": 1.0 - dist_global,
        "moe_opt_twostage": 1.0 - dist_twostage,
    }

    metrics_path = args.outdir / "optimized_moe_metrics.tsv"
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
        for name, pred in preds.items():
            for split, mask in masks.items():
                out = {"model": name, "split": split}
                out.update(metrics(y_ani[mask], pred[mask]))
                writer.writerow(out)

    params_path = args.outdir / "optimized_moe_params.tsv"
    with params_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("block", "index", "value"))
        for i, value in enumerate(global_fit.p):
            writer.writerow(("global_p", i, value))
        for i, value in enumerate(n95_fit.p):
            writer.writerow(("n95_p", i, value))
        for i, value in enumerate(global_fit.coef):
            writer.writerow(("global_coef", i, value))
        writer.writerow(("global_alpha", 0, global_fit.alpha))
        for i, value in enumerate(n95_fit.coef):
            writer.writerow(("n95_coef", i, value))
        writer.writerow(("n95_alpha", 0, n95_fit.alpha))

    write_c_arrays(args.outdir / "optimized_moe_model_ani_arrays.hfrag", global_fit, n95_fit)
    print(f"metrics -> {metrics_path}")
    print(f"params -> {params_path}")
    print(f"C arrays -> {args.outdir / 'optimized_moe_model_ani_arrays.hfrag'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
