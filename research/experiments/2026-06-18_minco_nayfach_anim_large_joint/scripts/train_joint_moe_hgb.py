#!/usr/bin/env python3
"""Train drop-in HGB models on newly optimized MoE predictions.

This script is the joint-training bridge for the large Nayfach calibration run:
`optimize_moe_model.py` fits the raw CtxMoE counts model, then this script
computes that optimized raw ANI for every feature row and trains the C-compatible
11-feature HGB layer on the same train/test split.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import train_test_split


SCRIPT_DIR = Path(__file__).resolve().parent
PILOT_SCRIPTS = SCRIPT_DIR.parents[1] / "2026-06-18_minco_nayfach_anim_recalibration_pilot" / "scripts"
sys.path.insert(0, str(PILOT_SCRIPTS))

from minco_nayfach_calibration_pilot import export_refaf_hgb_header  # noqa: E402
from optimize_moe_model import MoeFit, predict_dist  # noqa: E402


DROPIN_FEATURES = (
    "raw_ani",
    "ref_af",
    "log1p_xny",
    "diff_obj_rate",
    "diff_section_rate",
    "mut2_rate",
    "ref_aaf_ani",
    "raw_minus_ref_aaf",
    "ref_af_lt_0_2",
    "ref_af_0_2_to_0_5",
    "ref_af_ge_0_5",
)


@dataclass(frozen=True)
class Candidate:
    name: str
    max_iter: int
    learning_rate: float
    max_leaf_nodes: int
    l2_regularization: float


HGB_CANDIDATES = (
    Candidate("hgb_moe_160", 160, 0.05, 31, 0.01),
    Candidate("hgb_moe_240", 240, 0.04, 31, 0.01),
)


def as_float(row: Dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def clamp01(value: float) -> float:
    if not np.isfinite(value):
        return value
    return min(1.0, max(0.0, value))


def safe_div(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0.0:
        return float("nan")
    return num / den


def read_feature_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row.get("status") == "ok"]


def read_moe_params(path: Path) -> Tuple[MoeFit, MoeFit]:
    blocks: Dict[str, Dict[int, float]] = {
        "global_p": {},
        "n95_p": {},
        "global_coef": {},
        "n95_coef": {},
        "global_alpha": {},
        "n95_alpha": {},
    }
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            block = row["block"]
            if block not in blocks:
                continue
            blocks[block][int(row["index"])] = float(row["value"])

    def arr(block: str, n: int) -> np.ndarray:
        values = blocks[block]
        missing = [i for i in range(n) if i not in values]
        if missing:
            raise ValueError(f"{path}: block {block} missing indices {missing[:5]}")
        return np.array([values[i] for i in range(n)], dtype=float)

    global_fit = MoeFit(
        p=arr("global_p", 12),
        coef=arr("global_coef", 17),
        alpha=blocks["global_alpha"].get(0, float("nan")),
    )
    n95_fit = MoeFit(
        p=arr("n95_p", 12),
        coef=arr("n95_coef", 17),
        alpha=blocks["n95_alpha"].get(0, float("nan")),
    )
    return global_fit, n95_fit


def counts_arrays(rows: Sequence[Dict[str, str]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xny = np.array([as_float(row, "minco_XnY_ctx") for row in rows], dtype=float)
    nsection = np.array([as_float(row, "minco_N_diff_obj_section") for row in rows], dtype=float)
    nmut = np.array([as_float(row, "minco_N_mut2_ctx") for row in rows], dtype=float)
    ndiff = np.array([as_float(row, "minco_N_diff_obj") for row in rows], dtype=float)
    return xny, nsection, nmut, ndiff


def optimized_moe_predictions(rows: Sequence[Dict[str, str]], global_fit: MoeFit, n95_fit: MoeFit) -> np.ndarray:
    xny, nsection, nmut, ndiff = counts_arrays(rows)
    dist_global = np.clip(predict_dist(global_fit, xny, nsection, nmut, ndiff), 0.0, 1.0)
    dist_twostage = dist_global.copy()
    gate = dist_global < 0.05
    if np.any(gate):
        dist_twostage[gate] = np.clip(
            predict_dist(n95_fit, xny[gate], nsection[gate], nmut[gate], ndiff[gate]),
            0.0,
            1.0,
        )
    return np.clip(1.0 - dist_twostage, 0.0, 1.0)


def dropin_matrix(rows: Sequence[Dict[str, str]], raw_ani: np.ndarray, fold_scale: float) -> np.ndarray:
    values: List[List[float]] = []
    for row, raw in zip(rows, raw_ani):
        raw = clamp01(float(raw))
        ref_af = as_float(row, "minco_Ref_align_fraction")
        if not np.isfinite(ref_af) or ref_af < 1e-12:
            ref_af = 1e-12
        ref_af = min(1.0, ref_af)
        xny = as_float(row, "minco_XnY_ctx")
        nmut = as_float(row, "minco_N_mut2_ctx")
        ndiff = as_float(row, "minco_N_diff_obj")
        nsection = as_float(row, "minco_N_diff_obj_section")
        ref_aaf_ani = clamp01(1.0 + math.log(ref_af) / 11.0)
        values.append(
            [
                raw,
                ref_af,
                math.log1p(xny * fold_scale) if np.isfinite(xny) and xny >= 0 else float("nan"),
                safe_div(ndiff, xny),
                safe_div(nsection, xny),
                safe_div(nmut, xny),
                ref_aaf_ani,
                raw - ref_aaf_ani if np.isfinite(raw) else float("nan"),
                1.0 if ref_af < 0.2 else 0.0,
                1.0 if 0.2 <= ref_af < 0.5 else 0.0,
                1.0 if ref_af >= 0.5 else 0.0,
            ]
        )
    return np.array(values, dtype=float)


def pearson_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2:
        return float("nan")
    y0 = y_true - y_true.mean()
    y1 = y_pred - y_pred.mean()
    den = math.sqrt(float((y0 * y0).sum() * (y1 * y1).sum()))
    return float((y0 * y1).sum() / den) if den > 0 else float("nan")


def metric_row(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    y_pred = np.clip(y_pred, 0.0, 1.0)
    err = y_pred - y_true
    truth_pos = y_true >= 0.95
    pred_pos = y_pred >= 0.95
    return {
        "n": float(len(y_true)),
        "bias": float(err.mean()),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "max_abs": float(np.max(np.abs(err))) if len(err) else float("nan"),
        "pearson_r": pearson_r(y_true, y_pred),
        "tp95": float(np.sum(truth_pos & pred_pos)),
        "fp95": float(np.sum(~truth_pos & pred_pos)),
        "fn95": float(np.sum(truth_pos & ~pred_pos)),
        "tn95": float(np.sum(~truth_pos & ~pred_pos)),
    }


def cutoff_errors(row: Dict[str, float]) -> float:
    return float(row["fp95"] + row["fn95"])


def train_hgb(candidate: Candidate, X: np.ndarray, y: np.ndarray, train_idx: np.ndarray, seed: int) -> Tuple[HistGradientBoostingRegressor, np.ndarray]:
    model = HistGradientBoostingRegressor(
        max_iter=candidate.max_iter,
        learning_rate=candidate.learning_rate,
        max_leaf_nodes=candidate.max_leaf_nodes,
        l2_regularization=candidate.l2_regularization,
        random_state=seed,
    )
    model.fit(X[train_idx], y[train_idx])
    return model, np.clip(model.predict(X), 0.0, 1.0)


def write_metrics(
    path: Path,
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    bins: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    predictions: Dict[str, np.ndarray],
) -> None:
    indices = np.arange(len(rows))
    split_masks = {
        "train": np.isin(indices, train_idx),
        "test": np.isin(indices, test_idx),
        "all": np.ones(len(rows), dtype=bool),
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("model", "split", "n", "bias", "mae", "rmse", "max_abs", "pearson_r", "tp95", "fp95", "fn95", "tn95"),
            delimiter="\t",
        )
        writer.writeheader()
        for name, pred in predictions.items():
            for split, mask in split_masks.items():
                out = {"model": name, "split": split}
                out.update(metric_row(y[mask], pred[mask]))
                writer.writerow(out)


def write_metrics_by_bin(
    path: Path,
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    bins: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    predictions: Dict[str, np.ndarray],
) -> None:
    indices = np.arange(len(rows))
    split_indices = {"train": train_idx, "test": test_idx, "all": indices}
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("model", "split", "ani_bin", "n", "bias", "mae", "rmse", "max_abs", "pearson_r", "tp95", "fp95", "fn95", "tn95"),
            delimiter="\t",
        )
        writer.writeheader()
        for name, pred in predictions.items():
            for split, split_idx in split_indices.items():
                for label in sorted(set(bins)):
                    mask = np.isin(indices, split_idx) & (bins == label)
                    if not np.any(mask):
                        continue
                    out = {"model": name, "split": split, "ani_bin": label}
                    out.update(metric_row(y[mask], pred[mask]))
                    writer.writerow(out)


def write_predictions(
    path: Path,
    rows: Sequence[Dict[str, str]],
    y: np.ndarray,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    predictions: Dict[str, np.ndarray],
    selected_name: str,
) -> None:
    split_lookup = {int(i): "train" for i in train_idx}
    split_lookup.update({int(i): "test" for i in test_idx})
    fieldnames = [
        "pair_id",
        "split",
        "ani_bin",
        "anim_truth",
        *predictions.keys(),
        "hgb_refaf11_dropin",
        "selected_hgb_model",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows):
            out = {
                "pair_id": row["pair_id"],
                "split": split_lookup[idx],
                "ani_bin": row["ani_bin"],
                "anim_truth": y[idx],
                "hgb_refaf11_dropin": predictions[selected_name][idx],
                "selected_hgb_model": selected_name,
            }
            for name, pred in predictions.items():
                out[name] = pred[idx]
            writer.writerow(out)


def replace_generated_comment(path: Path) -> None:
    text = path.read_text()
    text = text.replace(
        "/* Generated by minco_nayfach_calibration_pilot.py.",
        "/* Generated by train_joint_moe_hgb.py.",
        1,
    )
    path.write_text(text)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--moe-params", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--fold-scale", type=float, default=1.0)
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_feature_rows(args.features)
    if len(rows) < 100:
        raise RuntimeError(f"too few ok feature rows: {len(rows)}")

    y = np.array([as_float(row, "anim_truth") for row in rows], dtype=float)
    bins = np.array([row["ani_bin"] for row in rows])
    indices = np.arange(len(rows))
    train_idx, test_idx = train_test_split(indices, test_size=0.30, random_state=args.seed, stratify=bins)

    global_fit, n95_fit = read_moe_params(args.moe_params)
    current_raw = np.array([as_float(row, "minco_ANI") for row in rows], dtype=float)
    optimized_raw = optimized_moe_predictions(rows, global_fit, n95_fit)

    predictions: Dict[str, np.ndarray] = {
        "current_raw_ctxmoe": np.clip(current_raw, 0.0, 1.0),
        "optimized_moe_twostage": optimized_raw,
    }

    X_current = dropin_matrix(rows, predictions["current_raw_ctxmoe"], args.fold_scale)
    current_candidate = Candidate("hgb_currentraw_160", 160, 0.05, 31, 0.01)
    _, current_hgb_pred = train_hgb(current_candidate, X_current, y, train_idx, args.seed)
    predictions[current_candidate.name] = current_hgb_pred

    X_optimized = dropin_matrix(rows, optimized_raw, args.fold_scale)
    hgb_models: Dict[str, HistGradientBoostingRegressor] = {}
    for offset, candidate in enumerate(HGB_CANDIDATES):
        model, pred = train_hgb(candidate, X_optimized, y, train_idx, args.seed + offset)
        hgb_models[candidate.name] = model
        predictions[candidate.name] = pred

    # Choose an install candidate by held-out MAE, then by 0.95 cutoff errors.
    test_scores = []
    for name in hgb_models:
        row = metric_row(y[test_idx], predictions[name][test_idx])
        test_scores.append((row["mae"], cutoff_errors(row), name))
    test_scores.sort()
    selected_name = test_scores[0][2]

    write_metrics(args.outdir / "joint_model_metrics.tsv", rows, y, bins, train_idx, test_idx, predictions)
    write_metrics_by_bin(args.outdir / "joint_model_metrics_by_bin.tsv", rows, y, bins, train_idx, test_idx, predictions)
    write_predictions(args.outdir / "joint_model_predictions.tsv", rows, y, train_idx, test_idx, predictions, selected_name)

    header_path = args.outdir / "model_refaf_hgb.joint_generated.h"
    export_refaf_hgb_header(
        hgb_models[selected_name],
        header_path,
        note=f"Joint large-run HGB on optimized MoE raw ANI. Selected {selected_name}; trained on {len(train_idx)} stratified Nayfach pairs, tested on {len(test_idx)}; sketch-size 10000.",
    )
    replace_generated_comment(header_path)

    with (args.outdir / "joint_selected_model.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("key", "value"))
        writer.writerow(("selected_hgb_model", selected_name))
        writer.writerow(("train_pairs", len(train_idx)))
        writer.writerow(("test_pairs", len(test_idx)))
        writer.writerow(("feature_rows", len(rows)))
        writer.writerow(("fold_scale", args.fold_scale))
        for mae, errors, name in test_scores:
            writer.writerow((f"{name}_test_mae", mae))
            writer.writerow((f"{name}_test_95cutoff_errors", errors))

    print(f"rows={len(rows)} train={len(train_idx)} test={len(test_idx)} selected={selected_name}")
    print(f"metrics -> {args.outdir / 'joint_model_metrics.tsv'}")
    print(f"predictions -> {args.outdir / 'joint_model_predictions.tsv'}")
    print(f"header -> {header_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
