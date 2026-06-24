#!/usr/bin/env python3
"""Calibrate MinCO abundance from selected-call feature rows.

This is an offline experiment. It imports cached MinCO result rows from the
previous threshold-search pipeline, keeps presence callsets fixed, and tests
whether abundance mass can be improved by learning a per-selected-species score.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import HuberRegressor, Ridge, TweedieRegressor
from sklearn.preprocessing import StandardScaler


warnings.filterwarnings("ignore", message="scipy._lib.messagestream.MessageStream size changed")
warnings.filterwarnings("ignore", category=ConvergenceWarning)

EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
ABN_PATH = ROOT / "research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py"
OUT = EXP / "results"
EPS = 1e-12


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


abn = load_module(ABN_PATH, "abundance_norm_prev")
prev = abn.prev


@dataclass(frozen=True)
class CallSetDef:
    name: str
    marker_combo: object | None
    ctxobj_combo: object | None
    add_abs_min: float = 0.0
    add_rel_min: float | None = None
    include_ctxobj_overlap: bool = False
    ctxobj_only: bool = False


@dataclass(frozen=True)
class ModelSpec:
    name: str
    estimator: str
    target_mode: str
    feature_set: str
    weight_mode: str


NUMERIC_ROW_COLS = [
    "XnY_ctx",
    "Raw_XnY_ctx",
    "ANI",
    "ANI_naive_calc",
    "ANI_AF_delta",
    "ANI_from_Reliable_ztp_af",
    "Reliable_Ref_breadth",
    "Reliable_Ref_mean_depth",
    "Reliable_Ref_hit_ctx",
    "Reliable_Ref_hit_mean_depth",
    "Reliable_Ref_hit_median_depth",
    "Reliable_Ref_hit_depth_variance",
    "Reliable_Ref_zip_af",
    "Reliable_ztp_af",
    "Reliable_depth_vmr",
    "Normalized_abundance_depth",
    "Normalized_effective_abundance_depth",
    "Effective_abundance_depth",
    "Relative_abundance_depth",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_variance",
    "Ref_depth_cv",
    "Ref_zip_af",
    "Ref_zero_fraction",
    "Ref_align_fraction",
    "Real_Ref_align_fraction",
    "blastn_Ref_align_fraction",
    "Fake_ctx_fraction",
    "Fake_ctx_prob_mean",
    "Fake_ctx_prob_weighted",
    "Unique_query_ctx",
    "Unique_query_ctx_hit",
    "Unique_ref_ctx_hit",
    "Blocks_with_ctx_match",
    "Reads_with_ctx_match",
    "Density_block_ctx",
    "ctx_marker_size",
    "N_diff_obj",
    "N_mut2_ctx",
    "Rejected_ctx",
    "Rejected_diff_ctx",
]

VALUE_FEATURES = [
    "marker_value",
    "add_value",
    "value_sum",
    "value_max",
    "value_marker_first",
    "marker_present",
    "add_present",
    "both_present",
    "marker_fraction",
    "add_fraction",
    "marker_norm",
    "add_norm",
    "sum_norm",
    "marker_rank_frac",
    "add_rank_frac",
    "sum_rank_frac",
    "sample_marker_total",
    "sample_add_total",
    "sample_sum_total",
    "sample_selected_count",
]


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except Exception:
        return default
    return out if math.isfinite(out) else default


def add_row_features(rec: dict[str, object], prefix: str, row: pd.Series | None) -> None:
    for col in NUMERIC_ROW_COLS:
        rec[f"{prefix}_{col}"] = finite_float(row[col]) if row is not None and col in row.index else 0.0


def candidate_records(rows: pd.DataFrame, mask: pd.Series, depths: pd.Series, source: str, scale: float, reason: str) -> list[dict[str, object]]:
    if rows.empty or not bool(mask.any()):
        return []
    selected = rows.loc[mask].copy()
    out: list[dict[str, object]] = []
    for idx, row in selected.iterrows():
        target = str(row.get("target_id", ""))
        value = finite_float(depths.loc[idx]) * scale
        if not target or value <= 0.0:
            continue
        rec: dict[str, object] = {
            "target_id": target,
            "candidate_value": value,
            "candidate_source": source,
            "candidate_reason": reason,
        }
        add_row_features(rec, "row", row)
        out.append(rec)
    return out


def select_best_candidate(candidates: list[dict[str, object]], target: str) -> dict[str, object] | None:
    best = None
    best_value = -1.0
    for rec in candidates:
        if rec["target_id"] == target and float(rec["candidate_value"]) > best_value:
            best = rec
            best_value = float(rec["candidate_value"])
    return best


def values_and_features_for_combo(sample, combo_obj, source: str) -> tuple[dict[str, float], dict[str, dict[str, object]]]:
    if sample.marker.empty:
        return {}, {}

    candidates: list[dict[str, object]] = []
    marker_depth = prev.effective_depth(sample.marker, combo_obj)
    marker_active = prev.active_mask(sample.marker, combo_obj)
    marker_values = prev.collect_values(sample.marker, marker_active, marker_depth, 1.0)
    candidates.extend(candidate_records(sample.marker, marker_active, marker_depth, source, 1.0, "active"))

    rescued = prev.rescue_values(sample.marker, marker_active, marker_depth, marker_values, combo_obj)
    if rescued:
        rescue_mask = sample.marker["target_id"].isin(rescued)
        candidates.extend(candidate_records(sample.marker, rescue_mask, marker_depth, source, 1.0, "rescue"))

    values = dict(marker_values)
    if combo_obj.fallback_marker_cutoff > 0 and not sample.full.empty:
        full_base = (
            (sample.full["ctx_marker_size"] < combo_obj.fallback_marker_cutoff)
            & (sample.full["XnY_ctx"] >= combo_obj.fallback_xny_min)
        )
        if bool(full_base.any()):
            full_rows = sample.full.loc[full_base].copy()
            full_depth = prev.effective_depth(full_rows, combo_obj)
            full_active = prev.active_mask(full_rows, combo_obj)
            full_values = prev.collect_values(full_rows, full_active, full_depth, combo_obj.fallback_scale)
            for target in list(full_values):
                if target in values:
                    del full_values[target]
            if full_values:
                kept_targets = set(full_values)
                keep_mask = full_active & full_rows["target_id"].isin(kept_targets)
                candidates.extend(
                    candidate_records(full_rows, keep_mask, full_depth, f"{source}_full_fallback", combo_obj.fallback_scale, "fallback")
                )
                prev.merge_values(values, full_values)

    features: dict[str, dict[str, object]] = {}
    for target, value in values.items():
        best = select_best_candidate(candidates, target)
        rec: dict[str, object] = {
            "target_id": target,
            "value": float(value),
            "source": source,
            "reason": "unknown" if best is None else str(best["candidate_reason"]),
        }
        if best is not None:
            for key, val in best.items():
                if key.startswith("row_"):
                    rec[key] = val
        else:
            for col in NUMERIC_ROW_COLS:
                rec[f"row_{col}"] = 0.0
        features[target] = rec
    return values, features


def callsets() -> list[CallSetDef]:
    return [
        CallSetDef("marker_l1_only", abn.MARKER_L1, None),
        CallSetDef("marker_f1_only", abn.MARKER_F1, None),
        CallSetDef("ctxobj_f1_only", None, abn.CTXOBJ_F1, ctxobj_only=True),
        CallSetDef("marker_l1_ctxobj_f1_all", abn.MARKER_L1, abn.CTXOBJ_F1),
        CallSetDef("marker_l1_ctxobj_l1_blend", abn.MARKER_L1, abn.CTXOBJ_L1, include_ctxobj_overlap=True),
        CallSetDef("marker_l1_ctxobj_f1_blend", abn.MARKER_L1, abn.CTXOBJ_F1, include_ctxobj_overlap=True),
        CallSetDef("marker_l1_ctxobj_cami_blend", abn.MARKER_L1, abn.CTXOBJ_CAMI, include_ctxobj_overlap=True),
        CallSetDef("marker_l1_ctxobj_f1_rel0002", abn.MARKER_L1, abn.CTXOBJ_F1, add_rel_min=0.002),
        CallSetDef("marker_f1_ctxobj_f1_all", abn.MARKER_F1, abn.CTXOBJ_F1),
    ]


def sample_callset_rows(marker_sample, ctx_sample, cs: CallSetDef) -> list[dict[str, object]]:
    if cs.ctxobj_only:
        marker_values: dict[str, float] = {}
        marker_features: dict[str, dict[str, object]] = {}
        add_values, add_features = values_and_features_for_combo(ctx_sample, cs.ctxobj_combo, "ctxobj")
    else:
        marker_values, marker_features = values_and_features_for_combo(marker_sample, cs.marker_combo, "marker")
        add_values, add_features = ({}, {})
        if cs.ctxobj_combo is not None:
            add_values, add_features = values_and_features_for_combo(ctx_sample, cs.ctxobj_combo, "ctxobj")

    marker_max = max(marker_values.values()) if marker_values else 0.0
    filtered_add: dict[str, float] = {}
    for target, value in add_values.items():
        if target in marker_values and not cs.include_ctxobj_overlap:
            continue
        if value < cs.add_abs_min:
            continue
        if cs.add_rel_min is not None and (marker_max <= 0.0 or value < cs.add_rel_min * marker_max):
            continue
        filtered_add[target] = value

    targets = sorted(set(marker_values) | set(filtered_add))
    rows: list[dict[str, object]] = []
    for target in targets:
        marker_value = float(marker_values.get(target, 0.0))
        add_value = float(filtered_add.get(target, 0.0))
        truth_value = float(marker_sample.truth.get(target, 0.0))
        rec: dict[str, object] = {
            "callset": cs.name,
            "dataset": marker_sample.dataset,
            "sample_id": marker_sample.sample_id,
            "target_id": target,
            "truth_abundance": truth_value,
            "is_truth": 1 if truth_value > 0.0 else 0,
            "marker_value": marker_value,
            "add_value": add_value,
            "value_sum": marker_value + add_value,
            "value_max": max(marker_value, add_value),
            "value_marker_first": marker_value if marker_value > 0.0 else add_value,
            "marker_present": 1 if marker_value > 0.0 else 0,
            "add_present": 1 if add_value > 0.0 else 0,
            "both_present": 1 if marker_value > 0.0 and add_value > 0.0 else 0,
        }
        denom = marker_value + add_value
        rec["marker_fraction"] = marker_value / denom if denom > 0.0 else 0.0
        rec["add_fraction"] = add_value / denom if denom > 0.0 else 0.0

        mf = marker_features.get(target)
        af = add_features.get(target) if target in filtered_add else None
        for col in NUMERIC_ROW_COLS:
            rec[f"marker_{col}"] = finite_float(mf.get(f"row_{col}", 0.0)) if mf else 0.0
            rec[f"add_{col}"] = finite_float(af.get(f"row_{col}", 0.0)) if af else 0.0
        rec["marker_reason_active"] = 1 if mf and mf.get("reason") == "active" else 0
        rec["marker_reason_rescue"] = 1 if mf and mf.get("reason") == "rescue" else 0
        rec["marker_reason_fallback"] = 1 if mf and "fallback" in str(mf.get("source", "")) else 0
        rec["add_reason_active"] = 1 if af and af.get("reason") == "active" else 0
        rec["add_reason_rescue"] = 1 if af and af.get("reason") == "rescue" else 0
        rec["dataset_is_mouse"] = 1 if marker_sample.dataset == "mouse_gtdb" else 0
        rec["dataset_is_cami3"] = 1 if marker_sample.dataset == "cami3_ncbi" else 0
        rows.append(rec)

    if not rows:
        return rows
    df = pd.DataFrame(rows)
    marker_total = float(df["marker_value"].sum())
    add_total = float(df["add_value"].sum())
    sum_total = float(df["value_sum"].sum())
    selected_count = float(len(df))
    df["sample_marker_total"] = marker_total
    df["sample_add_total"] = add_total
    df["sample_sum_total"] = sum_total
    df["sample_selected_count"] = selected_count
    df["marker_norm"] = df["marker_value"] / marker_total if marker_total > 0.0 else 0.0
    df["add_norm"] = df["add_value"] / add_total if add_total > 0.0 else 0.0
    df["sum_norm"] = df["value_sum"] / sum_total if sum_total > 0.0 else 0.0
    for col, out in [("marker_value", "marker_rank_frac"), ("add_value", "add_rank_frac"), ("value_sum", "sum_rank_frac")]:
        ranks = df[col].rank(method="average", ascending=False)
        df[out] = ranks / max(float(len(df)), 1.0)
    return df.to_dict("records")


def build_feature_table() -> tuple[pd.DataFrame, dict[tuple[str, int], dict[str, float]]]:
    marker_samples = prev.load_mouse_samples() + prev.load_cami3_samples()
    ctxobj_samples = prev.load_mouse_ctxobj_samples() + prev.load_cami3_ctxobj_samples()
    ctx_by_key = {(s.dataset, s.sample_id): s for s in ctxobj_samples}
    truth_by_key: dict[tuple[str, int], dict[str, float]] = {}
    rows: list[dict[str, object]] = []
    for sample in marker_samples:
        truth_by_key[(sample.dataset, int(sample.sample_id))] = dict(sample.truth)
        ctx_sample = ctx_by_key[(sample.dataset, sample.sample_id)]
        for cs in callsets():
            rows.extend(sample_callset_rows(sample, ctx_sample, cs))
    df = pd.DataFrame(rows)
    for col in df.columns:
        if col not in {"callset", "dataset", "sample_id", "target_id"}:
            df[col] = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return df, truth_by_key


def truth_to_frame(truth_by_key: dict[tuple[str, int], dict[str, float]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (dataset, sample_id), truth in sorted(truth_by_key.items()):
        for target_id, abundance in sorted(truth.items()):
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "target_id": target_id,
                    "truth_abundance": abundance,
                }
            )
    return pd.DataFrame(rows)


def truth_from_frame(rows: pd.DataFrame) -> dict[tuple[str, int], dict[str, float]]:
    truth_by_key: dict[tuple[str, int], dict[str, float]] = {}
    for (dataset, sample_id), sub in rows.groupby(["dataset", "sample_id"]):
        truth_by_key[(str(dataset), int(sample_id))] = dict(
            zip(sub["target_id"].astype(str), pd.to_numeric(sub["truth_abundance"], errors="coerce").fillna(0.0).astype(float))
        )
    return truth_by_key


def normalize_scores(scores: pd.Series, fallback: pd.Series) -> dict[str, float]:
    vals = pd.to_numeric(scores, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(float)
    vals = vals.clip(lower=0.0)
    if float(vals.sum()) <= 0.0:
        vals = pd.to_numeric(fallback, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(float)
        vals = vals.clip(lower=0.0)
    if float(vals.sum()) <= 0.0:
        vals = pd.Series(np.ones(len(vals), dtype=float), index=vals.index)
    vals = vals + EPS
    total = float(vals.sum())
    return {str(idx): float(val / total) for idx, val in vals.items()}


def score_sample_rows(
    rows: pd.DataFrame,
    score_col: str,
    method: str,
    truth_by_key: dict[tuple[str, int], dict[str, float]],
) -> dict[str, object]:
    rows = rows.copy()
    rows.index = rows["target_id"].astype(str)
    pred = normalize_scores(rows[score_col], rows["value_sum"])
    dataset = str(rows["dataset"].iloc[0])
    sample_id = int(rows["sample_id"].iloc[0])
    truth_dict = truth_by_key[(dataset, sample_id)]
    pred_set = set(rows["target_id"].astype(str))
    truth_set = set(truth_dict)
    tp = pred_set & truth_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set
    precision = len(tp) / len(pred_set) if pred_set else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    keys = sorted(truth_dict)
    y_true = [truth_dict[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    pearson = float(pd.Series(y_pred, dtype=float).corr(pd.Series(y_true, dtype=float), method="pearson")) if len(keys) > 1 else float("nan")
    return {
        "method": method,
        "callset": str(rows["callset"].iloc[0]),
        "dataset": dataset,
        "sample_id": sample_id,
        "truth_taxa": len(truth_set),
        "pred_taxa": len(pred_set),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) * 100.0,
        "pearson": pearson,
        "pred_sum_on_truth": sum(y_pred),
    }


def feature_columns(df: pd.DataFrame, feature_set: str) -> list[str]:
    base = list(VALUE_FEATURES)
    reason = [
        "marker_reason_active",
        "marker_reason_rescue",
        "marker_reason_fallback",
        "add_reason_active",
        "add_reason_rescue",
    ]
    if feature_set == "value":
        cols = base
    elif feature_set == "nodataset":
        cols = base + reason
        cols += [c for c in df.columns if c.startswith("marker_") or c.startswith("add_")]
        cols = [c for c in cols if c not in {"marker_present", "add_present"} or c in base]
        cols = [c for c in cols if not c.startswith("dataset_")]
    elif feature_set == "dataset":
        cols = base + reason + ["dataset_is_mouse", "dataset_is_cami3"]
        cols += [c for c in df.columns if c.startswith("marker_") or c.startswith("add_")]
    else:
        raise ValueError(feature_set)
    drop = {"marker_value", "add_value", "marker_fraction", "add_fraction"} if feature_set == "value_light" else set()
    ordered: list[str] = []
    for col in cols:
        if col in df.columns and col not in drop and col not in ordered:
            ordered.append(col)
    return ordered


def make_estimator(name: str):
    if name == "ridge":
        return Ridge(alpha=1.0, random_state=1), True
    if name == "ridge10":
        return Ridge(alpha=10.0, random_state=1), True
    if name == "huber":
        return HuberRegressor(alpha=0.001, epsilon=1.35, max_iter=500), True
    if name == "rf":
        return RandomForestRegressor(n_estimators=160, min_samples_leaf=4, random_state=1, n_jobs=1), False
    if name == "et":
        return ExtraTreesRegressor(n_estimators=200, min_samples_leaf=3, random_state=1, n_jobs=1), False
    if name == "hgb":
        return HistGradientBoostingRegressor(max_iter=180, learning_rate=0.04, l2_regularization=0.05, random_state=1), False
    if name == "tweedie":
        return TweedieRegressor(power=1.0, alpha=0.001, link="log", max_iter=1000), True
    raise ValueError(name)


def transform_target(y: pd.Series, mode: str) -> np.ndarray:
    arr = y.to_numpy(dtype=float)
    if mode == "raw":
        return arr
    if mode == "sqrt":
        return np.sqrt(np.maximum(arr, 0.0))
    if mode == "log100":
        return np.log1p(np.maximum(arr, 0.0) * 100.0)
    if mode == "log1000":
        return np.log1p(np.maximum(arr, 0.0) * 1000.0)
    raise ValueError(mode)


def inverse_target(pred: np.ndarray, mode: str) -> np.ndarray:
    pred = np.asarray(pred, dtype=float)
    if mode == "raw":
        return pred
    if mode == "sqrt":
        return pred**2
    if mode == "log100":
        return np.expm1(pred) / 100.0
    if mode == "log1000":
        return np.expm1(pred) / 1000.0
    raise ValueError(mode)


def sample_weights(y: pd.Series, mode: str) -> np.ndarray | None:
    arr = y.to_numpy(dtype=float)
    if mode == "none":
        return None
    if mode == "positive5":
        return np.where(arr > 0.0, 5.0, 1.0)
    if mode == "abundance":
        return 1.0 + 200.0 * arr
    raise ValueError(mode)


def model_specs(max_models: int | None = None) -> list[ModelSpec]:
    specs: list[ModelSpec] = []
    for feature_set in ["value", "nodataset", "dataset"]:
        for target_mode in ["raw", "sqrt", "log100", "log1000"]:
            for weight_mode in ["none", "positive5", "abundance"]:
                for est in ["ridge", "ridge10", "huber", "tweedie", "hgb", "rf", "et"]:
                    if est == "tweedie" and target_mode != "raw":
                        continue
                    if est in {"rf", "et"} and target_mode == "log1000":
                        continue
                    specs.append(ModelSpec(f"{est}_{feature_set}_{target_mode}_{weight_mode}", est, target_mode, feature_set, weight_mode))
    if max_models and len(specs) > max_models:
        ordered = sorted(specs, key=lambda s: s.name)
        idx = np.linspace(0, len(ordered) - 1, max_models, dtype=int)
        specs = [ordered[int(i)] for i in idx]
    return specs


def fit_predict(train: pd.DataFrame, test: pd.DataFrame, spec: ModelSpec) -> np.ndarray:
    cols = feature_columns(train, spec.feature_set)
    x_train = train[cols].to_numpy(dtype=float)
    x_test = test[cols].to_numpy(dtype=float)
    y_train = transform_target(train["truth_abundance"], spec.target_mode)
    weights = sample_weights(train["truth_abundance"], spec.weight_mode)
    model, scale = make_estimator(spec.estimator)
    if scale:
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x_train)
        x_test = scaler.transform(x_test)
    try:
        model.fit(x_train, y_train, sample_weight=weights)
    except TypeError:
        model.fit(x_train, y_train)
    pred = inverse_target(model.predict(x_test), spec.target_mode)
    pred = np.nan_to_num(pred, nan=0.0, posinf=0.0, neginf=0.0)
    return np.clip(pred, 0.0, None)


def summarize(rows: pd.DataFrame, method: str, callset: str, extra: dict[str, object] | None = None) -> dict[str, object]:
    rec: dict[str, object] = {"method": method, "callset": callset}
    if extra:
        rec.update(extra)
    for prefix, sub in [
        ("all", rows),
        ("mouse", rows.loc[rows["dataset"] == "mouse_gtdb"]),
        ("cami3", rows.loc[rows["dataset"] == "cami3_ncbi"]),
    ]:
        rec[f"{prefix}_F1"] = float(sub["F1"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_L1"] = float(sub["l1_pct_points"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_Pearson"] = float(sub["pearson"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_TP"] = float(sub["TP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FP"] = float(sub["FP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FN"] = float(sub["FN"].mean()) if len(sub) else float("nan")
    return rec


def score_raw_baselines(
    df: pd.DataFrame,
    callset: str,
    truth_by_key: dict[tuple[str, int], dict[str, float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for score_col in ["value_marker_first", "value_sum", "value_max", "sum_norm"]:
        per_sample = []
        for (_, _), sub in df.groupby(["dataset", "sample_id"]):
            rows = sub.copy()
            rows["score"] = rows[score_col]
            per_sample.append(score_sample_rows(rows, "score", f"raw_{score_col}", truth_by_key))
        sample_df = pd.DataFrame(per_sample)
        sample_rows.extend(per_sample)
        summary_rows.append(summarize(sample_df, f"raw_{score_col}", callset, {"eval_mode": "raw"}))
    return sample_rows, summary_rows


def score_loso_models(
    df: pd.DataFrame,
    callset: str,
    specs: list[ModelSpec],
    truth_by_key: dict[tuple[str, int], dict[str, float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_keys = sorted(df[["dataset", "sample_id"]].drop_duplicates().itertuples(index=False, name=None))
    all_sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for i, spec in enumerate(specs, start=1):
        per_sample: list[dict[str, object]] = []
        for dataset, sample_id in sample_keys:
            test_mask = (df["dataset"] == dataset) & (df["sample_id"] == sample_id)
            train = df.loc[~test_mask].copy()
            test = df.loc[test_mask].copy()
            pred = fit_predict(train, test, spec)
            rows = test.copy()
            rows["score"] = pred
            per_sample.append(score_sample_rows(rows, "score", spec.name, truth_by_key))
        sample_df = pd.DataFrame(per_sample)
        all_sample_rows.extend(per_sample)
        summary_rows.append(
            summarize(
                sample_df,
                spec.name,
                callset,
                {
                    "eval_mode": "leave_one_sample_out",
                    "estimator": spec.estimator,
                    "target_mode": spec.target_mode,
                    "feature_set": spec.feature_set,
                    "weight_mode": spec.weight_mode,
                },
            )
        )
        if i % 25 == 0:
            print(f"{callset}: scored {i}/{len(specs)} model specs", file=sys.stderr)
    return all_sample_rows, summary_rows


def load_sylph() -> pd.DataFrame:
    return pd.read_csv(ROOT / "research/experiments/2026-06-24_addback_abundance_normalization_search/results/sylph_baseline.tsv", sep="\t")


def write_final_comparison(summary: pd.DataFrame) -> None:
    sylph = load_sylph()
    rows = []
    s_all = sylph.loc[sylph["dataset"] == "all"].iloc[0]
    s_mouse = sylph.loc[sylph["dataset"] == "mouse_gtdb"].iloc[0]
    s_cami3 = sylph.loc[sylph["dataset"] == "cami3_ncbi"].iloc[0]
    rows.append(
        {
            "method": "Sylph",
            "callset": "baseline",
            "F1": s_all["F1"],
            "L1": s_all["l1_pct_points"],
            "mouse_F1": s_mouse["F1"],
            "mouse_L1": s_mouse["l1_pct_points"],
            "cami3_F1": s_cami3["F1"],
            "cami3_L1": s_cami3["l1_pct_points"],
            "note": "baseline",
        }
    )

    loso = summary.loc[summary["eval_mode"] == "leave_one_sample_out"].copy()
    raw = summary.loc[summary["eval_mode"] == "raw"].copy()
    candidates = [
        ("best_loso_L1", loso.sort_values(["all_L1", "all_F1"], ascending=[True, False]).head(1)),
        ("best_loso_F1_then_L1", loso.sort_values(["all_F1", "all_L1"], ascending=[False, True]).head(1)),
        ("best_raw_L1", raw.sort_values(["all_L1", "all_F1"], ascending=[True, False]).head(1)),
    ]
    for label, sub in candidates:
        if sub.empty:
            continue
        row = sub.iloc[0]
        rows.append(
            {
                "method": label,
                "callset": row["callset"],
                "F1": row["all_F1"],
                "L1": row["all_L1"],
                "mouse_F1": row["mouse_F1"],
                "mouse_L1": row["mouse_L1"],
                "cami3_F1": row["cami3_F1"],
                "cami3_L1": row["cami3_L1"],
                "note": row["method"],
            }
        )
    pd.DataFrame(rows).to_csv(OUT / "final_comparison.tsv", sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-models", type=int, default=90)
    parser.add_argument("--save-feature-table", action="store_true")
    parser.add_argument("--feature-table", type=Path)
    parser.add_argument("--truth-table", type=Path)
    args = parser.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    if args.feature_table and args.truth_table:
        features = pd.read_csv(args.feature_table, sep="\t")
        truth_by_key = truth_from_frame(pd.read_csv(args.truth_table, sep="\t"))
        for col in features.columns:
            if col not in {"callset", "dataset", "sample_id", "target_id"}:
                features[col] = pd.to_numeric(features[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    else:
        features, truth_by_key = build_feature_table()
    if args.save_feature_table:
        features.to_csv(OUT / "selected_call_features.tsv", sep="\t", index=False)
        truth_to_frame(truth_by_key).to_csv(OUT / "truth_profiles.tsv", sep="\t", index=False)

    specs = model_specs(args.max_models)
    sample_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for callset, sub in features.groupby("callset"):
        raw_samples, raw_summary = score_raw_baselines(sub.copy(), callset, truth_by_key)
        sample_rows.extend(raw_samples)
        summary_rows.extend(raw_summary)
        model_samples, model_summary = score_loso_models(sub.copy(), callset, specs, truth_by_key)
        sample_rows.extend(model_samples)
        summary_rows.extend(model_summary)

    sample_df = pd.DataFrame(sample_rows)
    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values(["all_L1", "all_F1"], ascending=[True, False])
    summary_df.to_csv(OUT / "calibration_rule_summary.tsv", sep="\t", index=False)
    summary_df.head(200).to_csv(OUT / "calibration_top200_by_l1.tsv", sep="\t", index=False)
    summary_df.sort_values(["all_F1", "all_L1"], ascending=[False, True]).head(200).to_csv(
        OUT / "calibration_top200_by_f1.tsv", sep="\t", index=False
    )
    sample_df.to_csv(OUT / "calibration_sample_metrics.tsv", sep="\t", index=False)
    write_final_comparison(summary_df)
    print(pd.read_csv(OUT / "final_comparison.tsv", sep="\t").to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
