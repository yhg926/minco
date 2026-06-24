#!/usr/bin/env python3
"""Multi-sample minco call-model calibration.

This script keeps direct minco row-level scoring separate from experimental
taxid-level model scoring. The direct baseline reproduces the CLI recipe:
XnY_ctx >= 10, ANI >= 0.94, Real_min_align_fraction >= 0.05.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
HELPER_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(HELPER_DIR))

from analyze_readwise_corrections import (  # noqa: E402
    load_minco,
    load_sylph,
    parse_gold_profile,
    parse_species_taxmap,
    score_taxids,
)


NUMERIC_COLS = [
    "ANI",
    "Distance",
    "XnY_ctx",
    "Qry_align_fraction",
    "Ref_align_fraction",
    "N_diff_obj",
    "N_diff_obj_section",
    "N_mut2_ctx",
    "Real_Qry_align_fraction",
    "Real_Ref_align_fraction",
    "Real_min_align_fraction",
    "Reads_with_ctx_match",
    "Read_match_fraction",
    "Unique_ref_ctx_hit",
    "Blocks_with_ctx_match",
    "Block_match_fraction",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_variance",
    "Ref_depth_cv",
    "Ref_zero_fraction",
    "Relative_abundance_depth",
    "Normalized_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
]

MODEL_COLS = [
    "u_ANI_max",
    "u_XnY_ctx_max",
    "u_Real_min_align_fraction_max",
    "u_Ref_breadth_max",
    "u_Ref_mean_depth_max",
    "u_Ref_hit_mean_depth_max",
    "u_Ref_depth_cv_min",
    "u_Normalized_abundance_depth_max",
    "u_Ref_zip_af_max",
    "u_Ref_zip_aaf_ani_max",
    "u_rows",
    "u_direct_call",
    "u_relaxed_call",
    "s_ANI_max",
    "s_XnY_ctx_max",
    "s_Real_min_align_fraction_max",
    "s_Ref_breadth_max",
    "s_Ref_mean_depth_max",
    "s_Ref_hit_mean_depth_max",
    "s_Ref_depth_cv_min",
    "s_Normalized_abundance_depth_max",
    "s_Ref_zip_af_max",
    "s_Ref_zip_aaf_ani_max",
    "s_rows",
    "s_direct_call",
    "s_relaxed_call",
    "split_unique_xny_ratio",
    "split_unique_breadth_ratio",
    "split_unique_ani_delta",
]


def _to_float_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def gold_scope(gold: pd.DataFrame, scope: str) -> pd.DataFrame:
    if scope == "all":
        return gold
    if scope == "bacteria":
        return gold.loc[gold["taxpathsn"].str.contains("Bacteria", regex=False)].copy()
    if scope == "virus":
        return gold.loc[gold["taxpathsn"].str.contains("Viruses", regex=False)].copy()
    raise ValueError(f"unsupported scope: {scope}")


def taxid_scope_map(taxmap: Mapping[str, Mapping[str, str]]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in taxmap.values():
        taxid = str(rec.get("taxid", ""))
        path = str(rec.get("taxpathsn", ""))
        if not taxid:
            continue
        if "Viruses" in path:
            out.setdefault(taxid, "virus")
        elif "Bacteria" in path:
            out.setdefault(taxid, "bacteria")
        elif "Archaea" in path:
            out.setdefault(taxid, "archaea")
        else:
            out.setdefault(taxid, "other")
    return out


def filter_scope_taxids(taxids: Iterable[str], scope: str, scope_by_taxid: Mapping[str, str]) -> Set[str]:
    values = {str(t) for t in taxids if str(t)}
    if scope == "all":
        return values
    return {t for t in values if scope_by_taxid.get(t, "other") == scope}


def score_prediction(pred_taxids: Iterable[str], gold_taxids: Set[str]) -> Dict[str, object]:
    pred = {str(t) for t in pred_taxids if str(t)}
    out: Dict[str, object] = {"pred_taxa": len(pred)}
    out.update(score_taxids(pred, gold_taxids))
    out["FP_plus_FN"] = out["FP"] + out["FN"]
    return out


def load_numeric_minco(path: Path, taxmap: Mapping[str, Mapping[str, str]]) -> pd.DataFrame:
    rows = load_minco(path, taxmap, 11.0)
    rows = rows.loc[rows["taxid"].astype(bool)].copy()
    for col in NUMERIC_COLS:
        rows[col] = _to_float_series(rows, col)
    return rows


def direct_mask(rows: pd.DataFrame) -> pd.Series:
    return (
        (_to_float_series(rows, "XnY_ctx") >= 10.0)
        & (_to_float_series(rows, "ANI") >= 0.94)
        & (_to_float_series(rows, "Real_min_align_fraction") >= 0.05)
    )


def relaxed_mask(rows: pd.DataFrame) -> pd.Series:
    return (
        (_to_float_series(rows, "XnY_ctx") >= 1.0)
        & (_to_float_series(rows, "ANI") >= 0.90)
        & (_to_float_series(rows, "Ref_breadth") >= 0.005)
    )


def _agg_one(rows: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if rows.empty:
        return pd.DataFrame(columns=["taxid"])
    work = rows.copy()
    work["_direct_call"] = direct_mask(work).astype(float)
    work["_relaxed_call"] = relaxed_mask(work).astype(float)
    grouped = work.groupby("taxid", as_index=False)
    out = pd.DataFrame({"taxid": sorted(work["taxid"].astype(str).unique())})
    out = out.merge(grouped.size().rename(columns={"size": f"{prefix}_rows"}), on="taxid", how="left")
    for flag in ["_direct_call", "_relaxed_call"]:
        flag_df = grouped[flag].max().rename(columns={flag: f"{prefix}{flag}"})
        out = out.merge(flag_df, on="taxid", how="left")
    for col in NUMERIC_COLS:
        if col not in work.columns:
            continue
        max_df = grouped[col].max().rename(columns={col: f"{prefix}_{col}_max"})
        out = out.merge(max_df, on="taxid", how="left")
        if col == "Ref_depth_cv":
            min_df = grouped[col].min().rename(columns={col: f"{prefix}_{col}_min"})
            out = out.merge(min_df, on="taxid", how="left")
    return out


def joined_features(unique_rows: pd.DataFrame, split_rows: pd.DataFrame) -> pd.DataFrame:
    u = _agg_one(unique_rows, "u")
    s = _agg_one(split_rows, "s")
    joined = u.merge(s, on="taxid", how="outer")
    for col in joined.columns:
        if col != "taxid":
            joined[col] = pd.to_numeric(joined[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    for col in MODEL_COLS:
        if col not in joined.columns:
            joined[col] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        joined["split_unique_xny_ratio"] = joined["s_XnY_ctx_max"] / joined["u_XnY_ctx_max"].replace(0.0, np.nan)
        joined["split_unique_breadth_ratio"] = joined["s_Ref_breadth_max"] / joined["u_Ref_breadth_max"].replace(0.0, np.nan)
        joined["split_unique_ani_delta"] = joined["s_ANI_max"] - joined["u_ANI_max"]
    for col in ["split_unique_xny_ratio", "split_unique_breadth_ratio", "split_unique_ani_delta"]:
        joined[col] = pd.to_numeric(joined[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return joined


def score_sylph(path: Path, taxmap: Mapping[str, Mapping[str, str]], scope: str, scope_by_taxid: Mapping[str, str], gold_taxids: Set[str]) -> Dict[str, object]:
    if not path or not path.exists():
        return {"pred_taxa": 0, **score_taxids(set(), gold_taxids), "FP_plus_FN": len(gold_taxids)}
    rows = load_sylph(path, taxmap)
    pred = filter_scope_taxids(rows.loc[rows["taxid"].astype(bool), "taxid"].astype(str), scope, scope_by_taxid)
    return score_prediction(pred, gold_taxids)


def tune_threshold(prob: np.ndarray, y: np.ndarray, sample: np.ndarray, gold_counts: Mapping[str, int]) -> float:
    best_t = 0.5
    best_key: Tuple[float, float, float] = (-1.0, -1e9, -1e9)
    for t in np.linspace(0.05, 0.95, 91):
        f1_values: List[float] = []
        fpfn_values: List[float] = []
        for sample_key in sorted(set(sample)):
            mask = sample == sample_key
            pred = prob[mask] >= t
            tp = int(np.sum(pred & (y[mask] == 1)))
            fp = int(np.sum(pred & (y[mask] == 0)))
            fn = int(gold_counts[sample_key] - tp)
            precision = tp / (tp + fp) if tp + fp else 0.0
            recall = tp / (tp + fn) if tp + fn else 0.0
            f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
            f1_values.append(f1)
            fpfn_values.append(fp + fn)
        key = (float(np.mean(f1_values)), -float(np.mean(fpfn_values)), -abs(t - 0.5))
        if key > best_key:
            best_key = key
            best_t = float(t)
    return best_t


def model_suite(random_state: int = 17):
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"sklearn is required for model calibration: {exc}")

    return {
        "logreg": make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            LogisticRegression(max_iter=1000, class_weight="balanced", random_state=random_state),
        ),
        "rf": make_pipeline(
            SimpleImputer(strategy="median"),
            RandomForestClassifier(
                n_estimators=300,
                min_samples_leaf=3,
                class_weight="balanced_subsample",
                n_jobs=8,
                random_state=random_state,
            ),
        ),
        "hgb": make_pipeline(
            SimpleImputer(strategy="median"),
            HistGradientBoostingClassifier(max_iter=220, learning_rate=0.04, max_leaf_nodes=15, random_state=random_state),
        ),
    }


def predict_probability(model, x: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x)[:, 1]
    pred = model.predict(x)
    return np.asarray(pred, dtype=float)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--taxmap", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    taxmap = parse_species_taxmap(args.taxmap)
    scope_by_taxid = taxid_scope_map(taxmap)

    sample_frames: List[pd.DataFrame] = []
    summary_records: List[Dict[str, object]] = []
    gold_counts: Dict[str, int] = {}

    for _, rec in manifest.iterrows():
        sample_key = str(rec["sample_key"])
        dataset = str(rec["dataset"])
        scope = str(rec["scope"])
        gold = gold_scope(parse_gold_profile(Path(rec["gold_profile"]), str(rec["gold_sample_id"])), scope)
        gold_taxids = set(gold["taxid"].astype(str))
        gold_counts[sample_key] = len(gold_taxids)

        unique_rows = load_numeric_minco(Path(rec["unique_unfiltered"]), taxmap)
        split_rows = load_numeric_minco(Path(rec["split_unfiltered"]), taxmap)
        unique_pred = filter_scope_taxids(unique_rows.loc[direct_mask(unique_rows), "taxid"].astype(str), scope, scope_by_taxid)
        split_pred = filter_scope_taxids(split_rows.loc[direct_mask(split_rows), "taxid"].astype(str), scope, scope_by_taxid)
        unique_relaxed_pred = filter_scope_taxids(unique_rows.loc[relaxed_mask(unique_rows), "taxid"].astype(str), scope, scope_by_taxid)

        for label, pred in [
            ("minco_unique_direct", unique_pred),
            ("minco_split_direct", split_pred),
            ("minco_unique_relaxed", unique_relaxed_pred),
        ]:
            row: Dict[str, object] = {"sample_key": sample_key, "dataset": dataset, "scope": scope, "method": label, "gold_taxa": len(gold_taxids)}
            row.update(score_prediction(pred, gold_taxids))
            summary_records.append(row)

        row = {"sample_key": sample_key, "dataset": dataset, "scope": scope, "method": "sylph_default", "gold_taxa": len(gold_taxids)}
        row.update(score_sylph(Path(rec["sylph_profile"]), taxmap, scope, scope_by_taxid, gold_taxids))
        summary_records.append(row)

        features = joined_features(unique_rows, split_rows)
        if scope != "all":
            features = features.loc[features["taxid"].astype(str).map(lambda t: scope_by_taxid.get(t, "other") == scope)].copy()
        features["sample_key"] = sample_key
        features["dataset"] = dataset
        features["scope"] = scope
        features["label"] = features["taxid"].astype(str).isin(gold_taxids).astype(int)
        sample_frames.append(features)
        features.to_csv(args.outdir / f"{sample_key}.joined_features.tsv", sep="\t", index=False)

    all_features = pd.concat(sample_frames, ignore_index=True, sort=False).fillna(0.0)
    all_features.to_csv(args.outdir / "all_samples.joined_features.tsv", sep="\t", index=False)

    models = model_suite()
    x_all = all_features[MODEL_COLS].to_numpy(dtype=float)
    y_all = all_features["label"].to_numpy(dtype=int)
    sample_all = all_features["sample_key"].astype(str).to_numpy()
    taxid_all = all_features["taxid"].astype(str).to_numpy()
    dataset_all = all_features["dataset"].astype(str).to_numpy()
    scope_all = all_features["scope"].astype(str).to_numpy()

    prediction_records: List[Dict[str, object]] = []
    for heldout in sorted(set(sample_all)):
        train_idx = sample_all != heldout
        test_idx = sample_all == heldout
        train_samples = sample_all[train_idx]
        sample_key = str(heldout)
        gold_taxa = int(gold_counts[sample_key])
        gold_subset = set(all_features.loc[all_features["sample_key"] == sample_key].loc[all_features["label"] == 1, "taxid"].astype(str))

        def record_probability_model(model_name: str, train_prob: np.ndarray, test_prob: np.ndarray) -> None:
            threshold = tune_threshold(train_prob, y_all[train_idx], train_samples, gold_counts)
            pred_taxids = set(taxid_all[test_idx][test_prob >= threshold])
            row = {
                "sample_key": sample_key,
                "dataset": str(dataset_all[test_idx][0]),
                "scope": str(scope_all[test_idx][0]),
                "method": f"loso_{model_name}",
                "gold_taxa": gold_taxa,
                "threshold": threshold,
                "train_samples": ",".join(sorted(set(train_samples))),
            }
            row.update(score_prediction(pred_taxids, gold_subset))
            # FNs must include gold taxa absent from minco candidates.
            row["FN"] = gold_taxa - row["TP"]
            precision = row["TP"] / (row["TP"] + row["FP"]) if row["TP"] + row["FP"] else 0.0
            recall = row["TP"] / (row["TP"] + row["FN"]) if row["TP"] + row["FN"] else 0.0
            row["precision"] = precision
            row["recall"] = recall
            row["F1"] = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
            row["FP_plus_FN"] = row["FP"] + row["FN"]
            summary_records.append(row)
            for taxid, prob in zip(taxid_all[test_idx], test_prob):
                prediction_records.append(
                    {
                        "sample_key": sample_key,
                        "model": model_name,
                        "threshold": threshold,
                        "taxid": taxid,
                        "probability": float(prob),
                        "predicted": int(prob >= threshold),
                        "label": int(taxid in gold_subset),
                    }
                )

        train_prob_by_model: Dict[str, np.ndarray] = {}
        test_prob_by_model: Dict[str, np.ndarray] = {}
        for model_name, model in models.items():
            model.fit(x_all[train_idx], y_all[train_idx])
            train_prob = predict_probability(model, x_all[train_idx])
            test_prob = predict_probability(model, x_all[test_idx])
            train_prob_by_model[model_name] = train_prob
            test_prob_by_model[model_name] = test_prob
            record_probability_model(model_name, train_prob, test_prob)

        if "rf" in train_prob_by_model and "hgb" in train_prob_by_model:
            record_probability_model(
                "rf_hgb_avg",
                0.5 * (train_prob_by_model["rf"] + train_prob_by_model["hgb"]),
                0.5 * (test_prob_by_model["rf"] + test_prob_by_model["hgb"]),
            )

    summary = pd.DataFrame(summary_records)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    pd.DataFrame(prediction_records).to_csv(args.outdir / "loso_predictions.tsv", sep="\t", index=False)

    mean_rows = []
    for method, group in summary.groupby("method"):
        mean_rows.append(
            {
                "method": method,
                "samples": int(len(group)),
                "mean_precision": float(pd.to_numeric(group["precision"], errors="coerce").mean()),
                "mean_recall": float(pd.to_numeric(group["recall"], errors="coerce").mean()),
                "mean_F1": float(pd.to_numeric(group["F1"], errors="coerce").mean()),
                "mean_FP": float(pd.to_numeric(group["FP"], errors="coerce").mean()),
                "mean_FN": float(pd.to_numeric(group["FN"], errors="coerce").mean()),
            }
        )
    mean_summary = pd.DataFrame(mean_rows).sort_values("mean_F1", ascending=False)
    mean_summary.to_csv(args.outdir / "mean_summary.tsv", sep="\t", index=False)

    rf_model = models.get("rf")
    if rf_model is not None:
        rf_model.fit(x_all, y_all)
        final_estimator = rf_model.steps[-1][1] if hasattr(rf_model, "steps") else rf_model
        if hasattr(final_estimator, "feature_importances_"):
            pd.DataFrame(
                {
                    "feature": MODEL_COLS,
                    "importance": final_estimator.feature_importances_,
                }
            ).sort_values("importance", ascending=False).to_csv(args.outdir / "rf_feature_importance.tsv", sep="\t", index=False)

    print(mean_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
