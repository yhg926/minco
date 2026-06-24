#!/usr/bin/env python3
"""Search deployable MinCO call rules across cached multi-domain features.

The priority is presence/absence F1. Abundance and ANI are evaluated in
separate scripts/notes because the cached S1000 joined features do not contain
the newer markerdb robust effective-depth and adjusted-ANI fields.
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
from typing import Iterable

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


def f1_from_counts(tp: int, fp: int, gold: int) -> tuple[int, int, float, float, float]:
    fn = max(gold - tp, 0)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return fn, precision, recall, f1, float(fp + fn)


def score_mask(
    df: pd.DataFrame,
    mask: np.ndarray,
    gold_counts: dict[str, int],
    *,
    method: str,
    detail_rows: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    work = df.loc[mask, ["sample_key", "dataset", "label"]].copy()
    sample_keys = sorted(df["sample_key"].astype(str).unique())
    rows = []
    for sample_key in sample_keys:
        sub = work.loc[work["sample_key"] == sample_key]
        full = df.loc[df["sample_key"] == sample_key]
        gold = int(gold_counts.get(sample_key, int(full["label"].sum())))
        tp = int((sub["label"] > 0).sum())
        pred = int(len(sub))
        fp = pred - tp
        fn, precision, recall, f1, fpfn = f1_from_counts(tp, fp, gold)
        dataset = str(full["dataset"].iloc[0]) if len(full) else ""
        row = {
            "method": method,
            "sample_key": sample_key,
            "dataset": dataset,
            "gold_taxa": gold,
            "pred_taxa": pred,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "F1": f1,
            "FP_plus_FN": fpfn,
        }
        rows.append(row)
        if detail_rows is not None:
            detail_rows.append(row)
    sample_df = pd.DataFrame(rows)
    by_dataset = sample_df.groupby("dataset", as_index=False).agg(
        dataset_F1=("F1", "mean"),
        dataset_precision=("precision", "mean"),
        dataset_recall=("recall", "mean"),
        dataset_FP=("FP", "mean"),
        dataset_FN=("FN", "mean"),
    )
    ds = {f"{r.dataset}_F1": r.dataset_F1 for r in by_dataset.itertuples(index=False)}
    return {
        "method": method,
        "samples": int(len(sample_df)),
        "mean_precision": float(sample_df["precision"].mean()),
        "mean_recall": float(sample_df["recall"].mean()),
        "mean_F1": float(sample_df["F1"].mean()),
        "mean_FP": float(sample_df["FP"].mean()),
        "mean_FN": float(sample_df["FN"].mean()),
        "mean_FP_plus_FN": float(sample_df["FP_plus_FN"].mean()),
        "min_dataset_F1": float(by_dataset["dataset_F1"].min()),
        **ds,
    }


def prepare_fast_context(df: pd.DataFrame, gold_counts: dict[str, int]) -> dict[str, object]:
    sample_order = sorted(df["sample_key"].astype(str).unique())
    sample_to_id = {sample: i for i, sample in enumerate(sample_order)}
    sample_ids = df["sample_key"].astype(str).map(sample_to_id).to_numpy(dtype=int)
    labels = (numeric(df, "label").to_numpy() > 0).astype(np.int8)
    datasets = {
        sample: str(df.loc[df["sample_key"].astype(str) == sample, "dataset"].iloc[0])
        for sample in sample_order
    }
    gold = np.array(
        [
            int(gold_counts.get(sample, int(labels[sample_ids == i].sum())))
            for i, sample in enumerate(sample_order)
        ],
        dtype=float,
    )
    return {
        "sample_order": sample_order,
        "sample_ids": sample_ids,
        "labels": labels,
        "datasets": datasets,
        "gold": gold,
    }


def score_mask_fast(ctx: dict[str, object], mask: np.ndarray, *, method: str) -> dict[str, object]:
    sample_order: list[str] = ctx["sample_order"]  # type: ignore[assignment]
    sample_ids: np.ndarray = ctx["sample_ids"]  # type: ignore[assignment]
    labels: np.ndarray = ctx["labels"]  # type: ignore[assignment]
    datasets: dict[str, str] = ctx["datasets"]  # type: ignore[assignment]
    gold: np.ndarray = ctx["gold"]  # type: ignore[assignment]
    n = len(sample_order)

    mask = np.asarray(mask, dtype=bool)
    pred = np.bincount(sample_ids[mask], minlength=n).astype(float)
    tp = np.bincount(sample_ids[mask & (labels > 0)], minlength=n).astype(float)
    fp = pred - tp
    fn = np.maximum(gold - tp, 0.0)
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp), where=(tp + fn) > 0)
    f1 = np.divide(
        2.0 * precision * recall,
        precision + recall,
        out=np.zeros_like(precision),
        where=(precision + recall) > 0,
    )

    sample_df = pd.DataFrame(
        {
            "sample_key": sample_order,
            "dataset": [datasets[s] for s in sample_order],
            "precision": precision,
            "recall": recall,
            "F1": f1,
            "FP": fp,
            "FN": fn,
            "FP_plus_FN": fp + fn,
        }
    )
    by_dataset = sample_df.groupby("dataset", as_index=False).agg(
        dataset_F1=("F1", "mean"),
        dataset_precision=("precision", "mean"),
        dataset_recall=("recall", "mean"),
        dataset_FP=("FP", "mean"),
        dataset_FN=("FN", "mean"),
    )
    ds = {f"{r.dataset}_F1": r.dataset_F1 for r in by_dataset.itertuples(index=False)}
    return {
        "method": method,
        "samples": int(n),
        "mean_precision": float(precision.mean()),
        "mean_recall": float(recall.mean()),
        "mean_F1": float(f1.mean()),
        "mean_FP": float(fp.mean()),
        "mean_FN": float(fn.mean()),
        "mean_FP_plus_FN": float((fp + fn).mean()),
        "min_dataset_F1": float(by_dataset["dataset_F1"].min()),
        **ds,
    }


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
    out = []
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
            "source": "baseline",
        }
        for row in by_dataset.itertuples(index=False):
            rec[f"{row.dataset}_F1"] = row.dataset_F1
        out.append(rec)
    return pd.DataFrame(out)


def candidate_masks(df: pd.DataFrame) -> Iterable[tuple[str, np.ndarray]]:
    u_direct = numeric(df, "u_direct_call").to_numpy() > 0
    s_direct = numeric(df, "s_direct_call").to_numpy() > 0
    u_relaxed = numeric(df, "u_relaxed_call").to_numpy() > 0
    s_relaxed = numeric(df, "s_relaxed_call").to_numpy() > 0
    either_direct = u_direct | s_direct
    either_relaxed = u_relaxed | s_relaxed

    max_ani = np.maximum(numeric(df, "u_ANI_max"), numeric(df, "s_ANI_max")).to_numpy()
    max_xny = np.maximum(numeric(df, "u_XnY_ctx_max"), numeric(df, "s_XnY_ctx_max")).to_numpy()
    max_minaf = np.maximum(
        numeric(df, "u_Real_min_align_fraction_max"),
        numeric(df, "s_Real_min_align_fraction_max"),
    ).to_numpy()
    max_breadth = np.maximum(numeric(df, "u_Ref_breadth_max"), numeric(df, "s_Ref_breadth_max")).to_numpy()
    min_cv = np.minimum(
        numeric(df, "u_Ref_depth_cv_min").replace(0, np.nan),
        numeric(df, "s_Ref_depth_cv_min").replace(0, np.nan),
    ).fillna(0.0).to_numpy()
    max_norm_ab = np.maximum(
        numeric(df, "u_Normalized_abundance_depth_max"),
        numeric(df, "s_Normalized_abundance_depth_max"),
    ).to_numpy()
    max_zip_ani = np.maximum(numeric(df, "u_Ref_zip_aaf_ani_max"), numeric(df, "s_Ref_zip_aaf_ani_max")).to_numpy()
    max_zip_af = np.maximum(numeric(df, "u_Ref_zip_af_max"), numeric(df, "s_Ref_zip_af_max")).to_numpy()

    yield "rule_u_direct", u_direct
    yield "rule_s_direct", s_direct
    yield "rule_u_or_s_direct", either_direct
    yield "rule_u_direct_or_s_direct_lowcv20", either_direct & (min_cv <= 20.0)
    yield "rule_u_or_s_relaxed", either_relaxed

    for ani, xny, minaf in itertools.product(
        [0.90, 0.92, 0.94, 0.96, 0.98, 0.99, 0.995, 0.999],
        [1, 2, 3, 5, 10, 20, 50, 100],
        [0.001, 0.002, 0.005, 0.01, 0.02, 0.05],
    ):
        mask = (max_xny >= xny) & (max_ani >= ani) & (max_minaf >= minaf)
        yield f"rule_any_xny{xny}_ani{ani:g}_minaf{minaf:g}", mask

    for ani, xny, breadth in itertools.product(
        [0.90, 0.92, 0.94, 0.96, 0.98, 0.99, 0.995, 0.999],
        [1, 2, 3, 5, 10, 20, 50, 100],
        [0.001, 0.002, 0.005, 0.01, 0.02, 0.05],
    ):
        mask = (max_xny >= xny) & (max_ani >= ani) & (max_breadth >= breadth)
        yield f"rule_any_xny{xny}_ani{ani:g}_breadth{breadth:g}", mask

    # Conservative direct plus guarded rescue. These are closer to a deployable
    # universal rule than a free-form threshold surface.
    for ani, xny, minaf, cvmax in itertools.product(
        [0.92, 0.94, 0.96, 0.98, 0.99],
        [1, 2, 3, 5, 10, 20, 50],
        [0.001, 0.002, 0.005, 0.01, 0.02],
        [10.0, 20.0, 35.0, 50.0],
    ):
        rescue = (
            (~either_direct)
            & either_relaxed
            & (max_ani >= ani)
            & (max_xny >= xny)
            & (max_minaf >= minaf)
            & ((min_cv <= cvmax) | (max_norm_ab >= 1e-5))
        )
        yield f"rule_direct_plus_relaxed_ani{ani:g}_xny{xny}_minaf{minaf:g}_cv{cvmax:g}", either_direct | rescue

    for zani, xny, zaf in itertools.product(
        [0.90, 0.92, 0.94, 0.96, 0.98, 0.99],
        [1, 2, 3, 5, 10, 20],
        [0.005, 0.01, 0.02, 0.05, 0.10],
    ):
        rescue = (
            (~either_direct)
            & either_relaxed
            & (max_zip_ani >= zani)
            & (max_xny >= xny)
            & (max_zip_af >= zaf)
        )
        yield f"rule_direct_plus_ziprescue_zani{zani:g}_xny{xny}_zaf{zaf:g}", either_direct | rescue


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load_features()
    gold_counts = load_gold_counts()
    fast_ctx = prepare_fast_context(df, gold_counts)

    summary_rows = []
    seen: set[bytes] = set()
    for name, mask in candidate_masks(df):
        mask = np.asarray(mask, dtype=bool)
        key = np.packbits(mask).tobytes()
        if key in seen:
            continue
        seen.add(key)
        rec = score_mask_fast(fast_ctx, mask, method=name)
        rec["source"] = "candidate_rule"
        summary_rows.append(rec)

    rules = pd.DataFrame(summary_rows)
    baselines = load_baseline_summary()
    combined = pd.concat([rules, baselines], ignore_index=True, sort=False)
    sort_cols = ["mean_F1", "mean_FP_plus_FN", "min_dataset_F1"]
    combined = combined.sort_values(sort_cols, ascending=[False, True, False])
    rules = rules.sort_values(sort_cols, ascending=[False, True, False])

    rules.to_csv(OUT / "candidate_rule_summary.tsv", sep="\t", index=False)
    combined.to_csv(OUT / "combined_summary.tsv", sep="\t", index=False)
    detail_rows: list[dict[str, object]] = []
    for name, mask in candidate_masks(df):
        if name in {"rule_u_direct", "rule_s_direct", "rule_u_or_s_direct"}:
            score_mask(df, np.asarray(mask, dtype=bool), gold_counts, method=name, detail_rows=detail_rows)
    pd.DataFrame(detail_rows).to_csv(OUT / "selected_baseline_rule_sample_metrics.tsv", sep="\t", index=False)
    print(combined.head(40).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
