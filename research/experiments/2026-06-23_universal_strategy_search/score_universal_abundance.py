#!/usr/bin/env python3
"""Score abundance for saved universal RF/HGB MinCO calls.

The scorer avoids the lost accession taxmap by using cached joined-feature
taxids directly. It therefore scores MinCO abundance estimators; Sylph abundance
on this exact 12-sample panel needs the accession-to-taxid taxmap to be rebuilt.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
OUT = EXP / "results"
TRAIN_MANIFEST = Path("research/experiments/2026-06-21_minco_third_domain_plant_calibration/manifest_marine_toy0_2_plant0_2.tsv")
STRAIN_MANIFEST = Path("research/experiments/2026-06-21_minco_external_strain_holdout/manifest_strain0_2.tsv")
PANEL = Path("/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2")
TRAIN_FEATURES = PANEL / "train.joined_features.tsv"
TEST_FEATURES = PANEL / "test.joined_features.tsv"
TRAIN_PRED = Path("/mnt/new3T/minco_cami2_plant_20260621/calibration_marine_toy_plant_20260621/model_noleak_rules_ensemble/loso_predictions.tsv")
STRAIN_PRED = PANEL / "model_predictions.tsv"

DROP_GOLD_NAMES = {"unidentified", "unidentified plasmid", "unidentified virus"}


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def parse_gold_profile(path: Path, sample_id: str, scope: str) -> pd.DataFrame:
    rows = []
    target = str(sample_id).strip()
    active = False
    seen = False
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("@SampleID:"):
                sid = line.split(":", 1)[1].strip()
                if seen and sid != target:
                    break
                active = sid == target
                seen = seen or active
                continue
            if not active or not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] == "TAXID":
                continue
            taxid, rank, _taxpath, taxpathsn, pct = fields[:5]
            if rank != "species":
                continue
            if scope == "bacteria" and "Bacteria" not in taxpathsn:
                continue
            try:
                pct_f = float(pct)
            except ValueError:
                continue
            if pct_f <= 0.0:
                continue
            species_name = taxpathsn.split("|")[-1].strip()
            if species_name.lower() in DROP_GOLD_NAMES:
                continue
            rows.append({"taxid": str(taxid), "gold_percentage_raw": pct_f, "species_name": species_name})
    if not rows:
        raise RuntimeError(f"no scoped species rows found for sample {sample_id!r} in {path}")
    df = pd.DataFrame(rows).groupby("taxid", as_index=False).agg(gold_percentage_raw=("gold_percentage_raw", "sum"))
    total = float(df["gold_percentage_raw"].sum())
    df["gold_percentage"] = df["gold_percentage_raw"] / total * 100.0 if total > 0 else 0.0
    return df


def load_manifest() -> pd.DataFrame:
    frames = [pd.read_csv(TRAIN_MANIFEST, sep="\t", keep_default_na=False), pd.read_csv(STRAIN_MANIFEST, sep="\t", keep_default_na=False)]
    manifest = pd.concat(frames, ignore_index=True, sort=False)
    manifest["sample_key"] = manifest["sample_key"].astype(str)
    return manifest


def load_features() -> pd.DataFrame:
    df = pd.concat([pd.read_csv(TRAIN_FEATURES, sep="\t"), pd.read_csv(TEST_FEATURES, sep="\t")], ignore_index=True, sort=False)
    df["sample_key"] = df["sample_key"].astype(str)
    df["taxid"] = df["taxid"].astype(str)
    df["direct"] = (numeric(df, "u_direct_call").to_numpy() > 0) | (numeric(df, "s_direct_call").to_numpy() > 0)
    df["u_norm"] = numeric(df, "u_Normalized_abundance_depth_max")
    df["s_norm"] = numeric(df, "s_Normalized_abundance_depth_max")
    df["max_norm"] = np.maximum(df["u_norm"].to_numpy(), df["s_norm"].to_numpy())
    df["mean_nonzero_norm"] = np.where(
        (df["u_norm"] > 0) & (df["s_norm"] > 0),
        0.5 * (df["u_norm"] + df["s_norm"]),
        np.maximum(df["u_norm"], df["s_norm"]),
    )
    return df


def load_rf_hgb_predictions() -> pd.DataFrame:
    train = pd.read_csv(TRAIN_PRED, sep="\t")
    train = train.loc[train["model"] == "rf_hgb_avg", ["sample_key", "taxid", "probability", "predicted", "label"]].copy()
    strain = pd.read_csv(STRAIN_PRED, sep="\t")
    strain = strain.loc[
        strain["method"] == "train9_rf_hgb_avg",
        ["sample_key", "taxid", "probability", "predicted", "label"],
    ].copy()
    out = pd.concat([train, strain], ignore_index=True, sort=False)
    out["sample_key"] = out["sample_key"].astype(str)
    out["taxid"] = out["taxid"].astype(str)
    return out


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy()
    yr = pd.Series(y).rank(method="average").to_numpy()
    return pearson(xr, yr)


def abundance_metrics(pred: pd.DataFrame, truth: pd.DataFrame) -> dict[str, float]:
    pred = pred[["taxid", "pred_percentage"]].copy()
    truth = truth[["taxid", "gold_percentage"]].copy()
    merged = truth.merge(pred, on="taxid", how="outer").fillna(0.0)
    x = merged["pred_percentage"].to_numpy(dtype=float)
    y = merged["gold_percentage"].to_numpy(dtype=float)
    return {
        "L1": float(np.abs(x - y).sum()),
        "MAE": float(np.abs(x - y).mean()) if len(x) else float("nan"),
        "Pearson": pearson(x, y),
        "Spearman": spearman(x, y),
    }


def build_predicted_abundance(rows: pd.DataFrame, abundance_col: str) -> pd.DataFrame:
    selected = rows.loc[rows["selected"]].copy()
    if selected.empty:
        return pd.DataFrame(columns=["taxid", "pred_percentage"])
    selected["raw_abundance"] = numeric(selected, abundance_col)
    if abundance_col == "prob_x_max_norm":
        selected["raw_abundance"] = numeric(selected, "probability") * numeric(selected, "max_norm")
    total = float(selected["raw_abundance"].sum())
    if total <= 0.0:
        selected["raw_abundance"] = numeric(selected, "probability")
        total = float(selected["raw_abundance"].sum())
    if total <= 0.0:
        selected["raw_abundance"] = 1.0
        total = float(selected["raw_abundance"].sum())
    out = selected.groupby("taxid", as_index=False).agg(raw_abundance=("raw_abundance", "sum"))
    out["pred_percentage"] = out["raw_abundance"] / float(out["raw_abundance"].sum()) * 100.0
    return out[["taxid", "pred_percentage"]]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest()
    features = load_features()
    probs = load_rf_hgb_predictions()
    df = probs.merge(
        features[
            [
                "sample_key",
                "taxid",
                "dataset",
                "scope",
                "direct",
                "u_norm",
                "s_norm",
                "max_norm",
                "mean_nonzero_norm",
            ]
        ],
        on=["sample_key", "taxid"],
        how="left",
    )
    df["dataset"] = df["dataset"].fillna("")
    df["scope"] = df["scope"].fillna("")
    df["direct"] = df["direct"].fillna(False).astype(bool)

    selectors = {
        "rf_hgb_native_threshold": df["predicted"].to_numpy(dtype=int) > 0,
        "rf_hgb_t035": numeric(df, "probability").to_numpy() >= 0.35,
        "rf_hgb_t035_direct_or_prob": (numeric(df, "probability").to_numpy() >= 0.35) | df["direct"].to_numpy(dtype=bool),
        "u_or_s_direct": df["direct"].to_numpy(dtype=bool),
    }
    abundance_cols = ["max_norm", "s_norm", "u_norm", "mean_nonzero_norm", "prob_x_max_norm"]

    sample_rows = []
    for _, rec in manifest.iterrows():
        sample_key = str(rec["sample_key"])
        truth = parse_gold_profile(Path(rec["gold_profile"]), str(rec["gold_sample_id"]), str(rec["scope"]))
        sample_df = df.loc[df["sample_key"] == sample_key].copy()
        for selector_name, selector in selectors.items():
            sample_selector = selector[df["sample_key"].to_numpy() == sample_key]
            for abundance_col in abundance_cols:
                work = sample_df.copy()
                work["selected"] = sample_selector
                pred = build_predicted_abundance(work, abundance_col)
                metrics = abundance_metrics(pred, truth)
                sample_rows.append(
                    {
                        "method": f"{selector_name}_{abundance_col}",
                        "selector": selector_name,
                        "abundance_col": abundance_col,
                        "sample_key": sample_key,
                        "dataset": str(rec["dataset"]),
                        "scope": str(rec["scope"]),
                        "pred_taxa": int(len(pred)),
                        "gold_taxa": int(len(truth)),
                        **metrics,
                    }
                )

    sample = pd.DataFrame(sample_rows)
    sample.to_csv(OUT / "universal_abundance_sample_metrics.tsv", sep="\t", index=False)
    summary_rows = []
    for method, group in sample.groupby("method"):
        by_dataset = group.groupby("dataset", as_index=False).agg(dataset_L1=("L1", "mean"))
        rec = {
            "method": method,
            "samples": int(len(group)),
            "mean_L1": float(group["L1"].mean()),
            "mean_MAE": float(group["MAE"].mean()),
            "mean_Pearson": float(group["Pearson"].mean()),
            "mean_Spearman": float(group["Spearman"].mean()),
            "max_dataset_L1": float(by_dataset["dataset_L1"].max()),
        }
        for row in by_dataset.itertuples(index=False):
            rec[f"{row.dataset}_L1"] = float(row.dataset_L1)
        summary_rows.append(rec)
    summary = pd.DataFrame(summary_rows).sort_values(["mean_L1", "max_dataset_L1"], ascending=[True, True])
    summary.to_csv(OUT / "universal_abundance_summary.tsv", sep="\t", index=False)
    print(summary.head(40).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
