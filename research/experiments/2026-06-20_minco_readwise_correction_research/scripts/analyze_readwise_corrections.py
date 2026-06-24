#!/usr/bin/env python3
"""Evaluate retrospective minco readwise correction ideas on CAMI marine sample0.

This script intentionally keeps large joined tables in /tmp and writes compact
summaries into the experiment directory.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+|[A-Z]{1,4}_[0-9]+\.[0-9]+)")
DROP_GOLD_NAMES = {"unidentified", "unidentified plasmid", "unidentified virus"}
FEATURE_COLS = [
    "ANI",
    "raw_aaf_ani",
    "aaf_ani_zip",
    "aaf_ani_zinb",
    "log_support",
    "Ref_breadth",
    "Real_min_align_fraction",
    "diff_rate",
    "section_rate",
    "mut2_rate",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_cv",
    "Ref_zero_fraction",
    "Read_match_fraction",
    "Block_match_fraction",
    "Relative_abundance_depth",
]


def extract_accession(text: object) -> str:
    if not isinstance(text, str):
        return ""
    m = ACC_RE.search(text)
    return m.group(1) if m else ""


def bounded_ani(x: np.ndarray | float) -> np.ndarray | float:
    return np.clip(x, 0.0, 1.0)


def aaf_ani_from_af(af: np.ndarray, ctx_k: float) -> np.ndarray:
    af = np.asarray(af, dtype=float)
    out = np.zeros_like(af)
    mask = af > 0.0
    out[mask] = 1.0 + np.log(np.clip(af[mask], 1e-300, 1.0)) / ctx_k
    return bounded_ani(out)


def parse_gold_profile(path: Path, sample_id: str) -> pd.DataFrame:
    rows = []
    active = False
    seen_sample = False
    target_sample_id = sample_id.strip()
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("@SampleID:"):
                sid = line.split(":", 1)[1].strip()
                if seen_sample and sid != target_sample_id:
                    break
                active = sid == target_sample_id
                seen_sample = seen_sample or active
                continue
            if not active or not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] == "TAXID":
                continue
            taxid, rank, taxpath, taxpathsn, pct = fields[:5]
            if rank != "species":
                continue
            try:
                pct_f = float(pct)
            except ValueError:
                continue
            if pct_f <= 0.0:
                continue
            name = taxpathsn.split("|")[-1].strip()
            if name.lower() in DROP_GOLD_NAMES:
                continue
            rows.append(
                {
                    "taxid": taxid,
                    "rank": rank,
                    "taxpath": taxpath,
                    "taxpathsn": taxpathsn,
                    "species_name": name,
                    "gold_percentage": pct_f,
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(f"no species rows found for {sample_id} in {path}")
    return df.drop_duplicates("taxid", keep="first").reset_index(drop=True)


def parse_species_taxmap(path: Path) -> Dict[str, Dict[str, str]]:
    taxmap: Dict[str, Dict[str, str]] = {}
    with path.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            key, taxid, rank, taxpath, taxpathsn = fields[:5]
            if key.lower() in {"key", "ref_key", "accession"}:
                continue
            if rank != "species":
                continue
            rec = {
                "taxid": taxid,
                "rank": rank,
                "taxpath": taxpath,
                "taxpathsn": taxpathsn,
                "species_name": taxpathsn.split("|")[-1] if taxpathsn else "",
            }
            taxmap.setdefault(key, rec)
            acc = extract_accession(key)
            if acc:
                taxmap.setdefault(acc, rec)
    return taxmap


def map_tax_record(acc: str, taxmap: Mapping[str, Mapping[str, str]]) -> Mapping[str, str]:
    if acc in taxmap:
        return taxmap[acc]
    if "." in acc:
        short = acc.rsplit(".", 1)[0]
        if short in taxmap:
            return taxmap[short]
    return {}


def load_sylph(path: Path, taxmap: Mapping[str, Mapping[str, str]]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["accession"] = df["Genome_file"].map(extract_accession)
    tax_records = [map_tax_record(acc, taxmap) for acc in df["accession"]]
    df["taxid"] = [r.get("taxid", "") for r in tax_records]
    df["species_name"] = [r.get("species_name", "") for r in tax_records]
    for col in ["Taxonomic_abundance", "Sequence_abundance", "Adjusted_ANI", "Eff_cov", "Naive_ANI"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_minco(path: Path, taxmap: Mapping[str, Mapping[str, str]], ctx_k: float) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["accession"] = df["Ref"].map(extract_accession)
    missing = df["accession"] == ""
    if missing.any() and "Ref_annotation" in df.columns:
        df.loc[missing, "accession"] = df.loc[missing, "Ref_annotation"].map(extract_accession)
    tax_records = [map_tax_record(acc, taxmap) for acc in df["accession"]]
    df["taxid"] = [r.get("taxid", "") for r in tax_records]
    df["species_name"] = [r.get("species_name", "") for r in tax_records]

    numeric_cols = [
        "ANI",
        "Distance",
        "XnY_ctx",
        "N_diff_obj",
        "N_diff_obj_section",
        "N_mut2_ctx",
        "Real_min_align_fraction",
        "Reads_with_ctx_match",
        "Total_reads",
        "Read_match_fraction",
        "Unique_query_ctx",
        "Unique_query_ctx_hit",
        "Unique_ref_ctx_hit",
        "Density_block_ctx",
        "Total_density_blocks",
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
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        support = df["XnY_ctx"].replace(0, np.nan)
        df["minco_dist"] = 1.0 - df["ANI"]
        df["raw_aaf_ani"] = aaf_ani_from_af(df["Ref_breadth"].fillna(0.0).to_numpy(), ctx_k)
        df["diff_rate"] = (df["N_diff_obj"] / support).fillna(0.0)
        df["section_rate"] = (df["N_diff_obj_section"] / support).fillna(0.0)
        df["mut2_rate"] = (df["N_mut2_ctx"] / support).fillna(0.0)
        df["log_support"] = np.log1p(df["XnY_ctx"].fillna(0.0))
        df["log_mean_depth"] = np.log1p(df["Ref_mean_depth"].fillna(0.0))
        df["log_hit_mean_depth"] = np.log1p(df["Ref_hit_mean_depth"].fillna(0.0))
    return df


def solve_zip_pi_vec(q: np.ndarray, mu: np.ndarray, grid_n: int = 512, chunk: int = 8192) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    mu = np.asarray(mu, dtype=float)
    out = q.copy()
    grid = np.linspace(0.0, 1.0, grid_n)
    for start in range(0, len(q), chunk):
        stop = min(start + chunk, len(q))
        qc = np.clip(q[start:stop], 0.0, 1.0)
        muc = np.maximum(mu[start:stop], 0.0)
        valid_rows = (qc > 0.0) & (muc > 0.0) & (qc < 1.0)
        if not valid_rows.any():
            out[start:stop] = np.where(qc >= 1.0, 1.0, qc)
            continue
        pi = qc[:, None] + (1.0 - qc[:, None]) * grid[None, :]
        pi = np.clip(pi, 1e-12, 1.0)
        lam = muc[:, None] / pi
        q_pred = pi * (1.0 - np.exp(-np.clip(lam, 0.0, 700.0)))
        err = np.square(q_pred - qc[:, None])
        best = np.argmin(err, axis=1)
        out[start:stop] = pi[np.arange(stop - start), best]
        out[start:stop] = np.where(valid_rows, out[start:stop], np.where(qc >= 1.0, 1.0, qc))
    return np.clip(out, 0.0, 1.0)


def solve_zinb_pi_vec(
    q: np.ndarray, mu: np.ndarray, var: np.ndarray, grid_n: int = 512, chunk: int = 4096
) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    mu = np.asarray(mu, dtype=float)
    var = np.asarray(var, dtype=float)
    out = q.copy()
    grid = np.linspace(0.0, 1.0, grid_n)
    for start in range(0, len(q), chunk):
        stop = min(start + chunk, len(q))
        qc = np.clip(q[start:stop], 0.0, 1.0)
        muc = np.maximum(mu[start:stop], 0.0)
        varc = np.maximum(var[start:stop], 0.0)
        valid_rows = (qc > 0.0) & (muc > 0.0) & (qc < 1.0) & (varc > muc)
        if not valid_rows.any():
            out[start:stop] = np.where(qc >= 1.0, 1.0, qc)
            continue
        pi = qc[:, None] + (1.0 - qc[:, None]) * grid[None, :]
        pi = np.clip(pi, 1e-12, 1.0)
        inv_k = ((varc - muc)[:, None] * pi / np.square(np.maximum(muc, 1e-12))[:, None]) - (
            1.0 - pi
        )
        valid = inv_k > 1e-9
        k = np.where(valid, 1.0 / np.clip(inv_k, 1e-9, np.inf), np.nan)
        mean_present = muc[:, None] / pi
        log_p0 = k * (np.log(k) - np.log(k + mean_present))
        p0 = np.exp(np.clip(log_p0, -745.0, 0.0))
        q_pred = pi * (1.0 - p0)
        err = np.square(q_pred - qc[:, None])
        err[~valid] = np.inf
        all_bad = ~np.isfinite(err).any(axis=1)
        best = np.argmin(err, axis=1)
        cand = pi[np.arange(stop - start), best]
        cand[all_bad] = qc[all_bad]
        out[start:stop] = np.where(valid_rows, cand, np.where(qc >= 1.0, 1.0, qc))
    return np.clip(out, 0.0, 1.0)


def score_taxids(pred_taxids: Iterable[str], gold_taxids: Set[str]) -> Dict[str, float]:
    pred = {str(x) for x in pred_taxids if str(x)}
    tp = len(pred & gold_taxids)
    fp = len(pred - gold_taxids)
    fn = len(gold_taxids - pred)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "pred_taxa": len(pred),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
    }


def abundance_corr(rows: pd.DataFrame, gold: pd.DataFrame) -> Tuple[float, float]:
    if rows.empty:
        return (float("nan"), float("nan"))
    work = rows.copy()
    work["Relative_abundance_depth"] = pd.to_numeric(work["Relative_abundance_depth"], errors="coerce").fillna(0.0)
    pred = work.groupby("taxid", as_index=False)["Relative_abundance_depth"].sum()
    total = pred["Relative_abundance_depth"].sum()
    if total > 0.0:
        pred["pred_pct"] = 100.0 * pred["Relative_abundance_depth"] / total
    else:
        pred["pred_pct"] = 0.0
    merged = pred.merge(gold[["taxid", "gold_percentage"]], on="taxid", how="inner")
    if len(merged) < 2:
        return (float("nan"), float("nan"))
    pearson = float(np.corrcoef(merged["pred_pct"], merged["gold_percentage"])[0, 1])
    mae = float(np.mean(np.abs(merged["pred_pct"] - merged["gold_percentage"])))
    return pearson, mae


def eval_strategy(
    name: str,
    rows: pd.DataFrame,
    mask: pd.Series,
    gold: pd.DataFrame,
    sylph_only_tp: Set[str],
    note: str = "",
) -> Dict[str, object]:
    selected = rows.loc[mask & rows["taxid"].astype(bool)].copy()
    gold_taxids = set(gold["taxid"].astype(str))
    stats = score_taxids(selected["taxid"], gold_taxids)
    pearson, abundance_mae = abundance_corr(selected, gold)
    pred = set(selected["taxid"].astype(str))
    stats.update(
        {
            "strategy": name,
            "selected_rows": int(len(selected)),
            "recovered_sylph_only_tp": len(pred & sylph_only_tp),
            "tp_abundance_pearson": pearson,
            "tp_abundance_mae_pct": abundance_mae,
            "note": note,
        }
    )
    return stats


def grid_best(
    rows: pd.DataFrame,
    metric: str,
    gold: pd.DataFrame,
    sylph_only_tp: Set[str],
    label: str,
    support: int = 10,
    breadth: float = 0.05,
    cuts: Sequence[float] = tuple(np.round(np.arange(0.86, 0.991, 0.005), 3)),
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    records = []
    for cut in cuts:
        mask = (rows["XnY_ctx"] >= support) & (rows["Ref_breadth"] >= breadth) & (rows[metric] >= cut)
        rec = eval_strategy(
            f"{label}_cut_{cut:.3f}",
            rows,
            mask,
            gold,
            sylph_only_tp,
            note=f"support>={support}; Ref_breadth>={breadth}; {metric}>={cut:.3f}",
        )
        rec["cut"] = cut
        rec["metric"] = metric
        records.append(rec)
    best = min(records, key=lambda r: (r["FP"] + r["FN"], -r["F1"], -r["TP"]))
    best = dict(best)
    best["strategy"] = f"{label}_best_{metric}"
    return best, records


def representative_called_rows(rows: pd.DataFrame, mask: pd.Series, metric: str) -> pd.DataFrame:
    selected = rows.loc[mask & rows["taxid"].astype(bool)].copy()
    if selected.empty:
        return selected
    selected["_rank_metric"] = pd.to_numeric(selected[metric], errors="coerce").fillna(-1.0)
    selected = selected.sort_values(["taxid", "_rank_metric", "XnY_ctx"], ascending=[True, False, False])
    return selected.groupby("taxid", as_index=False).head(1).drop(columns=["_rank_metric"])


def add_feature_strata_summary(
    rows: pd.DataFrame,
    strategy_defs: Sequence[Tuple[str, pd.Series, str]],
    gold_taxids: Set[str],
    out_path: Path,
) -> None:
    feature_cols = [
        "ANI",
        "raw_aaf_ani",
        "aaf_ani_zip",
        "aaf_ani_zinb",
        "aaf_ani_sylph_effcov_oracle",
        "hgb_sylphcal_ani",
        "XnY_ctx",
        "Ref_breadth",
        "Ref_mean_depth",
        "Ref_hit_mean_depth",
        "Ref_depth_cv",
        "diff_rate",
        "section_rate",
        "mut2_rate",
        "Relative_abundance_depth",
    ]
    out_records = []
    for name, mask, metric in strategy_defs:
        reps = representative_called_rows(rows, mask, metric)
        if reps.empty:
            continue
        reps["category"] = np.where(reps["taxid"].astype(str).isin(gold_taxids), "TP", "FP")
        for category, sub in reps.groupby("category"):
            rec: Dict[str, object] = {"strategy": name, "category": category, "taxa": int(len(sub))}
            for col in feature_cols:
                vals = pd.to_numeric(sub[col], errors="coerce").dropna().to_numpy(dtype=float)
                if len(vals) == 0:
                    continue
                rec[f"{col}_median"] = float(np.median(vals))
                rec[f"{col}_p10"] = float(np.quantile(vals, 0.10))
                rec[f"{col}_p90"] = float(np.quantile(vals, 0.90))
            out_records.append(rec)
    pd.DataFrame(out_records).to_csv(out_path, sep="\t", index=False)


def add_sylph_columns(minco: pd.DataFrame, sylph: pd.DataFrame) -> pd.DataFrame:
    cols = ["accession", "Adjusted_ANI", "Eff_cov", "Taxonomic_abundance", "Naive_ANI"]
    s = sylph[cols].dropna(subset=["accession"]).drop_duplicates("accession", keep="first")
    s = s.rename(
        columns={
            "Adjusted_ANI": "sylph_adjusted_ani_pct",
            "Eff_cov": "sylph_eff_cov",
            "Taxonomic_abundance": "sylph_tax_abundance_pct",
            "Naive_ANI": "sylph_naive_ani_pct",
        }
    )
    return minco.merge(s, on="accession", how="left")


def run_sklearn_models(rows: pd.DataFrame, gold: pd.DataFrame, outdir: Path) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []
    usable = rows.loc[rows["taxid"].astype(bool)].copy()
    usable[FEATURE_COLS] = usable[FEATURE_COLS].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    gold_taxids = set(gold["taxid"].astype(str))
    usable["is_gold_taxid"] = usable["taxid"].astype(str).isin(gold_taxids).astype(int)

    try:
        from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
        from sklearn.linear_model import Ridge
        from sklearn.metrics import (
            average_precision_score,
            mean_absolute_error,
            r2_score,
            roc_auc_score,
        )
        from sklearn.model_selection import KFold, StratifiedShuffleSplit, cross_val_predict
    except Exception as exc:  # pragma: no cover - environment diagnostic
        records.append({"analysis": "sklearn_unavailable", "note": repr(exc)})
        return records

    X = usable[FEATURE_COLS].to_numpy(dtype=float)
    y = usable["is_gold_taxid"].to_numpy(dtype=int)
    if len(np.unique(y)) == 2:
        splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=11)
        train_idx, test_idx = next(splitter.split(X, y))
        clf = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.06, max_leaf_nodes=31, random_state=11)
        clf.fit(X[train_idx], y[train_idx])
        prob = clf.predict_proba(X[test_idx])[:, 1]
        clf_all = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.06, max_leaf_nodes=31, random_state=17)
        clf_all.fit(X, y)
        usable_out = usable[["accession", "taxid", "species_name", "ANI", "XnY_ctx", "Ref_breadth"]].copy()
        usable_out["gold_classifier_probability"] = clf_all.predict_proba(X)[:, 1]
        usable_out.to_csv(outdir / "row_gold_classifier_allrows_probability.tsv", sep="\t", index=False)
        records.append(
            {
                "analysis": "row_gold_classifier_holdout",
                "rows": int(len(y)),
                "positive_rows": int(y.sum()),
                "test_rows": int(len(test_idx)),
                "roc_auc": float(roc_auc_score(y[test_idx], prob)),
                "average_precision": float(average_precision_score(y[test_idx], prob)),
                "note": "diagnostic only; row labels come from one CAMI sample and are not an independent training set",
            }
        )

    matched = usable.loc[usable["sylph_adjusted_ani_pct"].notna()].copy()
    if len(matched) >= 20:
        Xm = matched[FEATURE_COLS].to_numpy(dtype=float)
        target = matched["sylph_adjusted_ani_pct"].to_numpy(dtype=float) / 100.0
        raw = matched["ANI"].to_numpy(dtype=float)
        kfold = KFold(n_splits=min(5, len(matched)), shuffle=True, random_state=13)
        ridge = Ridge(alpha=1.0)
        hgb = HistGradientBoostingRegressor(max_iter=120, learning_rate=0.05, max_leaf_nodes=15, random_state=13)
        pred_ridge = cross_val_predict(ridge, Xm, target, cv=kfold)
        pred_hgb = cross_val_predict(hgb, Xm, target, cv=kfold)
        regression_out = matched[
            ["accession", "taxid", "species_name", "ANI", "raw_aaf_ani", "aaf_ani_zinb", "sylph_adjusted_ani_pct"]
        ].copy()
        regression_out["ridge_pred_sylph_ani"] = pred_ridge
        regression_out["hgb_pred_sylph_ani"] = pred_hgb
        regression_out.to_csv(outdir / "sylph_matched_distance_cv_predictions.tsv", sep="\t", index=False)
        for name, pred in [
            ("raw_minco_ani_vs_sylph", raw),
            ("ridge_feature_dist_calibration_cv", pred_ridge),
            ("hgb_feature_dist_calibration_cv", pred_hgb),
        ]:
            records.append(
                {
                    "analysis": name,
                    "rows": int(len(matched)),
                    "MAE": float(mean_absolute_error(target, pred)),
                    "RMSE": float(math.sqrt(np.mean(np.square(target - pred)))),
                    "pearson": float(np.corrcoef(target, pred)[0, 1]) if len(target) > 1 else float("nan"),
                    "r2": float(r2_score(target, pred)),
                    "note": "target is Sylph adjusted ANI, used here as pseudo-truth for sequencing-error/distance inflation",
                }
            )
    return records


def fit_sylph_calibrated_ani(rows: pd.DataFrame, outdir: Path) -> Optional[np.ndarray]:
    try:
        from sklearn.ensemble import HistGradientBoostingRegressor
    except Exception:
        return None
    matched = rows.loc[rows["taxid"].astype(bool) & rows["sylph_adjusted_ani_pct"].notna()].copy()
    if len(matched) < 20:
        return None
    model = HistGradientBoostingRegressor(max_iter=160, learning_rate=0.05, max_leaf_nodes=15, random_state=23)
    X_train = matched[FEATURE_COLS].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)
    y_train = matched["sylph_adjusted_ani_pct"].to_numpy(dtype=float) / 100.0
    model.fit(X_train, y_train)
    X_all = rows[FEATURE_COLS].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)
    pred = np.clip(model.predict(X_all), 0.0, 1.0)
    out = rows[["accession", "taxid", "species_name", "ANI", "XnY_ctx", "Ref_breadth"]].copy()
    out["hgb_sylphcal_ani"] = pred
    out.to_csv(outdir / "hgb_sylphcal_ani_allrows.tsv", sep="\t", index=False)
    return pred


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--minco", type=Path, default=Path("/tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_unfiltered.tsv"))
    ap.add_argument("--sylph", type=Path, default=Path("/tmp/sylph_marine_sample0/profile.tsv"))
    ap.add_argument("--gold", type=Path, default=Path("/tmp/gs_marine_short.profile"))
    ap.add_argument("--taxmap", type=Path, default=Path("/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv"))
    ap.add_argument("--sample-id", default="marmgCAMI2_short_read_sample_0")
    ap.add_argument("--ctx-k", type=float, default=11.0)
    ap.add_argument("--outdir", type=Path, default=Path("/tmp/minco_readwise_correction_20260620"))
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gold = parse_gold_profile(args.gold, args.sample_id)
    taxmap = parse_species_taxmap(args.taxmap)
    sylph = load_sylph(args.sylph, taxmap)
    minco = load_minco(args.minco, taxmap, args.ctx_k)
    minco = add_sylph_columns(minco, sylph)
    minco["is_gold_taxid"] = minco["taxid"].astype(str).isin(set(gold["taxid"].astype(str)))

    q = minco["Ref_breadth"].fillna(0.0).to_numpy(dtype=float)
    mu = minco["Ref_mean_depth"].fillna(0.0).to_numpy(dtype=float)
    var = minco["Ref_depth_variance"].fillna(0.0).to_numpy(dtype=float)
    minco["af_zip"] = solve_zip_pi_vec(q, mu)
    minco["af_zinb"] = solve_zinb_pi_vec(q, mu, var)
    minco["aaf_ani_zip"] = aaf_ani_from_af(minco["af_zip"].to_numpy(), args.ctx_k)
    minco["aaf_ani_zinb"] = aaf_ani_from_af(minco["af_zinb"].to_numpy(), args.ctx_k)
    sylph_eff = minco["sylph_eff_cov"].fillna(0.0).to_numpy(dtype=float)
    p_seen = 1.0 - np.exp(-np.clip(sylph_eff, 0.0, 700.0))
    af_oracle = np.divide(q, p_seen, out=np.zeros_like(q), where=p_seen > 0.0)
    minco["af_sylph_effcov_oracle"] = np.clip(af_oracle, 0.0, 1.0)
    minco["aaf_ani_sylph_effcov_oracle"] = aaf_ani_from_af(
        minco["af_sylph_effcov_oracle"].to_numpy(), args.ctx_k
    )
    sylphcal_pred = fit_sylph_calibrated_ani(minco, args.outdir)
    if sylphcal_pred is not None:
        minco["hgb_sylphcal_ani"] = sylphcal_pred
    else:
        minco["hgb_sylphcal_ani"] = np.nan

    joined_path = args.outdir / "joined_s10000_rows.tsv"
    minco.to_csv(joined_path, sep="\t", index=False)
    gold.to_csv(args.outdir / "gold_species_sample0.tsv", sep="\t", index=False)
    sylph.to_csv(args.outdir / "sylph_with_taxids.tsv", sep="\t", index=False)

    gold_taxids = set(gold["taxid"].astype(str))
    sylph_pred = set(sylph.loc[sylph["taxid"].astype(bool), "taxid"].astype(str))
    baseline_mask = (minco["XnY_ctx"] >= 10) & (minco["Ref_breadth"] >= 0.05) & (minco["ANI"] >= 0.92)
    baseline_tp = set(minco.loc[baseline_mask, "taxid"].astype(str)) & gold_taxids
    sylph_only_tp = (sylph_pred & gold_taxids) - baseline_tp

    strategy_records: List[Dict[str, object]] = []
    strategy_defs: List[Tuple[str, pd.Series, str]] = []
    strategy_records.append(
        eval_strategy(
            "sylph_gtdb_r226_default",
            sylph.rename(columns={"Adjusted_ANI": "ANI"}).assign(
                XnY_ctx=1,
                Ref_breadth=1.0,
                Relative_abundance_depth=sylph["Taxonomic_abundance"].fillna(0.0),
            ),
            sylph["taxid"].astype(bool),
            gold,
            sylph_only_tp,
            note="Sylph profile rows mapped through the same taxmap; not a minco strategy",
        )
    )
    strategy_defs.append(("minco_s10000_opt_combo_f0.05_n0.92_t10", baseline_mask, "ANI"))
    strategy_records.append(
        eval_strategy(
            "minco_s10000_opt_combo_f0.05_n0.92_t10",
            minco,
            baseline_mask,
            gold,
            sylph_only_tp,
            note="previous optimized combined species+genus filter",
        )
    )
    strategy_records.append(
        eval_strategy(
            "minco_s10000_species_f0.05_n0.94_t10",
            minco,
            (minco["XnY_ctx"] >= 10) & (minco["Ref_breadth"] >= 0.05) & (minco["ANI"] >= 0.94),
            gold,
            sylph_only_tp,
            note="previous species-only threshold",
        )
    )
    strategy_defs.append(
        (
            "minco_s10000_species_f0.05_n0.94_t10",
            (minco["XnY_ctx"] >= 10) & (minco["Ref_breadth"] >= 0.05) & (minco["ANI"] >= 0.94),
            "ANI",
        )
    )
    grid_records: List[Dict[str, object]] = []
    for metric, label in [
        ("ANI", "minco_raw_ani"),
        ("raw_aaf_ani", "raw_ctx_aaf"),
        ("aaf_ani_zip", "zip_coverage_ctx_aaf"),
        ("aaf_ani_zinb", "zinb_coverage_ctx_aaf"),
        ("aaf_ani_sylph_effcov_oracle", "oracle_sylph_effcov_ctx_aaf"),
        ("hgb_sylphcal_ani", "hgb_sylphcal_dist"),
    ]:
        if minco[metric].notna().sum() == 0:
            continue
        best, records = grid_best(minco, metric, gold, sylph_only_tp, label)
        strategy_records.append(best)
        grid_records.extend(records)
        cut = float(best["cut"])
        strategy_defs.append(
            (
                str(best["strategy"]),
                (minco["XnY_ctx"] >= 10) & (minco["Ref_breadth"] >= 0.05) & (minco[metric] >= cut),
                metric,
            )
        )

    pd.DataFrame(strategy_records).to_csv(args.outdir / "strategy_summary.tsv", sep="\t", index=False)
    pd.DataFrame(grid_records).to_csv(args.outdir / "strategy_grid.tsv", sep="\t", index=False)
    add_feature_strata_summary(
        minco, strategy_defs, gold_taxids, args.outdir / "feature_strata_summary.tsv"
    )

    model_records = run_sklearn_models(minco, gold, args.outdir)
    pd.DataFrame(model_records).to_csv(args.outdir / "model_summary.tsv", sep="\t", index=False)

    # Detail for Sylph true positives not called by the current minco S10000 threshold.
    best_rows = (
        minco.loc[minco["taxid"].astype(str).isin(sylph_only_tp)]
        .sort_values(["taxid", "ANI", "XnY_ctx"], ascending=[True, False, False])
        .groupby("taxid", as_index=False)
        .head(1)
    )
    cols = [
        "taxid",
        "species_name",
        "accession",
        "ANI",
        "raw_aaf_ani",
        "aaf_ani_zip",
        "aaf_ani_zinb",
        "aaf_ani_sylph_effcov_oracle",
        "hgb_sylphcal_ani",
        "XnY_ctx",
        "Ref_breadth",
        "Ref_mean_depth",
        "Ref_hit_mean_depth",
        "Ref_depth_variance",
        "Ref_depth_cv",
        "N_diff_obj",
        "N_diff_obj_section",
        "N_mut2_ctx",
        "sylph_adjusted_ani_pct",
        "sylph_eff_cov",
        "sylph_tax_abundance_pct",
    ]
    best_rows[cols].to_csv(args.outdir / "sylph_only_tp_minco_correction_features.tsv", sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {"key": "minco_rows", "value": len(minco)},
            {"key": "minco_mapped_species_rows", "value": int(minco["taxid"].astype(bool).sum())},
            {"key": "gold_species_count", "value": len(gold_taxids)},
            {"key": "sylph_mapped_species_predictions", "value": len(sylph_pred)},
            {"key": "sylph_only_tp_not_minco_baseline", "value": len(sylph_only_tp)},
            {"key": "joined_table", "value": str(joined_path)},
        ]
    )
    manifest.to_csv(args.outdir / "manifest.tsv", sep="\t", index=False)

    print(f"wrote {args.outdir}")
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
