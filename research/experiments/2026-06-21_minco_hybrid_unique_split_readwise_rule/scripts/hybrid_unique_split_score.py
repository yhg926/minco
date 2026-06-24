#!/usr/bin/env python3
"""Evaluate hybrid minco readwise rules from unique and split candidate tables."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
ANALYSIS_SCRIPT_DIR = (
    REPO_ROOT
    / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
)
sys.path.insert(0, str(ANALYSIS_SCRIPT_DIR))

from analyze_readwise_corrections import (  # noqa: E402
    load_minco,
    load_sylph,
    parse_gold_profile,
    parse_species_taxmap,
    score_taxids,
)


NUMERIC_COLS = [
    "ANI",
    "Real_min_align_fraction",
    "XnY_ctx",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_cv",
    "Relative_abundance_depth",
    "Normalized_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
]


def _as_num(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(0.0)


def gold_scope(gold: pd.DataFrame, scope: str) -> pd.DataFrame:
    if scope == "all":
        return gold
    if scope == "bacteria":
        return gold.loc[gold["taxpathsn"].str.contains("Bacteria", regex=False)].copy()
    if scope == "virus":
        return gold.loc[gold["taxpathsn"].str.contains("Viruses", regex=False)].copy()
    raise ValueError(f"unsupported scope: {scope}")


def build_taxid_scope_map(taxmap: Mapping[str, Mapping[str, str]]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in taxmap.values():
        taxid = str(rec.get("taxid", ""))
        path = rec.get("taxpathsn", "")
        if not taxid:
            continue
        if "Viruses" in path:
            out.setdefault(taxid, "virus")
        elif "Bacteria" in path:
            out.setdefault(taxid, "bacteria")
        else:
            out.setdefault(taxid, "other")
    return out


def aggregate_minco_by_taxid(path: Path, taxmap: Mapping[str, Mapping[str, str]], prefix: str) -> pd.DataFrame:
    rows = load_minco(path, taxmap, 11.0)
    rows = rows.loc[rows["taxid"].astype(bool)].copy()
    if rows.empty:
        return pd.DataFrame(columns=["taxid"])
    for col in NUMERIC_COLS:
        rows[col] = _as_num(rows, col)
    row_counts = rows.groupby("taxid").size().rename(f"{prefix}_rows").reset_index()
    sort_cols = [c for c in ["ANI", "Ref_breadth", "XnY_ctx", "Ref_mean_depth"] if c in rows.columns]
    rows = rows.sort_values(sort_cols, ascending=[False] * len(sort_cols))
    out = rows.groupby("taxid", as_index=False).head(1)
    keep_cols = ["taxid"] + [col for col in NUMERIC_COLS if col in out.columns]
    out = out[keep_cols].merge(row_counts, on="taxid", how="left")
    for col in NUMERIC_COLS:
        if col in out.columns:
            out = out.rename(columns={col: f"{prefix}_{col}"})
    return out


def join_unique_split(unique_path: Path, split_path: Path, taxmap: Mapping[str, Mapping[str, str]]) -> pd.DataFrame:
    unique = aggregate_minco_by_taxid(unique_path, taxmap, "u")
    split = aggregate_minco_by_taxid(split_path, taxmap, "s")
    joined = unique.merge(split, on="taxid", how="outer")
    for col in joined.columns:
        if col == "taxid":
            continue
        joined[col] = pd.to_numeric(joined[col], errors="coerce").fillna(0.0)
    joined["u_present"] = joined.get("u_rows", 0.0) > 0
    joined["s_present"] = joined.get("s_rows", 0.0) > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        joined["split_unique_xny_ratio"] = joined["s_XnY_ctx"] / joined["u_XnY_ctx"].replace(0, np.nan)
        joined["split_unique_breadth_ratio"] = joined["s_Ref_breadth"] / joined["u_Ref_breadth"].replace(0, np.nan)
        joined["split_unique_ani_delta"] = joined["s_ANI"] - joined["u_ANI"]
    for col in ["split_unique_xny_ratio", "split_unique_breadth_ratio", "split_unique_ani_delta"]:
        joined[col] = pd.to_numeric(joined[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return joined


def score_prediction(pred_taxids: Iterable[str], gold_taxids: Set[str]) -> Dict[str, object]:
    pred = set(str(x) for x in pred_taxids if str(x))
    out: Dict[str, object] = {"pred_taxa": len(pred)}
    out.update(score_taxids(pred, gold_taxids))
    out["FP_plus_FN"] = out["FP"] + out["FN"]
    return out


def mask_unique_current(df: pd.DataFrame) -> pd.Series:
    return (
        (df["u_XnY_ctx"] >= 10)
        & (df["u_ANI"] >= 0.94)
        & (df["u_Real_min_align_fraction"] >= 0.05)
    )


def mask_split_current(df: pd.DataFrame) -> pd.Series:
    return (
        (df["s_XnY_ctx"] >= 10)
        & (df["s_ANI"] >= 0.94)
        & (df["s_Real_min_align_fraction"] >= 0.05)
    )


def mask_unique_relaxed(df: pd.DataFrame) -> pd.Series:
    return (df["u_XnY_ctx"] >= 1) & (df["u_ANI"] >= 0.90) & (df["u_Ref_breadth"] >= 0.005)


def mask_split_relaxed(df: pd.DataFrame) -> pd.Series:
    return (df["s_XnY_ctx"] >= 1) & (df["s_ANI"] >= 0.90) & (df["s_Ref_breadth"] >= 0.005)


def hybrid_mask(df: pd.DataFrame, p: Mapping[str, float]) -> pd.Series:
    primary = (
        (df["u_XnY_ctx"] >= p["u_xny"])
        & (df["u_ANI"] >= p["u_ani"])
        & (df["u_Ref_breadth"] >= p["u_breadth"])
    )
    rescue = (
        (df["s_XnY_ctx"] >= p["s_xny"])
        & (df["s_ANI"] >= p["s_ani"])
        & (df["s_Ref_breadth"] >= p["s_breadth"])
        & (df["u_XnY_ctx"] >= p["u_confirm_xny"])
        & (df["u_Ref_breadth"] >= p["u_confirm_breadth"])
        & (df["split_unique_xny_ratio"] <= p["max_xny_ratio"])
        & (df["split_unique_breadth_ratio"] <= p["max_breadth_ratio"])
    )
    return primary | rescue


def grid_parameters() -> List[Dict[str, float]]:
    params: List[Dict[str, float]] = []
    for u_xny in [1, 10]:
        for u_ani in [0.90, 0.94, 0.96]:
            for u_breadth in [0.005, 0.03]:
                for s_xny in [25, 100]:
                    for s_ani in [0.92, 0.94, 0.96]:
                        for s_breadth in [0.01, 0.05]:
                            for u_confirm_xny in [1, 10]:
                                for max_xny_ratio in [2.0, 5.0, 9999.0]:
                                    params.append(
                                        {
                                            "u_xny": float(u_xny),
                                            "u_ani": float(u_ani),
                                            "u_breadth": float(u_breadth),
                                            "s_xny": float(s_xny),
                                            "s_ani": float(s_ani),
                                            "s_breadth": float(s_breadth),
                                            "u_confirm_xny": float(u_confirm_xny),
                                            "u_confirm_breadth": 0.0,
                                            "max_xny_ratio": float(max_xny_ratio),
                                            "max_breadth_ratio": 9999.0,
                                        }
                                    )
    return params


def evaluate_dataset(
    joined: pd.DataFrame,
    gold_taxids: Set[str],
    dataset: str,
    scope: str,
    params: Sequence[Mapping[str, float]] | None = None,
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    baselines = {
        "unique_current": mask_unique_current(joined),
        "split_current": mask_split_current(joined),
        "unique_relaxed": mask_unique_relaxed(joined),
        "split_relaxed": mask_split_relaxed(joined),
    }
    for label, mask in baselines.items():
        rec = {"dataset": dataset, "scope": scope, "label": label}
        rec.update(score_prediction(joined.loc[mask, "taxid"], gold_taxids))
        records.append(rec)
    if params is not None:
        for i, p in enumerate(params):
            mask = hybrid_mask(joined, p)
            rec = {"dataset": dataset, "scope": scope, "label": "hybrid", "param_id": i}
            rec.update(p)
            rec.update(score_prediction(joined.loc[mask, "taxid"], gold_taxids))
            records.append(rec)
    return pd.DataFrame(records)


def add_sylph_score(
    records: List[Dict[str, object]],
    dataset: str,
    scope: str,
    sylph_path: Path | None,
    taxmap: Mapping[str, Mapping[str, str]],
    gold_taxids: Set[str],
    taxid_scope: Mapping[str, str],
) -> None:
    if not sylph_path or not sylph_path.exists():
        return
    rows = load_sylph(sylph_path, taxmap)
    rows = rows.loc[rows["taxid"].astype(bool)].copy()
    if scope != "all":
        rows = rows.loc[rows["taxid"].astype(str).map(lambda x: taxid_scope.get(x, "other") == scope)]
    rec = {"dataset": dataset, "scope": scope, "label": "sylph_default"}
    rec.update(score_prediction(rows["taxid"], gold_taxids))
    records.append(rec)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxmap", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--toy-gold", type=Path, required=True)
    ap.add_argument("--toy-sample-id", default="0")
    ap.add_argument("--toy-unique", type=Path, required=True)
    ap.add_argument("--toy-split", type=Path, required=True)
    ap.add_argument("--toy-sylph", type=Path)
    ap.add_argument("--marine-gold", type=Path, required=True)
    ap.add_argument("--marine-sample-id", default="marmgCAMI2_short_read_sample_0")
    ap.add_argument("--marine-unique", type=Path, required=True)
    ap.add_argument("--marine-split", type=Path, required=True)
    ap.add_argument("--marine-sylph", type=Path)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    taxmap = parse_species_taxmap(args.taxmap)
    taxid_scope = build_taxid_scope_map(taxmap)
    datasets = {
        "toy_gtdb_bacteria": {
            "gold": parse_gold_profile(args.toy_gold, args.toy_sample_id),
            "unique": args.toy_unique,
            "split": args.toy_split,
            "sylph": args.toy_sylph,
            "scope": "bacteria",
        },
        "marine_gtdb_species": {
            "gold": parse_gold_profile(args.marine_gold, args.marine_sample_id),
            "unique": args.marine_unique,
            "split": args.marine_split,
            "sylph": args.marine_sylph,
            "scope": "all",
        },
    }

    joined_by_dataset: Dict[str, pd.DataFrame] = {}
    gold_by_dataset: Dict[str, Set[str]] = {}
    for name, cfg in datasets.items():
        joined = join_unique_split(cfg["unique"], cfg["split"], taxmap)
        scope = str(cfg["scope"])
        if scope != "all":
            joined = joined.loc[joined["taxid"].astype(str).map(lambda x: taxid_scope.get(x, "other") == scope)].copy()
        gold = gold_scope(cfg["gold"], scope)
        joined_by_dataset[name] = joined
        gold_by_dataset[name] = set(gold["taxid"].astype(str))
        joined.to_csv(args.outdir / f"{name}.joined_features.tsv", sep="\t", index=False)

    params = grid_parameters()
    grid_frames = [
        evaluate_dataset(joined_by_dataset[name], gold_by_dataset[name], name, str(datasets[name]["scope"]), params)
        for name in datasets
    ]
    grid = pd.concat(grid_frames, ignore_index=True, sort=False)
    grid.to_csv(args.outdir / "hybrid_grid_all.tsv", sep="\t", index=False)

    hybrid_grid = grid.loc[grid["label"] == "hybrid"].copy()
    mean_grid = (
        hybrid_grid.groupby("param_id", as_index=False)
        .agg(
            mean_F1=("F1", "mean"),
            mean_FP_plus_FN=("FP_plus_FN", "mean"),
            min_F1=("F1", "min"),
            sum_TP=("TP", "sum"),
            sum_FP=("FP", "sum"),
            sum_FN=("FN", "sum"),
        )
        .sort_values(["mean_F1", "mean_FP_plus_FN", "min_F1"], ascending=[False, True, False])
    )
    best_param_id = int(mean_grid.iloc[0]["param_id"])
    best_params = params[best_param_id]

    final_records: List[Dict[str, object]] = []
    for name, cfg in datasets.items():
        joined = joined_by_dataset[name]
        gold_taxids = gold_by_dataset[name]
        scope = str(cfg["scope"])
        baseline = evaluate_dataset(joined, gold_taxids, name, scope, None)
        final_records.extend(baseline.to_dict("records"))
        mask = hybrid_mask(joined, best_params)
        rec = {"dataset": name, "scope": scope, "label": "hybrid_joint_best", "param_id": best_param_id}
        rec.update(best_params)
        rec.update(score_prediction(joined.loc[mask, "taxid"], gold_taxids))
        final_records.append(rec)
        add_sylph_score(final_records, name, scope, cfg.get("sylph"), taxmap, gold_taxids, taxid_scope)

    summary = pd.DataFrame(final_records)
    for col in ["precision", "recall", "F1"]:
        if col in summary.columns:
            summary[col] = summary[col].astype(float)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    mean_grid.head(50).to_csv(args.outdir / "hybrid_grid_best50.tsv", sep="\t", index=False)
    pd.DataFrame([{"param_id": best_param_id, **best_params}]).to_csv(
        args.outdir / "hybrid_best_params.tsv", sep="\t", index=False
    )
    print(summary.to_string(index=False))
    print("\nBest hybrid params:")
    print(pd.DataFrame([{"param_id": best_param_id, **best_params}]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
