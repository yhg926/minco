#!/usr/bin/env python3
"""Sweep fixed-call MinCO abundance formulas on CAMI3 source-readmap truth."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, Mapping

import numpy as np
import pandas as pd

import score_cami3_gtdb_source_readmap as source_score
import score_cami3_gtdb_taxid_transfer as taxid_score


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def safe_div(num: pd.Series, den: pd.Series, floor: float = 1e-6) -> pd.Series:
    return numeric_value(num) / np.maximum(numeric_value(den), floor)


def numeric_value(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)


def truth_for_sample(sample_id: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "cami3_gtdb_source_readmap_truth.tsv", sep="\t")
    sub = truth.loc[truth["sample"].astype(int) == sample_id].copy()
    if sub.empty:
        raise SystemExit(f"missing source-readmap truth for sample {sample_id}; run score_cami3_gtdb_source_readmap.py first")
    return sub[["gtdb_species", "truth_abundance"]]


def map_minco_calls(
    sample_id: int,
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
    by_accession,
    by_core,
    taxmap,
) -> pd.DataFrame:
    path = source_score.PROFILE_PATHS[sample_id]["minco"]
    raw = pd.read_csv(path, sep="\t")
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[call].copy()
    best_ref_species, _diag = source_score.best_raw_ref_species_by_taxid(
        source_score.RAW_TABLES[sample_id],
        taxmap,
        by_accession,
        by_core,
    )

    mapped_species = []
    map_methods = []
    for row in selected.itertuples(index=False):
        taxid = str(getattr(row, "taxid", ""))
        species_name = source_score.normalize_name(getattr(row, "species_name", ""))
        gtdb_name = ""
        method = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                gtdb_name = source_score.gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if gtdb_name:
                    method = col
                    break
        if not gtdb_name:
            gtdb_name = best_ref_species.get(taxid, "")
            if gtdb_name:
                method = "raw_best_ref"
        if not gtdb_name:
            gtdb_name = taxid_to_species.get(taxid, "")
            if gtdb_name:
                method = "unique_taxid"
        if not gtdb_name and species_name:
            gtdb_name = name_to_species.get(species_name, "")
            if gtdb_name:
                method = "unique_name"
        mapped_species.append(gtdb_name)
        map_methods.append(method or "unmapped")
    selected["gtdb_species"] = mapped_species
    selected["gtdb_mapping_method"] = map_methods
    return selected.loc[selected["gtdb_species"].astype(bool)].copy()


def candidate_raws(df: pd.DataFrame) -> dict[str, np.ndarray]:
    s_mean = numeric(df, "s_Ref_mean_depth_max")
    s_hit = numeric(df, "s_Ref_hit_mean_depth_max")
    s_breadth = numeric(df, "s_Ref_breadth_max")
    s_zip = numeric(df, "s_Ref_zip_af_max")
    s_norm = numeric(df, "s_Normalized_abundance_depth_max")
    s_rel = numeric(df, "s_Relative_abundance_depth_max")
    u_mean = numeric(df, "u_Ref_mean_depth_max")
    prob = numeric(df, "calibrated_probability", 1.0)

    out: dict[str, np.ndarray] = {
        "current_calibrated_abundance": numeric_value(numeric(df, "calibrated_abundance")),
        "current_calibrated_raw": numeric_value(numeric(df, "calibrated_abundance_raw")),
        "split_mean_depth": numeric_value(s_mean),
        "split_hit_mean_x_breadth": numeric_value(s_hit * s_breadth),
        "split_normalized_depth": numeric_value(s_norm),
        "split_relative_depth": numeric_value(s_rel),
        "unique_mean_depth": numeric_value(u_mean),
        "split_mean_prob": numeric_value(s_mean * prob),
        "split_mean_prob2": numeric_value(s_mean * prob * prob),
    }
    for power in [0.25, 0.5, 0.75, 1.0, 1.25]:
        out[f"split_mean_over_zip_p{power:g}"] = numeric_value(s_mean / np.maximum(numeric_value(s_zip), 1e-6) ** power)
    for floor in [0.01, 0.05, 0.10]:
        clipped = np.maximum(numeric_value(s_zip), floor)
        out[f"split_mean_over_zip_floor{floor:g}"] = numeric_value(s_mean / clipped)
    return {name: np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0) for name, values in out.items()}


def score_abundance(sample_id: int, method: str, calls: pd.DataFrame, raw: np.ndarray) -> dict[str, object]:
    work = calls[["gtdb_species"]].copy()
    work["raw"] = np.where(raw > 0.0, raw, 0.0)
    grouped = work.groupby("gtdb_species", as_index=False)["raw"].sum()
    total = float(grouped["raw"].sum())
    grouped["pred_abundance"] = grouped["raw"] / total if total > 0.0 else 0.0
    grouped["pred_ani"] = 0.0
    grouped = grouped.rename(columns={"raw": "pred_abundance_raw"})
    return taxid_score.score_prediction(sample_id, method, grouped, truth_for_sample(sample_id), {})


def summarize(rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    out = []
    for method, sub in df.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        out.append(
            {
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique()))),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_L1_truth_only_pp": sub["L1_truth_only_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
            }
        )
    return pd.DataFrame(out).sort_values(["mean_L1_union_pp", "method"])


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    wgs_to_species, taxid_to_species, name_to_species, _map_diag = source_score.build_transfer_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(source_score.TAXMAP)

    detail_rows = []
    map_rows = []
    for sample_id in sorted(source_score.READ_MAPPINGS):
        calls = map_minco_calls(sample_id, taxid_to_species, name_to_species, by_accession, by_core, taxmap)
        map_rows.append(
            {
                "sample": sample_id,
                "called_rows_mapped": len(calls),
                "called_species_mapped": calls["gtdb_species"].nunique(),
                "mapping_methods": ";".join(
                    f"{k}:{v}" for k, v in calls["gtdb_mapping_method"].value_counts().sort_index().items()
                ),
            }
        )
        for method, raw in candidate_raws(calls).items():
            detail_rows.append(score_abundance(sample_id, method, calls, raw))

    detail = pd.DataFrame(detail_rows)
    summary = summarize(detail)
    detail.to_csv(RESULTS / "cami3_source_readmap_abundance_sweep_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "cami3_source_readmap_abundance_sweep_summary.tsv", sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(RESULTS / "cami3_source_readmap_abundance_sweep_mapping.tsv", sep="\t", index=False)
    print(summary.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
