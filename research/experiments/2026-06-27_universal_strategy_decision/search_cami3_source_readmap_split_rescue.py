#!/usr/bin/env python3
"""Search CAMI3 source-readmap high-depth split rescues.

This is an exploratory diagnostic. It uses CAMI3 source-readmap GTDB truth to
test whether high-depth split candidates below the calibrated probability gate
can repair high-abundance false negatives. Results must not be promoted without
cross-panel validation.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import score_cami3_gtdb_source_readmap as source_score
import score_cami3_gtdb_taxid_transfer as taxid_score
import sweep_cami3_source_readmap_abundance as abundance_sweep


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def genus(species: str) -> str:
    clean = str(species).replace("s__", "").strip()
    parts = clean.split()
    if not parts:
        return ""
    return parts[0]


def topn_by_genus(df: pd.DataFrame, mask: np.ndarray, score: np.ndarray, topn: int) -> np.ndarray:
    if topn <= 0:
        return mask
    idx = np.flatnonzero(mask)
    keep = np.zeros(len(df), dtype=bool)
    if len(idx) == 0:
        return keep
    work = pd.DataFrame(
        {
            "idx": idx,
            "genus": [genus(df.iloc[i]["gtdb_species"]) for i in idx],
            "score": score[idx],
        }
    ).sort_values(["genus", "score"], ascending=[True, False])
    selected = work.groupby("genus", as_index=False).head(topn)
    keep[selected["idx"].to_numpy(dtype=int)] = True
    return keep


def load_candidates(sample_id: int, taxid_to_species, name_to_species, by_accession, by_core, taxmap) -> pd.DataFrame:
    raw = pd.read_csv(source_score.PROFILE_PATHS[sample_id]["minco"], sep="\t")
    best_ref_species, _diag = source_score.best_raw_ref_species_by_taxid(
        source_score.RAW_TABLES[sample_id],
        taxmap,
        by_accession,
        by_core,
    )
    species = []
    method = []
    for row in raw.itertuples(index=False):
        taxid = str(getattr(row, "taxid", ""))
        species_name = source_score.normalize_name(getattr(row, "species_name", ""))
        gtdb_name = best_ref_species.get(taxid, "")
        map_method = "raw_best_ref" if gtdb_name else ""
        if not gtdb_name:
            gtdb_name = taxid_to_species.get(taxid, "")
            map_method = "unique_taxid" if gtdb_name else ""
        if not gtdb_name and species_name:
            gtdb_name = name_to_species.get(species_name, "")
            map_method = "unique_name" if gtdb_name else ""
        species.append(gtdb_name)
        method.append(map_method or "unmapped")
    raw["gtdb_species"] = species
    raw["gtdb_mapping_method"] = method
    return raw.loc[raw["gtdb_species"].astype(bool)].copy()


def score_sample(sample_id: int, method: str, df: pd.DataFrame, call_mask: np.ndarray, abundance_power: float) -> dict[str, object]:
    called = df.loc[call_mask].copy()
    s_mean = numeric(called, "s_Ref_mean_depth_max")
    s_zip = numeric(called, "s_Ref_zip_af_max")
    raw = s_mean.to_numpy(dtype=float) / np.maximum(s_zip.to_numpy(dtype=float), 1e-6) ** abundance_power
    pred = pd.DataFrame({"gtdb_species": called["gtdb_species"].astype(str), "raw": np.nan_to_num(raw)})
    pred["raw"] = np.where(pred["raw"] > 0.0, pred["raw"], 0.0)
    pred = pred.groupby("gtdb_species", as_index=False)["raw"].sum()
    total = float(pred["raw"].sum())
    pred["pred_abundance"] = pred["raw"] / total if total > 0.0 else 0.0
    pred["pred_ani"] = 0.0
    pred = pred.rename(columns={"raw": "pred_abundance_raw"})
    return taxid_score.score_prediction(sample_id, method, pred, abundance_sweep.truth_for_sample(sample_id), {})


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, sub in detail.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique()))),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_added": sub["added"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(["mean_F1", "mean_L1_union_pp"], ascending=[False, True])


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    _wgs, taxid_to_species, name_to_species, _diag = source_score.build_transfer_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(source_score.TAXMAP)

    candidates = {
        sample_id: load_candidates(sample_id, taxid_to_species, name_to_species, by_accession, by_core, taxmap)
        for sample_id in sorted(source_score.READ_MAPPINGS)
    }

    detail_rows = []
    added_rows = []
    configs = [
        ("baseline_current_p1", None),
        ("baseline_split_p025", None),
    ]
    for p in [0.0, 0.02, 0.05, 0.10, 0.20, 0.25]:
        for sx in [300.0, 500.0, 650.0]:
            for sani in [0.93, 0.95]:
                for saf in [0.30, 0.50]:
                    for breadth in [0.20, 0.30, 0.35, 0.45]:
                        for depth in [1.0, 2.0, 5.0, 10.0]:
                            for topn in [0, 1, 2]:
                                name = f"split_rescue_p{p:g}_sx{sx:g}_ani{sani:g}_af{saf:g}_b{breadth:g}_d{depth:g}_top{topn}"
                                configs.append(
                                    (
                                        name,
                                        {
                                            "p": p,
                                            "sx": sx,
                                            "sani": sani,
                                            "saf": saf,
                                            "breadth": breadth,
                                            "depth": depth,
                                            "topn": topn,
                                        },
                                    )
                                )

    for name, cfg in configs:
        abundance_power = 1.0 if name == "baseline_current_p1" else 0.25
        for sample_id, df in candidates.items():
            base = df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()
            if cfg is None:
                mask = base.copy()
                added = np.zeros(len(df), dtype=bool)
            else:
                prob = numeric(df, "calibrated_probability").to_numpy(dtype=float)
                sx = numeric(df, "s_XnY_ctx_max").to_numpy(dtype=float)
                sani = numeric(df, "s_ANI_max").to_numpy(dtype=float)
                saf = numeric(df, "s_Real_min_align_fraction_max").to_numpy(dtype=float)
                breadth = numeric(df, "s_Ref_breadth_max").to_numpy(dtype=float)
                depth = numeric(df, "s_Ref_mean_depth_max").to_numpy(dtype=float)
                raw_rescue = (
                    (~base)
                    & (prob >= cfg["p"])
                    & (sx >= cfg["sx"])
                    & (sani >= cfg["sani"])
                    & (saf >= cfg["saf"])
                    & (breadth >= cfg["breadth"])
                    & (depth >= cfg["depth"])
                )
                score = depth * np.maximum(breadth, 1e-6) * np.maximum(sani, 0.0) * np.maximum(sx, 1.0)
                added = topn_by_genus(df, raw_rescue, score, int(cfg["topn"]))
                mask = base | added
            row = score_sample(sample_id, name, df, mask, abundance_power)
            row["added"] = int(np.sum(added))
            detail_rows.append(row)
            if np.any(added):
                for rec in df.loc[added].itertuples(index=False):
                    added_rows.append(
                        {
                            "method": name,
                            "sample": sample_id,
                            "gtdb_species": getattr(rec, "gtdb_species"),
                            "prob": getattr(rec, "calibrated_probability", ""),
                            "s_XnY_ctx_max": getattr(rec, "s_XnY_ctx_max", ""),
                            "s_ANI_max": getattr(rec, "s_ANI_max", ""),
                            "s_Ref_breadth_max": getattr(rec, "s_Ref_breadth_max", ""),
                            "s_Ref_mean_depth_max": getattr(rec, "s_Ref_mean_depth_max", ""),
                        }
                    )

    detail = pd.DataFrame(detail_rows)
    summary = summarize(detail)
    detail.to_csv(RESULTS / "cami3_source_readmap_split_rescue_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "cami3_source_readmap_split_rescue_summary.tsv", sep="\t", index=False)
    pd.DataFrame(added_rows).to_csv(RESULTS / "cami3_source_readmap_split_rescue_added.tsv", sep="\t", index=False)
    print(summary.head(30).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
