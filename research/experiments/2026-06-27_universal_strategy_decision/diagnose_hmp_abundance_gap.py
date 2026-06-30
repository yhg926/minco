#!/usr/bin/env python3
"""Diagnose the same-release HMP r232 MinCO/Sylph abundance gap."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

import score_cami3_gtdb_taxid_transfer as taxid_score
import score_hmp_gtdb_source_abundance as hmp_score


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

MINCO_PROFILES = {
    6: Path("/tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv"),
    11: Path("/tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv"),
}


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def safe_div(num: pd.Series, den: pd.Series, power: float = 1.0) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        out = numeric_value(num) / np.maximum(numeric_value(den), 1e-6) ** power
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def numeric_value(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)


def map_minco_calls(
    sample_id: int,
    path: Path,
    taxmap,
    by_accession,
    by_core,
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[call].copy()
    best_ref_species, best_ref_diag = hmp_score.best_raw_ref_species_by_taxid(
        {
            "unique": hmp_score.SAMPLES[sample_id]["unique"],
            "split": hmp_score.SAMPLES[sample_id]["split"],
        },
        taxmap,
        by_accession,
        by_core,
    )

    mapped_species: list[str] = []
    mapping_methods: list[str] = []
    for row in selected.itertuples(index=False):
        species = ""
        method = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species, lookup_method = hmp_score.gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if species:
                    method = f"{col}:{lookup_method}"
                    break
        if not species:
            species = best_ref_species.get(str(getattr(row, "taxid", "")), "")
            if species:
                method = "raw_best_ref"
        mapped_species.append(species)
        mapping_methods.append(method or "unmapped")

    selected["gtdb_species"] = mapped_species
    selected["gtdb_mapping_method"] = mapping_methods
    mapped = selected.loc[selected["gtdb_species"].astype(bool)].copy()
    diag = {
        "sample": sample_id,
        "profile": str(path),
        "called_rows": int(len(selected)),
        "mapped_called_rows": int(len(mapped)),
        "mapped_called_species": int(mapped["gtdb_species"].nunique()),
        "mapping_methods": ";".join(
            f"{k}:{v}" for k, v in mapped["gtdb_mapping_method"].value_counts().sort_index().items()
        ),
        **best_ref_diag,
    }
    return mapped, diag


def genus_from_gtdb_species(value: object) -> str:
    text = str(value or "")
    if text.startswith("s__"):
        text = text[3:]
    text = text.strip()
    return text.split()[0] if text else ""


def genus_reallocated_raw(calls: pd.DataFrame, base_raw: np.ndarray, quality: np.ndarray) -> np.ndarray:
    """Preserve each genus' total raw mass while reweighting species inside it."""

    work = calls[["gtdb_species"]].copy()
    work["genus"] = work["gtdb_species"].map(genus_from_gtdb_species)
    work["base_raw"] = np.where(np.asarray(base_raw, dtype=float) > 0.0, base_raw, 0.0)
    work["quality"] = np.where(np.asarray(quality, dtype=float) > 0.0, quality, 0.0)
    out = work["base_raw"].to_numpy(dtype=float).copy()
    for _genus, idx in work.groupby("genus").groups.items():
        idx_array = np.asarray(list(idx), dtype=int)
        if len(idx_array) <= 1:
            continue
        total = float(work.loc[idx_array, "base_raw"].sum())
        adjusted = work.loc[idx_array, "base_raw"].to_numpy(dtype=float) * work.loc[idx_array, "quality"].to_numpy(dtype=float)
        adjusted_sum = float(np.sum(adjusted))
        if total <= 0.0 or adjusted_sum <= 0.0:
            continue
        out[idx_array] = adjusted / adjusted_sum * total
    return out


def candidate_raws(df: pd.DataFrame) -> dict[str, np.ndarray]:
    s_mean = numeric(df, "s_Ref_mean_depth_max")
    s_hit = numeric(df, "s_Ref_hit_mean_depth_max")
    s_breadth = numeric(df, "s_Ref_breadth_max")
    s_zip = numeric(df, "s_Ref_zip_af_max")
    s_norm = numeric(df, "s_Normalized_abundance_depth_max")
    s_rel = numeric(df, "s_Relative_abundance_depth_max")
    u_mean = numeric(df, "u_Ref_mean_depth_max")
    u_hit = numeric(df, "u_Ref_hit_mean_depth_max")
    u_breadth = numeric(df, "u_Ref_breadth_max")
    u_zip = numeric(df, "u_Ref_zip_af_max")
    u_norm = numeric(df, "u_Normalized_abundance_depth_max")
    u_rel = numeric(df, "u_Relative_abundance_depth_max")
    prob = numeric(df, "calibrated_probability", 1.0)
    s_xny = numeric(df, "s_XnY_ctx_max")
    s_mut = numeric(df, "s_N_mut2_ctx_max") / np.maximum(s_xny, 1.0)
    s_diff = numeric(df, "s_N_diff_obj_max") / np.maximum(s_xny, 1.0)

    s_zip1 = safe_div(s_mean, s_zip, 1.0)
    u_zip1 = safe_div(u_mean, u_zip, 1.0)
    current_raw = numeric_value(numeric(df, "calibrated_abundance_raw"))
    xny_quality = np.clip(numeric_value(s_xny) / 1000.0, 0.0, 1.0)
    breadth_quality = np.clip(numeric_value(s_breadth), 0.0, 1.0)
    mutdiff_light_quality = np.exp(-numeric_value(s_mut) - 0.5 * numeric_value(s_diff))
    mutdiff_strong_quality = np.exp(-4.0 * numeric_value(s_mut) - 2.0 * numeric_value(s_diff))
    out: dict[str, np.ndarray] = {
        "minco_current_calibrated_abundance": numeric_value(numeric(df, "calibrated_abundance")),
        "minco_current_calibrated_raw": current_raw,
        "minco_split_mean_depth": numeric_value(s_mean),
        "minco_unique_mean_depth": numeric_value(u_mean),
        "minco_split_hit_mean_x_breadth": numeric_value(s_hit * s_breadth),
        "minco_unique_hit_mean_x_breadth": numeric_value(u_hit * u_breadth),
        "minco_split_normalized_depth": numeric_value(s_norm),
        "minco_split_relative_depth": numeric_value(s_rel),
        "minco_unique_normalized_depth": numeric_value(u_norm),
        "minco_unique_relative_depth": numeric_value(u_rel),
        "minco_max_split_unique_zip_p1": np.maximum(s_zip1, u_zip1),
        "minco_split_unique_zip_geom_p1": np.sqrt(np.maximum(s_zip1, 0.0) * np.maximum(u_zip1, 0.0)),
        "minco_probability_x_current_raw": current_raw * numeric_value(prob),
        "minco_quality_xny_current_raw": current_raw * xny_quality,
        "minco_quality_breadth_current_raw": current_raw * breadth_quality,
        "minco_quality_xny_breadth_current_raw": current_raw * xny_quality * breadth_quality,
        "minco_quality_mutdiff_light_current_raw": current_raw * mutdiff_light_quality,
        "minco_quality_mutdiff_strong_current_raw": current_raw * mutdiff_strong_quality,
    }
    genus_reallocation_qualities = {
        "current_raw_power_p0.75": np.power(np.maximum(current_raw, 0.0), -0.25),
        "current_raw_power_p1.25": np.power(np.maximum(current_raw, 0.0), 0.25),
        "xny": xny_quality,
        "breadth": breadth_quality,
        "xny_breadth": xny_quality * breadth_quality,
        "mutdiff_light": mutdiff_light_quality,
        "mutdiff_strong": mutdiff_strong_quality,
    }
    for name, quality in genus_reallocation_qualities.items():
        out[f"minco_genus_realloc_{name}"] = genus_reallocated_raw(df, current_raw, quality)
    for power in [1.25, 1.5]:
        powered_raw = np.power(np.maximum(current_raw, 0.0), power)
        out[f"minco_raw_p{power:g}_genus_realloc_mutdiff_light"] = genus_reallocated_raw(
            df,
            powered_raw,
            mutdiff_light_quality,
        )
        out[f"minco_raw_p{power:g}_genus_realloc_mutdiff_strong"] = genus_reallocated_raw(
            df,
            powered_raw,
            mutdiff_strong_quality,
        )
        out[f"minco_raw_p{power:g}_genus_realloc_xny"] = genus_reallocated_raw(
            df,
            powered_raw,
            xny_quality,
        )
    for power in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]:
        out[f"minco_split_mean_over_zip_p{power:g}"] = safe_div(s_mean, s_zip, power)
        out[f"minco_unique_mean_over_zip_p{power:g}"] = safe_div(u_mean, u_zip, power)
        out[f"minco_current_raw_power_p{power:g}"] = np.power(np.maximum(current_raw, 0.0), power)
    return {name: np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0) for name, values in out.items()}


def prediction_from_raw(calls: pd.DataFrame, raw_values: np.ndarray) -> pd.DataFrame:
    work = calls[["gtdb_species"]].copy()
    work["raw"] = np.where(np.asarray(raw_values, dtype=float) > 0.0, raw_values, 0.0)
    if "reported_ani" in calls.columns:
        work["ani"] = numeric(calls, "reported_ani")
    else:
        work["ani"] = 0.0
    grouped = work.groupby("gtdb_species", as_index=False).agg(
        pred_abundance_raw=("raw", "sum"),
        pred_ani=("ani", "max"),
    )
    total = float(grouped["pred_abundance_raw"].sum())
    grouped["pred_abundance"] = grouped["pred_abundance_raw"] / total if total > 0.0 else 0.0
    return grouped


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, sub in scores.groupby("method"):
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
                "mean_L1_truth_only_pp": sub["L1_truth_only_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(["mean_L1_union_pp", "method"])


def species_rows_for_sample(
    sample_id: int,
    truth_df: pd.DataFrame,
    predictions: Mapping[str, pd.DataFrame],
    keep_methods: list[str],
) -> list[dict[str, object]]:
    truth_abund = dict(zip(truth_df["gtdb_species"].astype(str), numeric(truth_df, "truth_abundance")))
    species = set(truth_abund)
    for method in keep_methods:
        pred = predictions[method]
        species.update(pred["gtdb_species"].astype(str))

    rows = []
    for name in sorted(species):
        rec: dict[str, object] = {
            "sample": sample_id,
            "gtdb_species": name,
            "truth_abundance": truth_abund.get(name, 0.0),
        }
        for method in keep_methods:
            pred = predictions[method]
            pred_abund = dict(zip(pred["gtdb_species"].astype(str), numeric(pred, "pred_abundance")))
            value = pred_abund.get(name, 0.0)
            rec[method] = value
            rec[f"{method}_abs_error_pp"] = abs(value - rec["truth_abundance"]) * 100.0
        rows.append(rec)
    return rows


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    for path in MINCO_PROFILES.values():
        if not path.exists():
            raise SystemExit(f"missing MinCO profile: {path}")

    by_accession, by_core = hmp_score.truth.load_gtdb_metadata(hmp_score.truth.GTDB_METADATA)
    source_accessions = hmp_score.load_source_accessions(hmp_score.TRUTH_ROOT / "genome_to_id.tsv")
    taxmap = hmp_score.parse_species_taxmap(hmp_score.TAXMAP)

    score_rows = []
    mapping_rows = []
    species_rows = []
    top_error_rows = []
    for sample_id in sorted(hmp_score.SAMPLES):
        truth_df, _quality, _source = hmp_score.build_source_truth(sample_id, source_accessions, by_accession, by_core)
        calls, mapping_diag = map_minco_calls(sample_id, MINCO_PROFILES[sample_id], taxmap, by_accession, by_core)
        mapping_rows.append(mapping_diag)

        predictions: dict[str, pd.DataFrame] = {}
        for method, raw in candidate_raws(calls).items():
            pred = prediction_from_raw(calls, raw)
            predictions[method] = pred
            score_rows.append(taxid_score.score_prediction(sample_id, method, pred, truth_df, {}))

        sylph_profile = hmp_score.R232_RUN / f"sylph_sample{sample_id}/profile.tsv"
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(sylph_profile, by_accession, by_core)
        predictions["sylph_r232_taxonomic_abundance"] = sylph_pred
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_r232_taxonomic_abundance",
                sylph_pred,
                truth_df,
                {"pred_rows_called": sylph_extra.get("pred_rows_called", 0)},
            )
        )

        keep_methods = [
            "minco_current_calibrated_abundance",
            "minco_max_split_unique_zip_p1",
            "minco_split_mean_depth",
            "minco_split_mean_over_zip_p0.25",
            "minco_current_raw_power_p0.5",
            "sylph_r232_taxonomic_abundance",
        ]
        sample_species_rows = species_rows_for_sample(sample_id, truth_df, predictions, keep_methods)
        species_rows.extend(sample_species_rows)
        for method in ["minco_current_calibrated_abundance", "sylph_r232_taxonomic_abundance"]:
            key = f"{method}_abs_error_pp"
            for row in sorted(sample_species_rows, key=lambda rec: rec[key], reverse=True)[:20]:
                top_error_rows.append(
                    {
                        "sample": sample_id,
                        "method": method,
                        "gtdb_species": row["gtdb_species"],
                        "truth_abundance": row["truth_abundance"],
                        "pred_abundance": row[method],
                        "abs_error_pp": row[key],
                    }
                )

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    mapping = pd.DataFrame(mapping_rows)
    species = pd.DataFrame(species_rows)
    top_errors = pd.DataFrame(top_error_rows).sort_values(["sample", "method", "abs_error_pp"], ascending=[True, True, False])

    scores.to_csv(RESULTS / "hmp_r232_abundance_gap_variant_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "hmp_r232_abundance_gap_variant_summary.tsv", sep="\t", index=False)
    mapping.to_csv(RESULTS / "hmp_r232_abundance_gap_mapping.tsv", sep="\t", index=False)
    species.to_csv(RESULTS / "hmp_r232_abundance_gap_species.tsv", sep="\t", index=False)
    top_errors.to_csv(RESULTS / "hmp_r232_abundance_gap_top_errors.tsv", sep="\t", index=False)
    print(summary.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
