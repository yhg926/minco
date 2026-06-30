#!/usr/bin/env python3
"""Sweep fixed-call abundance formulas across cached GTDB panels.

This is a post-hoc diagnostic: calls are held fixed and only predicted
abundance mass is recomputed. A variant is promotable only if it improves the
current abundance default across panels without relying on panel-specific
selection.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

import decompose_abundance_errors as decomp
import score_cami3_gtdb_source_readmap as cami3
import score_cami3_gtdb_taxid_transfer as taxid_score
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def arr(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(dtype=float)


def safe_div(num: pd.Series, den: pd.Series, power: float = 1.0) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        out = arr(num) / np.maximum(arr(den), 1e-6) ** power
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def genus_from_species(value: object) -> str:
    text = str(value or "")
    if text.startswith("s__"):
        text = text[3:]
    text = text.strip()
    return text.split()[0] if text else ""


def genus_reallocated_raw(calls: pd.DataFrame, base_raw: np.ndarray, quality: np.ndarray) -> np.ndarray:
    """Preserve each genus' raw mass while reweighting called rows inside it."""

    genus_source = "profile_species_name" if "profile_species_name" in calls.columns else "gtdb_species"
    work = calls[[genus_source]].copy()
    work["genus"] = work[genus_source].map(genus_from_species)
    work["base_raw"] = np.maximum(np.asarray(base_raw, dtype=float), 0.0)
    work["quality"] = np.maximum(np.asarray(quality, dtype=float), 0.0)
    out = work["base_raw"].to_numpy(dtype=float).copy()
    for _genus, idx in work.groupby("genus").groups.items():
        idx_array = np.asarray(list(idx), dtype=int)
        if len(idx_array) <= 1:
            continue
        total = float(work.loc[idx_array, "base_raw"].sum())
        adjusted = (
            work.loc[idx_array, "base_raw"].to_numpy(dtype=float)
            * work.loc[idx_array, "quality"].to_numpy(dtype=float)
        )
        adjusted_sum = float(adjusted.sum())
        if total <= 0.0 or adjusted_sum <= 0.0:
            continue
        out[idx_array] = adjusted / adjusted_sum * total
    return out


def candidate_raws(calls: pd.DataFrame) -> dict[str, np.ndarray]:
    s_mean = numeric(calls, "s_Ref_mean_depth_max")
    s_hit = numeric(calls, "s_Ref_hit_mean_depth_max")
    s_breadth = numeric(calls, "s_Ref_breadth_max")
    s_zip = numeric(calls, "s_Ref_zip_af_max")
    s_norm = numeric(calls, "s_Normalized_abundance_depth_max")
    s_rel = numeric(calls, "s_Relative_abundance_depth_max")
    s_reads = numeric(calls, "s_Reads_with_ctx_match_max")
    s_read_frac = numeric(calls, "s_Read_match_fraction_max")
    s_blocks = numeric(calls, "s_Blocks_with_ctx_match_max")
    s_block_frac = numeric(calls, "s_Block_match_fraction_max")
    s_xny = numeric(calls, "s_XnY_ctx_max")
    s_mut = numeric(calls, "s_N_mut2_ctx_max") / np.maximum(s_xny, 1.0)
    s_diff = numeric(calls, "s_N_diff_obj_max") / np.maximum(s_xny, 1.0)
    u_mean = numeric(calls, "u_Ref_mean_depth_max")
    u_hit = numeric(calls, "u_Ref_hit_mean_depth_max")
    u_breadth = numeric(calls, "u_Ref_breadth_max")
    u_zip = numeric(calls, "u_Ref_zip_af_max")
    u_norm = numeric(calls, "u_Normalized_abundance_depth_max")
    u_rel = numeric(calls, "u_Relative_abundance_depth_max")
    u_reads = numeric(calls, "u_Reads_with_ctx_match_max")
    u_read_frac = numeric(calls, "u_Read_match_fraction_max")
    u_blocks = numeric(calls, "u_Blocks_with_ctx_match_max")
    u_block_frac = numeric(calls, "u_Block_match_fraction_max")
    prob = numeric(calls, "calibrated_probability", 1.0)
    current_raw = arr(numeric(calls, "calibrated_abundance_raw"))

    out: dict[str, np.ndarray] = {
        "current_calibrated_abundance": arr(numeric(calls, "calibrated_abundance")),
        "current_calibrated_raw": current_raw,
        "split_mean_depth": arr(s_mean),
        "unique_mean_depth": arr(u_mean),
        "split_hit_mean_x_breadth": arr(s_hit * s_breadth),
        "unique_hit_mean_x_breadth": arr(u_hit * u_breadth),
        "split_normalized_depth": arr(s_norm),
        "unique_normalized_depth": arr(u_norm),
        "split_relative_depth": arr(s_rel),
        "unique_relative_depth": arr(u_rel),
        "split_reads_with_ctx": arr(s_reads),
        "unique_reads_with_ctx": arr(u_reads),
        "max_split_unique_reads_with_ctx": np.maximum(arr(s_reads), arr(u_reads)),
        "mean_split_unique_reads_with_ctx": (arr(s_reads) + arr(u_reads)) / 2.0,
        "split_read_match_fraction": arr(s_read_frac),
        "unique_read_match_fraction": arr(u_read_frac),
        "max_split_unique_read_match_fraction": np.maximum(arr(s_read_frac), arr(u_read_frac)),
        "split_blocks_with_ctx": arr(s_blocks),
        "unique_blocks_with_ctx": arr(u_blocks),
        "max_split_unique_blocks_with_ctx": np.maximum(arr(s_blocks), arr(u_blocks)),
        "mean_split_unique_blocks_with_ctx": (arr(s_blocks) + arr(u_blocks)) / 2.0,
        "split_block_match_fraction": arr(s_block_frac),
        "unique_block_match_fraction": arr(u_block_frac),
        "max_split_unique_block_match_fraction": np.maximum(arr(s_block_frac), arr(u_block_frac)),
        "current_raw_x_probability": current_raw * arr(prob),
        "current_raw_x_probability2": current_raw * arr(prob) * arr(prob),
    }

    for power in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
        split = safe_div(s_mean, s_zip, power)
        unique = safe_div(u_mean, u_zip, power)
        split_reads = safe_div(s_reads, s_zip, power)
        unique_reads = safe_div(u_reads, u_zip, power)
        out[f"split_mean_over_zip_p{power:g}"] = split
        out[f"unique_mean_over_zip_p{power:g}"] = unique
        out[f"max_split_unique_zip_p{power:g}"] = np.maximum(split, unique)
        out[f"mean_split_unique_zip_p{power:g}"] = (split + unique) / 2.0
        out[f"split_reads_over_zip_p{power:g}"] = split_reads
        out[f"unique_reads_over_zip_p{power:g}"] = unique_reads
        out[f"max_split_unique_reads_zip_p{power:g}"] = np.maximum(split_reads, unique_reads)
        out[f"mean_split_unique_reads_zip_p{power:g}"] = (split_reads + unique_reads) / 2.0
        out[f"current_raw_power_p{power:g}"] = np.power(np.maximum(current_raw, 0.0), power)

    xny_quality = np.clip(arr(s_xny) / 1000.0, 0.0, 1.0)
    breadth_quality = np.clip(arr(s_breadth), 0.0, 1.0)
    prob_quality = np.clip(arr(prob), 0.0, 1.0)
    mutdiff_light = np.exp(-arr(s_mut) - 0.5 * arr(s_diff))
    mutdiff_strong = np.exp(-4.0 * arr(s_mut) - 2.0 * arr(s_diff))
    qualities = {
        "xny": xny_quality,
        "breadth": breadth_quality,
        "probability": prob_quality,
        "xny_breadth": xny_quality * breadth_quality,
        "prob_breadth": prob_quality * breadth_quality,
        "mutdiff_light": mutdiff_light,
        "mutdiff_strong": mutdiff_strong,
    }
    for name, quality in qualities.items():
        out[f"current_raw_x_{name}"] = current_raw * quality
        out[f"genus_realloc_{name}"] = genus_reallocated_raw(calls, current_raw, quality)
    for power in [0.75, 1.25, 1.5]:
        powered = np.power(np.maximum(current_raw, 0.0), power)
        out[f"raw_p{power:g}_genus_realloc_xny"] = genus_reallocated_raw(calls, powered, xny_quality)
        out[f"raw_p{power:g}_genus_realloc_breadth"] = genus_reallocated_raw(calls, powered, breadth_quality)
        out[f"raw_p{power:g}_genus_realloc_mutdiff_light"] = genus_reallocated_raw(calls, powered, mutdiff_light)

    blend_targets = [
        "genus_realloc_mutdiff_light",
        "genus_realloc_xny",
        "unique_mean_over_zip_p1",
        "mean_split_unique_zip_p1",
        "mean_split_unique_reads_with_ctx",
        "mean_split_unique_blocks_with_ctx",
        "unique_reads_over_zip_p1",
    ]
    for target in blend_targets:
        target_raw = np.asarray(out.get(target, np.zeros_like(current_raw)), dtype=float)
        for alpha in [0.1, 0.25, 0.5, 0.75]:
            out[f"blend_current_{target}_a{alpha:g}"] = (
                (1.0 - alpha) * current_raw + alpha * target_raw
            )

    return {name: np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0) for name, value in out.items()}


def collapse_prediction(calls: pd.DataFrame, raw_values: np.ndarray, collapse: str) -> pd.DataFrame:
    work = calls[["gtdb_species"]].copy()
    work["raw"] = np.maximum(np.asarray(raw_values, dtype=float), 0.0)
    work = work.loc[work["gtdb_species"].fillna("").astype(str).astype(bool)].copy()
    if collapse == "max":
        grouped = work.groupby("gtdb_species", as_index=False)["raw"].max()
    elif collapse == "sum":
        grouped = work.groupby("gtdb_species", as_index=False)["raw"].sum()
    else:
        raise ValueError(f"unsupported collapse rule: {collapse}")
    total = float(grouped["raw"].sum())
    grouped["pred_abundance"] = grouped["raw"] / total if total > 0.0 else 0.0
    grouped["pred_ani"] = 0.0
    return grouped.rename(columns={"raw": "pred_abundance_raw"})


def toy_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / f"mouse{sample}_gtdb_species_profile.tsv", sep="\t")
    return truth.rename(columns={"relative_abundance": "truth_abundance"})[
        ["gtdb_species", "truth_abundance"]
    ]


def load_toy_calls(sample: int) -> tuple[pd.DataFrame, str, Path]:
    path = Path(f"/tmp/minco_current_code_toymouse_refresh_20260627/sample{sample}_current.tsv")
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    calls = raw.loc[call].copy()
    calls["profile_species_name"] = calls["species_name"].fillna("").astype(str)
    calls["gtdb_species"] = calls["species_name"].fillna("").astype(str)
    return calls.loc[calls["gtdb_species"].astype(bool)].copy(), "max", path


def load_hmp_calls(
    sample: int,
    taxmap: Mapping[str, Mapping[str, str]],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, str, Path]:
    alt_root = Path("/mnt/new3T/minco_release_holdouts_20260628")
    raw_tables = {
        18: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin18_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        21: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin21_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        13: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin13_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin13_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        17: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin17_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin17_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        16: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin16_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin16_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        15: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin15_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin15_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        4: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin4_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin4_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        23: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin23_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin23_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        20: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin20_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin20_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        7: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin7_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin7_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        9: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin9_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin9_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        10: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin10_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin10_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        19: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin19_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin19_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        14: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin14_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin14_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        24: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin24_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin24_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        25: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin25_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin25_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        3: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin3_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin3_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
        1: {
            "unique": alt_root
            / "minco_current_code_hmp_airskin1_20260627/work/minco.best_diff_unique.unfiltered.tsv",
            "split": alt_root
            / "minco_current_code_hmp_airskin1_20260627/work/minco.best_diff_split.unfiltered.tsv",
        },
    }
    path_map = {
        0: "/tmp/minco_current_code_hmp_airskin0_20260627/minco_sample0_current_default.tsv",
        1: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/minco_sample1_current_default.tsv",
        3: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/minco_sample3_current_default.tsv",
        4: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin4_20260627/minco_sample4_current_default.tsv",
        5: "/tmp/minco_current_code_hmp_airskin5_20260627/minco_sample5_current_default.tsv",
        6: "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv",
        7: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin7_20260627/minco_sample7_current_default.tsv",
        9: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin9_20260627/minco_sample9_current_default.tsv",
        10: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin10_20260627/minco_sample10_current_default.tsv",
        11: "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv",
        13: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/minco_sample13_current_default.tsv",
        14: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin14_20260627/minco_sample14_current_default.tsv",
        15: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin15_20260627/minco_sample15_current_default.tsv",
        16: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/minco_sample16_current_default.tsv",
        17: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/minco_sample17_current_default.tsv",
        18: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/minco_sample18_current_default.tsv",
        19: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin19_20260627/minco_sample19_current_default.tsv",
        20: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin20_20260627/minco_sample20_current_default.tsv",
        21: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv",
        22: "/tmp/minco_current_code_hmp_airskin22_20260627/minco_sample22_current_default.tsv",
        23: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin23_20260627/minco_sample23_current_default.tsv",
        24: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/minco_sample24_current_default.tsv",
        25: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/minco_sample25_current_default.tsv",
        28: "/tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv",
    }
    if sample in path_map:
        path = Path(path_map[sample])
    else:
        hmp.ensure_sample_record(int(sample))
        path = Path(hmp.SAMPLES[int(sample)]["minco"])
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    calls = raw.loc[call].copy()
    calls["profile_species_name"] = calls["species_name"].fillna("").astype(str)
    hmp.ensure_sample_record(int(sample))
    selected_raw_tables = raw_tables.get(
        sample,
        {"unique": hmp.SAMPLES[sample]["unique"], "split": hmp.SAMPLES[sample]["split"]},
    )
    best_ref_species, _diag = hmp.best_raw_ref_species_by_taxid(
        selected_raw_tables,
        taxmap,
        by_accession,
        by_core,
    )
    species: list[str] = []
    for row in calls.itertuples(index=False):
        mapped = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                mapped, _method = hmp.gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if mapped:
                    break
        if not mapped:
            mapped = best_ref_species.get(str(getattr(row, "taxid", "")), "")
        species.append(mapped)
    calls["gtdb_species"] = species
    return calls.loc[calls["gtdb_species"].astype(bool)].copy(), "sum", path


def load_hmp_gastro_calls(
    sample: int,
    taxmap: Mapping[str, Mapping[str, str]],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, str, Path]:
    raw_default_paths = {
        0: Path("/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv"),
        6: Path("/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv"),
    }
    path = raw_default_paths.get(int(sample), hmp_gastro.SAMPLES[int(sample)]["minco"])
    if not path.exists():
        path = hmp_gastro.SAMPLES[int(sample)]["minco"]
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    calls = raw.loc[call].copy()
    calls["profile_species_name"] = calls["species_name"].fillna("").astype(str)
    best_ref_species, _diag = hmp_gastro.best_raw_ref_species_by_taxid(
        {
            "unique": hmp_gastro.SAMPLES[int(sample)]["unique"],
            "split": hmp_gastro.SAMPLES[int(sample)]["split"],
        },
        taxmap,
        by_accession,
        by_core,
    )
    species: list[str] = []
    for row in calls.itertuples(index=False):
        mapped = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                mapped, _method = hmp_gastro.gtdb_from_accession(
                    getattr(row, col, ""),
                    by_accession,
                    by_core,
                )
                if mapped:
                    break
        if not mapped:
            mapped = best_ref_species.get(str(getattr(row, "taxid", "")), "")
        species.append(mapped)
    calls["gtdb_species"] = species
    return calls.loc[calls["gtdb_species"].astype(bool)].copy(), "sum", path


def load_cami3_calls(
    sample: int,
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
    taxmap: Mapping[str, Mapping[str, str]],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, str, Path]:
    path = Path(
        {
            0: "/tmp/cami3_toy_human_gut_20260626/run/sample0_universal_autoexact.tsv",
            1: "/tmp/cami3_toy_human_gut_20260626/run/sample1_universal_autoexact.tsv",
            2: "/tmp/cami3_toy_human_gut_20260626/run/sample2_universal_autoexact_tail_p025.tsv",
        }[sample]
    )
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    calls = raw.loc[call].copy()
    calls["profile_species_name"] = calls["species_name"].fillna("").astype(str)
    best_ref_species, _diag = cami3.best_raw_ref_species_by_taxid(
        cami3.RAW_TABLES[sample],
        taxmap,
        by_accession,
        by_core,
    )
    species: list[str] = []
    for row in calls.itertuples(index=False):
        mapped = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                mapped = cami3.gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if mapped:
                    break
        taxid = str(getattr(row, "taxid", ""))
        species_name = cami3.normalize_name(getattr(row, "species_name", ""))
        if not mapped:
            mapped = best_ref_species.get(taxid, "")
        if not mapped:
            mapped = taxid_to_species.get(taxid, "")
        if not mapped and species_name:
            mapped = name_to_species.get(species_name, "")
        species.append(mapped)
    calls["gtdb_species"] = species
    return calls.loc[calls["gtdb_species"].astype(bool)].copy(), "sum", path


def hmp_truth(sample: int) -> pd.DataFrame:
    single_sample = RESULTS / f"hmp_airskin{sample}_r232_source_abundance_truth.tsv"
    if single_sample.exists():
        truth = pd.read_csv(single_sample, sep="\t")
    else:
        truth = pd.read_csv(
            RESULTS / "hmp_current_refresh_r232_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28_truth.tsv",
            sep="\t",
        )
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def hmp_gastro_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "hmp_gastrooral_r232_source_abundance_truth.tsv", sep="\t")
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def cami3_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "cami3_gtdb_source_readmap_truth.tsv", sep="\t")
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def summarize(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for (panel, method), sub in scores.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        official_l1_col = "L1_truth_only_pp" if panel == "cami2_toy_mouse_gut" else "L1_union_pp"
        official_pearson_col = "Pearson_truth_only" if panel == "cami2_toy_mouse_gut" else "Pearson_union"
        rows.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique(), key=str))),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_L1_truth_only_pp": sub["L1_truth_only_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
                "official_L1_pp": sub[official_l1_col].mean(),
                "official_Pearson": sub[official_pearson_col].mean(),
            }
        )
    panel_summary = pd.DataFrame(rows)

    current = panel_summary.loc[panel_summary["method"].eq("current_calibrated_abundance")][
        ["panel", "official_L1_pp", "official_Pearson"]
    ].rename(columns={"official_L1_pp": "current_official_L1_pp", "official_Pearson": "current_official_Pearson"})
    merged = panel_summary.merge(current, on="panel", how="left")
    merged["delta_current_L1_pp"] = merged["official_L1_pp"] - merged["current_official_L1_pp"]
    merged["delta_current_Pearson"] = merged["official_Pearson"] - merged["current_official_Pearson"]

    overall_rows: list[dict[str, object]] = []
    expected_panel_count = int(panel_summary["panel"].nunique())
    for method, sub in merged.groupby("method", sort=True):
        panel_count = int(sub["panel"].nunique())
        improved = int((sub["delta_current_L1_pp"] < -1e-9).sum())
        worsened = int((sub["delta_current_L1_pp"] > 1e-9).sum())
        max_worse = float(sub["delta_current_L1_pp"].max())
        mean_delta = float(sub["delta_current_L1_pp"].mean())
        mean_l1 = float(sub["official_L1_pp"].mean())
        mean_pearson = float(sub["official_Pearson"].mean())
        decision = "diagnostic"
        if method == "current_calibrated_abundance":
            decision = "current_default"
        elif panel_count == expected_panel_count and max_worse <= 0.0 and improved > 0:
            decision = "promotable_by_L1_if_validated_from_raw"
        elif mean_delta < 0.0 and worsened > 0:
            decision = "local_tradeoff_not_default"
        else:
            decision = "worse_or_equal"
        overall_rows.append(
            {
                "method": method,
                "panels": ",".join(sorted(sub["panel"].astype(str).unique())),
                "panel_count": panel_count,
                "mean_official_L1_pp": mean_l1,
                "mean_delta_current_L1_pp": mean_delta,
                "max_worse_current_L1_pp": max_worse,
                "improved_panel_count": improved,
                "worsened_panel_count": worsened,
                "mean_official_Pearson": mean_pearson,
                "mean_delta_current_Pearson": float(sub["delta_current_Pearson"].mean()),
                "decision": decision,
            }
        )
    overall = pd.DataFrame(overall_rows).sort_values(
        ["decision", "mean_delta_current_L1_pp", "max_worse_current_L1_pp", "method"]
    )
    return merged.sort_values(["panel", "official_L1_pp", "method"]), overall


def add_external_sylph(panel_summary: pd.DataFrame) -> pd.DataFrame:
    decomp_summary = pd.read_csv(RESULTS / "abundance_error_decomposition_summary.tsv", sep="\t")
    rows = []
    for row in decomp_summary.itertuples(index=False):
        method = str(getattr(row, "method"))
        if not method.startswith("sylph"):
            continue
        panel = str(getattr(row, "panel"))
        official_l1 = (
            float(getattr(row, "mean_truth_only_L1_pp"))
            if panel == "cami2_toy_mouse_gut"
            else float(getattr(row, "mean_union_L1_pp"))
        )
        official_pearson = (
            float(getattr(row, "mean_Pearson_truth_only"))
            if panel == "cami2_toy_mouse_gut"
            else float(getattr(row, "mean_Pearson_union"))
        )
        rows.append(
            {
                "panel": panel,
                "method": "sylph_external_baseline",
                "samples": getattr(row, "samples"),
                "mean_F1": getattr(row, "mean_F1"),
                "pooled_F1": getattr(row, "pooled_F1"),
                "pooled_TP": getattr(row, "pooled_TP"),
                "pooled_FP": getattr(row, "pooled_FP"),
                "pooled_FN": getattr(row, "pooled_FN"),
                "mean_L1_union_pp": getattr(row, "mean_union_L1_pp"),
                "mean_L1_truth_only_pp": getattr(row, "mean_truth_only_L1_pp"),
                "mean_Pearson_union": getattr(row, "mean_Pearson_union"),
                "mean_Pearson_truth_only": getattr(row, "mean_Pearson_truth_only"),
                "official_L1_pp": official_l1,
                "official_Pearson": official_pearson,
                "current_official_L1_pp": np.nan,
                "current_official_Pearson": np.nan,
                "delta_current_L1_pp": np.nan,
                "delta_current_Pearson": np.nan,
            }
        )
    return pd.concat([panel_summary, pd.DataFrame(rows)], ignore_index=True, sort=False)


def validation(panel_summary: pd.DataFrame) -> pd.DataFrame:
    decomp_summary = pd.read_csv(RESULTS / "abundance_error_decomposition_summary.tsv", sep="\t")
    rows = []
    for panel in sorted(panel_summary["panel"].unique()):
        ours = panel_summary.loc[
            panel_summary["panel"].eq(panel)
            & panel_summary["method"].eq("current_calibrated_abundance")
        ].iloc[0]
        theirs = decomp_summary.loc[
            decomp_summary["panel"].eq(panel)
            & ~decomp_summary["method"].astype(str).str.startswith("sylph")
        ].iloc[0]
        if panel == "cami2_toy_mouse_gut":
            official_l1_col = "mean_truth_only_L1_pp"
            official_pearson_col = "mean_Pearson_truth_only"
        else:
            official_l1_col = "mean_union_L1_pp"
            official_pearson_col = "mean_Pearson_union"
        rows.append(
            {
                "panel": panel,
                "L1_delta_pp": float(ours["official_L1_pp"]) - float(theirs[official_l1_col]),
                "Pearson_delta": float(ours["official_Pearson"]) - float(theirs[official_pearson_col]),
                "F1_delta": float(ours["mean_F1"]) - float(theirs["mean_F1"]),
                "source": "abundance_error_decomposition_summary.tsv",
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)

    toy_mod = decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    hmp_taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
    hmp_gastro_taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)
    _cami3_wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _map_diag = cami3.build_transfer_maps()
    cami3_taxmap = cami3.parse_species_taxmap(cami3.TAXMAP)

    sample_records = [
        ("cami2_toy_mouse_gut", 5),
        ("cami2_toy_mouse_gut", 6),
        ("cami2_toy_mouse_gut", 7),
        ("hmp_airskin_gtdb_source_abundance", 0),
        ("hmp_airskin_gtdb_source_abundance", 1),
        ("hmp_airskin_gtdb_source_abundance", 3),
        ("hmp_airskin_gtdb_source_abundance", 4),
        ("hmp_airskin_gtdb_source_abundance", 5),
        ("hmp_airskin_gtdb_source_abundance", 6),
        ("hmp_airskin_gtdb_source_abundance", 7),
        ("hmp_airskin_gtdb_source_abundance", 9),
        ("hmp_airskin_gtdb_source_abundance", 10),
        ("hmp_airskin_gtdb_source_abundance", 11),
        ("hmp_airskin_gtdb_source_abundance", 13),
        ("hmp_airskin_gtdb_source_abundance", 14),
        ("hmp_airskin_gtdb_source_abundance", 15),
        ("hmp_airskin_gtdb_source_abundance", 16),
        ("hmp_airskin_gtdb_source_abundance", 17),
        ("hmp_airskin_gtdb_source_abundance", 18),
        ("hmp_airskin_gtdb_source_abundance", 19),
        ("hmp_airskin_gtdb_source_abundance", 20),
        ("hmp_airskin_gtdb_source_abundance", 21),
        ("hmp_airskin_gtdb_source_abundance", 22),
        ("hmp_airskin_gtdb_source_abundance", 23),
        ("hmp_airskin_gtdb_source_abundance", 24),
        ("hmp_airskin_gtdb_source_abundance", 25),
        ("hmp_airskin_gtdb_source_abundance", 28),
        ("hmp_gastrooral_gtdb_source_abundance", 0),
        ("hmp_gastrooral_gtdb_source_abundance", 6),
        ("cami3_toy_human_gut_gtdb_source_readmap", 0),
        ("cami3_toy_human_gut_gtdb_source_readmap", 1),
        ("cami3_toy_human_gut_gtdb_source_readmap", 2),
    ]

    score_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    for panel, sample in sample_records:
        if panel == "cami2_toy_mouse_gut":
            truth = toy_truth(sample)
            calls, collapse, path = load_toy_calls(sample)
        elif panel == "hmp_airskin_gtdb_source_abundance":
            truth = hmp_truth(sample)
            calls, collapse, path = load_hmp_calls(sample, hmp_taxmap, by_accession, by_core)
        elif panel == "hmp_gastrooral_gtdb_source_abundance":
            truth = hmp_gastro_truth(sample)
            calls, collapse, path = load_hmp_gastro_calls(
                sample,
                hmp_gastro_taxmap,
                by_accession,
                by_core,
            )
        elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
            truth = cami3_truth(sample)
            calls, collapse, path = load_cami3_calls(
                sample,
                cami3_taxid_to_species,
                cami3_name_to_species,
                cami3_taxmap,
                by_accession,
                by_core,
            )
        else:
            raise AssertionError(panel)

        map_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "profile": str(path),
                "collapse_rule": collapse,
                "called_rows_mapped": len(calls),
                "called_species_mapped": calls["gtdb_species"].nunique(),
            }
        )
        for method, raw_values in candidate_raws(calls).items():
            pred = collapse_prediction(calls, raw_values, collapse)
            row = taxid_score.score_prediction(sample, method, pred, truth, {})
            row["panel"] = panel
            row["collapse_rule"] = collapse
            score_rows.append(row)

    scores = pd.DataFrame(score_rows)
    panel_summary, overall = summarize(scores)
    validation_df = validation(panel_summary)
    panel_with_sylph = add_external_sylph(panel_summary)

    scores.to_csv(RESULTS / "cross_panel_abundance_variant_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(RESULTS / "cross_panel_abundance_variant_panel_summary.tsv", sep="\t", index=False)
    panel_with_sylph.to_csv(RESULTS / "cross_panel_abundance_variant_panel_summary_with_sylph.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / "cross_panel_abundance_variant_overall.tsv", sep="\t", index=False)
    validation_df.to_csv(RESULTS / "cross_panel_abundance_variant_validation.tsv", sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(RESULTS / "cross_panel_abundance_variant_mapping.tsv", sep="\t", index=False)

    print("VALIDATION")
    print(validation_df.to_string(index=False))
    print("\nTOP OVERALL")
    print(overall.head(25).to_string(index=False))
    print("\nPANEL TOP")
    print(panel_summary.groupby("panel", group_keys=False).head(8).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
