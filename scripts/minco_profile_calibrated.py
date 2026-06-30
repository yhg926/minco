#!/usr/bin/env python3
"""MinCO calibrated species profiler.

This wrapper builds the joined unique/split readwise feature table used by the
current RF/HGB call gate. It can run concurrent unique/split MinCO passes from
FASTQ input, run experimental same-stream sidecar modes, or consume precomputed
unfiltered unique/split readwise tables. The default strategy is the
F1-priority `universal-auto-exact` gate; use `--strategy probability` to
reproduce the legacy RF/HGB threshold-only output.
"""

from __future__ import annotations

import argparse
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import joblib
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
CAL_DIR = REPO_ROOT / "research/experiments/2026-06-21_minco_multisample_call_calibration/scripts"
HELPER_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(CAL_DIR))
sys.path.insert(0, str(HELPER_DIR))

from analyze_readwise_corrections import parse_species_taxmap  # noqa: E402
from calibrate_multisample_calls import (  # noqa: E402
    MODEL_COLS,
    joined_features,
    load_numeric_minco,
    model_suite,
    predict_probability,
    taxid_scope_map,
)


DEFAULT_THRESHOLD = 0.35
ZIP_POWER = 1.0
EXACT_SPLIT_TRIGGER = 0.10
EXACT_SPLIT_MEDIAN_UAF_GUARD = 0.35
EXACT_SPLIT_RAW_UNIQUE_RATIO_GUARD = 0.80
EXACT_SPLIT_LOW_EXTRA_MODE = "skip"
TAIL_RESCUE_MEDIAN_UAF = 0.35
TAIL_RESCUE_PROBABILITY_THRESHOLD = 0.25
LOW_EXTRA_SPLIT_RESCUE_MEDIAN_UAF = 0.45
LOW_EXTRA_SPLIT_RESCUE_PROBABILITY_THRESHOLD = 0.20
LOW_EXTRA_SPLIT_RESCUE_TOPN_PER_GENUS = 1
MODEL_CACHE_VERSION = 1
DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA = 0.0
ABUNDANCE_GENUS_XNY_BLEND_STRATEGIES = {"universal", "universal-auto-exact"}
ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF = "off"
ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002 = "guarded-genus-hit-breadth-a002"
ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230 = "guarded-genus-hit-breadth-a002-xny230"
ABUNDANCE_FEATURE_ALLOCATOR_ALPHA = 0.02
ABUNDANCE_FEATURE_ALLOCATOR_BASE_MULTI_GENUS_FRAC_MIN = 0.37501510201405147
ABUNDANCE_FEATURE_ALLOCATOR_S_XNY_MEDIAN_MIN = 230.3
ADAPTIVE_CALL_FILTER_SWITCH_OFF = "off"
ADAPTIVE_CALL_FILTER_SWITCH_LOPO_MIN_XNY25 = "lopo-min-xny25"
ADAPTIVE_CALL_FILTER_MAX_XNY_MEDIAN_THRESHOLD = 253.0
ADAPTIVE_CALL_FILTER_MIN_XNY_THRESHOLD = 25.0
CANDIDATE_RESCUE_SWITCH_OFF = "off"
CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70 = "emitted-ani90-xny100-br01-af70"
CANDIDATE_RESCUE_ANI_MIN = 0.90
CANDIDATE_RESCUE_XNY_MIN = 100.0
CANDIDATE_RESCUE_BREADTH_MIN = 0.01
CANDIDATE_RESCUE_REAL_AF_MIN = 0.70
CANDIDATE_SURFACE_SWITCH_OFF = "off"
CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70 = "accession-ani90-xny100-br01-af70"
CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI93_XNY650_BR20 = "accession-ani93-xny650-br20"
CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20 = (
    "accession-current-or-ani93-xny650-br20"
)
CANDIDATE_SURFACE_ANI_MIN = 0.90
CANDIDATE_SURFACE_XNY_MIN = 100.0
CANDIDATE_SURFACE_BREADTH_MIN = 0.01
CANDIDATE_SURFACE_REAL_AF_MIN = 0.70
CANDIDATE_SURFACE_CROSS_PANEL_ANI_MIN = 0.93
CANDIDATE_SURFACE_CROSS_PANEL_XNY_MIN = 650.0
CANDIDATE_SURFACE_CROSS_PANEL_BREADTH_MIN = 0.20
CANDIDATE_SURFACE_CROSS_PANEL_REAL_AF_MIN = 0.0
CANDIDATE_ABUNDANCE_POLICY_ZERO = "zero"
CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2 = "normalized-depth-alpha2"
CANDIDATE_ABUNDANCE_NORMALIZED_DEPTH_ALPHA2 = 2.0


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def validate_candidate_abundance_policy(policy: str) -> None:
    allowed = {
        CANDIDATE_ABUNDANCE_POLICY_ZERO,
        CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    }
    if policy not in allowed:
        raise ValueError(f"unsupported candidate abundance policy: {policy}")


def candidate_policy_alpha(policy: str) -> float:
    validate_candidate_abundance_policy(policy)
    if policy == CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2:
        return CANDIDATE_ABUNDANCE_NORMALIZED_DEPTH_ALPHA2
    return 0.0


def candidate_normalized_depth_mass(features: pd.DataFrame, policy: str) -> np.ndarray:
    alpha = candidate_policy_alpha(policy)
    if alpha <= 0.0:
        return np.zeros(len(features), dtype=float)
    split = numeric(features, "s_Normalized_abundance_depth_max").to_numpy(dtype=float)
    unique = numeric(features, "u_Normalized_abundance_depth_max").to_numpy(dtype=float)
    mass = alpha * np.maximum(split, unique)
    return np.nan_to_num(np.where(mass > 0.0, mass, 0.0), nan=0.0, posinf=0.0, neginf=0.0)


def candidate_abundance_policy_rule(policy: str) -> str:
    validate_candidate_abundance_policy(policy)
    if policy == CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2:
        return (
            "candidate normalized mass = 2*max(s,u Normalized_abundance_depth), "
            "converted to raw scale at output normalization"
        )
    return "candidate rows receive zero abundance mass"


def candidate_surface_thresholds(mode: str) -> tuple[float, float, float, float]:
    if mode in {
        CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20,
    }:
        return (
            CANDIDATE_SURFACE_ANI_MIN,
            CANDIDATE_SURFACE_XNY_MIN,
            CANDIDATE_SURFACE_BREADTH_MIN,
            CANDIDATE_SURFACE_REAL_AF_MIN
            if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
            else CANDIDATE_SURFACE_CROSS_PANEL_REAL_AF_MIN,
        )
    if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI93_XNY650_BR20:
        return (
            CANDIDATE_SURFACE_CROSS_PANEL_ANI_MIN,
            CANDIDATE_SURFACE_CROSS_PANEL_XNY_MIN,
            CANDIDATE_SURFACE_CROSS_PANEL_BREADTH_MIN,
            CANDIDATE_SURFACE_CROSS_PANEL_REAL_AF_MIN,
        )
    raise ValueError(f"unsupported candidate surface switch: {mode}")


def candidate_surface_mask(raw: pd.DataFrame, mode: str) -> pd.Series:
    current_mask = (
        (raw["ANI"] >= CANDIDATE_SURFACE_ANI_MIN)
        & (raw["XnY_ctx"] >= CANDIDATE_SURFACE_XNY_MIN)
        & (raw["Ref_breadth"] >= CANDIDATE_SURFACE_BREADTH_MIN)
        & (raw["Real_min_align_fraction"] >= CANDIDATE_SURFACE_REAL_AF_MIN)
    )
    strict_mask = (
        (raw["ANI"] >= CANDIDATE_SURFACE_CROSS_PANEL_ANI_MIN)
        & (raw["XnY_ctx"] >= CANDIDATE_SURFACE_CROSS_PANEL_XNY_MIN)
        & (raw["Ref_breadth"] >= CANDIDATE_SURFACE_CROSS_PANEL_BREADTH_MIN)
        & (raw["Real_min_align_fraction"] >= CANDIDATE_SURFACE_CROSS_PANEL_REAL_AF_MIN)
    )
    if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70:
        return current_mask
    if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI93_XNY650_BR20:
        return strict_mask
    if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20:
        return current_mask | strict_mask
    raise ValueError(f"unsupported candidate surface switch: {mode}")


def candidate_surface_rule_text(
    ani_min: float,
    xny_min: float,
    breadth_min: float,
    real_af_min: float,
    abundance_policy: str,
) -> str:
    real_af_clause = (
        f", and real AF>={real_af_min:.2f}"
        if real_af_min > 0.0
        else ", with no real-AF gate"
    )
    return (
        "append one accession-level row per uncalled species with "
        f"ANI>={ani_min:.2f}, XnY>={xny_min:g}, breadth>={breadth_min:.2f}"
        f"{real_af_clause}; {candidate_abundance_policy_rule(abundance_policy)}"
    )


def candidate_surface_rule_text_for_mode(mode: str, abundance_policy: str) -> str:
    if mode == CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20:
        return (
            "append one accession-level row per uncalled species passing either "
            "current surface rule ANI>=0.90, XnY>=100, breadth>=0.01, real AF>=0.70 "
            "or raw-retention rule ANI>=0.93, XnY>=650, breadth>=0.20 with no real-AF gate; "
            f"{candidate_abundance_policy_rule(abundance_policy)}"
        )
    ani_min, xny_min, breadth_min, real_af_min = candidate_surface_thresholds(mode)
    return candidate_surface_rule_text(
        ani_min,
        xny_min,
        breadth_min,
        real_af_min,
        abundance_policy,
    )


def eprint(*parts: object) -> None:
    print(*parts, file=sys.stderr)


def scope_keep(taxid: str, scope: str, scope_by_taxid: Mapping[str, str]) -> bool:
    if scope == "all":
        return True
    value = scope_by_taxid.get(str(taxid), "other")
    if scope == "prokaryote":
        return value in {"bacteria", "archaea"}
    return value == scope


def filter_scope(df: pd.DataFrame, scope: str, scope_by_taxid: Mapping[str, str]) -> pd.DataFrame:
    if scope == "all" or "taxid" not in df.columns:
        return df.copy()
    mask = df["taxid"].astype(str).map(lambda t: scope_keep(t, scope, scope_by_taxid))
    return df.loc[mask].copy()


def genus_from_name(name: object) -> str:
    if not isinstance(name, str) or not name.strip():
        return ""
    clean = name.replace("[", "").replace("]", "").strip()
    parts = clean.split()
    if not parts:
        return ""
    if parts[0] == "Candidatus" and len(parts) > 1:
        return parts[1]
    return parts[0]


def topn_by_genus(
    features: pd.DataFrame,
    names: Mapping[str, str],
    mask: np.ndarray,
    score: np.ndarray,
    topn: int,
) -> np.ndarray:
    idx = np.flatnonzero(mask)
    keep = np.zeros(len(features), dtype=bool)
    if len(idx) == 0 or topn <= 0:
        return keep
    taxids = features["taxid"].astype(str).to_numpy()
    work = pd.DataFrame(
        {
            "idx": idx,
            "genus": [genus_from_name(names.get(str(taxids[i]), "")) for i in idx],
            "score": score[idx],
        }
    ).sort_values(["genus", "score"], ascending=[True, False])
    selected = work.groupby("genus", as_index=False).head(topn)
    keep[selected["idx"].to_numpy(dtype=int)] = True
    return keep


def direct_unique_taxids(unique_rows: pd.DataFrame) -> set[str]:
    mask = (
        (numeric(unique_rows, "XnY_ctx") >= 10.0)
        & (numeric(unique_rows, "ANI") >= 0.95)
        & (numeric(unique_rows, "Real_min_align_fraction") >= 0.05)
    )
    return {str(t) for t in unique_rows.loc[mask, "taxid"].astype(str) if str(t)}


def best_ref_metadata(rows: pd.DataFrame, prefix: str) -> pd.DataFrame:
    if rows.empty or "taxid" not in rows.columns:
        return pd.DataFrame(columns=["taxid"])
    work = rows.copy()
    for col in ["accession", "Ref", "Ref_annotation"]:
        if col not in work.columns:
            work[col] = ""
    sort_cols = [
        "taxid",
        "XnY_ctx",
        "Real_min_align_fraction",
        "ANI",
        "Ref_breadth",
        "Ref_mean_depth",
    ]
    for col in sort_cols:
        if col != "taxid" and col not in work.columns:
            work[col] = 0.0
    work = work.sort_values(
        sort_cols,
        ascending=[True, False, False, False, False, False],
        kind="mergesort",
    )
    best = work.drop_duplicates("taxid", keep="first")
    out = best[["taxid", "accession", "Ref", "Ref_annotation"]].copy()
    out = out.rename(
        columns={
            "accession": f"{prefix}_best_accession",
            "Ref": f"{prefix}_best_ref",
            "Ref_annotation": f"{prefix}_best_ref_annotation",
        }
    )
    return out


def panel_abundance_raw(features: pd.DataFrame) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        split_raw = numeric(features, "s_Ref_mean_depth_max").to_numpy(dtype=float) / np.maximum(
            numeric(features, "s_Ref_zip_af_max").to_numpy(dtype=float), 1e-6
        ) ** ZIP_POWER
        unique_raw = numeric(features, "u_Ref_mean_depth_max").to_numpy(dtype=float) / np.maximum(
            numeric(features, "u_Ref_zip_af_max").to_numpy(dtype=float), 1e-6
        ) ** ZIP_POWER
    raw = np.maximum(split_raw, unique_raw)
    if "tail_rescue_added" in features.columns:
        tail_added = features["tail_rescue_added"].astype(bool).to_numpy()
        raw = np.asarray(raw, dtype=float)
        raw[tail_added] = numeric(features, "s_Ref_mean_depth_max").to_numpy(dtype=float)[tail_added]
    return np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)


def min_positive_support(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    out = np.maximum(left, right)
    both = (left > 0.0) & (right > 0.0)
    out[both] = np.minimum(left[both], right[both])
    return out


def apply_adaptive_call_filter_switch(
    features: pd.DataFrame,
    call_mask: np.ndarray,
    mode: str,
) -> tuple[np.ndarray, dict[str, object]]:
    """Apply an experimental output-derived post-call filter.

    The current candidate comes from the cross-panel adaptive switch audit. It
    is intentionally off by default and records enough metadata to make wrapper
    replays auditable.
    """

    call = np.asarray(call_mask, dtype=bool).copy()
    details: dict[str, object] = {
        "adaptive_call_filter_switch": mode,
        "adaptive_call_filter_applied": False,
        "adaptive_call_filter_removed_n": 0,
        "adaptive_call_filter_input_call_n": int(call.sum()),
        "adaptive_call_filter_output_call_n": int(call.sum()),
        "adaptive_call_filter_rule": "",
        "adaptive_call_filter_max_xny_median": 0.0,
        "adaptive_call_filter_max_xny_median_threshold": 0.0,
        "adaptive_call_filter_min_xny_threshold": 0.0,
    }
    if mode == ADAPTIVE_CALL_FILTER_SWITCH_OFF:
        return call, details
    if mode != ADAPTIVE_CALL_FILTER_SWITCH_LOPO_MIN_XNY25:
        raise ValueError(f"unsupported adaptive call-filter switch: {mode}")
    if not np.any(call):
        details.update(
            {
                "adaptive_call_filter_rule": (
                    "if median(max(s_XnY_ctx,u_XnY_ctx)) among current calls >= 253, "
                    "keep only called rows with min_positive(s_XnY_ctx,u_XnY_ctx) >= 25"
                ),
                "adaptive_call_filter_max_xny_median_threshold": ADAPTIVE_CALL_FILTER_MAX_XNY_MEDIAN_THRESHOLD,
                "adaptive_call_filter_min_xny_threshold": ADAPTIVE_CALL_FILTER_MIN_XNY_THRESHOLD,
            }
        )
        return call, details

    s_xny = numeric(features, "s_XnY_ctx_max").to_numpy(dtype=float)
    u_xny = numeric(features, "u_XnY_ctx_max").to_numpy(dtype=float)
    max_xny = np.maximum(s_xny, u_xny)
    min_xny = min_positive_support(s_xny, u_xny)
    median_max_xny = float(np.median(max_xny[call]))
    switch_applies = median_max_xny >= ADAPTIVE_CALL_FILTER_MAX_XNY_MEDIAN_THRESHOLD
    if switch_applies:
        filtered = call & (min_xny >= ADAPTIVE_CALL_FILTER_MIN_XNY_THRESHOLD)
    else:
        filtered = call
    details.update(
        {
            "adaptive_call_filter_applied": switch_applies,
            "adaptive_call_filter_removed_n": int(call.sum() - filtered.sum()),
            "adaptive_call_filter_output_call_n": int(filtered.sum()),
            "adaptive_call_filter_rule": (
                "if median(max(s_XnY_ctx,u_XnY_ctx)) among current calls >= 253, "
                "keep only called rows with min_positive(s_XnY_ctx,u_XnY_ctx) >= 25"
            ),
            "adaptive_call_filter_max_xny_median": median_max_xny,
            "adaptive_call_filter_max_xny_median_threshold": ADAPTIVE_CALL_FILTER_MAX_XNY_MEDIAN_THRESHOLD,
            "adaptive_call_filter_min_xny_threshold": ADAPTIVE_CALL_FILTER_MIN_XNY_THRESHOLD,
        }
    )
    return filtered, details


def apply_candidate_rescue_switch(
    features: pd.DataFrame,
    call_mask: np.ndarray,
    abundance_raw: np.ndarray,
    mode: str,
    abundance_policy: str = CANDIDATE_ABUNDANCE_POLICY_ZERO,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Apply an experimental rescue for strong uncalled emitted candidates.

    The support columns intentionally mirror the emitted-profile candidate
    rescue audit: max split/unique ZIP-AAF ANI, XnY, breadth, and real AF. Raw
    MinCO ANI is not used here because it can saturate in readwise profiles.
    """

    validate_candidate_abundance_policy(abundance_policy)
    call = np.asarray(call_mask, dtype=bool).copy()
    abundance = np.nan_to_num(
        np.asarray(abundance_raw, dtype=float).copy(),
        nan=0.0,
        posinf=0.0,
        neginf=0.0,
    )
    added = np.zeros(len(features), dtype=bool)
    details: dict[str, object] = {
        "candidate_rescue_switch": mode,
        "candidate_rescue_applied": False,
        "candidate_rescue_added_n": 0,
        "candidate_rescue_input_call_n": int(call.sum()),
        "candidate_rescue_output_call_n": int(call.sum()),
        "candidate_rescue_zero_mass": False,
        "candidate_rescue_rule": "",
        "candidate_rescue_ani_min": 0.0,
        "candidate_rescue_xny_min": 0.0,
        "candidate_rescue_breadth_min": 0.0,
        "candidate_rescue_real_af_min": 0.0,
        "candidate_rescue_ani_source_rule": "",
        "candidate_abundance_policy": abundance_policy,
        "candidate_abundance_policy_alpha": candidate_policy_alpha(abundance_policy),
        "candidate_abundance_policy_rule": candidate_abundance_policy_rule(abundance_policy),
    }
    if mode == CANDIDATE_RESCUE_SWITCH_OFF:
        features["candidate_rescue_added"] = added
        features["candidate_abundance_norm_mass"] = np.zeros(len(features), dtype=float)
        return call, abundance, details
    if mode != CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70:
        raise ValueError(f"unsupported candidate rescue switch: {mode}")

    support_ani = np.maximum(
        numeric(features, "s_Ref_zip_aaf_ani_max").to_numpy(dtype=float),
        numeric(features, "u_Ref_zip_aaf_ani_max").to_numpy(dtype=float),
    )
    support_xny = np.maximum(
        numeric(features, "s_XnY_ctx_max").to_numpy(dtype=float),
        numeric(features, "u_XnY_ctx_max").to_numpy(dtype=float),
    )
    support_breadth = np.maximum(
        numeric(features, "s_Ref_breadth_max").to_numpy(dtype=float),
        numeric(features, "u_Ref_breadth_max").to_numpy(dtype=float),
    )
    support_real_af = np.maximum(
        numeric(features, "s_Real_min_align_fraction_max").to_numpy(dtype=float),
        numeric(features, "u_Real_min_align_fraction_max").to_numpy(dtype=float),
    )
    added = (
        (~call)
        & (support_ani >= CANDIDATE_RESCUE_ANI_MIN)
        & (support_xny >= CANDIDATE_RESCUE_XNY_MIN)
        & (support_breadth >= CANDIDATE_RESCUE_BREADTH_MIN)
        & (support_real_af >= CANDIDATE_RESCUE_REAL_AF_MIN)
    )
    call = call | added
    candidate_mass = candidate_normalized_depth_mass(features, abundance_policy)
    abundance[added] = 0.0
    features["candidate_rescue_added"] = added
    features["candidate_abundance_norm_mass"] = np.where(added, candidate_mass, 0.0)
    added_norm_mass = candidate_mass[added]
    zero_mass = not bool(np.any(added_norm_mass > 0.0))
    details.update(
        {
            "candidate_rescue_applied": bool(np.any(added)),
            "candidate_rescue_added_n": int(np.sum(added)),
            "candidate_rescue_output_call_n": int(call.sum()),
            "candidate_rescue_zero_mass": zero_mass,
            "candidate_rescue_rule": (
                "add uncalled emitted-profile candidates with "
                "max(split,unique) ZIP-AAF ANI>=0.90, XnY>=100, "
                f"breadth>=0.01, and real AF>=0.70; {candidate_abundance_policy_rule(abundance_policy)}"
            ),
            "candidate_rescue_ani_min": CANDIDATE_RESCUE_ANI_MIN,
            "candidate_rescue_xny_min": CANDIDATE_RESCUE_XNY_MIN,
            "candidate_rescue_breadth_min": CANDIDATE_RESCUE_BREADTH_MIN,
            "candidate_rescue_real_af_min": CANDIDATE_RESCUE_REAL_AF_MIN,
            "candidate_rescue_ani_source_rule": "max(s_Ref_zip_aaf_ani_max,u_Ref_zip_aaf_ani_max)",
        }
    )
    return call, abundance, details


def tax_record_for_accession(
    accession: object,
    taxmap: Mapping[str, Mapping[str, str]],
) -> Mapping[str, str]:
    acc = str(accession or "")
    if acc in taxmap:
        return taxmap[acc]
    if "." in acc:
        short = acc.rsplit(".", 1)[0]
        if short in taxmap:
            return taxmap[short]
    return {}


def species_for_accession(
    accession: object,
    taxmap: Mapping[str, Mapping[str, str]],
    fallback: str = "",
) -> str:
    rec = tax_record_for_accession(accession, taxmap)
    return str(rec.get("species_name", "") or fallback or "")


def called_accession_species(
    features: pd.DataFrame,
    call_mask: np.ndarray,
    taxmap: Mapping[str, Mapping[str, str]],
    names: Mapping[str, str],
) -> set[str]:
    called = np.asarray(call_mask, dtype=bool)
    out: set[str] = set()
    if not len(features):
        return out
    for row in features.loc[called].itertuples(index=False):
        species = ""
        for col in ["s_best_accession", "u_best_accession"]:
            accession = str(getattr(row, col, "") or "")
            if accession and accession.lower() != "nan":
                species = species_for_accession(accession, taxmap)
                if species:
                    break
        if not species:
            taxid = str(getattr(row, "taxid", "") or "")
            species = names.get(taxid, "")
        if species:
            out.add(species)
    return out


def called_taxid_species(
    features: pd.DataFrame,
    call_mask: np.ndarray,
    names: Mapping[str, str],
) -> set[str]:
    called = np.asarray(call_mask, dtype=bool)
    out: set[str] = set()
    if not len(features) or "taxid" not in features.columns:
        return out
    for taxid in features.loc[called, "taxid"].astype(str):
        species = names.get(str(taxid), "")
        if species:
            out.add(species)
    return out


def candidate_surface_has_numeric_candidates(
    unique_rows: pd.DataFrame,
    split_rows: pd.DataFrame,
    mode: str,
) -> bool:
    if mode == CANDIDATE_SURFACE_SWITCH_OFF:
        return False
    for rows in [unique_rows, split_rows]:
        if rows.empty:
            continue
        work = pd.DataFrame(index=rows.index)
        for col in ["ANI", "XnY_ctx", "Ref_breadth", "Real_min_align_fraction"]:
            work[col] = numeric(rows, col)
        if bool(candidate_surface_mask(work, mode).any()):
            return True
    return False


def build_accession_candidate_surface(
    unique_rows: pd.DataFrame,
    split_rows: pd.DataFrame,
    features: pd.DataFrame,
    call_mask: np.ndarray,
    mode: str,
    taxmap: Mapping[str, Mapping[str, str]],
    names: Mapping[str, str],
    abundance_policy: str = CANDIDATE_ABUNDANCE_POLICY_ZERO,
    max_called_species: int = 0,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Build off-by-default accession-level candidate rows.

    The normal calibrated model is taxid-collapsed. This experimental surface
    keeps a small number of strong accession-level candidates that would
    otherwise be hidden behind an already-emitted taxid row. It is output-only
    and does not change the fitted model or base calls.
    """

    validate_candidate_abundance_policy(abundance_policy)
    details: dict[str, object] = {
        "candidate_surface_switch": mode,
        "candidate_surface_added": False,
        "candidate_surface_applied": False,
        "candidate_surface_added_n": 0,
        "candidate_surface_zero_mass": False,
        "candidate_surface_rule": "",
        "candidate_surface_ani_min": 0.0,
        "candidate_surface_xny_min": 0.0,
        "candidate_surface_breadth_min": 0.0,
        "candidate_surface_real_af_min": 0.0,
        "candidate_surface_source_rule": "",
        "candidate_surface_called_species_n": 0,
        "candidate_surface_max_called_species": int(max_called_species),
        "candidate_surface_guard_blocked": False,
    }
    if mode == CANDIDATE_SURFACE_SWITCH_OFF:
        return pd.DataFrame(), details
    ani_min, xny_min, breadth_min, real_af_min = candidate_surface_thresholds(mode)
    rule_text = candidate_surface_rule_text_for_mode(mode, abundance_policy)

    frames: list[pd.DataFrame] = []
    for raw_mode, rows in [("unique", unique_rows), ("split", split_rows)]:
        if rows.empty:
            continue
        work = rows.copy()
        work["candidate_surface_raw_mode"] = raw_mode
        frames.append(work)
    if not frames:
        details.update(
            {
                "candidate_surface_rule": rule_text,
                "candidate_surface_ani_min": ani_min,
                "candidate_surface_xny_min": xny_min,
                "candidate_surface_breadth_min": breadth_min,
                "candidate_surface_real_af_min": real_af_min,
                "candidate_surface_source_rule": "raw unique/split accession rows",
            }
        )
        return pd.DataFrame(), details

    raw = pd.concat(frames, ignore_index=True, sort=False)
    for col in [
        "accession",
        "taxid",
        "species_name",
        "Ref",
        "Ref_annotation",
        "candidate_surface_raw_mode",
    ]:
        if col not in raw.columns:
            raw[col] = ""
    for col in ["ANI", "XnY_ctx", "Ref_breadth", "Real_min_align_fraction"]:
        raw[col] = numeric(raw, col)
    called_species = called_taxid_species(features, call_mask, names)
    called_species_n = len(called_species)
    if max_called_species > 0 and called_species_n > max_called_species:
        details.update(
            {
                "candidate_surface_rule": (
                    f"{rule_text}; skipped when called species count > {int(max_called_species)}"
                ),
                "candidate_surface_ani_min": ani_min,
                "candidate_surface_xny_min": xny_min,
                "candidate_surface_breadth_min": breadth_min,
                "candidate_surface_real_af_min": real_af_min,
                "candidate_surface_source_rule": "raw unique/split accession rows",
                "candidate_surface_called_species_n": called_species_n,
                "candidate_surface_max_called_species": int(max_called_species),
                "candidate_surface_guard_blocked": True,
            }
        )
        return pd.DataFrame(), details

    numeric_mask = candidate_surface_mask(raw, mode)
    if not bool(numeric_mask.any()):
        details.update(
            {
                "candidate_surface_rule": rule_text,
                "candidate_surface_ani_min": ani_min,
                "candidate_surface_xny_min": xny_min,
                "candidate_surface_breadth_min": breadth_min,
                "candidate_surface_real_af_min": real_af_min,
                "candidate_surface_source_rule": "raw unique/split accession rows",
                "candidate_surface_called_species_n": called_species_n,
                "candidate_surface_max_called_species": int(max_called_species),
                "candidate_surface_guard_blocked": False,
            }
        )
        return pd.DataFrame(), details

    raw = raw.loc[numeric_mask].copy()
    raw["candidate_surface_species_name"] = raw.apply(
        lambda row: species_for_accession(
            row.get("accession", ""),
            taxmap,
            str(row.get("species_name", "") or ""),
        ),
        axis=1,
    )

    called_species = called_accession_species(features, call_mask, taxmap, names)
    called_species_n = len(called_species)
    if max_called_species > 0 and called_species_n > max_called_species:
        details.update(
            {
                "candidate_surface_rule": (
                    f"{rule_text}; skipped when called species count > {int(max_called_species)}"
                ),
                "candidate_surface_ani_min": ani_min,
                "candidate_surface_xny_min": xny_min,
                "candidate_surface_breadth_min": breadth_min,
                "candidate_surface_real_af_min": real_af_min,
                "candidate_surface_source_rule": "raw unique/split accession rows",
                "candidate_surface_called_species_n": called_species_n,
                "candidate_surface_max_called_species": int(max_called_species),
                "candidate_surface_guard_blocked": True,
            }
        )
        return pd.DataFrame(), details

    mask = (
        raw["candidate_surface_species_name"].astype(str).astype(bool)
        & ~raw["candidate_surface_species_name"].astype(str).isin(called_species)
    )
    candidates = raw.loc[mask].copy()
    if not candidates.empty:
        candidates = candidates.sort_values(
            [
                "candidate_surface_species_name",
                "XnY_ctx",
                "Real_min_align_fraction",
                "ANI",
                "Ref_breadth",
            ],
            ascending=[True, False, False, False, False],
            kind="mergesort",
        ).drop_duplicates("candidate_surface_species_name", keep="first")
    alpha = candidate_policy_alpha(abundance_policy)
    if alpha > 0.0 and not candidates.empty:
        candidate_mass = alpha * numeric(candidates, "Normalized_abundance_depth")
        surface_zero_mass = not bool(np.any(candidate_mass.to_numpy(dtype=float) > 0.0))
    else:
        surface_zero_mass = True if not candidates.empty else False
    details.update(
        {
            "candidate_surface_applied": bool(len(candidates)),
            "candidate_surface_added_n": int(len(candidates)),
            "candidate_surface_zero_mass": surface_zero_mass,
            "candidate_surface_rule": rule_text,
            "candidate_surface_ani_min": ani_min,
            "candidate_surface_xny_min": xny_min,
            "candidate_surface_breadth_min": breadth_min,
            "candidate_surface_real_af_min": real_af_min,
            "candidate_surface_source_rule": "raw unique/split accession rows",
            "candidate_surface_called_species_n": called_species_n,
            "candidate_surface_max_called_species": int(max_called_species),
            "candidate_surface_guard_blocked": False,
        }
    )
    return candidates, details


def ensure_candidate_surface_output_columns(work: pd.DataFrame) -> pd.DataFrame:
    needed = [
        "candidate_surface_added",
        "candidate_surface_accession",
        "candidate_surface_species_name",
        "candidate_surface_taxid",
        "candidate_surface_raw_mode",
        "candidate_surface_zero_mass",
        "candidate_abundance_policy",
        "candidate_abundance_policy_alpha",
        "candidate_abundance_policy_rule",
        "candidate_abundance_norm_mass",
        "u_best_accession",
        "u_best_ref",
        "u_best_ref_annotation",
        "s_best_accession",
        "s_best_ref",
        "s_best_ref_annotation",
    ]
    metric_suffixes = [
        "rows",
        "direct_call",
        "relaxed_call",
        "ANI_max",
        "XnY_ctx_max",
        "Real_min_align_fraction_max",
        "Ref_breadth_max",
        "Ref_mean_depth_max",
        "Ref_hit_mean_depth_max",
        "Ref_depth_cv_min",
        "Normalized_abundance_depth_max",
        "Ref_zip_af_max",
        "Ref_zip_aaf_ani_max",
    ]
    needed.extend(f"{prefix}_{suffix}" for prefix in ["u", "s"] for suffix in metric_suffixes)
    for col in needed:
        if col not in work.columns:
            work[col] = False if col in {"candidate_surface_added", "candidate_surface_zero_mass"} else ""
    return work


def append_accession_candidate_surface_rows(
    work: pd.DataFrame,
    candidates: pd.DataFrame,
    abundance_policy: str = CANDIDATE_ABUNDANCE_POLICY_ZERO,
    candidate_mass_scale: float = 1.0,
) -> pd.DataFrame:
    validate_candidate_abundance_policy(abundance_policy)
    mass_scale = finite(candidate_mass_scale, 1.0)
    if mass_scale <= 0.0:
        mass_scale = 1.0
    work = ensure_candidate_surface_output_columns(work)
    if candidates.empty:
        return work
    columns = list(work.columns)
    added_rows: list[dict[str, object]] = []
    for raw in candidates.itertuples(index=False):
        raw_mode = str(getattr(raw, "candidate_surface_raw_mode", "") or "")
        prefix = "u" if raw_mode == "unique" else "s"
        accession = str(getattr(raw, "accession", "") or "")
        species_name = str(getattr(raw, "candidate_surface_species_name", "") or "")
        norm_mass = 0.0
        raw_mass = 0.0
        if abundance_policy == CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2:
            norm_mass = CANDIDATE_ABUNDANCE_NORMALIZED_DEPTH_ALPHA2 * finite(
                getattr(raw, "Normalized_abundance_depth", 0.0)
            )
            norm_mass = norm_mass if norm_mass > 0.0 else 0.0
            raw_mass = norm_mass * mass_scale
        row = {col: "" for col in columns}
        row.update(
            {
                "taxid": str(getattr(raw, "taxid", "") or species_name),
                "species_name": species_name,
                "calibrated_call": True,
                "calibrated_probability": 0.0,
                "calibrated_abundance": 0.0,
                "calibrated_abundance_raw": raw_mass,
                "reported_ani": finite(getattr(raw, "ANI", 0.0)),
                "reported_ani_source": "candidate_surface_ANI",
                "raw_unique_fallback": False,
                "candidate_surface_added": True,
                "candidate_surface_accession": accession,
                "candidate_surface_species_name": species_name,
                "candidate_surface_taxid": str(getattr(raw, "taxid", "") or ""),
                "candidate_surface_raw_mode": raw_mode,
                "candidate_surface_zero_mass": raw_mass <= 0.0,
                "candidate_abundance_policy": abundance_policy,
                "candidate_abundance_policy_alpha": candidate_policy_alpha(abundance_policy),
                "candidate_abundance_policy_rule": candidate_abundance_policy_rule(abundance_policy),
                "candidate_abundance_norm_mass": norm_mass,
            }
        )
        for col in [
            "profile_strategy",
            "calibrated_threshold",
            "train_pool",
            "scope",
            "abundance_rule",
            "abundance_genus_xny_blend_alpha",
            "abundance_genus_xny_quality",
            "candidate_surface_switch",
            "candidate_surface_applied",
            "candidate_surface_added_n",
            "candidate_surface_rule",
            "candidate_surface_ani_min",
            "candidate_surface_xny_min",
            "candidate_surface_breadth_min",
            "candidate_surface_real_af_min",
            "candidate_surface_source_rule",
            "candidate_surface_called_species_n",
            "candidate_surface_max_called_species",
            "candidate_surface_guard_blocked",
            "candidate_abundance_policy",
            "candidate_abundance_policy_alpha",
            "candidate_abundance_policy_rule",
        ]:
            if col in work.columns:
                values = work[col].dropna()
                if len(values):
                    value = values.iloc[0]
                    if not (isinstance(value, str) and value == ""):
                        row[col] = value
        row[f"{prefix}_rows"] = 1
        row[f"{prefix}_direct_call"] = 1
        row[f"{prefix}_relaxed_call"] = 1
        row[f"{prefix}_ANI_max"] = finite(getattr(raw, "ANI", 0.0))
        row[f"{prefix}_XnY_ctx_max"] = finite(getattr(raw, "XnY_ctx", 0.0))
        row[f"{prefix}_Real_min_align_fraction_max"] = finite(
            getattr(raw, "Real_min_align_fraction", 0.0)
        )
        row[f"{prefix}_Ref_breadth_max"] = finite(getattr(raw, "Ref_breadth", 0.0))
        row[f"{prefix}_Ref_mean_depth_max"] = finite(getattr(raw, "Ref_mean_depth", 0.0))
        row[f"{prefix}_Ref_hit_mean_depth_max"] = finite(
            getattr(raw, "Ref_hit_mean_depth", 0.0)
        )
        row[f"{prefix}_Ref_depth_cv_min"] = finite(getattr(raw, "Ref_depth_cv", 0.0))
        row[f"{prefix}_Normalized_abundance_depth_max"] = finite(
            getattr(raw, "Normalized_abundance_depth", 0.0)
        )
        row[f"{prefix}_Ref_zip_af_max"] = finite(getattr(raw, "Ref_zip_af", 0.0))
        row[f"{prefix}_Ref_zip_aaf_ani_max"] = finite(getattr(raw, "Ref_zip_aaf_ani", 0.0))
        row[f"{prefix}_best_accession"] = accession
        row[f"{prefix}_best_ref"] = getattr(raw, "Ref", "")
        row[f"{prefix}_best_ref_annotation"] = getattr(raw, "Ref_annotation", "")
        added_rows.append(row)

    return pd.concat(
        [work, pd.DataFrame(added_rows, columns=columns)],
        ignore_index=True,
        sort=False,
    )


def genus_xny_blended_abundance_raw(
    features: pd.DataFrame,
    call_mask: np.ndarray,
    abundance_raw: np.ndarray,
    alpha: float = DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA,
) -> np.ndarray:
    """Blend current raw abundance with within-genus XnY-weighted mass.

    The adjustment is intentionally conservative: only called species
    participate, each genus keeps the same total raw abundance, and `alpha`
    controls how much of the mass is reallocated by split XnY support.
    """

    base = np.nan_to_num(np.asarray(abundance_raw, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    base = np.where(base > 0.0, base, 0.0)
    adjusted = base.copy()
    called = np.asarray(call_mask, dtype=bool)
    idx = np.flatnonzero(called & (base > 0.0))
    if len(idx) <= 1 or alpha <= 0.0:
        return adjusted

    species_names = features.get("species_name")
    if species_names is None:
        return adjusted
    quality = np.clip(numeric(features, "s_XnY_ctx_max").to_numpy(dtype=float) / 1000.0, 0.0, 1.0)
    work = pd.DataFrame(
        {
            "idx": idx,
            "genus": [genus_from_name(species_names.iloc[i]) for i in idx],
            "base": base[idx],
            "quality": quality[idx],
        }
    )
    for _genus, sub in work.groupby("genus", sort=False):
        group_idx = sub["idx"].to_numpy(dtype=int)
        if len(group_idx) <= 1:
            continue
        total = float(base[group_idx].sum())
        weighted = base[group_idx] * quality[group_idx]
        weighted_sum = float(weighted.sum())
        if total <= 0.0 or weighted_sum <= 0.0:
            continue
        reallocated = weighted / weighted_sum * total
        adjusted[group_idx] = (1.0 - alpha) * base[group_idx] + alpha * reallocated
    return np.nan_to_num(adjusted, nan=0.0, posinf=0.0, neginf=0.0)


def boolean_column(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    if pd.api.types.is_bool_dtype(df[col]):
        return df[col].fillna(False).astype(bool)
    return df[col].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def feature_allocator_species_labels(
    work: pd.DataFrame,
    taxmap: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> pd.Series:
    fallback = (
        work["species_name"].fillna("").astype(str)
        if "species_name" in work.columns
        else pd.Series([""] * len(work), index=work.index, dtype=str)
    )
    if not taxmap:
        return fallback
    labels: list[str] = []
    for raw in work.itertuples(index=False):
        label = ""
        for col in ["s_best_accession", "u_best_accession", "candidate_surface_accession", "taxid"]:
            if not hasattr(raw, col):
                continue
            key = str(getattr(raw, col, "") or "")
            if not key:
                continue
            rec = taxmap.get(key)
            if rec:
                label = str(rec.get("species_name", "") or "")
                if label:
                    break
        labels.append(label)
    mapped = pd.Series(labels, index=work.index, dtype=str)
    return mapped.where(mapped.astype(bool), fallback)


def base_multi_genus_mass_fraction(
    work: pd.DataFrame,
    call_mask: np.ndarray,
    abundance_raw: np.ndarray,
    taxmap: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> float:
    base = np.nan_to_num(np.asarray(abundance_raw, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    base = np.where(base > 0.0, base, 0.0)
    called = np.asarray(call_mask, dtype=bool)
    candidate_added = boolean_column(work, "candidate_rescue_added") | boolean_column(
        work, "candidate_surface_added"
    )
    eligible = called & (~candidate_added.to_numpy(dtype=bool)) & (base > 0.0)
    if not bool(eligible.any()):
        return 0.0
    species_names = feature_allocator_species_labels(work, taxmap).to_numpy()
    genus = np.asarray([genus_from_name(name) for name in species_names], dtype=object)
    idx = np.flatnonzero(eligible)
    total = float(base[idx].sum())
    if total <= 0.0:
        return 0.0
    counts = pd.Series(genus[idx]).value_counts()
    multi = {str(key) for key, value in counts.items() if int(value) > 1 and str(key)}
    if not multi:
        return 0.0
    multi_mask = np.asarray([str(value) in multi for value in genus], dtype=bool)
    return float(base[eligible & multi_mask].sum() / total)


def guarded_feature_allocator_abundance_raw(
    work: pd.DataFrame,
    call_mask: np.ndarray,
    abundance_raw: np.ndarray,
    mode: str,
    taxmap: Optional[Mapping[str, Mapping[str, str]]] = None,
) -> tuple[np.ndarray, dict[str, object]]:
    base = np.nan_to_num(np.asarray(abundance_raw, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    base = np.where(base > 0.0, base, 0.0)
    details: dict[str, object] = {
        "abundance_feature_allocator_switch": mode,
        "abundance_feature_allocator_applied": False,
        "abundance_feature_allocator_guard_passed": False,
        "abundance_feature_allocator_method": "",
        "abundance_feature_allocator_alpha": 0.0,
        "abundance_feature_allocator_target": "",
        "abundance_feature_allocator_guard_feature": "base_mass_multi_genus_frac",
        "abundance_feature_allocator_guard_threshold": ABUNDANCE_FEATURE_ALLOCATOR_BASE_MULTI_GENUS_FRAC_MIN,
        "abundance_feature_allocator_base_mass_multi_genus_frac": 0.0,
        "abundance_feature_allocator_s_xny_median": 0.0,
        "abundance_feature_allocator_s_xny_median_threshold": 0.0,
        "abundance_feature_allocator_adjusted_rows_n": 0,
        "abundance_feature_allocator_rule": "off",
    }
    if mode == ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF:
        return base, details
    if mode not in {
        ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002,
        ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230,
    }:
        raise ValueError(f"unsupported abundance feature allocator switch: {mode}")

    called = np.asarray(call_mask, dtype=bool)
    candidate_added = boolean_column(work, "candidate_rescue_added") | boolean_column(
        work, "candidate_surface_added"
    )
    eligible = called & (~candidate_added.to_numpy(dtype=bool)) & (base > 0.0)
    base_multi_frac = base_multi_genus_mass_fraction(work, called, base, taxmap)
    s_xny_median = 0.0
    s_xny_guard_passed = True
    if mode == ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230:
        called_s_xny = numeric(work.loc[called], "s_XnY_ctx_max")
        s_xny_median = float(called_s_xny.median()) if len(called_s_xny) else 0.0
        s_xny_guard_passed = s_xny_median >= ABUNDANCE_FEATURE_ALLOCATOR_S_XNY_MEDIAN_MIN
    guard_passed = (
        base_multi_frac >= ABUNDANCE_FEATURE_ALLOCATOR_BASE_MULTI_GENUS_FRAC_MIN
        and s_xny_guard_passed
    )
    guard_feature = "base_mass_multi_genus_frac"
    guard_threshold: object = ABUNDANCE_FEATURE_ALLOCATOR_BASE_MULTI_GENUS_FRAC_MIN
    method = "genus_raw_hit_breadth_a0.02"
    rule = (
        "if base_mass_multi_genus_frac>=0.375015 then blend 0.02 of "
        "each multi-row called genus toward raw*hit-depth*breadth; "
        "candidate-added row mass is preserved"
    )
    if mode == ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230:
        guard_feature = "base_mass_multi_genus_frac;s_xny_median"
        guard_threshold = (
            f"{ABUNDANCE_FEATURE_ALLOCATOR_BASE_MULTI_GENUS_FRAC_MIN};"
            f"{ABUNDANCE_FEATURE_ALLOCATOR_S_XNY_MEDIAN_MIN}"
        )
        method = "genus_raw_hit_breadth_a0.02_xny230_guard"
        rule = (
            "if base_mass_multi_genus_frac>=0.375015 and s_xny_median>=230.3 "
            "then blend 0.02 of each multi-row called genus toward "
            "raw*hit-depth*breadth; candidate-added row mass is preserved"
        )
    details.update(
        {
            "abundance_feature_allocator_guard_passed": guard_passed,
            "abundance_feature_allocator_method": method,
            "abundance_feature_allocator_alpha": ABUNDANCE_FEATURE_ALLOCATOR_ALPHA,
            "abundance_feature_allocator_target": "raw*max(s,u Ref_hit_mean_depth*Ref_breadth)",
            "abundance_feature_allocator_guard_feature": guard_feature,
            "abundance_feature_allocator_guard_threshold": guard_threshold,
            "abundance_feature_allocator_base_mass_multi_genus_frac": base_multi_frac,
            "abundance_feature_allocator_s_xny_median": s_xny_median,
            "abundance_feature_allocator_s_xny_median_threshold": (
                ABUNDANCE_FEATURE_ALLOCATOR_S_XNY_MEDIAN_MIN
                if mode == ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230
                else 0.0
            ),
            "abundance_feature_allocator_rule": rule,
        }
    )
    if not guard_passed or not bool(eligible.any()):
        return base, details

    species_names = feature_allocator_species_labels(work, taxmap).to_numpy()
    genus = np.asarray([genus_from_name(name) for name in species_names], dtype=object)
    split_quality = numeric(work, "s_Ref_hit_mean_depth_max").to_numpy(dtype=float) * numeric(
        work, "s_Ref_breadth_max"
    ).to_numpy(dtype=float)
    unique_quality = numeric(work, "u_Ref_hit_mean_depth_max").to_numpy(dtype=float) * numeric(
        work, "u_Ref_breadth_max"
    ).to_numpy(dtype=float)
    quality = np.nan_to_num(np.maximum(split_quality, unique_quality), nan=0.0, posinf=0.0, neginf=0.0)
    adjusted = base.copy()
    adjusted_rows: set[int] = set()
    for value in sorted(set(genus[eligible])):
        group = eligible & (genus == value)
        group_idx = np.flatnonzero(group)
        if len(group_idx) <= 1:
            continue
        total = float(base[group_idx].sum())
        weighted = base[group_idx] * quality[group_idx]
        weighted_sum = float(weighted.sum())
        if total <= 0.0 or weighted_sum <= 0.0:
            continue
        target = weighted / weighted_sum * total
        adjusted[group_idx] = (
            (1.0 - ABUNDANCE_FEATURE_ALLOCATOR_ALPHA) * base[group_idx]
            + ABUNDANCE_FEATURE_ALLOCATOR_ALPHA * target
        )
        adjusted_rows.update(int(i) for i in group_idx)
    details["abundance_feature_allocator_applied"] = bool(adjusted_rows)
    details["abundance_feature_allocator_adjusted_rows_n"] = len(adjusted_rows)
    return np.nan_to_num(adjusted, nan=0.0, posinf=0.0, neginf=0.0), details


def normalized_split_abundance_raw(features: pd.DataFrame) -> np.ndarray:
    raw = numeric(features, "s_Normalized_abundance_depth_max").to_numpy(dtype=float)
    return np.where(raw > 0.0, raw, 0.0)


def reported_ani_values(features: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return diagnostic ANI values and their source columns.

    The MinCO raw/emitted ANI field can saturate near 1.0 in readwise
    metagenome profiles. The ZIP-AAF ANI column is the continuous diagnostic
    signal supported by the cached source/ref ANI benchmarks.
    """

    split_ani = numeric(features, "s_Ref_zip_aaf_ani_max").to_numpy(dtype=float)
    unique_ani = numeric(features, "u_Ref_zip_aaf_ani_max").to_numpy(dtype=float)
    if "raw_unique_fallback" in features.columns:
        use_unique = features["raw_unique_fallback"].astype(bool).to_numpy() & (unique_ani > 0.0)
    else:
        use_unique = np.zeros(len(features), dtype=bool)
    values = np.where(use_unique, unique_ani, split_ani)
    source = np.where(use_unique, "u_Ref_zip_aaf_ani_max", "s_Ref_zip_aaf_ani_max")
    missing_split = (values <= 0.0) & (unique_ani > 0.0)
    values = np.where(missing_split, unique_ani, values)
    source = np.where(missing_split, "u_Ref_zip_aaf_ani_max", source)
    return np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0), source


def exact_split_guard_passes(details: Mapping[str, object]) -> bool:
    median_uaf = float(details.get("joined_base_median_uaf", 0.0))
    raw_ratio = float(details.get("raw_unique_to_base_ratio", 0.0))
    return bool(
        median_uaf < EXACT_SPLIT_MEDIAN_UAF_GUARD
        or raw_ratio >= EXACT_SPLIT_RAW_UNIQUE_RATIO_GUARD
    )


def adaptive_sub95_mask(
    features: pd.DataFrame,
    names: Mapping[str, str],
    probability_threshold: float,
) -> tuple[np.ndarray, dict[str, object]]:
    prob = numeric(features, "calibrated_probability").to_numpy(dtype=float)
    tail_prob_gate = prob >= TAIL_RESCUE_PROBABILITY_THRESHOLD
    u_x = numeric(features, "u_XnY_ctx_max").to_numpy(dtype=float)
    u_ani = numeric(features, "u_ANI_max").to_numpy(dtype=float)
    u_af = numeric(features, "u_Real_min_align_fraction_max").to_numpy(dtype=float)
    s_x = numeric(features, "s_XnY_ctx_max").to_numpy(dtype=float)
    s_ani = numeric(features, "s_ANI_max").to_numpy(dtype=float)
    s_af = numeric(features, "s_Real_min_align_fraction_max").to_numpy(dtype=float)
    s_b = numeric(features, "s_Ref_breadth_max").to_numpy(dtype=float)
    s_hit = numeric(features, "s_Ref_hit_mean_depth_max").to_numpy(dtype=float)
    s_norm = numeric(features, "s_Normalized_abundance_depth_max").to_numpy(dtype=float)

    base = (u_x >= 10.0) & (u_ani >= 0.95) & (u_af >= 0.05)
    prob_gate = prob >= probability_threshold
    split_strict = (s_x >= 50.0) & (s_ani >= 0.95) & (s_af >= 0.05) & (s_b >= 0.05)
    split_near_unique = split_strict & (u_x >= 50.0) & (u_ani >= 0.93)
    split_score = s_x * np.maximum(s_b, 1e-6) * np.maximum(s_ani, 0.0) * np.maximum(prob, 0.01)

    rescue = base | topn_by_genus(features, names, prob_gate & split_strict & (~base), split_score, 2)
    strict_plain = base | split_near_unique
    extra = prob_gate & (~base)
    base_mass = float(np.sum(s_norm[base]))
    extra_mass = float(np.sum(s_norm[extra]))
    mass_ratio = extra_mass / (base_mass + 1e-12)
    count_ratio = float(np.sum(prob_gate)) / max(float(np.sum(base)), 1.0)
    median_base_uaf = float(np.median(u_af[base])) if np.sum(base) else 0.0

    tail_rescue_added = np.zeros(len(features), dtype=bool)
    low_extra_split_rescue_added = np.zeros(len(features), dtype=bool)
    if mass_ratio < 0.005:
        adaptive = base.copy()
        mode = "base"
        if median_base_uaf < LOW_EXTRA_SPLIT_RESCUE_MEDIAN_UAF:
            low_extra = (
                (prob >= LOW_EXTRA_SPLIT_RESCUE_PROBABILITY_THRESHOLD)
                & split_strict
                & (~adaptive)
            )
            low_extra_split_rescue_added = topn_by_genus(
                features,
                names,
                low_extra,
                split_score,
                LOW_EXTRA_SPLIT_RESCUE_TOPN_PER_GENUS,
            )
            adaptive = adaptive | low_extra_split_rescue_added
            if np.any(low_extra_split_rescue_added):
                mode = "base_low_extra_split_rescue"
    elif mass_ratio >= 0.20 and count_ratio >= 1.20:
        if median_base_uaf < TAIL_RESCUE_MEDIAN_UAF:
            adaptive = tail_prob_gate.copy()
            mode = "probability_high_extra_low_uaf_tail025"
            tail_rescue_added = tail_prob_gate & (~rescue)
        else:
            adaptive = rescue.copy()
            mode = "rescue"
    elif int(np.sum(prob_gate)) <= int(np.sum(base)) and mass_ratio < 0.03:
        adaptive = strict_plain.copy()
        mode = "strict_plain"
    else:
        adaptive = prob_gate.copy()
        mode = "probability"

    near = (
        (s_x >= 300.0)
        & (s_ani >= 0.948)
        & (s_ani < 0.9505)
        & (s_af >= 0.05)
        & (s_b >= 0.20)
        & (s_hit >= 1.4)
        & (prob >= 0.20)
        & (u_ani >= 0.92)
        & (~adaptive)
    )
    near_score = prob * np.maximum(s_x, 0.0) * np.maximum(s_b, 1e-6) * np.maximum(s_hit, 1.0)
    final = adaptive | topn_by_genus(features, names, near, near_score, 1)
    features["tail_rescue_added"] = tail_rescue_added
    features["low_extra_split_rescue_added"] = low_extra_split_rescue_added
    details = {
        "adaptive_mode": mode,
        "joined_base_n": int(np.sum(base)),
        "probability_gate_n": int(np.sum(prob_gate)),
        "probability_extra_n": int(np.sum(extra)),
        "probability_extra_mass_ratio": mass_ratio,
        "probability_count_ratio": count_ratio,
        "joined_base_median_uaf": median_base_uaf,
        "tail_rescue_added_n": int(np.sum(tail_rescue_added)),
        "tail_rescue_median_uaf_threshold": TAIL_RESCUE_MEDIAN_UAF,
        "tail_rescue_probability_threshold": TAIL_RESCUE_PROBABILITY_THRESHOLD,
        "tail_rescue_abundance_rule": "s_Ref_mean_depth_for_tail_rescue_added_else_max_split_unique_mean_over_zip_af",
        "low_extra_split_rescue_added_n": int(np.sum(low_extra_split_rescue_added)),
        "low_extra_split_rescue_median_uaf_threshold": LOW_EXTRA_SPLIT_RESCUE_MEDIAN_UAF,
        "low_extra_split_rescue_probability_threshold": LOW_EXTRA_SPLIT_RESCUE_PROBABILITY_THRESHOLD,
        "low_extra_split_rescue_topn_per_genus": LOW_EXTRA_SPLIT_RESCUE_TOPN_PER_GENUS,
        "low_extra_split_rescue_rule": (
            "if probability_extra_mass_ratio<0.005 and joined_base_median_uaf<0.45, "
            "add top1/genus candidates with P>=0.20 and strict split evidence"
        ),
    }
    return final, details


def strategy_call_and_abundance(
    features: pd.DataFrame,
    unique_rows: pd.DataFrame,
    names: Mapping[str, str],
    strategy: str,
    probability_threshold: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    if strategy == "probability":
        call = numeric(features, "calibrated_probability").to_numpy(dtype=float) >= probability_threshold
        return call, normalized_split_abundance_raw(features), {"profile_strategy": strategy}

    incumbent, details = adaptive_sub95_mask(features, names, probability_threshold)
    if strategy == "adaptive-sub95":
        details["profile_strategy"] = strategy
        details["raw_unique_fallback"] = False
        return incumbent, panel_abundance_raw(features), details

    if strategy != "universal":
        raise ValueError(f"unsupported profile strategy: {strategy}")

    raw_taxids = direct_unique_taxids(unique_rows)
    raw_mask = features["taxid"].astype(str).isin(raw_taxids).to_numpy(dtype=bool)
    base_n = int(details["joined_base_n"])
    raw_ratio = float(len(raw_taxids)) / max(float(base_n), 1.0)
    use_raw = bool(
        base_n >= 15
        and float(details["joined_base_median_uaf"]) >= 0.50
        and float(details["probability_extra_mass_ratio"]) <= 0.03
        and raw_ratio >= 0.80
    )
    details.update(
        {
            "profile_strategy": strategy,
            "raw_unique_fallback": use_raw,
            "raw_unique95_n": len(raw_taxids),
            "raw_unique_to_base_ratio": raw_ratio,
            "raw_unique_fallback_rule": (
                "joined_base_n>=15;joined_base_median_uaf>=0.50;"
                "probability_extra_mass_ratio<=0.03;raw_unique95_n/joined_base_n>=0.80"
            ),
        }
    )
    return (raw_mask if use_raw else incumbent), panel_abundance_raw(features), details


def train_tables_from_dir(root: Path, pool: str) -> List[Path]:
    if pool == "train9":
        names = ["train.joined_features.tsv"]
    elif pool == "train12":
        names = ["train.joined_features.tsv", "test.joined_features.tsv"]
    else:
        raise ValueError(f"unsupported train pool: {pool}")

    paths = [root / name for name in names if (root / name).exists()]
    if not paths:
        raise FileNotFoundError(f"no joined feature training tables found in {root}")
    if pool == "train12" and len(paths) == 1:
        eprint(f"warning: {root}/test.joined_features.tsv is absent; using train.joined_features.tsv only")
    return paths


def load_training_features(
    train_dir: Optional[Path],
    train_tables: Sequence[Path],
    train_pool: str,
    scope: str,
    scope_by_taxid: Mapping[str, str],
    filter_training_scope: bool,
) -> pd.DataFrame:
    paths: List[Path] = []
    if train_dir is not None:
        paths.extend(train_tables_from_dir(train_dir, train_pool))
    paths.extend(train_tables)
    if not paths:
        raise SystemExit("--train-features or --train-table is required")

    frames = []
    for path in paths:
        eprint(f"loading training features: {path}")
        frame = pd.read_csv(path, sep="\t")
        frames.append(frame)
    train = pd.concat(frames, ignore_index=True, sort=False).fillna(0.0)
    if "label" not in train.columns:
        raise SystemExit("training feature table must contain a label column")
    if filter_training_scope:
        train = filter_scope(train, scope, scope_by_taxid)
    for col in MODEL_COLS:
        if col not in train.columns:
            train[col] = 0.0

    y = pd.to_numeric(train["label"], errors="coerce").fillna(0).to_numpy(dtype=int)
    if len(set(y.tolist())) < 2:
        raise SystemExit("training features must contain both positive and negative labels after scope filtering")
    return train


def fit_rf_hgb(train: pd.DataFrame):
    models = model_suite()
    x_train = train[MODEL_COLS].to_numpy(dtype=float)
    y_train = pd.to_numeric(train["label"], errors="coerce").fillna(0).to_numpy(dtype=int)
    fitted = []
    for name in ["rf", "hgb"]:
        model = models[name]
        eprint(f"fitting {name} on {len(train):,} candidate rows")
        model.fit(x_train, y_train)
        fitted.append(model)
    return fitted


def model_cache_scope_tag(scope: str, filter_training_scope: bool) -> str:
    return scope if filter_training_scope else "unfiltered"


def default_model_cache_paths(
    train_dir: Optional[Path],
    train_pool: str,
    scope: str,
    filter_training_scope: bool,
) -> list[Path]:
    if train_dir is None:
        return []
    tag = model_cache_scope_tag(scope, filter_training_scope)
    return [
        train_dir / f"minco_profile_rf_hgb.{train_pool}.{tag}.joblib",
        train_dir / f"minco_profile_rf_hgb.{train_pool}.joblib",
        train_dir / "minco_profile_rf_hgb.joblib",
    ]


def discover_model_cache(args: argparse.Namespace) -> Optional[Path]:
    explicit = getattr(args, "model_cache", None)
    if explicit is not None:
        if not explicit.is_file():
            raise SystemExit(f"--model-cache does not exist or is not a file: {explicit}")
        return explicit
    if getattr(args, "train_table", None):
        return None
    for path in default_model_cache_paths(
        getattr(args, "train_features", None),
        args.train_pool,
        args.scope,
        args.filter_training_scope,
    ):
        if path.is_file():
            return path
    return None


def validate_model_cache_bundle(
    bundle: object,
    path: Path,
    train_pool: str,
    scope: str,
    filter_training_scope: bool,
) -> Sequence[object]:
    if not isinstance(bundle, dict):
        raise SystemExit(f"{path}: model cache must contain a dictionary bundle")
    models = bundle.get("models")
    if not isinstance(models, (list, tuple)) or len(models) != 2:
        raise SystemExit(f"{path}: model cache must contain two fitted models")
    model_cols = list(bundle.get("model_cols", []))
    if model_cols != list(MODEL_COLS):
        raise SystemExit(f"{path}: model cache MODEL_COLS mismatch")
    if int(bundle.get("version", -1)) != MODEL_CACHE_VERSION:
        raise SystemExit(f"{path}: unsupported model cache version")
    expected = {
        "train_pool": train_pool,
        "scope": scope,
        "filter_training_scope": bool(filter_training_scope),
    }
    metadata = bundle.get("metadata", {})
    if not isinstance(metadata, dict):
        raise SystemExit(f"{path}: model cache metadata must be a dictionary")
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise SystemExit(
                f"{path}: model cache {key}={metadata.get(key)!r} does not match expected {value!r}"
            )
    return models


def load_model_cache(
    path: Path,
    train_pool: str,
    scope: str,
    filter_training_scope: bool,
) -> Sequence[object]:
    eprint(f"loading fitted RF/HGB model cache: {path}")
    bundle = joblib.load(path)
    return validate_model_cache_bundle(bundle, path, train_pool, scope, filter_training_scope)


def write_model_cache(
    path: Path,
    models: Sequence[object],
    train_pool: str,
    scope: str,
    filter_training_scope: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bundle = {
        "version": MODEL_CACHE_VERSION,
        "model_cols": list(MODEL_COLS),
        "models": list(models),
        "metadata": {
            "train_pool": train_pool,
            "scope": scope,
            "filter_training_scope": bool(filter_training_scope),
        },
    }
    joblib.dump(bundle, path)
    eprint(f"wrote fitted RF/HGB model cache: {path}")


def load_or_fit_models(
    args: argparse.Namespace,
    scope_by_taxid: Mapping[str, str],
) -> Sequence[object]:
    cache_path = discover_model_cache(args)
    if cache_path is not None:
        return load_model_cache(
            cache_path,
            args.train_pool,
            args.scope,
            args.filter_training_scope,
        )

    train = load_training_features(
        args.train_features,
        args.train_table,
        args.train_pool,
        args.scope,
        scope_by_taxid,
        args.filter_training_scope,
    )
    models = fit_rf_hgb(train)
    if getattr(args, "write_model_cache", None) is not None:
        write_model_cache(
            args.write_model_cache,
            models,
            args.train_pool,
            args.scope,
            args.filter_training_scope,
        )
    return models


def predict_rf_hgb(models: Sequence[object], features: pd.DataFrame) -> np.ndarray:
    for col in MODEL_COLS:
        if col not in features.columns:
            features[col] = 0.0
    x = features[MODEL_COLS].to_numpy(dtype=float)
    probs = [predict_probability(model, x) for model in models]
    return 0.5 * (probs[0] + probs[1])


def run_minco_pass(
    minco: str,
    ref: Path,
    reads: Path,
    out: Path,
    assign_mode: str,
    threads: int,
    pipecmd: str,
    density_block_ctx: Optional[int] = None,
    unique_sidecar_out: Optional[Path] = None,
    exact_split_sidecar_out: Optional[Path] = None,
) -> None:
    cmd = minco_pass_command(
        minco,
        ref,
        reads,
        out,
        assign_mode,
        threads,
        pipecmd,
        density_block_ctx,
        unique_sidecar_out,
        exact_split_sidecar_out,
    )
    log_path = out.with_suffix(out.suffix + ".log")
    eprint("running:", " ".join(cmd))
    with log_path.open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=log, check=True)


def minco_pass_command(
    minco: str,
    ref: Path,
    reads: Path,
    out: Path,
    assign_mode: str,
    threads: int,
    pipecmd: str,
    density_block_ctx: Optional[int] = None,
    unique_sidecar_out: Optional[Path] = None,
    exact_split_sidecar_out: Optional[Path] = None,
) -> list[str]:
    cmd = [
        minco,
        "ani",
        "-p",
        str(threads),
        "-r",
        str(ref),
        "--qraw",
        str(reads),
        "--query-density",
        "ref",
        "--abundance-est",
        "depth",
        "--readwise-profile-only",
        "--readwise-assign",
        assign_mode,
        "--readwise-ani",
        "zip-aaf",
        "-m0",
        "-f0",
        "-n0",
        "-t0",
        "-o",
        str(out),
    ]
    if density_block_ctx is not None:
        cmd.extend(["--density-block-ctx", str(density_block_ctx)])
    if unique_sidecar_out is not None:
        cmd.extend(["--readwise-unique-out", str(unique_sidecar_out)])
    if exact_split_sidecar_out is not None:
        cmd.extend(["--readwise-exact-split-out", str(exact_split_sidecar_out)])
    if pipecmd:
        cmd.extend(["--pipecmd", pipecmd])
    return cmd


def split_initial_pass_threads(total_threads: int) -> tuple[int, int]:
    if total_threads < 2:
        return total_threads, total_threads
    unique_threads = max(1, total_threads // 2)
    split_threads = max(1, total_threads - unique_threads)
    return unique_threads, split_threads


def run_minco_passes_parallel(
    minco: str,
    ref: Path,
    reads: Path,
    unique_out: Path,
    split_out: Path,
    threads: int,
    pipecmd: str,
    exact_split_sidecar_out: Optional[Path] = None,
) -> None:
    unique_threads, split_threads = split_initial_pass_threads(threads)
    if threads < 2:
        run_minco_pass(minco, ref, reads, unique_out, "best-diff-unique", threads, pipecmd)
        run_minco_pass(
            minco,
            ref,
            reads,
            split_out,
            "best-diff-split",
            threads,
            pipecmd,
            exact_split_sidecar_out=exact_split_sidecar_out,
        )
        return

    jobs = [
        (unique_out, "best-diff-unique", unique_threads, None),
        (split_out, "best-diff-split", split_threads, exact_split_sidecar_out),
    ]
    running: list[tuple[subprocess.Popen[bytes], object, list[str]]] = []
    try:
        for out, assign_mode, pass_threads, exact_sidecar in jobs:
            cmd = minco_pass_command(
                minco,
                ref,
                reads,
                out,
                assign_mode,
                pass_threads,
                pipecmd,
                exact_split_sidecar_out=exact_sidecar,
            )
            log = out.with_suffix(out.suffix + ".log").open("w")
            eprint(
                "running concurrent:",
                " ".join(cmd),
                f"(thread budget {pass_threads}/{threads})",
            )
            proc = subprocess.Popen(cmd, stdout=log, stderr=log)
            running.append((proc, log, cmd))

        failures: list[tuple[int, list[str]]] = []
        for proc, log, cmd in running:
            rc = proc.wait()
            log.close()
            if rc != 0:
                failures.append((rc, cmd))
        if failures:
            rc, cmd = failures[0]
            raise subprocess.CalledProcessError(rc, cmd)
    except BaseException:
        for proc, _log, _cmd in running:
            if proc.poll() is None:
                proc.terminate()
        for proc, log, _cmd in running:
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
            try:
                log.close()
            except Exception:
                pass
        raise


def use_same_stream_initial_passes(args: argparse.Namespace) -> bool:
    explicit = getattr(args, "same_stream_readwise_passes", None)
    if explicit is not None:
        return bool(explicit)
    if not getattr(args, "parallel_readwise_passes", True):
        return False
    return int(getattr(args, "threads", 1)) < 2


def use_same_stream_exact_split(args: argparse.Namespace) -> bool:
    explicit = getattr(args, "same_stream_exact_split", None)
    return bool(explicit) if explicit is not None else False


def run_initial_minco_passes(
    args: argparse.Namespace,
    unique_out: Path,
    split_out: Path,
    exact_split_sidecar_out: Optional[Path] = None,
) -> None:
    if use_same_stream_initial_passes(args):
        kwargs = {"unique_sidecar_out": unique_out}
        if exact_split_sidecar_out is not None:
            kwargs["exact_split_sidecar_out"] = exact_split_sidecar_out
        run_minco_pass(
            args.minco,
            args.ref,
            args.reads,
            split_out,
            "best-diff-split",
            args.threads,
            args.pipecmd,
            **kwargs,
        )
        return

    if args.parallel_readwise_passes:
        kwargs = {}
        if exact_split_sidecar_out is not None:
            kwargs["exact_split_sidecar_out"] = exact_split_sidecar_out
        run_minco_passes_parallel(
            args.minco,
            args.ref,
            args.reads,
            unique_out,
            split_out,
            args.threads,
            args.pipecmd,
            **kwargs,
        )
    else:
        run_minco_pass(args.minco, args.ref, args.reads, unique_out, "best-diff-unique", args.threads, args.pipecmd)
        kwargs = {}
        if exact_split_sidecar_out is not None:
            kwargs["exact_split_sidecar_out"] = exact_split_sidecar_out
        run_minco_pass(
            args.minco,
            args.ref,
            args.reads,
            split_out,
            "best-diff-split",
            args.threads,
            args.pipecmd,
            **kwargs,
        )


def ensure_feature_tables(args: argparse.Namespace) -> tuple[Path, Path, Optional[Path], Path]:
    if args.unique_table and args.split_table:
        workdir = args.workdir or args.out.parent
        workdir.mkdir(parents=True, exist_ok=True)
        return args.unique_table, args.split_table, None, workdir
    if args.unique_table or args.split_table:
        raise SystemExit("--unique-table and --split-table must be provided together")
    if not args.ref or not args.reads:
        raise SystemExit("provide either --unique-table/--split-table or both --ref and --reads")

    if args.workdir:
        workdir = args.workdir
        workdir.mkdir(parents=True, exist_ok=True)
        cleanup_dir: Optional[Path] = None
    else:
        cleanup_dir = Path(tempfile.mkdtemp(prefix="minco_calibrated_", dir="/tmp"))
        workdir = cleanup_dir

    unique_out = workdir / "minco.best_diff_unique.unfiltered.tsv"
    split_out = workdir / "minco.best_diff_split.unfiltered.tsv"
    exact_split_sidecar_out = (
        workdir / "minco.best_diff_split.exact.unfiltered.tsv"
        if args.strategy == "universal-auto-exact"
        and use_same_stream_exact_split(args)
        else None
    )
    if exact_split_sidecar_out is not None:
        run_initial_minco_passes(args, unique_out, split_out, exact_split_sidecar_out)
    else:
        run_initial_minco_passes(args, unique_out, split_out)
    return unique_out, split_out, cleanup_dir, workdir


def load_joined_predicted_features(
    unique_path: Path,
    split_path: Path,
    taxmap: Mapping[str, Mapping[str, str]],
    scope: str,
    scope_by_taxid: Mapping[str, str],
    models: Sequence[object],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    unique_rows = load_numeric_minco(unique_path, taxmap)
    split_rows = load_numeric_minco(split_path, taxmap)
    unique_rows = filter_scope(unique_rows, scope, scope_by_taxid)
    split_rows = filter_scope(split_rows, scope, scope_by_taxid)
    eprint(f"loaded candidates: unique rows={len(unique_rows):,}, split rows={len(split_rows):,}")
    features = joined_features(unique_rows, split_rows)
    for meta in [
        best_ref_metadata(unique_rows, "u"),
        best_ref_metadata(split_rows, "s"),
    ]:
        if len(meta.columns) > 1:
            features = features.merge(meta, on="taxid", how="left")
    features = filter_scope(features, scope, scope_by_taxid)
    eprint(f"joined taxid candidates: {len(features):,}")
    features["calibrated_probability"] = predict_rf_hgb(models, features)
    return unique_rows, split_rows, features


def taxid_name_map(taxmap: Mapping[str, Mapping[str, str]]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in taxmap.values():
        taxid = str(rec.get("taxid", ""))
        name = str(rec.get("species_name", ""))
        if taxid and name and taxid not in out:
            out[taxid] = name
    return out


def taxmap_has_gtdb_species_labels(taxmap: Mapping[str, Mapping[str, str]]) -> bool:
    for rec in taxmap.values():
        species = str(rec.get("species_name", "") or rec.get("taxid", ""))
        if species.startswith("s__"):
            return True
    return False


def write_output(
    out_path: Path,
    features: pd.DataFrame,
    threshold: float,
    train_pool: str,
    scope: str,
    names: Mapping[str, str],
    report_all: bool,
    strategy: str,
    call_mask: np.ndarray,
    abundance_raw: np.ndarray,
    strategy_details: Mapping[str, object],
    abundance_genus_xny_blend_alpha: float = DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA,
    abundance_feature_allocator_switch: str = ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF,
    abundance_feature_allocator_taxmap: Optional[Mapping[str, Mapping[str, str]]] = None,
    candidate_surface_rows: Optional[pd.DataFrame] = None,
    candidate_abundance_policy: str = CANDIDATE_ABUNDANCE_POLICY_ZERO,
) -> None:
    validate_candidate_abundance_policy(candidate_abundance_policy)
    work = features.copy()
    work["taxid"] = work["taxid"].astype(str)
    work["species_name"] = work["taxid"].map(lambda t: names.get(t, ""))
    work["calibrated_call"] = np.asarray(call_mask, dtype=bool)
    work["profile_strategy"] = strategy
    for key, value in strategy_details.items():
        if key == "profile_strategy":
            continue
        work[key] = value
    if "raw_unique_fallback" not in work.columns:
        work["raw_unique_fallback"] = False
    work["calibrated_threshold"] = threshold
    work["train_pool"] = train_pool
    work["scope"] = scope

    abundance = np.asarray(abundance_raw, dtype=float)
    abundance = np.nan_to_num(abundance, nan=0.0, posinf=0.0, neginf=0.0)
    abundance = np.where(abundance > 0.0, abundance, 0.0)
    abundance_blend_applied = (
        strategy in ABUNDANCE_GENUS_XNY_BLEND_STRATEGIES
        and abundance_genus_xny_blend_alpha > 0.0
    )
    if abundance_blend_applied:
        abundance = genus_xny_blended_abundance_raw(
            work,
            call_mask,
            abundance,
            abundance_genus_xny_blend_alpha,
        )
    work["abundance_rule"] = (
        "base_raw_then_within_genus_xny_blend"
        if abundance_blend_applied
        else "base_raw"
    )
    work["abundance_genus_xny_blend_alpha"] = (
        abundance_genus_xny_blend_alpha if abundance_blend_applied else 0.0
    )
    work["abundance_genus_xny_quality"] = (
        "clip(s_XnY_ctx_max/1000,0,1)" if abundance_blend_applied else ""
    )
    work["calibrated_abundance_raw"] = abundance
    work["calibrated_abundance"] = 0.0

    reported_ani, reported_ani_source = reported_ani_values(work)
    work["reported_ani"] = reported_ani
    work["reported_ani_source"] = reported_ani_source

    call_bool = work["calibrated_call"].astype(bool)
    if "candidate_rescue_added" in work.columns:
        rescue_added = work["candidate_rescue_added"].astype(bool)
    else:
        rescue_added = pd.Series(False, index=work.index)
    base_call_bool = call_bool & ~rescue_added
    base_called_mass = float(numeric(work, "calibrated_abundance_raw").loc[base_call_bool].sum())
    if base_called_mass <= 0.0:
        base_called_mass = 1.0
    if candidate_abundance_policy == CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2:
        rescue_norm_mass = candidate_normalized_depth_mass(work, candidate_abundance_policy)
        if "candidate_abundance_norm_mass" not in work.columns:
            work["candidate_abundance_norm_mass"] = 0.0
        work.loc[rescue_added, "candidate_abundance_norm_mass"] = rescue_norm_mass[rescue_added.to_numpy()]
        work.loc[rescue_added, "calibrated_abundance_raw"] = (
            work.loc[rescue_added, "candidate_abundance_norm_mass"].astype(float) * base_called_mass
        )

    work = append_accession_candidate_surface_rows(
        work,
        candidate_surface_rows if candidate_surface_rows is not None else pd.DataFrame(),
        candidate_abundance_policy,
        base_called_mass,
    )

    call_bool = work["calibrated_call"].astype(bool)
    allocated_raw, allocator_details = guarded_feature_allocator_abundance_raw(
        work,
        call_bool.to_numpy(dtype=bool),
        numeric(work, "calibrated_abundance_raw").to_numpy(dtype=float),
        abundance_feature_allocator_switch,
        abundance_feature_allocator_taxmap,
    )
    for key, value in allocator_details.items():
        work[key] = value
    if bool(allocator_details.get("abundance_feature_allocator_applied", False)):
        work["abundance_rule"] = (
            work["abundance_rule"].astype(str) + "_then_guarded_feature_allocator"
        )
    work["calibrated_abundance_raw"] = allocated_raw
    raw_abundance = numeric(work, "calibrated_abundance_raw")
    called_mass = float(raw_abundance.loc[call_bool].sum())
    if called_mass > 0.0:
        work["calibrated_abundance"] = np.where(
            call_bool,
            raw_abundance / called_mass,
            0.0,
        )
    else:
        work["calibrated_abundance"] = 0.0

    if not report_all:
        work = work.loc[work["calibrated_call"]].copy()

    first_cols = [
        "taxid",
        "species_name",
        "calibrated_call",
        "calibrated_probability",
        "calibrated_threshold",
        "calibrated_abundance",
        "calibrated_abundance_raw",
        "reported_ani",
        "reported_ani_source",
        "profile_strategy",
        "raw_unique_fallback",
        "train_pool",
        "scope",
    ]
    remaining = [c for c in work.columns if c not in first_cols]
    work = work[first_cols + remaining]
    work = work.sort_values(
        ["calibrated_call", "calibrated_probability", "calibrated_abundance_raw"],
        ascending=[False, False, False],
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    work.to_csv(out_path, sep="\t", index=False)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Run the MinCO joined unique/split RF-HGB calibrated profiler.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    io = ap.add_argument_group("input")
    io.add_argument("-r", "--ref", type=Path, help="Reference MinCO sketch; required when FASTQ is supplied.")
    io.add_argument("--reads", "--qraw", dest="reads", type=Path, help="Raw FASTA/FASTQ input.")
    io.add_argument("--unique-table", type=Path, help="Precomputed best-diff-unique unfiltered readwise table.")
    io.add_argument("--split-table", type=Path, help="Precomputed best-diff-split unfiltered readwise table.")
    io.add_argument("--taxmap", type=Path, required=True, help="Species-level reference taxmap TSV.")
    io.add_argument("-o", "--out", type=Path, required=True, help="Output calibrated species table.")

    run = ap.add_argument_group("minco execution")
    run.add_argument("--minco", default=str(REPO_ROOT / "bin/minco"), help="MinCO executable.")
    run.add_argument("-p", "--threads", type=int, default=1, help="Total thread budget for raw-read MinCO table generation. Auto scheduling uses same-stream unique sidecar at -p1 and concurrent unique/split passes at -p>=2.")
    run.add_argument("--pipecmd", default="", help="Optional MinCO --pipecmd for streaming reads.")
    run.add_argument("--workdir", type=Path, help="Directory for generated unique/split tables. Defaults to /tmp.")
    run.set_defaults(same_stream_readwise_passes=None)
    run.set_defaults(same_stream_exact_split=None)
    run.set_defaults(parallel_readwise_passes=True)
    run.add_argument("--same-stream-readwise-passes", dest="same_stream_readwise_passes", action="store_true", default=argparse.SUPPRESS, help="Force one split pass and write the unique table with --readwise-unique-out.")
    run.add_argument("--legacy-dual-readwise-passes", dest="same_stream_readwise_passes", action="store_false", default=argparse.SUPPRESS, help="Force the dual-pass initial table generation path.")
    run.add_argument("--sequential-readwise-passes", dest="parallel_readwise_passes", action="store_false", default=argparse.SUPPRESS, help="Run initial best-diff-unique and best-diff-split MinCO passes sequentially. This disables auto same-stream selection.")
    run.add_argument("--same-stream-exact-split", dest="same_stream_exact_split", action="store_true", default=argparse.SUPPRESS, help="During raw-read universal-auto-exact runs, write the exact per-read split table as a sidecar of the initial split pass so a later auto-exact trigger can avoid a third MinCO pass.")
    run.add_argument("--legacy-exact-split-rerun", dest="same_stream_exact_split", action="store_false", default=argparse.SUPPRESS, help="Disable eager exact-split sidecar generation and keep the legacy third-pass exact rerun path.")

    model = ap.add_argument_group("calibration")
    model.add_argument("--train-features", type=Path, help="Directory with train.joined_features.tsv and optional test.joined_features.tsv.")
    model.add_argument("--train-table", action="append", type=Path, default=[], help="Additional joined feature table; may be repeated.")
    model.add_argument("--model-cache", type=Path, help="Load a fitted RF/HGB model bundle instead of fitting from training features. If absent, a compatible cache is auto-discovered in --train-features when present.")
    model.add_argument("--write-model-cache", type=Path, help="After fitting from training features, write a fitted RF/HGB model bundle for future runs.")
    model.add_argument("--train-pool", choices=["train12", "train9"], default="train12", help="train12 includes strainmadness test.joined_features.tsv when present; train9 excludes it for strain holdout testing.")
    model.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help="RF/HGB average probability call threshold.")
    model.add_argument(
        "--strategy",
        choices=["probability", "adaptive-sub95", "universal", "universal-auto-exact"],
        default="universal-auto-exact",
        help="Final call strategy. Default: universal-auto-exact, the F1-priority guarded universal gate with optional exact split rerun. Use probability for legacy RF/HGB threshold-only output, adaptive-sub95 for the older adaptive guard, or universal to skip the exact split rerun.",
    )
    model.add_argument("--exact-split-trigger", type=float, default=EXACT_SPLIT_TRIGGER, help="For --strategy universal-auto-exact, rerun split with --density-block-ctx 0 when block-mode probability_extra_mass_ratio is at or below this value and the high-uAF/raw-unique guard passes.")
    model.add_argument("--exact-split-low-extra-mode", choices=["skip", "allow"], default=EXACT_SPLIT_LOW_EXTRA_MODE, help="For --strategy universal-auto-exact, skip the exact split rerun when the block-mode low-extra split rescue already added candidates. Use allow to keep the older exact behavior.")
    model.add_argument("--abundance-genus-xny-blend-alpha", type=float, default=DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA, help="Experimental abundance-only option for universal strategies: blend this fraction of each called genus' current raw abundance toward split-XnY-weighted within-genus mass. The default 0 keeps the selected default abundance unchanged.")
    model.add_argument(
        "--abundance-feature-allocator-switch",
        choices=[
            ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF,
            ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002,
            ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230,
        ],
        default=ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF,
        help="Experimental abundance-only guarded feature allocator. Default off keeps the selected default abundance unchanged. guarded-genus-hit-breadth-a002 applies the cached candidate rule when its output-derived guard passes; guarded-genus-hit-breadth-a002-xny230 adds the refined s_XnY median guard from the combined cached audit.",
    )
    model.add_argument("--adaptive-call-filter-switch", choices=[ADAPTIVE_CALL_FILTER_SWITCH_OFF, ADAPTIVE_CALL_FILTER_SWITCH_LOPO_MIN_XNY25], default=ADAPTIVE_CALL_FILTER_SWITCH_OFF, help="Experimental F1-priority post-call filter switch. Default off keeps the selected default unchanged. lopo-min-xny25 applies the current candidate rule from the adaptive call-filter audit.")
    model.add_argument("--candidate-rescue-switch", choices=[CANDIDATE_RESCUE_SWITCH_OFF, CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70], default=CANDIDATE_RESCUE_SWITCH_OFF, help="Experimental F1-priority rescue switch for strong uncalled emitted-profile candidates. Default off keeps the selected default unchanged.")
    model.add_argument(
        "--candidate-surface-switch",
        choices=[
            CANDIDATE_SURFACE_SWITCH_OFF,
            CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
            CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI93_XNY650_BR20,
            CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20,
        ],
        default=CANDIDATE_SURFACE_SWITCH_OFF,
        help=(
            "Experimental output-only accession-level candidate surface. "
            "Default off keeps the selected default unchanged."
        ),
    )
    model.add_argument("--candidate-abundance-policy", choices=[CANDIDATE_ABUNDANCE_POLICY_ZERO, CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2], default=CANDIDATE_ABUNDANCE_POLICY_ZERO, help="Experimental abundance policy for rows added by candidate rescue/surface switches. Default zero preserves the selected default abundance.")
    model.add_argument("--candidate-surface-taxmap", type=Path, help="Optional accession-level taxmap used only to label experimental candidate-surface rows. The normal --taxmap still drives model features.")
    model.add_argument("--candidate-surface-max-called-species", type=int, default=0, help="Experimental guard for candidate-surface rows. When >0, skip accession-level surface additions if the pre-surface called species count exceeds this value. Default 0 disables the guard.")
    model.add_argument("--scope", choices=["all", "prokaryote", "bacteria", "archaea", "virus"], default="all", help="Taxonomic scope for output candidates.")
    model.add_argument("--filter-training-scope", action="store_true", help="Also restrict training rows to --scope. By default the RF/HGB gate uses the full universal training pool.")
    model.add_argument("--report-all", action="store_true", help="Write all candidates instead of called species only.")
    return ap.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    if not (0.0 <= args.threshold <= 1.0):
        raise SystemExit("--threshold must be in 0..1")
    if args.exact_split_trigger < 0.0:
        raise SystemExit("--exact-split-trigger must be >= 0")
    if not (0.0 <= args.abundance_genus_xny_blend_alpha <= 1.0):
        raise SystemExit("--abundance-genus-xny-blend-alpha must be in 0..1")

    taxmap = parse_species_taxmap(args.taxmap)
    scope_by_taxid = taxid_scope_map(taxmap)
    names = taxid_name_map(taxmap)
    candidate_surface_switch = args.candidate_surface_switch
    candidate_surface_taxmap_path = str(args.candidate_surface_taxmap or "")
    lazy_candidate_surface_taxmap: Optional[Path] = None
    if args.candidate_surface_taxmap:
        lazy_candidate_surface_taxmap = args.candidate_surface_taxmap
        candidate_surface_taxmap = {}
        candidate_surface_names = names
    elif taxmap_has_gtdb_species_labels(taxmap):
        candidate_surface_taxmap = taxmap
        candidate_surface_names = names
        candidate_surface_taxmap_path = str(args.taxmap)
    else:
        candidate_surface_taxmap = {}
        candidate_surface_names = {}
        if candidate_surface_switch != CANDIDATE_SURFACE_SWITCH_OFF:
            eprint(
                "warning: candidate surface disabled because --candidate-surface-taxmap "
                "is missing and --taxmap lacks GTDB-style species labels"
            )
            candidate_surface_switch = CANDIDATE_SURFACE_SWITCH_OFF

    unique_path, split_path, cleanup_dir, workdir = ensure_feature_tables(args)

    models = load_or_fit_models(args, scope_by_taxid)

    unique_rows, split_rows, features = load_joined_predicted_features(
        unique_path,
        split_path,
        taxmap,
        args.scope,
        scope_by_taxid,
        models,
    )
    strategy_for_calls = "universal" if args.strategy == "universal-auto-exact" else args.strategy
    auto_exact_details: dict[str, object] = {}
    if args.strategy == "universal-auto-exact":
        _prelim_call, _prelim_abundance, prelim_details = strategy_call_and_abundance(
            features,
            unique_rows,
            names,
            "universal",
            args.threshold,
        )
        p_extra = float(prelim_details.get("probability_extra_mass_ratio", 1.0))
        p_extra_passed = p_extra <= args.exact_split_trigger
        exact_guard_passed = exact_split_guard_passes(prelim_details)
        low_extra_added_n = int(prelim_details.get("low_extra_split_rescue_added_n", 0))
        low_extra_skip_active = args.exact_split_low_extra_mode == "skip" and low_extra_added_n > 0
        low_extra_passed = not low_extra_skip_active
        use_exact = p_extra_passed and exact_guard_passed and low_extra_passed
        auto_exact_details = {
            "auto_exact_split_trigger": args.exact_split_trigger,
            "auto_exact_split_block_p_extra_mass_ratio": p_extra,
            "auto_exact_split_p_extra_gate_passed": p_extra_passed,
            "auto_exact_split_guard_passed": exact_guard_passed,
            "auto_exact_split_low_extra_mode": args.exact_split_low_extra_mode,
            "auto_exact_split_low_extra_added_n": low_extra_added_n,
            "auto_exact_split_low_extra_gate_passed": low_extra_passed,
            "auto_exact_split_guard_rule": (
                "probability_extra_mass_ratio<=trigger and "
                f"(joined_base_median_uaf<{EXACT_SPLIT_MEDIAN_UAF_GUARD:.2f} "
                f"or raw_unique_to_base_ratio>={EXACT_SPLIT_RAW_UNIQUE_RATIO_GUARD:.2f}) "
                "and not(block low-extra split rescue active) unless "
                "--exact-split-low-extra-mode allow"
            ),
            "auto_exact_split_requested": use_exact,
            "auto_exact_split_used": False,
            "auto_exact_split_path": "",
            "auto_exact_split_sidecar_requested": use_same_stream_exact_split(args),
            "auto_exact_split_source": "",
        }
        if p_extra_passed and not exact_guard_passed:
            auto_exact_details["auto_exact_split_unavailable_reason"] = "guard_high_uaf_low_raw_unique_support"
        elif p_extra_passed and exact_guard_passed and not low_extra_passed:
            auto_exact_details["auto_exact_split_unavailable_reason"] = "block_low_extra_split_rescue_already_active"
        if use_exact:
            if not args.ref or not args.reads:
                eprint("warning: universal-auto-exact trigger fired but --ref/--reads are unavailable; using supplied split table")
                auto_exact_details["auto_exact_split_unavailable_reason"] = "missing_ref_or_reads"
            else:
                exact_split_path = workdir / "minco.best_diff_split.exact.unfiltered.tsv"
                if not exact_split_path.exists():
                    run_minco_pass(
                        args.minco,
                        args.ref,
                        args.reads,
                        exact_split_path,
                        "best-diff-split",
                        args.threads,
                        args.pipecmd,
                        density_block_ctx=0,
                    )
                    auto_exact_details["auto_exact_split_source"] = "rerun"
                else:
                    auto_exact_details["auto_exact_split_source"] = "sidecar"
                split_path = exact_split_path
                unique_rows, split_rows, features = load_joined_predicted_features(
                    unique_path,
                    split_path,
                    taxmap,
                    args.scope,
                    scope_by_taxid,
                    models,
                )
                auto_exact_details["auto_exact_split_used"] = True
                auto_exact_details["auto_exact_split_path"] = str(exact_split_path)

    call_mask, abundance_raw, strategy_details = strategy_call_and_abundance(
        features,
        unique_rows,
        names,
        strategy_for_calls,
        args.threshold,
    )
    strategy_details.update(auto_exact_details)
    call_mask, call_filter_details = apply_adaptive_call_filter_switch(
        features,
        call_mask,
        args.adaptive_call_filter_switch,
    )
    strategy_details.update(call_filter_details)
    call_mask, abundance_raw, candidate_rescue_details = apply_candidate_rescue_switch(
        features,
        call_mask,
        abundance_raw,
        args.candidate_rescue_switch,
        args.candidate_abundance_policy,
    )
    strategy_details.update(candidate_rescue_details)
    if (
        lazy_candidate_surface_taxmap is not None
        and candidate_surface_switch != CANDIDATE_SURFACE_SWITCH_OFF
        and candidate_surface_has_numeric_candidates(unique_rows, split_rows, candidate_surface_switch)
    ):
        candidate_surface_taxmap = parse_species_taxmap(lazy_candidate_surface_taxmap)
        candidate_surface_names = taxid_name_map(candidate_surface_taxmap)
    candidate_surface_rows, candidate_surface_details = build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        call_mask,
        candidate_surface_switch,
        candidate_surface_taxmap,
        candidate_surface_names,
        args.candidate_abundance_policy,
        args.candidate_surface_max_called_species,
    )
    candidate_surface_details["candidate_surface_taxmap"] = candidate_surface_taxmap_path
    strategy_details.update(candidate_surface_details)
    write_output(
        args.out,
        features,
        args.threshold,
        args.train_pool,
        args.scope,
        names,
        args.report_all,
        args.strategy,
        call_mask,
        abundance_raw,
        strategy_details,
        args.abundance_genus_xny_blend_alpha,
        args.abundance_feature_allocator_switch,
        candidate_surface_taxmap,
        candidate_surface_rows,
        args.candidate_abundance_policy,
    )
    eprint(f"wrote {args.out}")
    if cleanup_dir:
        eprint(f"generated tables kept in {cleanup_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
