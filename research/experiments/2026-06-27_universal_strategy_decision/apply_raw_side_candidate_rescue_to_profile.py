#!/usr/bin/env python3
"""Append accession-level raw-side rescue candidates to a MinCO profile.

This is an experimental postprocessor for the universal-strategy decision note.
It preserves distinct GTDB species/accessions that are lost when the calibrated
wrapper emits one row per NCBI taxid. Added candidates carry zero abundance
mass, so the output tests call-set recovery only.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import sweep_hmp_raw_candidate_rescue as raw_rescue


ANI_MIN = 0.90
XNY_MIN = 100.0
BREADTH_MIN = 0.01
REAL_AF_MIN = 0.70
SWITCH_NAME = "raw-side-ani90-xny100-br01-af70-zero-mass"

RAW_COLS = {
    "gtdb_species",
    "gtdb_mapping_method",
    "accession",
    "Ref",
    "Ref_annotation",
    "raw_mode",
    "raw_path",
    "ANI",
    "XnY_ctx",
    "Ref_breadth",
    "Real_min_align_fraction",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_cv",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def metadata_for_accession(accession: object, mapper: raw_audit.RawMapper) -> Mapping[str, str]:
    return mapper.by_accession.get(str(accession or ""), {})


def current_called_species(profile: pd.DataFrame, panel: str, mapper: raw_audit.RawMapper) -> set[str]:
    if "calibrated_call" not in profile.columns:
        return set()
    called = profile.loc[bool_series(profile["calibrated_call"])].copy()
    out: set[str] = set()
    for row in called.itertuples(index=False):
        for col in ["s_best_accession", "u_best_accession"]:
            accession = str(getattr(row, col, "") or "")
            if not accession or accession.lower() == "nan":
                continue
            species = raw_rescue.gtdb_from_accession(panel, accession, mapper)
            if species:
                out.add(species)
                break
    return out


def load_rescue_candidates(raw_best_path: Path, called_species: set[str]) -> pd.DataFrame:
    raw = pd.read_csv(raw_best_path, sep="\t", usecols=lambda col: col in RAW_COLS, low_memory=False)
    for col in ["ANI", "XnY_ctx", "Ref_breadth", "Real_min_align_fraction"]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
    mask = (
        raw["gtdb_species"].astype(str).astype(bool)
        & ~raw["gtdb_species"].astype(str).isin(called_species)
        & (raw["ANI"] >= ANI_MIN)
        & (raw["XnY_ctx"] >= XNY_MIN)
        & (raw["Ref_breadth"] >= BREADTH_MIN)
        & (raw["Real_min_align_fraction"] >= REAL_AF_MIN)
    )
    candidates = raw.loc[mask].copy()
    if candidates.empty:
        return candidates
    candidates = candidates.sort_values(
        ["gtdb_species", "XnY_ctx", "Real_min_align_fraction", "ANI", "Ref_breadth"],
        ascending=[True, False, False, False, False],
        kind="mergesort",
    )
    return candidates.drop_duplicates("gtdb_species", keep="first")


def blank_output_row(columns: list[str]) -> dict[str, object]:
    return {col: "" for col in columns}


def make_added_row(
    columns: list[str],
    raw_row: object,
    panel: str,
    sample: int,
    mapper: raw_audit.RawMapper,
) -> dict[str, object]:
    accession = str(getattr(raw_row, "accession", "") or "")
    gtdb_species = str(getattr(raw_row, "gtdb_species", "") or "")
    meta = metadata_for_accession(accession, mapper)
    taxid = str(meta.get("ncbi_species_taxid") or meta.get("ncbi_taxid") or gtdb_species)
    row = blank_output_row(columns)
    row.update(
        {
            "taxid": taxid,
            "species_name": gtdb_species,
            "calibrated_call": True,
            "calibrated_probability": 0.0,
            "calibrated_threshold": 0.0,
            "calibrated_abundance": 0.0,
            "calibrated_abundance_raw": 0.0,
            "reported_ani": finite(getattr(raw_row, "ANI", 0.0)),
            "reported_ani_source": "raw_side_ANI",
            "profile_strategy": "raw_side_candidate_rescue",
            "raw_unique_fallback": False,
            "scope": "bacteria",
            "s_rows": 1,
            "s_direct_call": 1,
            "s_relaxed_call": 1,
            "s_ANI_max": finite(getattr(raw_row, "ANI", 0.0)),
            "s_XnY_ctx_max": finite(getattr(raw_row, "XnY_ctx", 0.0)),
            "s_Real_min_align_fraction_max": finite(
                getattr(raw_row, "Real_min_align_fraction", 0.0)
            ),
            "s_Ref_breadth_max": finite(getattr(raw_row, "Ref_breadth", 0.0)),
            "s_Ref_mean_depth_max": finite(getattr(raw_row, "Ref_mean_depth", 0.0)),
            "s_Ref_hit_mean_depth_max": finite(getattr(raw_row, "Ref_hit_mean_depth", 0.0)),
            "s_Ref_depth_cv_min": finite(getattr(raw_row, "Ref_depth_cv", 0.0)),
            "s_Ref_zip_af_max": finite(getattr(raw_row, "Ref_zip_af", 0.0)),
            "s_Ref_zip_aaf_ani_max": finite(getattr(raw_row, "Ref_zip_aaf_ani", 0.0)),
            "s_best_accession": accession,
            "s_best_ref": getattr(raw_row, "Ref", ""),
            "s_best_ref_annotation": getattr(raw_row, "Ref_annotation", ""),
            "candidate_rescue_added": False,
            "raw_side_rescue_added": True,
            "raw_side_rescue_switch": SWITCH_NAME,
            "raw_side_rescue_panel": panel,
            "raw_side_rescue_sample": sample,
            "raw_side_rescue_gtdb_species": gtdb_species,
            "raw_side_rescue_accession": accession,
            "raw_side_rescue_mapping_method": getattr(raw_row, "gtdb_mapping_method", ""),
            "raw_side_rescue_zero_mass": True,
            "raw_side_rescue_ani_min": ANI_MIN,
            "raw_side_rescue_xny_min": XNY_MIN,
            "raw_side_rescue_breadth_min": BREADTH_MIN,
            "raw_side_rescue_real_af_min": REAL_AF_MIN,
        }
    )
    return row


def apply_raw_side_rescue(
    profile_path: Path,
    raw_best_path: Path,
    panel: str,
    sample: int,
    out_path: Path,
    mapper: raw_audit.RawMapper | None = None,
) -> dict[str, object]:
    mapper = mapper or raw_audit.RawMapper()
    profile = pd.read_csv(profile_path, sep="\t", low_memory=False)
    for col in [
        "raw_side_rescue_added",
        "raw_side_rescue_switch",
        "raw_side_rescue_panel",
        "raw_side_rescue_sample",
        "raw_side_rescue_gtdb_species",
        "raw_side_rescue_accession",
        "raw_side_rescue_mapping_method",
        "raw_side_rescue_zero_mass",
        "raw_side_rescue_ani_min",
        "raw_side_rescue_xny_min",
        "raw_side_rescue_breadth_min",
        "raw_side_rescue_real_af_min",
        "raw_side_rescue_added_n",
    ]:
        if col not in profile.columns:
            profile[col] = False if col in {"raw_side_rescue_added", "raw_side_rescue_zero_mass"} else ""
    profile["raw_side_rescue_added"] = False
    called_species = current_called_species(profile, panel, mapper)
    candidates = load_rescue_candidates(raw_best_path, called_species)
    columns = list(profile.columns)
    added_rows = [
        make_added_row(columns, row, panel, sample, mapper)
        for row in candidates.itertuples(index=False)
    ]
    profile["raw_side_rescue_switch"] = SWITCH_NAME
    profile["raw_side_rescue_zero_mass"] = True
    profile["raw_side_rescue_ani_min"] = ANI_MIN
    profile["raw_side_rescue_xny_min"] = XNY_MIN
    profile["raw_side_rescue_breadth_min"] = BREADTH_MIN
    profile["raw_side_rescue_real_af_min"] = REAL_AF_MIN
    profile["raw_side_rescue_added_n"] = len(added_rows)
    if added_rows:
        out = pd.concat(
            [profile, pd.DataFrame(added_rows, columns=columns)],
            ignore_index=True,
            sort=False,
        )
    else:
        out = profile.copy()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return {
        "profile": str(profile_path),
        "raw_best": str(raw_best_path),
        "out": str(out_path),
        "panel": panel,
        "sample": sample,
        "input_rows": int(len(profile)),
        "added_rows": int(len(added_rows)),
        "output_rows": int(len(out)),
        "called_species_before": int(len(called_species)),
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", type=Path, required=True)
    ap.add_argument("--raw-best", type=Path, required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--sample", type=int, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    result = apply_raw_side_rescue(
        args.profile,
        args.raw_best,
        args.panel,
        args.sample,
        args.out,
    )
    print(pd.DataFrame([result]).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
