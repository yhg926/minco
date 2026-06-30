#!/usr/bin/env python3
"""Diagnose sample8 extra calls using the official HMP GTDB mapping policy."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
REPO_ROOT = NOTE_DIR.parents[2]

sys.path.insert(0, str(NOTE_DIR))

import score_hmp_gtdb_source_abundance as hmp_score  # noqa: E402


SAMPLE = 8
PROFILE = Path("/tmp/minco_current_code_hmp_airskin8_20260627/minco_sample8_current_default.tsv")
UNIQUE = Path("/tmp/minco_current_code_hmp_airskin8_20260627/work/minco.best_diff_unique.unfiltered.tsv")
SPLIT = Path("/tmp/minco_current_code_hmp_airskin8_20260627/work/minco.best_diff_split.unfiltered.tsv")


FEATURES = [
    "calibrated_probability",
    "calibrated_abundance",
    "reported_ani",
    "u_XnY_ctx_max",
    "s_XnY_ctx_max",
    "u_Ref_breadth_max",
    "s_Ref_breadth_max",
    "u_Real_min_align_fraction_max",
    "s_Real_min_align_fraction_max",
    "u_Ref_zip_af_max",
    "s_Ref_zip_af_max",
    "u_Ref_hit_mean_depth_max",
    "s_Ref_hit_mean_depth_max",
    "split_unique_xny_ratio",
    "split_unique_breadth_ratio",
    "probability_extra_mass_ratio",
    "joined_base_median_uaf",
]

FLAGS = [
    "tail_rescue_added",
    "low_extra_split_rescue_added",
    "candidate_rescue_added",
    "raw_unique_fallback",
    "auto_exact_split_used",
]


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def clean_numeric(value: object) -> float:
    out = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(out) or math.isinf(float(out)):
        return 0.0
    return float(out)


def map_profile_row(
    row: pd.Series,
    best_ref_species_by_taxid: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[str, str]:
    for col in ["s_best_accession", "u_best_accession"]:
        species, method = hmp_score.gtdb_from_accession(str(row.get(col, "")), by_accession, by_core)
        if species:
            return species, f"best_ref_{col}:{method}"
    species = best_ref_species_by_taxid.get(str(row.get("taxid", "")), "")
    if species:
        return species, "raw_taxid_fallback"
    return "", "unmapped"


def summarize_features(rows: pd.DataFrame) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for status, part in rows.groupby("call_status", dropna=False):
        record: dict[str, object] = {
            "call_status": status,
            "rows": int(len(part)),
            "species": int(part["gtdb_species"].nunique()),
            "abundance_sum": float(part["calibrated_abundance"].sum()),
        }
        for col in FEATURES:
            if col not in part.columns:
                continue
            vals = pd.to_numeric(part[col], errors="coerce").replace([math.inf, -math.inf], pd.NA).dropna()
            if vals.empty:
                continue
            record[f"{col}_median"] = float(vals.median())
            record[f"{col}_min"] = float(vals.min())
            record[f"{col}_max"] = float(vals.max())
        for col in FLAGS:
            if col in part.columns:
                record[f"{col}_true"] = int(part[col].map(as_bool).sum())
        out.append(record)
    return out


def filter_mask(calls: pd.DataFrame, method: str) -> pd.Series:
    prob = pd.to_numeric(calls.get("calibrated_probability", 0.0), errors="coerce").fillna(0.0)
    breadth = pd.to_numeric(calls.get("s_Ref_breadth_max", 0.0), errors="coerce").fillna(0.0)
    xny = pd.to_numeric(calls.get("s_XnY_ctx_max", 0.0), errors="coerce").fillna(0.0)
    realaf = pd.to_numeric(calls.get("s_Real_min_align_fraction_max", 0.0), errors="coerce").fillna(0.0)
    ani = pd.to_numeric(calls.get("reported_ani", 0.0), errors="coerce").fillna(0.0)
    if method == "current_calls":
        return pd.Series(True, index=calls.index)
    if method == "min_xny_ge_25":
        return xny >= 25
    if method == "min_breadth_ge_0.05":
        return breadth >= 0.05
    if method == "drop_prob_lt_0.5_and_xny_lt_500":
        return ~((prob < 0.5) & (xny < 500))
    if method == "drop_prob_lt_0.7_and_breadth_lt_0.2":
        return ~((prob < 0.7) & (breadth < 0.2))
    if method == "drop_prob_lt_0.5_and_breadth_lt_0.05":
        return ~((prob < 0.5) & (breadth < 0.05))
    if method == "drop_prob_lt_0.5_and_realaf_lt_0.3":
        return ~((prob < 0.5) & (realaf < 0.3))
    if method == "drop_ani_lt_0.95_and_breadth_lt_0.05":
        return ~((ani < 0.95) & (breadth < 0.05))
    raise ValueError(f"unknown filter method: {method}")


def score_filtered_calls(sample_id: int, calls: pd.DataFrame, truth_df: pd.DataFrame) -> pd.DataFrame:
    methods = [
        "current_calls",
        "min_xny_ge_25",
        "min_breadth_ge_0.05",
        "drop_prob_lt_0.5_and_xny_lt_500",
        "drop_prob_lt_0.7_and_breadth_lt_0.2",
        "drop_prob_lt_0.5_and_breadth_lt_0.05",
        "drop_prob_lt_0.5_and_realaf_lt_0.3",
        "drop_ani_lt_0.95_and_breadth_lt_0.05",
    ]
    rows = []
    for method in methods:
        mask = filter_mask(calls, method)
        kept = calls.loc[mask & calls["gtdb_species"].astype(bool)].copy()
        pred = hmp_score.taxid_score.collapse_prediction(
            [
                (
                    str(row.gtdb_species),
                    clean_numeric(row.calibrated_abundance),
                    clean_numeric(row.reported_ani),
                )
                for row in kept.itertuples(index=False)
            ]
        )
        score = hmp_score.taxid_score.score_prediction(
            sample_id,
            method,
            pred,
            truth_df,
            {
                "kept_rows": int(len(kept)),
                "removed_rows": int(len(calls) - len(kept)),
                "removed_tp_rows": int((calls.loc[~mask, "call_status"] == "TP").sum()),
                "removed_fp_rows": int((calls.loc[~mask, "call_status"] == "FP").sum()),
            },
        )
        rows.append(score)
    return pd.DataFrame(rows)


def main() -> int:
    required = [PROFILE, UNIQUE, SPLIT, hmp_score.TRUTH_ROOT / "genome_to_id.tsv", hmp_score.TAXMAP]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required inputs:\n" + "\n".join(map(str, missing)))

    RESULTS.mkdir(parents=True, exist_ok=True)
    hmp_score.ensure_sample_record(SAMPLE)
    by_accession, by_core = hmp_score.truth.load_gtdb_metadata(hmp_score.truth.GTDB_METADATA)
    source_accessions = hmp_score.load_source_accessions(hmp_score.TRUTH_ROOT / "genome_to_id.tsv")
    truth_df, quality, _sources = hmp_score.build_source_truth(SAMPLE, source_accessions, by_accession, by_core)
    truth_species = set(truth_df["gtdb_species"].astype(str))

    taxmap = hmp_score.parse_species_taxmap(hmp_score.TAXMAP)
    best_raw_species, raw_stats = hmp_score.best_raw_ref_species_by_taxid(
        {"unique": UNIQUE, "split": SPLIT},
        taxmap,
        by_accession,
        by_core,
    )

    profile = pd.read_csv(PROFILE, sep="\t", low_memory=False)
    calls = profile.loc[profile["calibrated_call"].map(as_bool)].copy()
    mapped_species = []
    mapping_methods = []
    for _, row in calls.iterrows():
        species, method = map_profile_row(row, best_raw_species, by_accession, by_core)
        mapped_species.append(species)
        mapping_methods.append(method)
    calls["gtdb_species"] = mapped_species
    calls["gtdb_mapping_method"] = mapping_methods
    calls["call_status"] = calls["gtdb_species"].map(lambda value: "TP" if value in truth_species else "FP")
    calls.loc[calls["gtdb_species"].eq(""), "call_status"] = "unmapped"

    for col in FEATURES:
        if col in calls.columns:
            calls[col] = pd.to_numeric(calls[col], errors="coerce").replace([math.inf, -math.inf], pd.NA)

    truth_only = truth_df.loc[~truth_df["gtdb_species"].isin(set(calls["gtdb_species"]))].copy()
    truth_only["call_status"] = "FN"

    selected_cols = [
        "taxid",
        "species_name",
        "gtdb_species",
        "call_status",
        "gtdb_mapping_method",
        "calibrated_probability",
        "calibrated_abundance",
        "reported_ani",
        "u_XnY_ctx_max",
        "s_XnY_ctx_max",
        "u_Ref_breadth_max",
        "s_Ref_breadth_max",
        "u_Real_min_align_fraction_max",
        "s_Real_min_align_fraction_max",
        "u_Ref_zip_af_max",
        "s_Ref_zip_af_max",
        "u_Ref_hit_mean_depth_max",
        "s_Ref_hit_mean_depth_max",
        "tail_rescue_added",
        "low_extra_split_rescue_added",
        "candidate_rescue_added",
        "raw_unique_fallback",
        "auto_exact_split_used",
        "s_best_accession",
        "u_best_accession",
    ]
    selected_cols = [col for col in selected_cols if col in calls.columns]
    calls[selected_cols].sort_values(
        ["call_status", "calibrated_abundance", "s_XnY_ctx_max"],
        ascending=[True, False, False],
        na_position="last",
    ).to_csv(RESULTS / "hmp_airskin8_minco_call_status.tsv", sep="\t", index=False)

    truth_only.sort_values("truth_abundance", ascending=False).to_csv(
        RESULTS / "hmp_airskin8_minco_missing_truth.tsv",
        sep="\t",
        index=False,
    )

    summary_rows = summarize_features(calls)
    summary_rows.append(
        {
            "call_status": "FN",
            "rows": int(len(truth_only)),
            "species": int(truth_only["gtdb_species"].nunique()),
            "truth_abundance_sum": float(truth_only["truth_abundance"].sum()) if len(truth_only) else 0.0,
        }
    )
    summary_rows.append(
        {
            "call_status": "mapping_quality",
            "rows": int(quality["source_genomes_positive"]),
            "species": int(quality["truth_gtdb_species"]),
            "truth_mass_mapped_pct": float(quality["truth_mass_mapped_pct"]),
            **raw_stats,
        }
    )
    pd.DataFrame(summary_rows).to_csv(
        RESULTS / "hmp_airskin8_minco_call_status_summary.tsv",
        sep="\t",
        index=False,
    )
    score_filtered_calls(SAMPLE, calls, truth_df).to_csv(
        RESULTS / "hmp_airskin8_minco_filter_diagnostic.tsv",
        sep="\t",
        index=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
