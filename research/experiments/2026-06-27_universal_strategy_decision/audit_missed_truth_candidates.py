#!/usr/bin/env python3
"""Audit whether missed truth species are visible as uncalled MinCO candidates.

The fixed-call abundance oracle showed that abundance allocation alone is not
enough on all panels. This diagnostic asks the next practical question: for the
truth species missed by the current MinCO calls, did the calibrated profile
contain any mapped candidate row that could potentially be rescued by a safer
gate, or was the species absent from the candidate surface altogether?
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import decompose_abundance_errors as decomp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
HIGH_TRUTH_ABUNDANCE = 0.01

PROFILE_COLS = {
    "taxid",
    "species_name",
    "calibrated_call",
    "calibrated_probability",
    "calibrated_threshold",
    "calibrated_abundance",
    "calibrated_abundance_raw",
    "reported_ani",
    "reported_ani_source",
    "u_XnY_ctx_max",
    "u_ANI_max",
    "u_Real_min_align_fraction_max",
    "u_Ref_breadth_max",
    "u_Ref_mean_depth_max",
    "u_Normalized_abundance_depth_max",
    "u_Ref_zip_af_max",
    "u_Ref_zip_aaf_ani_max",
    "u_best_accession",
    "s_XnY_ctx_max",
    "s_ANI_max",
    "s_Real_min_align_fraction_max",
    "s_Ref_breadth_max",
    "s_Ref_mean_depth_max",
    "s_Normalized_abundance_depth_max",
    "s_Ref_zip_af_max",
    "s_Ref_zip_aaf_ani_max",
    "s_best_accession",
    "tail_rescue_added",
    "low_extra_split_rescue_added",
    "adaptive_mode",
}


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize_name(value: object) -> str:
    text = str(value or "").strip()
    if text.startswith("s__"):
        return text
    return text


def call_mask(df: pd.DataFrame) -> pd.Series:
    return df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})


def read_profile(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=lambda col: col in PROFILE_COLS, low_memory=False)


class CandidateMapper:
    def __init__(self) -> None:
        toy_mod = decomp.load_toy_scorer()
        self.by_accession, self.by_core = toy_mod.score.truth.load_gtdb_metadata(
            toy_mod.score.truth.GTDB_METADATA
        )
        self.hmp_taxmap = decomp.hmp.parse_species_taxmap(decomp.hmp.TAXMAP)
        self.hmp_gastro_taxmap = decomp.hmp_gastro.parse_species_taxmap(decomp.hmp_gastro.TAXMAP)
        (
            _cami3_wgs_to_species,
            self.cami3_taxid_to_species,
            self.cami3_name_to_species,
            _map_diag,
        ) = decomp.cami3.build_transfer_maps()
        self.cami3_taxmap = decomp.cami3.parse_species_taxmap(decomp.cami3.TAXMAP)

    def hmp_best_ref_species(self, sample: int) -> dict[str, str]:
        decomp.hmp.ensure_sample_record(sample)
        hmp_paths = decomp.hmp.SAMPLES[sample]
        best, _diag = decomp.hmp.best_raw_ref_species_by_taxid(
            {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
            self.hmp_taxmap,
            self.by_accession,
            self.by_core,
        )
        return best

    def hmp_gastro_best_ref_species(self, sample: int) -> dict[str, str]:
        hmp_paths = decomp.hmp_gastro.SAMPLES[sample]
        best, _diag = decomp.hmp_gastro.best_raw_ref_species_by_taxid(
            {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
            self.hmp_gastro_taxmap,
            self.by_accession,
            self.by_core,
        )
        return best

    def cami3_best_ref_species(self, sample: int) -> dict[str, str]:
        best, _diag = decomp.cami3.best_raw_ref_species_by_taxid(
            decomp.cami3.RAW_TABLES[sample],
            self.cami3_taxmap,
            self.by_accession,
            self.by_core,
        )
        return best

    def accession_to_species(self, panel: str, accession: object) -> str:
        accession_text = str(accession or "")
        if not accession_text:
            return ""
        if panel == "cami3_toy_human_gut_gtdb_source_readmap":
            return decomp.cami3.gtdb_from_accession(
                accession_text,
                self.by_accession,
                self.by_core,
            )
        if panel == "hmp_gastrooral_gtdb_source_abundance":
            species, _method = decomp.hmp_gastro.gtdb_from_accession(
                accession_text,
                self.by_accession,
                self.by_core,
            )
            return species
        if panel == "hmp_airskin_gtdb_source_abundance":
            species, _method = decomp.hmp.gtdb_from_accession(
                accession_text,
                self.by_accession,
                self.by_core,
            )
            return species
        return ""

    def map_profile(self, panel: str, sample: int, path: Path) -> pd.DataFrame:
        raw = read_profile(path)
        raw["called"] = call_mask(raw)
        raw["mapped_gtdb_species"] = ""
        raw["mapping_method"] = ""

        if panel == "cami2_toy_mouse_gut":
            raw["mapped_gtdb_species"] = raw["species_name"].map(normalize_name)
            raw["mapping_method"] = "species_name"
            return raw

        best_ref_species: Mapping[str, str] = {}
        if panel == "hmp_airskin_gtdb_source_abundance":
            best_ref_species = self.hmp_best_ref_species(sample)
        elif panel == "hmp_gastrooral_gtdb_source_abundance":
            best_ref_species = self.hmp_gastro_best_ref_species(sample)
        elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
            best_ref_species = self.cami3_best_ref_species(sample)

        mapped_species: list[str] = []
        mapping_method: list[str] = []
        for row in raw.itertuples(index=False):
            species = ""
            method = ""
            for col in ["s_best_accession", "u_best_accession"]:
                if hasattr(row, col):
                    species = self.accession_to_species(panel, getattr(row, col))
                    if species:
                        method = col
                        break
            taxid = str(getattr(row, "taxid", ""))
            if not species and taxid in best_ref_species:
                species = best_ref_species[taxid]
                method = "best_raw_ref_taxid"
            if not species and panel == "cami3_toy_human_gut_gtdb_source_readmap":
                species = self.cami3_taxid_to_species.get(taxid, "")
                if species:
                    method = "taxid"
            if not species and panel == "cami3_toy_human_gut_gtdb_source_readmap":
                species_name = normalize_name(getattr(row, "species_name", ""))
                species = self.cami3_name_to_species.get(species_name, "")
                if species:
                    method = "species_name"
            mapped_species.append(species)
            mapping_method.append(method)
        raw["mapped_gtdb_species"] = mapped_species
        raw["mapping_method"] = mapping_method
        return raw


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([0.0] * len(df), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([math.inf, -math.inf], pd.NA).fillna(0.0)


def candidate_best_rows(candidates: pd.DataFrame) -> pd.DataFrame:
    work = candidates.loc[candidates["mapped_gtdb_species"].astype(str).astype(bool)].copy()
    if work.empty:
        return pd.DataFrame()
    for col in [
        "calibrated_probability",
        "calibrated_abundance",
        "calibrated_abundance_raw",
        "reported_ani",
        "s_XnY_ctx_max",
        "u_XnY_ctx_max",
        "s_Ref_breadth_max",
        "u_Ref_breadth_max",
        "s_Real_min_align_fraction_max",
        "u_Real_min_align_fraction_max",
        "s_Ref_mean_depth_max",
        "u_Ref_mean_depth_max",
        "s_Ref_zip_af_max",
        "u_Ref_zip_af_max",
        "s_Ref_zip_aaf_ani_max",
        "u_Ref_zip_aaf_ani_max",
    ]:
        work[col] = numeric(work, col)
    work["best_XnY"] = work[["s_XnY_ctx_max", "u_XnY_ctx_max"]].max(axis=1)
    work["best_breadth"] = work[["s_Ref_breadth_max", "u_Ref_breadth_max"]].max(axis=1)
    work["best_real_min_af"] = work[
        ["s_Real_min_align_fraction_max", "u_Real_min_align_fraction_max"]
    ].max(axis=1)
    work["best_depth"] = work[["s_Ref_mean_depth_max", "u_Ref_mean_depth_max"]].max(axis=1)
    work["best_zip_af"] = work[["s_Ref_zip_af_max", "u_Ref_zip_af_max"]].max(axis=1)
    work["best_zip_aaf_ani"] = work[
        ["s_Ref_zip_aaf_ani_max", "u_Ref_zip_aaf_ani_max", "reported_ani"]
    ].max(axis=1)
    work = work.sort_values(
        [
            "mapped_gtdb_species",
            "called",
            "calibrated_probability",
            "best_XnY",
            "best_real_min_af",
            "best_breadth",
            "best_depth",
        ],
        ascending=[True, False, False, False, False, False, False],
        kind="mergesort",
    )
    return work.drop_duplicates("mapped_gtdb_species", keep="first")


def called_species(candidates: pd.DataFrame) -> set[str]:
    return set(
        candidates.loc[candidates["called"] & candidates["mapped_gtdb_species"].astype(str).astype(bool),
                       "mapped_gtdb_species"].astype(str)
    )


def audit_sample(
    panel: str,
    sample: int,
    truth: Mapping[str, float],
    candidates: pd.DataFrame,
    profile_path: Path,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    truth_norm = decomp.normalize(truth)
    calls = called_species(candidates)
    best = candidate_best_rows(candidates)
    best_by_species = (
        {str(row.mapped_gtdb_species): row for row in best.itertuples(index=False)}
        if not best.empty
        else {}
    )
    fn_species = sorted(set(truth_norm) - calls, key=lambda species: truth_norm[species], reverse=True)
    rows: list[dict[str, object]] = []
    for species in fn_species:
        rec = best_by_species.get(species)
        candidate_present = rec is not None
        row = {
            "panel": panel,
            "sample": sample,
            "gtdb_species": species,
            "truth_abundance": truth_norm[species],
            "truth_abundance_pct": truth_norm[species] * 100.0,
            "high_truth_abundance": truth_norm[species] >= HIGH_TRUTH_ABUNDANCE,
            "candidate_present": candidate_present,
            "profile_path": str(profile_path),
        }
        if rec is not None:
            row.update(
                {
                    "candidate_called": bool(rec.called),
                    "candidate_probability": finite_float(getattr(rec, "calibrated_probability", 0.0)),
                    "candidate_threshold": finite_float(getattr(rec, "calibrated_threshold", 0.0)),
                    "candidate_abundance": finite_float(getattr(rec, "calibrated_abundance", 0.0)),
                    "candidate_abundance_raw": finite_float(getattr(rec, "calibrated_abundance_raw", 0.0)),
                    "candidate_best_XnY": finite_float(getattr(rec, "best_XnY", 0.0)),
                    "candidate_best_breadth": finite_float(getattr(rec, "best_breadth", 0.0)),
                    "candidate_best_real_min_af": finite_float(getattr(rec, "best_real_min_af", 0.0)),
                    "candidate_best_depth": finite_float(getattr(rec, "best_depth", 0.0)),
                    "candidate_best_zip_af": finite_float(getattr(rec, "best_zip_af", 0.0)),
                    "candidate_best_zip_aaf_ani": finite_float(getattr(rec, "best_zip_aaf_ani", 0.0)),
                    "mapping_method": str(getattr(rec, "mapping_method", "")),
                    "taxid": str(getattr(rec, "taxid", "")),
                    "species_name": str(getattr(rec, "species_name", "")),
                }
            )
        else:
            row.update(
                {
                    "candidate_called": False,
                    "candidate_probability": 0.0,
                    "candidate_threshold": 0.0,
                    "candidate_abundance": 0.0,
                    "candidate_abundance_raw": 0.0,
                    "candidate_best_XnY": 0.0,
                    "candidate_best_breadth": 0.0,
                    "candidate_best_real_min_af": 0.0,
                    "candidate_best_depth": 0.0,
                    "candidate_best_zip_af": 0.0,
                    "candidate_best_zip_aaf_ani": 0.0,
                    "mapping_method": "",
                    "taxid": "",
                    "species_name": "",
                }
            )
        rows.append(row)

    fn_mass = sum(truth_norm[species] for species in fn_species)
    visible_mass = sum(truth_norm[row["gtdb_species"]] for row in rows if row["candidate_present"])
    high_rows = [row for row in rows if row["high_truth_abundance"]]
    uncalled_candidate_rows = int(len(candidates) - candidates["called"].sum())
    summary = {
        "panel": panel,
        "sample": sample,
        "truth_species": len(truth_norm),
        "called_truth_species": len(set(truth_norm) & calls),
        "FN": len(fn_species),
        "FN_truth_mass_pct": fn_mass * 100.0,
        "FN_candidate_visible": sum(1 for row in rows if row["candidate_present"]),
        "FN_candidate_absent": sum(1 for row in rows if not row["candidate_present"]),
        "FN_candidate_visible_truth_mass_pct": visible_mass * 100.0,
        "high_truth_FN": len(high_rows),
        "high_truth_FN_mass_pct": sum(row["truth_abundance"] for row in high_rows) * 100.0,
        "high_truth_FN_candidate_visible": sum(1 for row in high_rows if row["candidate_present"]),
        "high_truth_FN_candidate_absent": sum(1 for row in high_rows if not row["candidate_present"]),
        "high_truth_threshold_pct": HIGH_TRUTH_ABUNDANCE * 100.0,
        "candidate_rows": int(len(candidates)),
        "candidate_rows_mapped": int(candidates["mapped_gtdb_species"].astype(str).astype(bool).sum()),
        "called_rows": int(candidates["called"].sum()),
        "uncalled_candidate_rows": uncalled_candidate_rows,
        "full_emitted_candidate_surface": uncalled_candidate_rows > 0,
        "profile_path": str(profile_path),
    }
    return rows, summary


def summarize_panel(sample_summaries: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, sub in sample_summaries.groupby("panel", sort=True):
        rows.append(
            {
                "panel": panel,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "total_FN": int(sub["FN"].sum()),
                "mean_FN_truth_mass_pct": sub["FN_truth_mass_pct"].mean(),
                "mean_FN_candidate_visible_truth_mass_pct": sub[
                    "FN_candidate_visible_truth_mass_pct"
                ].mean(),
                "total_high_truth_FN": int(sub["high_truth_FN"].sum()),
                "total_high_truth_FN_candidate_visible": int(
                    sub["high_truth_FN_candidate_visible"].sum()
                ),
                "total_high_truth_FN_candidate_absent": int(
                    sub["high_truth_FN_candidate_absent"].sum()
                ),
                "mean_high_truth_FN_mass_pct": sub["high_truth_FN_mass_pct"].mean(),
                "mean_candidate_rows_mapped": sub["candidate_rows_mapped"].mean(),
                "samples_with_uncalled_candidate_rows": int(
                    sub["full_emitted_candidate_surface"].astype(bool).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def build_audit(panel_summary: pd.DataFrame) -> pd.DataFrame:
    total_high = int(panel_summary["total_high_truth_FN"].sum()) if not panel_summary.empty else 0
    visible_high = (
        int(panel_summary["total_high_truth_FN_candidate_visible"].sum())
        if not panel_summary.empty
        else 0
    )
    absent_high = (
        int(panel_summary["total_high_truth_FN_candidate_absent"].sum())
        if not panel_summary.empty
        else 0
    )
    rows = [
        {
            "metric": "evaluated_panels",
            "value": int(len(panel_summary)),
            "evidence": ",".join(panel_summary["panel"].astype(str)) if not panel_summary.empty else "",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "high_truth_fn_threshold_pct",
            "value": HIGH_TRUTH_ABUNDANCE * 100.0,
            "evidence": "truth species with abundance >= this threshold are high-abundance missed calls",
            "decision": "diagnostic_parameter",
        },
        {
            "metric": "high_truth_fn_candidate_visibility",
            "value": f"{visible_high}/{total_high}",
            "evidence": f"absent={absent_high}; visible={visible_high}",
            "decision": (
                "rescue_gate_possible_for_many_high_abundance_misses"
                if visible_high > absent_high
                else "many_high_abundance_misses_not_visible_in_emitted_profile_surface"
            ),
        },
        {
            "metric": "samples_with_uncalled_candidate_rows",
            "value": int(panel_summary["samples_with_uncalled_candidate_rows"].sum())
            if not panel_summary.empty
            else 0,
            "evidence": "profiles with candidate_rows > called_rows; emitted-only profiles need raw-table follow-up",
            "decision": "candidate_visibility_is_profile_surface_limited",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "truth-aware missed-call audit over emitted profile rows; does not define an output-only rescue rule",
            "decision": "use_to_design_next_call_recovery_strategy",
        },
    ]
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mapper = CandidateMapper()
    fn_rows: list[dict[str, object]] = []
    sample_summaries: list[dict[str, object]] = []
    for panel, cfg in decomp.PANELS.items():
        for sample in cfg["samples"]:
            sample_id = int(sample)
            truth_path = Path(cfg["truth"](sample_id))
            profile_path = Path(cfg["minco"](sample_id))
            truth = decomp.load_truth(truth_path, sample_id, str(cfg["truth_abundance_col"]))
            candidates = mapper.map_profile(panel, sample_id, profile_path)
            rows, summary = audit_sample(panel, sample_id, truth, candidates, profile_path)
            fn_rows.extend(rows)
            sample_summaries.append(summary)

    fn_df = pd.DataFrame(fn_rows)
    sample_df = pd.DataFrame(sample_summaries)
    panel_df = summarize_panel(sample_df)
    audit_df = build_audit(panel_df)
    top_df = (
        fn_df.sort_values(["truth_abundance", "panel", "sample"], ascending=[False, True, True])
        .head(100)
        .copy()
        if not fn_df.empty
        else pd.DataFrame()
    )
    fn_df.to_csv(RESULTS / "missed_truth_candidate_species.tsv", sep="\t", index=False)
    top_df.to_csv(RESULTS / "missed_truth_candidate_top.tsv", sep="\t", index=False)
    sample_df.to_csv(RESULTS / "missed_truth_candidate_sample_summary.tsv", sep="\t", index=False)
    panel_df.to_csv(RESULTS / "missed_truth_candidate_panel_summary.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / "missed_truth_candidate_audit.tsv", sep="\t", index=False)
    print(panel_df.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
