#!/usr/bin/env python3
"""Check raw MinCO tables for high-abundance HMP misses absent from profiles.

`audit_missed_truth_candidates.py` found that many HMP false negatives were not
visible in emitted profile rows. This follow-up checks the underlying raw
unique/split tables for those high-abundance missed GTDB species before
deciding whether the next algorithmic target is a profile rescue rule or
candidate-generation/indexing changes.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
MISSED = RESULTS / "missed_truth_candidate_species.tsv"
HIGH_TRUTH_ABUNDANCE = 0.01

RAW_NUMERIC_COLS = [
    "ANI",
    "Distance",
    "XnY_ctx",
    "Real_min_align_fraction",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_cv",
    "Ref_zero_fraction",
    "Read_match_fraction",
    "Block_match_fraction",
    "Relative_abundance_depth",
    "Normalized_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
]


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def raw_paths_from_profile(panel: str, sample: int, profile_path: Path) -> dict[str, Path]:
    parent = profile_path.parent
    sample_work = parent / f"sample{sample}_work"
    candidates = [
        {
            "unique": parent / "work/minco.best_diff_unique.unfiltered.tsv",
            "split": parent / "work/minco.best_diff_split.unfiltered.tsv",
            "split_exact": parent / "work/minco.best_diff_split.exact.unfiltered.tsv",
        },
        {
            "unique": sample_work / "minco.best_diff_unique.unfiltered.tsv",
            "split": sample_work / "minco.best_diff_split.unfiltered.tsv",
            "split_exact": sample_work / "minco.best_diff_split.exact.unfiltered.tsv",
        },
    ]
    for paths in candidates:
        if paths["unique"].exists() or paths["split"].exists() or paths["split_exact"].exists():
            return {key: value for key, value in paths.items() if value.exists()}

    if panel == "hmp_airskin_gtdb_source_abundance":
        hmp.ensure_sample_record(sample)
        rec = hmp.SAMPLES[sample]
    elif panel == "hmp_gastrooral_gtdb_source_abundance":
        rec = hmp_gastro.SAMPLES[sample]
    else:
        return {}
    return {
        "unique": Path(rec["unique"]),
        "split": Path(rec["split"]),
    }


class RawMapper:
    def __init__(self) -> None:
        self.by_accession, self.by_core = hmp.truth.load_gtdb_metadata(hmp.truth.GTDB_METADATA)
        self.hmp_taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
        self.hmp_gastro_taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)

    def load_raw(self, panel: str, raw_paths: Mapping[str, Path]) -> pd.DataFrame:
        frames = []
        taxmap = (
            self.hmp_gastro_taxmap
            if panel == "hmp_gastrooral_gtdb_source_abundance"
            else self.hmp_taxmap
        )
        for mode, path in raw_paths.items():
            if not path.exists():
                continue
            df = hmp.load_minco(path, taxmap, 11.0)
            if df.empty:
                continue
            df["raw_mode"] = mode
            df["raw_path"] = str(path)
            frames.append(df)
        if not frames:
            return pd.DataFrame()
        raw = pd.concat(frames, ignore_index=True, sort=False)
        species = []
        mapping_methods = []
        for row in raw.itertuples(index=False):
            accession = str(getattr(row, "accession", ""))
            if panel == "hmp_gastrooral_gtdb_source_abundance":
                gtdb_species, method = hmp_gastro.gtdb_from_accession(
                    accession,
                    self.by_accession,
                    self.by_core,
                )
            else:
                gtdb_species, method = hmp.gtdb_from_accession(
                    accession,
                    self.by_accession,
                    self.by_core,
                )
            species.append(gtdb_species)
            mapping_methods.append(method)
        raw["gtdb_species"] = species
        raw["gtdb_mapping_method"] = mapping_methods
        for col in RAW_NUMERIC_COLS:
            if col in raw.columns:
                raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
            else:
                raw[col] = 0.0
        return raw


def best_raw_rows(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return raw
    work = raw.loc[raw["gtdb_species"].astype(str).astype(bool)].copy()
    if work.empty:
        return work
    work["raw_mode_rank"] = work["raw_mode"].map({"split_exact": 2, "split": 1, "unique": 0}).fillna(0)
    work = work.sort_values(
        [
            "gtdb_species",
            "raw_mode_rank",
            "XnY_ctx",
            "Real_min_align_fraction",
            "ANI",
            "Ref_breadth",
            "Ref_mean_depth",
        ],
        ascending=[True, False, False, False, False, False, False],
        kind="mergesort",
    )
    return work.drop_duplicates("gtdb_species", keep="first")


def audit_rows() -> pd.DataFrame:
    missed = pd.read_csv(MISSED, sep="\t")
    work = missed.loc[
        missed["panel"].isin(
            {
                "hmp_airskin_gtdb_source_abundance",
                "hmp_gastrooral_gtdb_source_abundance",
            }
        )
        & missed["high_truth_abundance"].map(truthy)
        & ~missed["candidate_present"].map(truthy)
    ].copy()
    return work


def build_audit() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapper = RawMapper()
    targets = audit_rows()
    raw_cache: dict[tuple[str, int, str], pd.DataFrame] = {}
    rows: list[dict[str, object]] = []
    for target in targets.itertuples(index=False):
        panel = str(target.panel)
        sample = int(target.sample)
        profile_path = Path(str(target.profile_path))
        key = (panel, sample, str(profile_path))
        if key not in raw_cache:
            paths = raw_paths_from_profile(panel, sample, profile_path)
            raw_cache[key] = best_raw_rows(mapper.load_raw(panel, paths))
        best = raw_cache[key]
        species = str(target.gtdb_species)
        hit = best.loc[best["gtdb_species"].astype(str).eq(species)] if not best.empty else pd.DataFrame()
        present = not hit.empty
        out = {
            "panel": panel,
            "sample": sample,
            "gtdb_species": species,
            "truth_abundance_pct": finite_float(target.truth_abundance_pct),
            "emitted_profile_candidate_present": False,
            "raw_candidate_present": present,
            "profile_path": str(profile_path),
        }
        if present:
            rec = hit.iloc[0]
            out.update(
                {
                    "raw_mode": rec.get("raw_mode", ""),
                    "raw_path": rec.get("raw_path", ""),
                    "raw_ref": rec.get("Ref", ""),
                    "raw_accession": rec.get("accession", ""),
                    "raw_taxid": rec.get("taxid", ""),
                    "raw_species_name": rec.get("species_name", ""),
                    "gtdb_mapping_method": rec.get("gtdb_mapping_method", ""),
                    "ANI": finite_float(rec.get("ANI", 0.0)),
                    "XnY_ctx": finite_float(rec.get("XnY_ctx", 0.0)),
                    "Real_min_align_fraction": finite_float(
                        rec.get("Real_min_align_fraction", 0.0)
                    ),
                    "Ref_breadth": finite_float(rec.get("Ref_breadth", 0.0)),
                    "Ref_mean_depth": finite_float(rec.get("Ref_mean_depth", 0.0)),
                    "Ref_zip_af": finite_float(rec.get("Ref_zip_af", 0.0)),
                    "Ref_zip_aaf_ani": finite_float(rec.get("Ref_zip_aaf_ani", 0.0)),
                    "Default_call": rec.get("Default_call", ""),
                    "Default_call_rule": rec.get("Default_call_rule", ""),
                }
            )
        else:
            out.update(
                {
                    "raw_mode": "",
                    "raw_path": "",
                    "raw_ref": "",
                    "raw_accession": "",
                    "raw_taxid": "",
                    "raw_species_name": "",
                    "gtdb_mapping_method": "",
                    "ANI": 0.0,
                    "XnY_ctx": 0.0,
                    "Real_min_align_fraction": 0.0,
                    "Ref_breadth": 0.0,
                    "Ref_mean_depth": 0.0,
                    "Ref_zip_af": 0.0,
                    "Ref_zip_aaf_ani": 0.0,
                    "Default_call": "",
                    "Default_call_rule": "",
                }
            )
        rows.append(out)
    detail = pd.DataFrame(rows)
    panel_rows = []
    for panel, sub in detail.groupby("panel", sort=True):
        visible = sub.loc[sub["raw_candidate_present"]]
        absent = sub.loc[~sub["raw_candidate_present"]]
        panel_rows.append(
            {
                "panel": panel,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "high_truth_emitted_absent": int(len(sub)),
                "raw_candidate_visible": int(len(visible)),
                "raw_candidate_absent": int(len(absent)),
                "raw_candidate_visible_truth_mass_pct": visible["truth_abundance_pct"].sum(),
                "raw_candidate_absent_truth_mass_pct": absent["truth_abundance_pct"].sum(),
                "mean_visible_XnY_ctx": visible["XnY_ctx"].mean() if not visible.empty else 0.0,
                "mean_visible_Ref_breadth": visible["Ref_breadth"].mean()
                if not visible.empty
                else 0.0,
                "mean_visible_ANI": visible["ANI"].mean() if not visible.empty else 0.0,
            }
        )
    panel_summary = pd.DataFrame(panel_rows)
    total = int(len(detail))
    visible_n = int(detail["raw_candidate_present"].sum()) if not detail.empty else 0
    absent_n = total - visible_n
    audit = pd.DataFrame(
        [
            {
                "metric": "high_truth_emitted_absent_targets",
                "value": total,
                "evidence": "HMP high-abundance false negatives absent from emitted profile rows",
                "decision": "diagnostic_scope",
            },
            {
                "metric": "raw_candidate_visibility",
                "value": f"{visible_n}/{total}",
                "evidence": f"raw_absent={absent_n}; raw_visible={visible_n}",
                "decision": (
                    "raw_tables_expose_many_profile_absent_misses"
                    if visible_n > absent_n
                    else "raw_tables_do_not_explain_most_profile_absent_misses"
                ),
            },
            {
                "metric": "promotion_decision",
                "value": "diagnostic_only_not_default",
                "evidence": "truth-aware HMP raw-table visibility audit; no output-only rule selected",
                "decision": "use_to_decide_between_profile_rescue_and_indexing_work",
            },
        ]
    )
    return detail, panel_summary, audit


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    detail, panel_summary, audit = build_audit()
    detail.to_csv(RESULTS / "hmp_missed_truth_raw_table_detail.tsv", sep="\t", index=False)
    panel_summary.to_csv(
        RESULTS / "hmp_missed_truth_raw_table_panel_summary.tsv",
        sep="\t",
        index=False,
    )
    audit.to_csv(RESULTS / "hmp_missed_truth_raw_table_audit.tsv", sep="\t", index=False)
    print(panel_summary.to_string(index=False))
    print("\nAUDIT")
    print(audit.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
