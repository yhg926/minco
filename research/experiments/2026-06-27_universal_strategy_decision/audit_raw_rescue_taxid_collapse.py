#!/usr/bin/env python3
"""Classify why raw-cache rescue candidates are absent from wrapper output.

The emitted-profile candidate rescue can only add rows already present in the
taxid-level joined profile. HMP raw-cache rescue found additional candidates by
accession-to-GTDB-species mapping. This audit checks whether those raw-side
candidates are lost because multiple GTDB species collapse to one NCBI taxid
row, because the taxid row exists but is uncalled/weak, or because the taxid is
absent from the emitted profile surface.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp
import sweep_hmp_raw_candidate_rescue as raw_rescue


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

ANI_MIN = 0.90
XNY_MIN = 100.0
BREADTH_MIN = 0.01
REAL_AF_MIN = 0.70
RULE_NAME = "raw_side_cross_rule_ani90_xny100_br01_af70"

RAW_COLS = {
    "gtdb_species",
    "gtdb_mapping_method",
    "accession",
    "Ref",
    "Ref_annotation",
    "raw_mode",
    "ANI",
    "XnY_ctx",
    "Ref_breadth",
    "Real_min_align_fraction",
    "Ref_mean_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
}

PROFILE_COLS = {
    "taxid",
    "species_name",
    "calibrated_call",
    "calibrated_probability",
    "calibrated_abundance",
    "reported_ani",
    "reported_ani_source",
    "u_XnY_ctx_max",
    "s_XnY_ctx_max",
    "u_Real_min_align_fraction_max",
    "s_Real_min_align_fraction_max",
    "u_Ref_breadth_max",
    "s_Ref_breadth_max",
    "u_Ref_zip_aaf_ani_max",
    "s_Ref_zip_aaf_ani_max",
    "u_best_accession",
    "s_best_accession",
}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def accession_species(panel: str, accession: object, mapper: raw_audit.RawMapper) -> str:
    return raw_rescue.gtdb_from_accession(panel, accession, mapper)


def metadata_for_accession(accession: object, mapper: raw_audit.RawMapper) -> Mapping[str, str]:
    acc = str(accession or "")
    return mapper.by_accession.get(acc, {})


def load_profile(path: Path, panel: str, mapper: raw_audit.RawMapper) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in PROFILE_COLS, low_memory=False)
    for col in PROFILE_COLS:
        if col not in df.columns:
            df[col] = ""
    df["taxid"] = df["taxid"].astype(str)
    df["calibrated_call_bool"] = df["calibrated_call"].map(truthy)
    for prefix in ["s", "u"]:
        df[f"{prefix}_best_gtdb_species"] = df[f"{prefix}_best_accession"].map(
            lambda acc: accession_species(panel, acc, mapper)
        )
    return df


def load_raw_candidates(panel: str, sample: int, called_species: set[str]) -> pd.DataFrame:
    cache = raw_rescue.cache_path(panel, sample)
    if not cache.exists():
        raise SystemExit(f"missing raw candidate cache {cache}; run sweep_hmp_raw_candidate_rescue.py")
    raw = pd.read_csv(cache, sep="\t", usecols=lambda col: col in RAW_COLS, low_memory=False)
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
    return raw.loc[mask].copy()


def collapse_reason(profile_hit: pd.DataFrame, raw_species: str, raw_accession: str) -> str:
    if profile_hit.empty:
        return "profile_taxid_absent"
    called = bool(profile_hit["calibrated_call_bool"].any())
    best_accessions = set(profile_hit["s_best_accession"].astype(str)) | set(
        profile_hit["u_best_accession"].astype(str)
    )
    best_species = set(profile_hit["s_best_gtdb_species"].astype(str)) | set(
        profile_hit["u_best_gtdb_species"].astype(str)
    )
    best_species.discard("")
    if raw_accession in best_accessions and raw_species in best_species:
        return "same_accession_profile_row_uncalled" if not called else "same_accession_called_mapping_gap"
    if called and best_species and raw_species not in best_species:
        return "taxid_collapsed_to_other_called_gtdb"
    if called:
        return "taxid_called_but_unmapped_to_raw_gtdb"
    if best_species and raw_species not in best_species:
        return "taxid_uncalled_best_other_gtdb"
    return "taxid_profile_row_uncalled"


def build_audit() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapper = raw_audit.RawMapper()
    detail_rows: list[dict[str, object]] = []
    for panel in raw_rescue.PANELS:
        cfg = decomp.PANELS[panel]
        for sample_raw in cfg["samples"]:
            sample = int(sample_raw)
            truth = decomp.load_truth(
                Path(cfg["truth"](sample)),
                sample,
                str(cfg["truth_abundance_col"]),
            )
            truth_norm = raw_rescue.normalize(truth)
            profile_path = Path(cfg["minco"](sample))
            current_pred = raw_rescue.load_current_pred_gtdb(
                panel,
                profile_path,
                str(cfg["minco_collapse"]),
                mapper,
            )
            called_species = set(current_pred)
            profile = load_profile(profile_path, panel, mapper)
            raw_candidates = load_raw_candidates(panel, sample, called_species)
            for row in raw_candidates.itertuples(index=False):
                raw_species = str(getattr(row, "gtdb_species", ""))
                raw_accession = str(getattr(row, "accession", ""))
                meta = metadata_for_accession(raw_accession, mapper)
                raw_taxid = str(
                    meta.get("ncbi_species_taxid")
                    or meta.get("ncbi_taxid")
                    or ""
                )
                profile_hit = (
                    profile.loc[profile["taxid"].astype(str).eq(raw_taxid)].copy()
                    if raw_taxid
                    else pd.DataFrame()
                )
                reason = collapse_reason(profile_hit, raw_species, raw_accession)
                first = profile_hit.iloc[0] if not profile_hit.empty else {}
                detail_rows.append(
                    {
                        "panel": panel,
                        "sample": sample,
                        "rule": RULE_NAME,
                        "gtdb_species": raw_species,
                        "truth_status": "TP" if raw_species in truth_norm else "FP",
                        "truth_abundance_pct": truth_norm.get(raw_species, 0.0) * 100.0,
                        "raw_accession": raw_accession,
                        "raw_ncbi_species_taxid": raw_taxid,
                        "raw_ncbi_species": meta.get("ncbi_species", ""),
                        "raw_gtdb_representative": meta.get("gtdb_representative", ""),
                        "raw_ANI": finite(getattr(row, "ANI", 0.0)),
                        "raw_XnY_ctx": finite(getattr(row, "XnY_ctx", 0.0)),
                        "raw_Ref_breadth": finite(getattr(row, "Ref_breadth", 0.0)),
                        "raw_Real_min_align_fraction": finite(
                            getattr(row, "Real_min_align_fraction", 0.0)
                        ),
                        "raw_Ref_mean_depth": finite(getattr(row, "Ref_mean_depth", 0.0)),
                        "profile_taxid_present": not profile_hit.empty,
                        "profile_taxid_called": bool(profile_hit["calibrated_call_bool"].any())
                        if not profile_hit.empty
                        else False,
                        "profile_taxid_species_name": first.get("species_name", "")
                        if not profile_hit.empty
                        else "",
                        "profile_calibrated_probability": finite(
                            first.get("calibrated_probability", 0.0)
                        )
                        if not profile_hit.empty
                        else 0.0,
                        "profile_s_best_accession": first.get("s_best_accession", "")
                        if not profile_hit.empty
                        else "",
                        "profile_u_best_accession": first.get("u_best_accession", "")
                        if not profile_hit.empty
                        else "",
                        "profile_s_best_gtdb_species": first.get("s_best_gtdb_species", "")
                        if not profile_hit.empty
                        else "",
                        "profile_u_best_gtdb_species": first.get("u_best_gtdb_species", "")
                        if not profile_hit.empty
                        else "",
                        "profile_s_XnY_ctx_max": finite(first.get("s_XnY_ctx_max", 0.0))
                        if not profile_hit.empty
                        else 0.0,
                        "profile_u_XnY_ctx_max": finite(first.get("u_XnY_ctx_max", 0.0))
                        if not profile_hit.empty
                        else 0.0,
                        "profile_s_Ref_breadth_max": finite(
                            first.get("s_Ref_breadth_max", 0.0)
                        )
                        if not profile_hit.empty
                        else 0.0,
                        "profile_u_Ref_breadth_max": finite(
                            first.get("u_Ref_breadth_max", 0.0)
                        )
                        if not profile_hit.empty
                        else 0.0,
                        "collapse_reason": reason,
                        "profile_path": str(profile_path),
                    }
                )
    detail = pd.DataFrame(detail_rows)
    if detail.empty:
        summary = pd.DataFrame()
        audit = pd.DataFrame(
            [
                {
                    "metric": "raw_rescue_candidates",
                    "value": 0,
                    "evidence": RULE_NAME,
                    "decision": "no_candidates",
                }
            ]
        )
        return detail, summary, audit

    summary = (
        detail.groupby(["collapse_reason", "truth_status"], as_index=False)
        .agg(
            candidates=("gtdb_species", "count"),
            samples=("sample", lambda values: ",".join(map(str, sorted(set(values))))),
            truth_abundance_pct=("truth_abundance_pct", "sum"),
        )
        .sort_values(["collapse_reason", "truth_status"])
    )
    tp = int(detail["truth_status"].eq("TP").sum())
    fp = int(detail["truth_status"].eq("FP").sum())
    reason_counts = detail["collapse_reason"].value_counts().to_dict()
    audit = pd.DataFrame(
        [
            {
                "metric": "raw_rescue_candidates",
                "value": int(len(detail)),
                "evidence": f"TP={tp};FP={fp};rule={RULE_NAME}",
                "decision": "diagnostic_scope",
            },
            {
                "metric": "taxid_collapsed_to_other_called_gtdb",
                "value": int(reason_counts.get("taxid_collapsed_to_other_called_gtdb", 0)),
                "evidence": "same NCBI taxid row is called but best accession maps to another GTDB species",
                "decision": "requires_accession_or_gtdb_species_level_candidate_surface",
            },
            {
                "metric": "profile_taxid_absent",
                "value": int(reason_counts.get("profile_taxid_absent", 0)),
                "evidence": "raw accession maps to a GTDB species whose NCBI taxid is absent from emitted profile rows",
                "decision": "requires_raw_side_candidate_surface",
            },
            {
                "metric": "taxid_profile_row_uncalled_or_other",
                "value": int(
                    len(detail)
                    - reason_counts.get("taxid_collapsed_to_other_called_gtdb", 0)
                    - reason_counts.get("profile_taxid_absent", 0)
                ),
                "evidence": "profile taxid exists but is uncalled or has nonmatching best accession",
                "decision": "could_be_handled_by_profile_or_accession_rescue",
            },
            {
                "metric": "promotion_decision",
                "value": "diagnostic_only_not_default",
                "evidence": "classification audit only; no default behavior changed",
                "decision": "use_to_design_raw_side_candidate_surface",
            },
        ]
    )
    return detail, summary, audit


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    detail, summary, audit = build_audit()
    detail.to_csv(RESULTS / "raw_rescue_taxid_collapse_detail.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "raw_rescue_taxid_collapse_summary.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "raw_rescue_taxid_collapse_audit.tsv", sep="\t", index=False)
    print(audit.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
