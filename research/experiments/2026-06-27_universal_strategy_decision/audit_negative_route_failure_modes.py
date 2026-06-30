#!/usr/bin/env python3
"""Audit why completed negative routes miss truth species.

This is a cached-only diagnostic. It compares selected-default calls with the
truth rows for the completed marine and CAMI3 extension routes, then checks the
underlying best-diff raw tables for each missed GTDB species. The goal is to
separate "not visible in raw evidence" from "visible before the final gate",
which points to different next algorithm work.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd


EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
RESULTS = EXP / "results"

TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import audit_cami3_source_readmap_recovery_inputs as cami3_inputs  # noqa: E402


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
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
    "Effective_abundance_depth",
    "Normalized_effective_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
    "Reliable_Ref_breadth",
    "Reliable_Ref_mean_depth",
    "Reliable_Ref_hit_ctx",
    "Reliable_Ref_hit_mean_depth",
    "Reliable_Ref_hit_median_depth",
    "Reliable_Ref_zip_af",
]

PROFILE_NUMERIC_COLS = [
    "calibrated_probability",
    "calibrated_abundance",
    "reported_ani",
    "s_XnY_ctx_max",
    "u_XnY_ctx_max",
    "s_Ref_breadth_max",
    "u_Ref_breadth_max",
    "s_Real_min_align_fraction_max",
    "u_Real_min_align_fraction_max",
    "s_Ref_hit_mean_depth_max",
    "u_Ref_hit_mean_depth_max",
]


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def extract_accession(value: object) -> str:
    match = ACC_RE.search(str(value or ""))
    return match.group(1) if match else ""


class GtdbMapper:
    def __init__(self) -> None:
        self.by_accession, self.by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)

    def species_from_accession(self, accession: object) -> tuple[str, str]:
        rec, method, _key = truth.lookup_accession(
            extract_accession(accession),
            self.by_accession,
            self.by_core,
        )
        if rec is None:
            return "", method
        return str(rec.get("gtdb_species", "")), method


def profile_calls(path: Path, mapper: GtdbMapper) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call_col = raw.get("calibrated_call")
    if call_col is None:
        return pd.DataFrame()
    calls = raw.loc[call_col.map(truthy)].copy()
    rows: list[dict[str, object]] = []
    for row in calls.itertuples(index=False):
        species = ""
        method = ""
        for col in ["candidate_surface_accession", "s_best_accession", "u_best_accession"]:
            value = getattr(row, col, "")
            species, method = mapper.species_from_accession(value)
            if species:
                break
        if not species:
            continue
        out = {
            "gtdb_species": species,
            "mapping_method": method,
            "profile_accession": extract_accession(
                getattr(row, "candidate_surface_accession", "")
                or getattr(row, "s_best_accession", "")
                or getattr(row, "u_best_accession", "")
            ),
        }
        for col in PROFILE_NUMERIC_COLS:
            out[col] = finite(getattr(row, col, 0.0))
        rows.append(out)
    if not rows:
        return pd.DataFrame(columns=["gtdb_species"])
    df = pd.DataFrame(rows)
    df["profile_best_xny"] = df[["s_XnY_ctx_max", "u_XnY_ctx_max"]].max(axis=1)
    df["profile_best_breadth"] = df[["s_Ref_breadth_max", "u_Ref_breadth_max"]].max(axis=1)
    df["profile_best_real_af"] = df[
        ["s_Real_min_align_fraction_max", "u_Real_min_align_fraction_max"]
    ].max(axis=1)
    df["profile_best_hit_depth"] = df[
        ["s_Ref_hit_mean_depth_max", "u_Ref_hit_mean_depth_max"]
    ].max(axis=1)
    df = df.sort_values(
        [
            "gtdb_species",
            "calibrated_probability",
            "profile_best_xny",
            "profile_best_breadth",
        ],
        ascending=[True, False, False, False],
        kind="mergesort",
    )
    return df.drop_duplicates("gtdb_species", keep="first")


def read_raw_table(path: Path, mode: str, mapper: GtdbMapper) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    cols = ["Qry", "Ref", "Ref_annotation", "Default_call", "Default_call_rule", *RAW_NUMERIC_COLS]
    raw = pd.read_csv(path, sep="\t", usecols=lambda col: col in set(cols), low_memory=False)
    if raw.empty:
        return raw
    raw["raw_mode"] = mode
    raw["raw_path"] = str(path)
    species: list[str] = []
    methods: list[str] = []
    accessions: list[str] = []
    for row in raw.itertuples(index=False):
        accession = extract_accession(getattr(row, "Ref", "")) or extract_accession(
            getattr(row, "Ref_annotation", "")
        )
        gtdb_species, method = mapper.species_from_accession(accession)
        species.append(gtdb_species)
        methods.append(method)
        accessions.append(accession)
    raw["gtdb_species"] = species
    raw["gtdb_mapping_method"] = methods
    raw["raw_accession"] = accessions
    for col in RAW_NUMERIC_COLS:
        if col in raw.columns:
            raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
        else:
            raw[col] = 0.0
    return raw.loc[raw["gtdb_species"].astype(str).ne("")].copy()


def best_raw_rows(paths: Mapping[str, Path], mapper: GtdbMapper) -> pd.DataFrame:
    frames = [read_raw_table(path, mode, mapper) for mode, path in paths.items()]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame()
    raw = pd.concat(frames, ignore_index=True, sort=False)
    raw["raw_mode_rank"] = raw["raw_mode"].map({"split_exact": 2, "split": 1, "unique": 0}).fillna(0)
    raw = raw.sort_values(
        [
            "gtdb_species",
            "raw_mode_rank",
            "XnY_ctx",
            "Real_min_align_fraction",
            "ANI",
            "Ref_breadth",
            "Ref_hit_mean_depth",
        ],
        ascending=[True, False, False, False, False, False, False],
        kind="mergesort",
    )
    return raw.drop_duplicates("gtdb_species", keep="first")


def route_configs() -> list[dict[str, object]]:
    marine_root = Path("/mnt/new3T/minco_marine_selected_default_20260630")
    return [
        {
            "route": "marine_setup_truth_profile_rescore",
            "samples": [3, 4, 5],
            "truth": RESULTS / "marine_setup_truth_profile_rescore_truth.tsv",
            "profile": lambda sample: marine_root / f"minco_sample{sample}_current_default.tsv",
            "raw_paths": lambda sample: {
                "unique": marine_root / f"sample{sample}_work/minco.best_diff_unique.unfiltered.tsv",
                "split": marine_root / f"sample{sample}_work/minco.best_diff_split.unfiltered.tsv",
                "split_exact": marine_root
                / f"sample{sample}_work/minco.best_diff_split.exact.unfiltered.tsv",
            },
        },
        {
            "route": "cami3_source_readmap_extension_after_recovery",
            "samples": [3, 4, 5],
            "truth": RESULTS / "cami3_source_readmap_extension_after_recovery_truth.tsv",
            "profile": lambda sample: Path(
                f"/tmp/minco_candidate_preset_replay_20260629/"
                f"cami3_toy_human_gut_gtdb_source_readmap.sample{sample}.candidate_preset.tsv"
            ),
            "raw_paths": lambda sample: {
                "unique": cami3_inputs.raw_unique(sample),
                "split": cami3_inputs.raw_split(sample),
            },
        },
    ]


def build_route_details(config: Mapping[str, object], mapper: GtdbMapper) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    truth_df = pd.read_csv(Path(str(config["truth"])), sep="\t")
    route = str(config["route"])
    detail_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for sample in config["samples"]:  # type: ignore[index]
        sample = int(sample)
        sample_truth = truth_df.loc[truth_df["sample"].astype(int).eq(sample)].copy()
        sample_truth = sample_truth.loc[sample_truth["gtdb_species"].astype(str).ne("")]
        truth_abundance = {
            str(row.gtdb_species): finite(row.truth_abundance)
            for row in sample_truth.itertuples(index=False)
        }
        profile_path = config["profile"](sample)  # type: ignore[operator]
        calls = profile_calls(Path(profile_path), mapper)
        called_species = set(calls["gtdb_species"].astype(str)) if not calls.empty else set()
        raw = best_raw_rows(config["raw_paths"](sample), mapper)  # type: ignore[operator]
        raw_by_species = (
            {str(row.gtdb_species): row for row in raw.itertuples(index=False)}
            if not raw.empty
            else {}
        )
        truth_species = set(truth_abundance)
        fn_species = sorted(truth_species - called_species, key=lambda sp: truth_abundance[sp], reverse=True)
        fp_species = sorted(called_species - truth_species)
        tp_species = truth_species & called_species
        visible_mass = 0.0
        visible_high = 0
        high_n = 0
        default_major_n = 0
        for species in fn_species:
            rec = raw_by_species.get(species)
            present = rec is not None
            abundance = truth_abundance[species]
            high = abundance >= HIGH_TRUTH_ABUNDANCE
            if high:
                high_n += 1
            if present:
                visible_mass += abundance
                if high:
                    visible_high += 1
                if str(getattr(rec, "Default_call", "")).lower() == "major":
                    default_major_n += 1
            row = {
                "route": route,
                "sample": sample,
                "status": "FN",
                "gtdb_species": species,
                "truth_abundance": abundance,
                "truth_abundance_pct": abundance * 100.0,
                "high_truth_abundance": high,
                "raw_candidate_present": present,
                "profile_path": str(profile_path),
            }
            if present:
                row.update(
                    {
                        "raw_mode": getattr(rec, "raw_mode", ""),
                        "raw_accession": getattr(rec, "raw_accession", ""),
                        "raw_path": getattr(rec, "raw_path", ""),
                        "gtdb_mapping_method": getattr(rec, "gtdb_mapping_method", ""),
                        "ANI": finite(getattr(rec, "ANI", 0.0)),
                        "XnY_ctx": finite(getattr(rec, "XnY_ctx", 0.0)),
                        "Real_min_align_fraction": finite(
                            getattr(rec, "Real_min_align_fraction", 0.0)
                        ),
                        "Ref_breadth": finite(getattr(rec, "Ref_breadth", 0.0)),
                        "Ref_hit_mean_depth": finite(getattr(rec, "Ref_hit_mean_depth", 0.0)),
                        "Ref_zip_af": finite(getattr(rec, "Ref_zip_af", 0.0)),
                        "Ref_zip_aaf_ani": finite(getattr(rec, "Ref_zip_aaf_ani", 0.0)),
                        "Default_call": getattr(rec, "Default_call", ""),
                        "Default_call_rule": getattr(rec, "Default_call_rule", ""),
                    }
                )
            else:
                row.update(
                    {
                        "raw_mode": "",
                        "raw_accession": "",
                        "raw_path": "",
                        "gtdb_mapping_method": "",
                        "ANI": 0.0,
                        "XnY_ctx": 0.0,
                        "Real_min_align_fraction": 0.0,
                        "Ref_breadth": 0.0,
                        "Ref_hit_mean_depth": 0.0,
                        "Ref_zip_af": 0.0,
                        "Ref_zip_aaf_ani": 0.0,
                        "Default_call": "",
                        "Default_call_rule": "",
                    }
                )
            detail_rows.append(row)
        sample_rows.append(
            {
                "route": route,
                "sample": sample,
                "truth_species": len(truth_species),
                "called_species": len(called_species),
                "TP": len(tp_species),
                "FP": len(fp_species),
                "FN": len(fn_species),
                "FN_truth_mass_pct": 100.0 * sum(truth_abundance[sp] for sp in fn_species),
                "FN_raw_visible": sum(1 for sp in fn_species if sp in raw_by_species),
                "FN_raw_absent": sum(1 for sp in fn_species if sp not in raw_by_species),
                "FN_raw_visible_truth_mass_pct": 100.0 * visible_mass,
                "high_truth_FN": high_n,
                "high_truth_FN_raw_visible": visible_high,
                "high_truth_threshold_pct": HIGH_TRUTH_ABUNDANCE * 100.0,
                "FN_raw_default_major": default_major_n,
                "raw_species_available": int(len(raw_by_species)),
                "profile_path": str(profile_path),
            }
        )
    return detail_rows, sample_rows


def summarize(sample_rows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for route, sub in sample_rows.groupby("route", sort=True):
        rows.append(
            {
                "route": route,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "truth_species": int(sub["truth_species"].sum()),
                "called_species": int(sub["called_species"].sum()),
                "TP": int(sub["TP"].sum()),
                "FP": int(sub["FP"].sum()),
                "FN": int(sub["FN"].sum()),
                "mean_FN_truth_mass_pct": sub["FN_truth_mass_pct"].mean(),
                "FN_raw_visible": int(sub["FN_raw_visible"].sum()),
                "FN_raw_absent": int(sub["FN_raw_absent"].sum()),
                "mean_FN_raw_visible_truth_mass_pct": sub[
                    "FN_raw_visible_truth_mass_pct"
                ].mean(),
                "high_truth_FN": int(sub["high_truth_FN"].sum()),
                "high_truth_FN_raw_visible": int(sub["high_truth_FN_raw_visible"].sum()),
                "FN_raw_default_major": int(sub["FN_raw_default_major"].sum()),
                "mean_raw_species_available": sub["raw_species_available"].mean(),
                "diagnosis": "raw_visible_gate_or_surface_problem"
                if int(sub["FN_raw_visible"].sum()) > int(sub["FN_raw_absent"].sum())
                else "raw_absence_or_mapping_problem",
            }
        )
    return pd.DataFrame(rows)


def hmp_context_rows() -> list[dict[str, object]]:
    path = RESULTS / "hmp_missed_truth_raw_table_panel_summary.tsv"
    if not path.is_file():
        return []
    rows = []
    for row in pd.read_csv(path, sep="\t").to_dict(orient="records"):
        rows.append(
            {
                "route": str(row.get("panel", "")),
                "samples": row.get("samples", ""),
                "high_truth_FN": int(finite(row.get("high_truth_emitted_absent", 0))),
                "high_truth_FN_raw_visible": int(finite(row.get("raw_candidate_visible", 0))),
                "high_truth_FN_raw_absent": int(finite(row.get("raw_candidate_absent", 0))),
                "raw_candidate_visible_truth_mass_pct": finite(
                    row.get("raw_candidate_visible_truth_mass_pct", 0.0)
                ),
                "mean_visible_XnY_ctx": finite(row.get("mean_visible_XnY_ctx", 0.0)),
                "mean_visible_Ref_breadth": finite(row.get("mean_visible_Ref_breadth", 0.0)),
                "mean_visible_ANI": finite(row.get("mean_visible_ANI", 0.0)),
                "source": str(path),
            }
        )
    return rows


def build_audit(summary: pd.DataFrame, hmp_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total_fn = int(summary["FN"].sum()) if not summary.empty else 0
    visible_fn = int(summary["FN_raw_visible"].sum()) if not summary.empty else 0
    absent_fn = int(summary["FN_raw_absent"].sum()) if not summary.empty else 0
    high_total = int(summary["high_truth_FN"].sum()) if not summary.empty else 0
    high_visible = int(summary["high_truth_FN_raw_visible"].sum()) if not summary.empty else 0
    hmp_total = sum(int(row["high_truth_FN"]) for row in hmp_rows)
    hmp_visible = sum(int(row["high_truth_FN_raw_visible"]) for row in hmp_rows)
    return [
        {
            "metric": "evaluated_routes",
            "value": ",".join(summary["route"].astype(str)) if not summary.empty else "",
            "evidence": "marine and CAMI3 selected-default raw-table visibility",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "fn_raw_visibility",
            "value": f"{visible_fn}/{total_fn}",
            "evidence": f"raw_absent={absent_fn}; raw_visible={visible_fn}",
            "decision": "most_misses_have_raw_evidence"
            if visible_fn > absent_fn
            else "many_misses_absent_from_raw_evidence",
        },
        {
            "metric": "high_truth_fn_raw_visibility",
            "value": f"{high_visible}/{high_total}",
            "evidence": f"threshold_pct={HIGH_TRUTH_ABUNDANCE * 100.0}",
            "decision": "high_abundance_misses_are_raw_visible"
            if high_total and high_visible == high_total
            else "some_high_abundance_misses_absent_or_not_tested",
        },
        {
            "metric": "hmp_context_raw_visibility",
            "value": f"{hmp_visible}/{hmp_total}",
            "evidence": "hmp_missed_truth_raw_table_panel_summary.tsv",
            "decision": "consistent_with_raw_visible_gate_or_surface_problem"
            if hmp_total and hmp_visible == hmp_total
            else "hmp_context_missing_or_mixed",
        },
        {
            "metric": "next_strategy_implication",
            "value": "in_pass_candidate_surface_or_gate_recovery",
            "evidence": "completed negative routes mostly retain missed truth in raw best-diff rows",
            "decision": "prioritize_raw_candidate_retention_over_more_output_only_threshold_sweeps",
        },
    ]


def main() -> int:
    mapper = GtdbMapper()
    detail_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    for config in route_configs():
        detail, samples = build_route_details(config, mapper)
        detail_rows.extend(detail)
        sample_rows.extend(samples)
    sample_df = pd.DataFrame(sample_rows)
    summary_df = summarize(sample_df)
    hmp_rows = hmp_context_rows()
    audit_rows = build_audit(summary_df, hmp_rows)
    top = (
        pd.DataFrame(detail_rows)
        .sort_values(["truth_abundance", "route", "sample"], ascending=[False, True, True])
        .head(100)
        .to_dict(orient="records")
        if detail_rows
        else []
    )

    write_tsv(
        RESULTS / "negative_route_failure_mode_fn_detail.tsv",
        detail_rows,
        [
            "route",
            "sample",
            "status",
            "gtdb_species",
            "truth_abundance",
            "truth_abundance_pct",
            "high_truth_abundance",
            "raw_candidate_present",
            "raw_mode",
            "raw_accession",
            "gtdb_mapping_method",
            "ANI",
            "XnY_ctx",
            "Real_min_align_fraction",
            "Ref_breadth",
            "Ref_hit_mean_depth",
            "Ref_zip_af",
            "Ref_zip_aaf_ani",
            "Default_call",
            "Default_call_rule",
            "profile_path",
            "raw_path",
        ],
    )
    write_tsv(
        RESULTS / "negative_route_failure_mode_fn_top.tsv",
        top,
        [
            "route",
            "sample",
            "status",
            "gtdb_species",
            "truth_abundance",
            "truth_abundance_pct",
            "high_truth_abundance",
            "raw_candidate_present",
            "raw_mode",
            "raw_accession",
            "gtdb_mapping_method",
            "ANI",
            "XnY_ctx",
            "Real_min_align_fraction",
            "Ref_breadth",
            "Ref_hit_mean_depth",
            "Ref_zip_af",
            "Ref_zip_aaf_ani",
            "Default_call",
            "Default_call_rule",
        ],
    )
    sample_df.to_csv(RESULTS / "negative_route_failure_mode_sample_summary.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / "negative_route_failure_mode_summary.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "negative_route_failure_mode_hmp_context.tsv",
        hmp_rows,
        [
            "route",
            "samples",
            "high_truth_FN",
            "high_truth_FN_raw_visible",
            "high_truth_FN_raw_absent",
            "raw_candidate_visible_truth_mass_pct",
            "mean_visible_XnY_ctx",
            "mean_visible_Ref_breadth",
            "mean_visible_ANI",
            "source",
        ],
    )
    write_tsv(
        RESULTS / "negative_route_failure_mode_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )
    print(summary_df.to_string(index=False))
    print("\nHMP CONTEXT")
    print(pd.DataFrame(hmp_rows).to_string(index=False) if hmp_rows else "no hmp context")
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
