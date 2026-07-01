#!/usr/bin/env python3
"""Offline sweep for the experimental exact-hit abundance sidecar trigger.

Calls are held fixed from the current MinCO profile tables. For samples with a
cached exact split raw table, this script simulates the abundance-only sidecar:
when the block extra-mass ratio is at least a threshold and the auto-exact
guards pass, called rows get abundance mass from exact split Ref_hit_mean_depth
matched by best reference accession. Samples without a cached exact table are
left at the current abundance and flagged in the trace output.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
DECISION_EXP = EXP.parents[0] / "2026-06-27_universal_strategy_decision"
sys.path.insert(0, str(DECISION_EXP))

import sweep_cross_panel_abundance_variants as sweep  # noqa: E402


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
THRESHOLDS = [math.inf, 0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def truthy(value: object, default: bool = False) -> bool:
    if pd.isna(value):
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    return default


def first_value(df: pd.DataFrame, col: str, default: object = "") -> object:
    if col not in df.columns or df.empty:
        return default
    return df[col].iloc[0]


def extract_accession(value: object) -> str:
    match = ACC_RE.search(str(value or ""))
    return match.group(1) if match else ""


def panel_abundance_raw(features: pd.DataFrame) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        split_raw = numeric(features, "s_Ref_mean_depth_max").to_numpy(dtype=float) / np.maximum(
            numeric(features, "s_Ref_zip_af_max").to_numpy(dtype=float),
            1e-6,
        )
        unique_raw = numeric(features, "u_Ref_mean_depth_max").to_numpy(dtype=float) / np.maximum(
            numeric(features, "u_Ref_zip_af_max").to_numpy(dtype=float),
            1e-6,
        )
    raw = np.maximum(split_raw, unique_raw)
    if "tail_rescue_added" in features.columns:
        tail_added = features["tail_rescue_added"].astype(bool).to_numpy()
        raw = np.asarray(raw, dtype=float)
        raw[tail_added] = numeric(features, "s_Ref_mean_depth_max").to_numpy(dtype=float)[tail_added]
    return np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)


def load_exact_by_accession(path: Path) -> Mapping[str, float]:
    raw = pd.read_csv(
        path,
        sep="\t",
        usecols=lambda col: col in {"Ref", "accession", "Ref_hit_mean_depth"},
        low_memory=False,
    )
    if "Ref_hit_mean_depth" not in raw.columns:
        return {}
    if "accession" in raw.columns:
        accessions = raw["accession"].astype(str)
    elif "Ref" in raw.columns:
        accessions = raw["Ref"].map(extract_accession)
    else:
        return {}
    work = pd.DataFrame(
        {
            "accession": accessions,
            "exact_hit": pd.to_numeric(raw["Ref_hit_mean_depth"], errors="coerce").fillna(0.0),
        }
    )
    work = work.loc[work["accession"].astype(str).astype(bool)]
    if work.empty:
        return {}
    return work.groupby("accession", sort=False)["exact_hit"].max().to_dict()


def discover_exact_path(panel: str, sample: int, profile_path: Path, calls: pd.DataFrame) -> Path | None:
    candidates: list[Path] = []
    for raw in calls.get("auto_exact_split_path", pd.Series(dtype=object)).dropna().astype(str).unique():
        if raw and raw.lower() != "nan":
            candidates.append(Path(raw))

    parent = profile_path.parent
    candidates.append(parent / "work/minco.best_diff_split.exact.unfiltered.tsv")
    if panel == "cami2_toy_mouse_gut":
        candidates.append(
            Path(f"/tmp/cami2_toymouse_current_default_20260626/run/sample{sample}_default/minco.best_diff_split.exact.unfiltered.tsv")
        )
    elif panel == "hmp_airskin_gtdb_source_abundance":
        candidates.extend(
            [
                Path(f"/tmp/minco_current_code_hmp_airskin{sample}_20260627/work/minco.best_diff_split.exact.unfiltered.tsv"),
                Path(f"/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin{sample}_20260627/work/minco.best_diff_split.exact.unfiltered.tsv"),
                Path(f"/tmp/minco_current_code_hmp_refresh_20260627/minco_sample{sample}_split_exact_unfiltered.tsv"),
            ]
        )
    elif panel == "hmp_gastrooral_gtdb_source_abundance":
        candidates.extend(
            [
                Path(f"/tmp/minco_current_code_hmp_gastrooral_20260627/sample{sample}_work/minco.best_diff_split.exact.unfiltered.tsv"),
                Path(f"/tmp/minco_current_code_hmp_gastrooral_sample{sample}_20260630/sample{sample}_work/minco.best_diff_split.exact.unfiltered.tsv"),
                Path(f"/tmp/minco_exactsplit_universal_20260626/hmp_gastrooral{sample}_s1000_gtdb_split_exact_unfiltered.tsv"),
            ]
        )
        if sample == 6:
            candidates.append(
                Path("/tmp/minco_exact_hit_default_20260630/sample6_forced_work/minco.best_diff_split.exact.unfiltered.tsv")
            )
    elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
        candidates.extend(
            [
                Path(f"/tmp/cami3_toy_human_gut_20260626/run/sample{sample}_exact_forced/minco.best_diff_split.exact.unfiltered.tsv"),
                Path(f"/tmp/cami3_toy_human_gut_20260626/run/sample{sample}_autoexact/minco.best_diff_split.exact.unfiltered.tsv"),
            ]
        )

    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.exists():
            return candidate
    return None


def trigger_state(calls: pd.DataFrame, threshold: float) -> tuple[bool, dict[str, object]]:
    p_extra = float(first_value(calls, "auto_exact_split_block_p_extra_mass_ratio", 0.0) or 0.0)
    guard = truthy(first_value(calls, "auto_exact_split_guard_passed", True), True)
    low_extra = truthy(first_value(calls, "auto_exact_split_low_extra_gate_passed", True), True)
    if math.isinf(threshold):
        requested = False
    else:
        requested = p_extra >= threshold and guard and low_extra
    return requested, {
        "p_extra": p_extra,
        "guard_passed": guard,
        "low_extra_gate_passed": low_extra,
        "threshold": "disabled" if math.isinf(threshold) else threshold,
        "trigger_requested": requested,
    }


def exact_hit_raw_values(
    calls: pd.DataFrame,
    exact_path: Path,
) -> tuple[np.ndarray, dict[str, object]]:
    exact_by_accession = load_exact_by_accession(exact_path)
    fallback = panel_abundance_raw(calls)
    if not exact_by_accession or "s_best_accession" not in calls.columns:
        return fallback, {
            "exact_accessions": len(exact_by_accession),
            "exact_matched_rows": 0,
            "exact_fallback_rows": len(calls),
        }
    exact = (
        calls["s_best_accession"]
        .astype(str)
        .map(exact_by_accession)
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    use_exact = exact > 0.0
    raw = np.where(use_exact, exact, fallback)
    return np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0), {
        "exact_accessions": len(exact_by_accession),
        "exact_matched_rows": int(use_exact.sum()),
        "exact_fallback_rows": int((~use_exact).sum()),
    }


def threshold_label(threshold: float) -> str:
    if math.isinf(threshold):
        return "current_no_sidecar"
    return f"exact_hit_abund_ge_{threshold:.2f}"


def sample_records() -> list[tuple[str, int]]:
    return [
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


def load_panel_contexts() -> dict[str, object]:
    toy_mod = sweep.decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    cami3_wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _diag = sweep.cami3.build_transfer_maps()
    return {
        "by_accession": by_accession,
        "by_core": by_core,
        "hmp_taxmap": sweep.hmp.parse_species_taxmap(sweep.hmp.TAXMAP),
        "hmp_gastro_taxmap": sweep.hmp_gastro.parse_species_taxmap(sweep.hmp_gastro.TAXMAP),
        "cami3_taxmap": sweep.cami3.parse_species_taxmap(sweep.cami3.TAXMAP),
        "cami3_taxid_to_species": cami3_taxid_to_species,
        "cami3_name_to_species": cami3_name_to_species,
        "cami3_wgs_to_species": cami3_wgs_to_species,
    }


def load_sample(panel: str, sample: int, ctx: Mapping[str, object]) -> tuple[pd.DataFrame, pd.DataFrame, str, Path]:
    by_accession = ctx["by_accession"]
    by_core = ctx["by_core"]
    if panel == "cami2_toy_mouse_gut":
        return (*sweep.load_toy_calls(sample), sweep.toy_truth(sample))  # type: ignore[misc]
    if panel == "hmp_airskin_gtdb_source_abundance":
        calls, collapse, path = sweep.load_hmp_calls(sample, ctx["hmp_taxmap"], by_accession, by_core)
        return calls, collapse, path, sweep.hmp_truth(sample)
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        calls, collapse, path = sweep.load_hmp_gastro_calls(
            sample,
            ctx["hmp_gastro_taxmap"],
            by_accession,
            by_core,
        )
        return calls, collapse, path, sweep.hmp_gastro_truth(sample)
    if panel == "cami3_toy_human_gut_gtdb_source_readmap":
        calls, collapse, path = sweep.load_cami3_calls(
            sample,
            ctx["cami3_taxid_to_species"],
            ctx["cami3_name_to_species"],
            ctx["cami3_taxmap"],
            by_accession,
            by_core,
        )
        return calls, collapse, path, sweep.cami3_truth(sample)
    raise ValueError(panel)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    ctx = load_panel_contexts()
    score_rows: list[dict[str, object]] = []
    trace_rows: list[dict[str, object]] = []
    mapping_rows: list[dict[str, object]] = []

    for panel, sample in sample_records():
        calls, collapse, profile_path, truth = load_sample(panel, sample, ctx)
        exact_path = discover_exact_path(panel, sample, profile_path, calls)
        mapping_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "profile": str(profile_path),
                "collapse_rule": collapse,
                "called_rows_mapped": len(calls),
                "called_species_mapped": calls["gtdb_species"].nunique(),
                "exact_path": str(exact_path or ""),
                "exact_available": bool(exact_path),
            }
        )

        baseline_methods = {
            "current_calibrated_abundance": numeric(calls, "calibrated_abundance").to_numpy(dtype=float),
            "current_calibrated_raw": numeric(calls, "calibrated_abundance_raw").to_numpy(dtype=float),
        }
        for method, raw_values in baseline_methods.items():
            pred = sweep.collapse_prediction(calls, raw_values, collapse)
            row = sweep.taxid_score.score_prediction(sample, method, pred, truth, {})
            row["panel"] = panel
            row["collapse_rule"] = collapse
            score_rows.append(row)

        exact_cache: tuple[np.ndarray, dict[str, object]] | None = None
        for threshold in THRESHOLDS:
            if math.isinf(threshold):
                continue
            requested, state = trigger_state(calls, threshold)
            if requested and exact_path:
                if exact_cache is None:
                    exact_cache = exact_hit_raw_values(calls, exact_path)
                raw_values, exact_info = exact_cache
                applied = True
            else:
                raw_values = numeric(calls, "calibrated_abundance_raw").to_numpy(dtype=float)
                exact_info = {
                    "exact_accessions": 0,
                    "exact_matched_rows": 0,
                    "exact_fallback_rows": 0,
                }
                applied = False
            method = threshold_label(threshold)
            pred = sweep.collapse_prediction(calls, raw_values, collapse)
            row = sweep.taxid_score.score_prediction(sample, method, pred, truth, {})
            row["panel"] = panel
            row["collapse_rule"] = collapse
            score_rows.append(row)
            trace_rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "method": method,
                    "profile": str(profile_path),
                    "exact_path": str(exact_path or ""),
                    "exact_available": bool(exact_path),
                    "exact_applied": applied,
                    **state,
                    **exact_info,
                }
            )

    scores = pd.DataFrame(score_rows)
    panel_summary, overall = sweep.summarize(scores)
    panel_with_sylph = sweep.add_external_sylph(panel_summary)

    scores.to_csv(RESULTS / "exact_hit_abundance_trigger_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(RESULTS / "exact_hit_abundance_trigger_panel_summary.tsv", sep="\t", index=False)
    panel_with_sylph.to_csv(RESULTS / "exact_hit_abundance_trigger_panel_summary_with_sylph.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / "exact_hit_abundance_trigger_overall.tsv", sep="\t", index=False)
    pd.DataFrame(trace_rows).to_csv(RESULTS / "exact_hit_abundance_trigger_trace.tsv", sep="\t", index=False)
    pd.DataFrame(mapping_rows).to_csv(RESULTS / "exact_hit_abundance_trigger_mapping.tsv", sep="\t", index=False)

    print("TOP OVERALL")
    print(overall.head(20).to_string(index=False))
    print("\nPANEL SUMMARY WITH SYLPH")
    show = panel_with_sylph.loc[
        panel_with_sylph["method"].isin(
            {
                "current_calibrated_abundance",
                "current_calibrated_raw",
                "exact_hit_abund_ge_0.30",
                "exact_hit_abund_ge_0.20",
                "exact_hit_abund_ge_0.10",
                "exact_hit_abund_ge_0.00",
                "sylph_external_baseline",
            }
        )
    ].copy()
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
