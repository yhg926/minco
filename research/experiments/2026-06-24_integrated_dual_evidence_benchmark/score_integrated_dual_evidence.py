#!/usr/bin/env python3
"""Score integrated full-S2000 readwise output with dual marker evidence."""

from __future__ import annotations

import math
import sys
import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
RESULTS = EXP_DIR / "results"
RUN_DIR = Path("/tmp/minco_dual_s2000_20260624")

HYBRID_EXP = ROOT / "research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback"
CAMI3_BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
CAMI3_RESCUE_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue"
MOUSE_ABUND_EXP = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"
MOUSE_MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
MOUSE_TRUTH_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
ABUND_CAL_EXP = ROOT / "research/experiments/2026-06-24_abundance_feature_calibration"
ABUND_ADD_EXP = ROOT / "research/experiments/2026-06-24_addback_abundance_normalization_search"

for path in [
    HYBRID_EXP,
    CAMI3_BASE_EXP,
    CAMI3_RESCUE_EXP,
    MOUSE_ABUND_EXP,
    MOUSE_MORE_EXP,
    MOUSE_TRUTH_EXP,
]:
    sys.path.insert(0, str(path))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import evaluate_abundance_estimators as mouse_ev  # noqa: E402
import score_cami3_abundance_rescue as cami3_rescue  # noqa: E402
import score_cami3_toygut as cami3_base  # noqa: E402
import score_hybrid_cami3 as hybrid_cami3  # noqa: E402
import score_hybrid_mouse as hybrid_mouse  # noqa: E402
import score_hybrid_mouse_baseline_plus as mouse_baseline_plus  # noqa: E402
import score_intragenus_winner_rescue as mouse_rescue  # noqa: E402
import score_more_toymouse_gtdb as mouse_more  # noqa: E402


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


abn = load_module("integrated_addback_abundance", ABUND_ADD_EXP / "search_abundance_normalization.py")


MARKER_MAP = {
    "Marker_XnY_ctx": "XnY_ctx",
    "Marker_N_diff_obj": "N_diff_obj",
    "Marker_N_diff_obj_section": "N_diff_obj_section",
    "Marker_N_mut2_ctx": "N_mut2_ctx",
    "Marker_Ref_breadth": "Ref_breadth",
    "Marker_Ref_mean_depth": "Ref_mean_depth",
    "Marker_Ref_hit_mean_depth": "Ref_hit_mean_depth",
    "Marker_Ref_depth_variance": "Ref_depth_variance",
    "Marker_Ref_depth_cv": "Ref_depth_cv",
    "Marker_Ref_zero_fraction": "Ref_zero_fraction",
    "Marker_Relative_abundance_depth": "Relative_abundance_depth",
    "Marker_Effective_abundance_depth": "Effective_abundance_depth",
    "Marker_Ref_zip_af": "Ref_zip_af",
    "Marker_Ref_zip_aaf_ani": "Ref_zip_aaf_ani",
    "Marker_Reliable_Ref_breadth": "Reliable_Ref_breadth",
    "Marker_Reliable_Ref_mean_depth": "Reliable_Ref_mean_depth",
    "Marker_Reliable_Ref_hit_ctx": "Reliable_Ref_hit_ctx",
    "Marker_Reliable_Ref_hit_mean_depth": "Reliable_Ref_hit_mean_depth",
    "Marker_Reliable_Ref_hit_median_depth": "Reliable_Ref_hit_median_depth",
    "Marker_Reliable_Ref_hit_depth_variance": "Reliable_Ref_hit_depth_variance",
    "Marker_Reliable_Ref_zip_af": "Reliable_Ref_zip_af",
}

MOUSE_READWISE = {
    0: RUN_DIR / "mouse_s0_dual.tsv",
    1: RUN_DIR / "mouse_s1_dual.tsv",
    2: RUN_DIR / "mouse_s2_dual.tsv",
}

CAMI3_READWISE = {
    0: RUN_DIR / "cami3_s0_dual.tsv",
    1: RUN_DIR / "cami3_s1_dual.tsv",
    2: RUN_DIR / "cami3_s2_dual.tsv",
}


def numeric(rows: pd.DataFrame, col: str) -> pd.Series:
    if col not in rows.columns:
        return pd.Series([0.0] * len(rows), index=rows.index)
    return pd.to_numeric(rows[col], errors="coerce").fillna(0.0)


def accession_from_rows(rows: pd.DataFrame) -> pd.DataFrame:
    rows = rows.copy()
    rows["accession"] = rows["Ref"].map(truth.extract_accession)
    if "Ref_annotation" in rows.columns:
        missing = rows["accession"] == ""
        rows.loc[missing, "accession"] = rows.loc[missing, "Ref_annotation"].map(truth.extract_accession)
    return rows


def marker_view(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    for src, dst in MARKER_MAP.items():
        out[dst] = numeric(out, src)
    out["ctx_marker_size"] = numeric(out, "ctx_marker_size").astype(int)
    out["source_mode"] = "ctx_marker"
    return out


def full_view(rows: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    out["ctx_marker_size"] = numeric(out, "ctx_marker_size").astype(int)
    out["source_mode"] = "full_fallback"
    return out


def add_mouse_metadata_and_gate(rows: pd.DataFrame, by_accession, by_core) -> pd.DataFrame:
    rows = accession_from_rows(rows)
    accession_map = {}
    for acc in rows["accession"].dropna().astype(str).unique():
        rec, method, key = truth.lookup_accession(acc, by_accession, by_core)
        accession_map[acc] = (
            rec["gtdb_species"] if rec else "",
            rec["gtdb_taxonomy"] if rec else "",
            method,
            key,
        )
    rows["gtdb_species"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[0])
    rows["gtdb_taxonomy"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[1])
    rows["gtdb_mapping_method"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[2])
    rows["gtdb_mapping_key"] = rows["accession"].map(lambda acc: accession_map.get(str(acc), ("", "", "", ""))[3])

    rows = truth.add_naive_ani(rows)
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_median_depth",
        "Reliable_Ref_hit_depth_variance",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
    ]:
        rows[col] = numeric(rows, col)

    rows["Reliable_ztp_af"] = [
        mouse_more.ztp_adjusted_af(float(b), float(m))
        for b, m in zip(rows["Reliable_Ref_breadth"], rows["Reliable_Ref_hit_mean_depth"])
    ]
    rows["ANI_from_Reliable_ztp_af"] = [
        1.0 + math.log(max(float(af), 1e-300)) / mouse_more.EFFECTIVE_CTX_LENGTH
        for af in rows["Reliable_ztp_af"]
    ]
    rows["ANI_AF_delta"] = rows["ANI_naive_calc"] - rows["ANI_from_Reliable_ztp_af"]
    rows["Reliable_depth_vmr"] = [
        (float(var) / float(mean)) if float(mean) > 0.0 else 0.0
        for mean, var in zip(rows["Reliable_Ref_hit_mean_depth"], rows["Reliable_Ref_hit_depth_variance"])
    ]
    rows["active_delta_trigger"] = (
        (rows["Reliable_Ref_hit_mean_depth"] > mouse_more.ACTIVE_MEAN_DEPTH_MIN)
        & (rows["Reliable_depth_vmr"] > mouse_more.ACTIVE_VMR_MIN)
    )
    rows["active_delta_pass"] = ~rows["active_delta_trigger"] | (rows["ANI_AF_delta"] < mouse_more.ACTIVE_DELTA_MAX)
    rows["active_gate_pass"] = (
        (rows["XnY_ctx"] >= mouse_more.ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > mouse_more.ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= mouse_more.ACTIVE_RELIABLE_ZTP_AF_FLOOR)
        & rows["active_delta_pass"]
        & rows["gtdb_species"].astype(bool)
    )

    median = numeric(rows, "Reliable_Ref_hit_median_depth")
    mean = numeric(rows, "Reliable_Ref_mean_depth")
    af = numeric(rows, "Reliable_Ref_zip_af").clip(lower=mouse_ev.EPS)
    rows["effective_depth"] = mean.where(median < mouse_rescue.MEDIAN_DEPTH_CUTOFF, median)
    rows.loc[median < mouse_rescue.MEDIAN_DEPTH_CUTOFF, "effective_depth"] = (
        mean.loc[median < mouse_rescue.MEDIAN_DEPTH_CUTOFF] /
        (af.loc[median < mouse_rescue.MEDIAN_DEPTH_CUTOFF] ** mouse_rescue.AF_EXPONENT)
    )
    rows.loc[mean <= 0.0, "effective_depth"] = 0.0
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(mouse_rescue.species_genus)
    return rows


def load_mouse_views(sample: int, by_accession, by_core) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = pd.read_csv(MOUSE_READWISE[sample], sep="\t")
    marker = add_mouse_metadata_and_gate(marker_view(rows), by_accession, by_core)
    full = add_mouse_metadata_and_gate(full_view(rows), by_accession, by_core)
    return marker, full


def add_cami3_metadata_and_gate(rows: pd.DataFrame, by_accession, by_core) -> pd.DataFrame:
    rows = accession_from_rows(rows)
    rows = cami3_base.add_metadata(rows, "accession", by_accession, by_core)
    return cami3_base.add_reliable_ztp_and_active_gate(rows)


def load_cami3_views(sample: int, by_accession, by_core) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = pd.read_csv(CAMI3_READWISE[sample], sep="\t")
    marker = add_cami3_metadata_and_gate(marker_view(rows), by_accession, by_core)
    full = add_cami3_metadata_and_gate(full_view(rows), by_accession, by_core)
    return marker, full


def load_integrated_addback_samples():
    marker_samples = []
    addback_samples = []

    mouse_by_accession, mouse_by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    for sample_id in hybrid_mouse.SAMPLES:
        marker, full = load_mouse_views(sample_id, mouse_by_accession, mouse_by_core)
        marker_prepared = abn.prev.prepare_rows(
            marker,
            dataset="mouse_gtdb",
            sample_id=sample_id,
            target_col="gtdb_species",
        )
        full_prepared = abn.prev.prepare_rows(
            full,
            dataset="mouse_gtdb",
            sample_id=sample_id,
            target_col="gtdb_species",
        )
        truth_profile = mouse_ev.load_truth_profile(sample_id)
        marker_samples.append(
            abn.prev.PreparedSample(
                dataset="mouse_gtdb",
                sample_id=sample_id,
                truth=truth_profile,
                marker=marker_prepared,
                full=full_prepared,
            )
        )
        addback_samples.append(
            abn.prev.PreparedSample(
                dataset="mouse_gtdb",
                sample_id=sample_id,
                truth=truth_profile,
                marker=full_prepared,
                full=pd.DataFrame(),
            )
        )

    cami_by_accession, cami_by_core = cami3_base.truth.load_gtdb_metadata(cami3_base.truth.GTDB_METADATA)
    for sample_id, paths in cami3_base.SAMPLES.items():
        marker, full = load_cami3_views(sample_id, cami_by_accession, cami_by_core)
        marker_prepared = abn.prev.prepare_rows(
            marker,
            dataset="cami3_ncbi",
            sample_id=sample_id,
            target_col="ncbi_species_taxid",
        )
        full_prepared = abn.prev.prepare_rows(
            full,
            dataset="cami3_ncbi",
            sample_id=sample_id,
            target_col="ncbi_species_taxid",
        )
        truth_rows = cami3_base.load_truth(paths["truth"])
        truth_profile = dict(
            zip(
                truth_rows["ncbi_species_taxid"].astype(str),
                truth_rows["truth_abundance_bacterial_norm"].astype(float),
            )
        )
        marker_samples.append(
            abn.prev.PreparedSample(
                dataset="cami3_ncbi",
                sample_id=sample_id,
                truth=truth_profile,
                marker=marker_prepared,
                full=full_prepared,
            )
        )
        addback_samples.append(
            abn.prev.PreparedSample(
                dataset="cami3_ncbi",
                sample_id=sample_id,
                truth=truth_profile,
                marker=full_prepared,
                full=pd.DataFrame(),
            )
        )

    return marker_samples, addback_samples


def filtered_add_values(marker_raw: dict[str, float], add_raw: dict[str, float], *, include_overlap: bool) -> dict[str, float]:
    out: dict[str, float] = {}
    for target, value in add_raw.items():
        if target in marker_raw and not include_overlap:
            continue
        if value > 0.0:
            out[target] = float(value)
    return out


def sum_values(marker_raw: dict[str, float], add_raw: dict[str, float]) -> dict[str, float]:
    values = dict(marker_raw)
    for target, value in add_raw.items():
        values[target] = values.get(target, 0.0) + float(value)
    return values


def em_proxy_values(
    marker_raw: dict[str, float],
    add_raw: dict[str, float],
    *,
    beta: float,
    seed_fraction: float,
    prior_power: float,
) -> dict[str, float]:
    if not add_raw:
        return dict(marker_raw)
    targets = sorted(set(marker_raw) | set(add_raw))
    priors = {}
    for target in targets:
        prior = float(marker_raw.get(target, 0.0)) + seed_fraction * float(add_raw.get(target, 0.0))
        if prior > 0.0:
            priors[target] = prior ** prior_power
    if not priors:
        priors = {target: float(value) ** prior_power for target, value in add_raw.items() if value > 0.0}
    prior_total = sum(priors.values())
    add_total = sum(float(v) for v in add_raw.values() if float(v) > 0.0)
    assigned = {
        target: add_total * value / prior_total
        for target, value in priors.items()
        if prior_total > 0.0
    }
    values = dict(marker_raw)
    for target, value in add_raw.items():
        values[target] = values.get(target, 0.0) + (1.0 - beta) * float(value)
    for target, value in assigned.items():
        values[target] = values.get(target, 0.0) + beta * float(value)
    return {target: value for target, value in values.items() if value > 0.0}


def score_integrated_addback() -> pd.DataFrame:
    marker_samples, addback_samples = load_integrated_addback_samples()
    addback_by_key = {(sample.dataset, int(sample.sample_id)): sample for sample in addback_samples}
    callsets = [
        ("integrated_value_sum_marker_l1_full_f1_blend", abn.MARKER_L1, abn.CTXOBJ_F1, "all"),
        ("integrated_value_sum_marker_l1_full_l1_blend", abn.MARKER_L1, abn.CTXOBJ_L1, "all"),
        ("integrated_value_sum_marker_l1_full_cami_blend", abn.MARKER_L1, abn.CTXOBJ_CAMI, "all"),
        ("integrated_value_sum_marker_l1_full_f1_marker_present", abn.MARKER_L1, abn.CTXOBJ_F1, "marker_present"),
        ("integrated_value_sum_marker_l1_full_l1_marker_present", abn.MARKER_L1, abn.CTXOBJ_L1, "marker_present"),
        ("integrated_value_sum_marker_l1_full_cami_marker_present", abn.MARKER_L1, abn.CTXOBJ_CAMI, "marker_present"),
    ]

    sample_rows: list[dict[str, object]] = []
    raw_by_callset: dict[tuple[str, str, int], tuple[dict[str, float], dict[str, float]]] = {}
    for method, marker_combo, add_combo, add_filter in callsets:
        for sample in marker_samples:
            add_sample = addback_by_key[(sample.dataset, int(sample.sample_id))]
            marker_raw = abn.raw_values_for_sample(sample, marker_combo)
            add_raw = filtered_add_values(
                marker_raw,
                abn.raw_values_for_sample(add_sample, add_combo),
                include_overlap=True,
            )
            if add_filter == "marker_present":
                add_raw = {target: value for target, value in add_raw.items() if target in marker_raw}
            raw_by_callset[(method, sample.dataset, int(sample.sample_id))] = (marker_raw, add_raw)
            sample_rows.append(
                abn.score_prediction(
                    sample,
                    abn.normalize(sum_values(marker_raw, add_raw)),
                    method,
                )
            )

    em_rows: list[dict[str, object]] = []
    for source_method, _, _, _ in callsets:
        for beta in [0.25, 0.5, 0.75, 1.0]:
            for seed_fraction in [0.0, 0.001, 0.01, 0.05, 0.1]:
                for prior_power in [0.5, 1.0, 2.0]:
                    method = (
                        f"{source_method}_em_proxy"
                        f"_b{beta:g}_s{seed_fraction:g}_p{prior_power:g}"
                    )
                    for sample in marker_samples:
                        marker_raw, add_raw = raw_by_callset[(source_method, sample.dataset, int(sample.sample_id))]
                        pred = abn.normalize(
                            em_proxy_values(
                                marker_raw,
                                add_raw,
                                beta=beta,
                                seed_fraction=seed_fraction,
                                prior_power=prior_power,
                            )
                        )
                        em_rows.append(abn.score_prediction(sample, pred, method))

    out = pd.DataFrame(sample_rows + em_rows)
    out.to_csv(RESULTS / "integrated_addback_sample_metrics.tsv", sep="\t", index=False)
    summary = summarize_addback(out)
    summary.to_csv(RESULTS / "integrated_addback_summary.tsv", sep="\t", index=False)
    summary.loc[summary["method"].str.contains("_em_proxy", regex=False)].head(50).to_csv(
        RESULTS / "integrated_em_proxy_top50.tsv",
        sep="\t",
        index=False,
    )
    return summary


def summarize_addback(sample_rows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, group in sample_rows.groupby("method", sort=False):
        rec = {"method": method}
        for prefix, sub in [
            ("all", group),
            ("mouse", group.loc[group["dataset"] == "mouse_gtdb"]),
            ("cami3", group.loc[group["dataset"] == "cami3_ncbi"]),
        ]:
            rec[f"{prefix}_F1"] = float(sub["F1"].mean()) if len(sub) else float("nan")
            rec[f"{prefix}_L1_pp"] = float(sub["l1_pct_points"].mean()) if len(sub) else float("nan")
            rec[f"{prefix}_Pearson"] = float(sub["pearson"].mean()) if len(sub) else float("nan")
            rec[f"{prefix}_TP"] = float(sub["TP"].mean()) if len(sub) else float("nan")
            rec[f"{prefix}_FP"] = float(sub["FP"].mean()) if len(sub) else float("nan")
            rec[f"{prefix}_FN"] = float(sub["FN"].mean()) if len(sub) else float("nan")
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["all_F1", "all_L1_pp"], ascending=[False, True])


def add_cami3_mode_columns(rows: pd.DataFrame, source_mode: str) -> pd.DataFrame:
    rows = cami3_rescue.add_effective_depth(rows)
    rows = rows.copy()
    rows["source_mode"] = source_mode
    rows["ctx_marker_size"] = numeric(rows, "ctx_marker_size").astype(int)
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(hybrid_cami3.species_genus)
    rows["pred_depth"] = cami3_base.numeric(rows, "effective_depth")
    rows.loc[rows["source_mode"] == "full_fallback", "pred_depth"] *= hybrid_cami3.FALLBACK_ABUNDANCE_SCALE
    return rows


def score_mouse() -> pd.DataFrame:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    sample_rows = []
    selected_rows = []
    for sample in hybrid_mouse.SAMPLES:
        marker, full = load_mouse_views(sample, by_accession, by_core)
        old_values, old_rescued = mouse_baseline_plus.values_scaled(marker, 1.0)
        sample_rows.append(
            mouse_baseline_plus.score_values(
                sample,
                "integrated_marker_only",
                old_values,
                old_rescued,
            )
        )
        for method, fallback_xny_min, scale in mouse_baseline_plus.STRATEGIES:
            rows = pd.concat(
                [
                    marker.copy(),
                    full.loc[
                        (full["ctx_marker_size"] < 200)
                        & (full["XnY_ctx"] >= fallback_xny_min)
                    ].copy(),
                ],
                ignore_index=True,
                sort=False,
            ).fillna(0)
            values, rescued = mouse_baseline_plus.values_scaled(rows, scale)
            sample_rows.append(score_mouse_values(sample, f"integrated_{method}", values, rescued))
            selected_rows.append(rows.loc[rows["active_gate_pass"]].assign(sample_id=sample, method=f"integrated_{method}"))

    sample_df = pd.DataFrame(sample_rows)
    sample_df.to_csv(RESULTS / "integrated_mouse_sample_metrics.tsv", sep="\t", index=False)
    if selected_rows:
        pd.concat(selected_rows, ignore_index=True).to_csv(
            RESULTS / "integrated_mouse_active_rows.tsv", sep="\t", index=False
        )
    mean_cols = [
        "TP",
        "FP",
        "FN",
        "missing_truth_mass",
        "fp_pred_mass",
        "pred_sum_on_truth",
        "pred_sum_all",
        "pearson",
        "spearman",
        "mae_pct_points",
        "l1_pct_points",
    ]
    mean_df = sample_df.groupby("method", as_index=False)[mean_cols].mean()
    mean_df["precision"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FP"])
    mean_df["recall"] = mean_df["TP"] / (mean_df["TP"] + mean_df["FN"])
    mean_df["F1"] = 2.0 * mean_df["precision"] * mean_df["recall"] / (mean_df["precision"] + mean_df["recall"])
    mean_df = mean_df.sort_values(["F1", "l1_pct_points"], ascending=[False, True])
    mean_df.to_csv(RESULTS / "integrated_mouse_mean_metrics.tsv", sep="\t", index=False)
    return mean_df


def score_mouse_values(sample: int, method: str, values: dict[str, float], rescued: list[str]) -> dict[str, object]:
    return mouse_rescue.score_values(sample, method, values, mouse_ev.load_truth_profile(sample), rescued)


def score_cami3() -> pd.DataFrame:
    by_accession, by_core = cami3_base.truth.load_gtdb_metadata(cami3_base.truth.GTDB_METADATA)
    score_rows = []
    abundance_rows = []
    detail_rows = []
    selected_rows = []
    rescue_rows = []
    for sample_id, paths in cami3_base.SAMPLES.items():
        truth_rows = cami3_base.load_truth(paths["truth"])
        marker, full = load_cami3_views(sample_id, by_accession, by_core)
        marker = add_cami3_mode_columns(marker, "ctx_marker")
        full = add_cami3_mode_columns(full, "full_fallback")

        methods = {"integrated_marker_only": marker}
        for method_name, fallback_xny_min in hybrid_cami3.HYBRID_STRATEGIES:
            methods[f"integrated_{method_name}"] = pd.concat(
                [
                    marker.copy(),
                    full.loc[
                        (full["ctx_marker_size"] < hybrid_cami3.MARKER_SIZE_CUTOFF)
                        & (full["XnY_ctx"] >= fallback_xny_min)
                    ].copy(),
                ],
                ignore_index=True,
                sort=False,
            ).fillna(0)

        for method, rows in methods.items():
            selected, rescued = hybrid_cami3.select_robust_rescue_fast(rows)
            for rec in rescued:
                rec["sample_id"] = sample_id
                rec["method"] = method
                rescue_rows.append(rec)
            score, abundance, details, selected_out = hybrid_cami3.score_selected(
                sample_id, method, selected, truth_rows
            )
            score_rows.append(score)
            abundance_rows.extend(abundance)
            detail_rows.extend(details)
            selected_rows.append(selected_out)

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    detail_df = pd.DataFrame(detail_rows)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()
    rescue_df = pd.DataFrame(rescue_rows)

    score_df.to_csv(RESULTS / "integrated_cami3_sample_presence.tsv", sep="\t", index=False)
    abundance_df.to_csv(RESULTS / "integrated_cami3_sample_abundance.tsv", sep="\t", index=False)
    detail_df.to_csv(RESULTS / "integrated_cami3_fp_fn_details.tsv", sep="\t", index=False)
    selected_df.to_csv(RESULTS / "integrated_cami3_selected_species.tsv", sep="\t", index=False)
    rescue_df.to_csv(RESULTS / "integrated_cami3_rescued_species.tsv", sep="\t", index=False)

    mean_presence = (
        score_df.groupby("method", as_index=False)[
            ["pred_species", "TP", "FP", "FN", "precision", "recall", "F1"]
        ].mean()
    )
    mean_abundance = (
        abundance_df.loc[abundance_df["renorm_pred"]]
        .groupby("method", as_index=False)[
            ["pred_sum_on_truth", "pred_sum_all", "pearson", "spearman", "mae_pct_points", "l1_pct_points"]
        ]
        .mean()
    )
    summary = mean_presence.merge(mean_abundance, on="method", how="outer")
    summary = summary.sort_values(["F1", "l1_pct_points"], ascending=[False, True])
    summary.to_csv(RESULTS / "integrated_cami3_summary.tsv", sep="\t", index=False)
    return summary


def load_cached_current_rows() -> pd.DataFrame:
    final_rows = pd.read_csv(ABUND_CAL_EXP / "results/final_comparison.tsv", sep="\t")
    top_l1 = pd.read_csv(ABUND_CAL_EXP / "results/calibration_top200_by_l1.tsv", sep="\t")
    top_f1 = pd.read_csv(ABUND_CAL_EXP / "results/calibration_top200_by_f1.tsv", sep="\t")
    rows = []
    sylph = final_rows.loc[final_rows["method"] == "Sylph"].iloc[0]
    rows.append(
        {
            "method": "Sylph baseline",
            "all_F1": float(sylph["F1"]),
            "all_L1_pp": float(sylph["L1"]),
            "mouse_F1": float(sylph["mouse_F1"]),
            "mouse_L1_pp": float(sylph["mouse_L1"]),
            "cami3_F1": float(sylph["cami3_F1"]),
            "cami3_L1_pp": float(sylph["cami3_L1"]),
        }
    )
    for table, label in [
        (top_f1, "MinCO current mixed F1-priority"),
        (top_l1, "MinCO current mixed L1-priority"),
    ]:
        row = table.loc[table["method"] == "raw_value_sum"].iloc[0]
        rows.append(
            {
                "method": label,
                "all_F1": float(row["all_F1"]),
                "all_L1_pp": float(row["all_L1"]),
                "mouse_F1": float(row["mouse_F1"]),
                "mouse_L1_pp": float(row["mouse_L1"]),
                "cami3_F1": float(row["cami3_F1"]),
                "cami3_L1_pp": float(row["cami3_L1"]),
            }
        )
    return pd.DataFrame(rows)


def combined_summary(mouse: pd.DataFrame, cami3: pd.DataFrame, addback: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in sorted(set(mouse["method"]) & set(cami3["method"])):
        m = mouse.loc[mouse["method"] == method].iloc[0]
        c = cami3.loc[cami3["method"] == method].iloc[0]
        rows.append(
            {
                "method": method,
                "all_F1": (float(m["F1"]) + float(c["F1"])) / 2.0,
                "all_L1_pp": (float(m["l1_pct_points"]) + float(c["l1_pct_points"])) / 2.0,
                "mouse_F1": float(m["F1"]),
                "mouse_L1_pp": float(m["l1_pct_points"]),
                "cami3_F1": float(c["F1"]),
                "cami3_L1_pp": float(c["l1_pct_points"]),
            }
        )
    addback_cols = [
        "method",
        "all_F1",
        "all_L1_pp",
        "mouse_F1",
        "mouse_L1_pp",
        "cami3_F1",
        "cami3_L1_pp",
    ]
    addback_best = pd.concat(
        [
            addback.loc[~addback["method"].str.contains("_em_proxy", regex=False)],
            addback.sort_values(["all_F1", "all_L1_pp"], ascending=[False, True])
            .loc[addback["method"].str.contains("_em_proxy", regex=False)]
            .head(3),
            addback.sort_values(["all_L1_pp", "all_F1"], ascending=[True, False])
            .loc[addback["method"].str.contains("_em_proxy", regex=False)]
            .head(3),
        ],
        ignore_index=True,
    )
    addback_best = addback_best.drop_duplicates("method")[addback_cols]
    out = pd.concat([load_cached_current_rows(), pd.DataFrame(rows), addback_best], ignore_index=True)
    out = out.sort_values(["all_F1", "all_L1_pp"], ascending=[False, True])
    out.to_csv(RESULTS / "integrated_f1_abundance_comparison.tsv", sep="\t", index=False)
    return out


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mouse = score_mouse()
    cami3 = score_cami3()
    addback = score_integrated_addback()
    combined = combined_summary(mouse, cami3, addback)
    combined.to_csv(EXP_DIR / "summary.tsv", sep="\t", index=False)
    print(combined.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
