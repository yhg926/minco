#!/usr/bin/env python3
"""Search MinCO gate and abundance thresholds across cached S2000 outputs.

This is an offline scorer. It does not rerun MinCO; it recombines cached
ctx-marker and full-sketch result rows, then scores presence F1, abundance L1,
and available source-to-reference ANI truth.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parent
RESULTS = EXP / "results"

MOUSE_HYBRID_EXP = ROOT / "research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback"
MOUSE_ABUND_EXP = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"
CAMI3_BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
CAMI3_RESCUE_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue"
MOUSE_ANI_PATH = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani/concrete_ani_check_coden11_coden15_vs_sylph.tsv"

sys.path.insert(0, str(MOUSE_HYBRID_EXP))
sys.path.insert(0, str(MOUSE_ABUND_EXP))
sys.path.insert(0, str(CAMI3_BASE_EXP))
sys.path.insert(0, str(CAMI3_RESCUE_EXP))

import score_hybrid_mouse as mouse_hybrid  # noqa: E402
import score_hybrid_cami3 as cami3_hybrid  # noqa: E402
import score_cami3_abundance_rescue as cami3_rescue  # noqa: E402
import score_cami3_toygut as cami3_base  # noqa: E402


PREFILTER_XNY = 5.0
EPS = 1e-12


@dataclass(frozen=True)
class Combo:
    combo_id: str
    ani_min: float
    xny_min: float
    af_mode: str
    af_value: float
    delta_mean_min: float
    delta_vmr_min: float
    delta_max: float
    median_cutoff: float
    abundance_af_exp: float
    fallback_marker_cutoff: int
    fallback_xny_min: float
    fallback_scale: float
    rescue_mode: str


@dataclass
class PreparedSample:
    dataset: str
    sample_id: int
    truth: dict[str, float]
    marker: pd.DataFrame
    full: pd.DataFrame


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)


def ztp_lambda_from_pos_mean(mean_pos: float) -> float:
    if not math.isfinite(mean_pos) or mean_pos <= 0.0:
        return 0.0
    if mean_pos <= 1.0 + 1e-10:
        return 1e-10
    lo = 1e-10
    hi = max(2.0, mean_pos * 2.0)

    def cond_mean(lam: float) -> float:
        return lam / (1.0 - math.exp(-lam))

    while cond_mean(hi) < mean_pos and hi < 1e6:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if cond_mean(mid) < mean_pos:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ztp_adjusted_af(breadth: float, mean_pos: float) -> float:
    if not math.isfinite(breadth) or breadth <= 0.0:
        return 0.0
    lam = ztp_lambda_from_pos_mean(mean_pos)
    p_nonzero = 1.0 - math.exp(-lam) if lam > 0.0 else 0.0
    if p_nonzero <= 0.0:
        return 1.0
    return min(1.0, breadth / p_nonzero)


def af_floor(combo: Combo) -> float:
    if combo.af_mode == "fixed":
        return combo.af_value
    return math.exp((combo.ani_min - 1.0) * combo.af_value)


def formula_ani_from_af(af: pd.Series, ctx_len: float) -> pd.Series:
    return 1.0 + np.log(np.maximum(af.to_numpy(dtype=float), 1e-300)) / ctx_len


def prepare_rows(rows: pd.DataFrame, *, dataset: str, sample_id: int, target_col: str) -> pd.DataFrame:
    rows = rows.copy()
    if "ANI_naive_calc" not in rows.columns:
        try:
            rows = mouse_hybrid.truth.add_naive_ani(rows)
        except Exception:
            rows["ANI_naive_calc"] = numeric(rows, "ANI")
    rows["dataset"] = dataset
    rows["sample_id"] = sample_id
    rows["target_id"] = rows[target_col].fillna("").astype(str) if target_col in rows.columns else ""

    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_mean_depth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_median_depth",
        "Reliable_Ref_hit_depth_variance",
        "Reliable_Ref_zip_af",
        "Normalized_abundance_depth",
        "ctx_marker_size",
    ]:
        rows[col] = numeric(rows, col)

    keep = (rows["XnY_ctx"] >= PREFILTER_XNY) & rows["target_id"].astype(bool)
    rows = rows.loc[keep].copy()
    rows["Reliable_ztp_af"] = [
        ztp_adjusted_af(float(b), float(m))
        for b, m in zip(rows["Reliable_Ref_breadth"], rows["Reliable_Ref_hit_mean_depth"])
    ]
    rows["Reliable_depth_vmr"] = np.divide(
        rows["Reliable_Ref_hit_depth_variance"].to_numpy(dtype=float),
        rows["Reliable_Ref_hit_mean_depth"].to_numpy(dtype=float),
        out=np.zeros(len(rows), dtype=float),
        where=rows["Reliable_Ref_hit_mean_depth"].to_numpy(dtype=float) > 0.0,
    )
    if "species_genus" not in rows.columns:
        rows["species_genus"] = rows.get("gtdb_species", "").astype(str).map(species_genus)
    return rows.reset_index(drop=True)


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def load_mouse_samples() -> list[PreparedSample]:
    ctx_marker_sizes = mouse_hybrid.load_psmp_sizes(mouse_hybrid.CTX_MARKER_PSM_PATH)
    by_accession, by_core = mouse_hybrid.truth.load_gtdb_metadata(mouse_hybrid.truth.GTDB_METADATA)
    samples: list[PreparedSample] = []
    for sample_id in mouse_hybrid.SAMPLES:
        marker, full = mouse_hybrid.load_sample_views(sample_id, by_accession, by_core, ctx_marker_sizes)
        samples.append(
            PreparedSample(
                dataset="mouse_gtdb",
                sample_id=sample_id,
                truth=mouse_hybrid.ev.load_truth_profile(sample_id),
                marker=prepare_rows(marker, dataset="mouse_gtdb", sample_id=sample_id, target_col="gtdb_species"),
                full=prepare_rows(full, dataset="mouse_gtdb", sample_id=sample_id, target_col="gtdb_species"),
            )
        )
    return samples


def load_mouse_ctxobj_samples() -> list[PreparedSample]:
    by_accession, by_core = mouse_hybrid.truth.load_gtdb_metadata(mouse_hybrid.truth.GTDB_METADATA)
    samples: list[PreparedSample] = []
    for sample_id in mouse_hybrid.SAMPLES:
        path = MOUSE_ABUND_EXP / f"toymouse_sample{sample_id}_ctxobjmarker_median_col.tsv"
        rows = mouse_hybrid.load_minco_active_rows_fast(path, by_accession, by_core)
        rows = mouse_hybrid.attach_mode(rows, "ctxobj_marker", {})
        samples.append(
            PreparedSample(
                dataset="mouse_gtdb",
                sample_id=sample_id,
                truth=mouse_hybrid.ev.load_truth_profile(sample_id),
                marker=prepare_rows(rows, dataset="mouse_gtdb", sample_id=sample_id, target_col="gtdb_species"),
                full=pd.DataFrame(),
            )
        )
    return samples


def load_cami3_samples() -> list[PreparedSample]:
    ctx_marker_sizes = cami3_hybrid.load_psmp_sizes(cami3_hybrid.CTX_MARKER_PSM_PATH)
    by_accession, by_core = cami3_base.truth.load_gtdb_metadata(cami3_base.truth.GTDB_METADATA)
    samples: list[PreparedSample] = []
    for sample_id, paths in cami3_base.SAMPLES.items():
        truth_rows = cami3_base.load_truth(paths["truth"])
        truth = dict(
            zip(
                truth_rows["ncbi_species_taxid"].astype(str),
                truth_rows["truth_abundance_bacterial_norm"].astype(float),
            )
        )
        marker = cami3_hybrid.add_hybrid_columns(
            cami3_base.load_minco(cami3_hybrid.marker_path(sample_id), by_accession, by_core),
            "ctx_marker",
            ctx_marker_sizes,
        )
        full = cami3_hybrid.add_hybrid_columns(
            cami3_base.load_minco(cami3_hybrid.full_path(sample_id), by_accession, by_core),
            "full_fallback",
            ctx_marker_sizes,
        )
        samples.append(
            PreparedSample(
                dataset="cami3_ncbi",
                sample_id=sample_id,
                truth=truth,
                marker=prepare_rows(marker, dataset="cami3_ncbi", sample_id=sample_id, target_col="ncbi_species_taxid"),
                full=prepare_rows(full, dataset="cami3_ncbi", sample_id=sample_id, target_col="ncbi_species_taxid"),
            )
        )
    return samples


def load_cami3_ctxobj_samples() -> list[PreparedSample]:
    by_accession, by_core = cami3_base.truth.load_gtdb_metadata(cami3_base.truth.GTDB_METADATA)
    samples: list[PreparedSample] = []
    for sample_id, paths in cami3_base.SAMPLES.items():
        truth_rows = cami3_base.load_truth(paths["truth"])
        truth = dict(
            zip(
                truth_rows["ncbi_species_taxid"].astype(str),
                truth_rows["truth_abundance_bacterial_norm"].astype(float),
            )
        )
        rows = cami3_base.load_minco(cami3_rescue.MINCO_OUTPUTS["ctxobj"][sample_id], by_accession, by_core)
        rows = cami3_rescue.add_effective_depth(rows)
        rows["source_mode"] = "ctxobj_marker"
        rows["ctx_marker_size"] = 999999
        samples.append(
            PreparedSample(
                dataset="cami3_ncbi",
                sample_id=sample_id,
                truth=truth,
                marker=prepare_rows(rows, dataset="cami3_ncbi", sample_id=sample_id, target_col="ncbi_species_taxid"),
                full=pd.DataFrame(),
            )
        )
    return samples


def effective_depth(rows: pd.DataFrame, combo: Combo) -> pd.Series:
    median = numeric(rows, "Reliable_Ref_hit_median_depth")
    mean = numeric(rows, "Reliable_Ref_mean_depth")
    zip_af = numeric(rows, "Reliable_Ref_zip_af").clip(lower=EPS)
    out = mean / (zip_af**combo.abundance_af_exp)
    out = out.where(median < combo.median_cutoff, median)
    return out.where(mean > 0.0, 0.0)


def active_mask(rows: pd.DataFrame, combo: Combo) -> pd.Series:
    reliable_ztp = numeric(rows, "Reliable_ztp_af")
    ani_from_af = formula_ani_from_af(reliable_ztp, combo.af_value if combo.af_mode == "formula" else 24.0)
    ani_delta = numeric(rows, "ANI_naive_calc") - ani_from_af
    delta_trigger = (
        (numeric(rows, "Reliable_Ref_hit_mean_depth") > combo.delta_mean_min)
        & (numeric(rows, "Reliable_depth_vmr") > combo.delta_vmr_min)
    )
    delta_pass = (~delta_trigger) | (ani_delta < combo.delta_max)
    return (
        (numeric(rows, "XnY_ctx") >= combo.xny_min)
        & (numeric(rows, "ANI_naive_calc") > combo.ani_min)
        & (reliable_ztp >= af_floor(combo))
        & delta_pass
        & rows["target_id"].astype(bool)
    )


def collect_values(rows: pd.DataFrame, mask: pd.Series, depths: pd.Series, scale: float) -> dict[str, float]:
    if rows.empty:
        return {}
    selected = rows.loc[mask, ["target_id"]].copy()
    if selected.empty:
        return {}
    selected["value"] = depths.loc[selected.index].to_numpy(dtype=float) * scale
    selected = selected.loc[selected["value"] > 0.0]
    if selected.empty:
        return {}
    return selected.groupby("target_id")["value"].max().to_dict()


def merge_values(dst: dict[str, float], src: dict[str, float]) -> None:
    for key, value in src.items():
        if value > dst.get(key, 0.0):
            dst[key] = float(value)


def rescue_values(
    rows: pd.DataFrame,
    active: pd.Series,
    depths: pd.Series,
    values: dict[str, float],
    combo: Combo,
) -> list[str]:
    if combo.rescue_mode == "off" or rows.empty:
        return []
    active_rows = rows.loc[active].copy()
    if active_rows.empty:
        return []
    active_rows["value"] = depths.loc[active_rows.index].to_numpy(dtype=float)
    active_genus_max = active_rows.loc[active_rows["value"] > 0.0].groupby("species_genus")["value"].max()
    if active_genus_max.empty:
        return []

    if combo.rescue_mode == "strict":
        rescue_ani_min, rescue_xny_min, rescue_effective_min, rescue_zip_max, rescue_ratio = 0.999, 100.0, 5.0, 0.25, 3.0
    elif combo.rescue_mode == "moderate":
        rescue_ani_min, rescue_xny_min, rescue_effective_min, rescue_zip_max, rescue_ratio = 0.997, 75.0, 4.0, 0.30, 2.5
    else:
        rescue_ani_min, rescue_xny_min, rescue_effective_min, rescue_zip_max, rescue_ratio = 0.995, 50.0, 3.0, 0.35, 2.0

    candidates = rows.loc[
        (~active)
        & (numeric(rows, "ANI_naive_calc") >= rescue_ani_min)
        & (numeric(rows, "XnY_ctx") >= rescue_xny_min)
        & (numeric(rows, "Reliable_Ref_zip_af") <= rescue_zip_max)
    ].copy()
    if candidates.empty:
        return []
    candidates["value"] = depths.loc[candidates.index].to_numpy(dtype=float)
    candidates["active_genus_value"] = candidates["species_genus"].map(active_genus_max).fillna(0.0)
    candidates = candidates.loc[
        (candidates["value"] >= rescue_effective_min)
        & (candidates["active_genus_value"] > 0.0)
        & (candidates["value"] >= rescue_ratio * candidates["active_genus_value"])
    ]
    if candidates.empty:
        return []
    winners = candidates.sort_values(["species_genus", "value", "XnY_ctx"]).groupby("species_genus", as_index=False).tail(1)
    rescued: list[str] = []
    for row in winners.itertuples(index=False):
        target = str(row.target_id)
        value = float(row.value)
        if value > values.get(target, 0.0):
            values[target] = value
            rescued.append(target)
    return rescued


def values_for_sample(sample: PreparedSample, combo: Combo) -> tuple[dict[str, float], list[str]]:
    marker_depth = effective_depth(sample.marker, combo)
    marker_active = active_mask(sample.marker, combo)
    values = collect_values(sample.marker, marker_active, marker_depth, 1.0)
    rescued = rescue_values(sample.marker, marker_active, marker_depth, values, combo)

    if combo.fallback_marker_cutoff > 0 and not sample.full.empty:
        full_base = (
            (sample.full["ctx_marker_size"] < combo.fallback_marker_cutoff)
            & (sample.full["XnY_ctx"] >= combo.fallback_xny_min)
        )
        if full_base.any():
            full_rows = sample.full.loc[full_base].copy()
            full_depth = effective_depth(full_rows, combo)
            full_active = active_mask(full_rows, combo)
            full_values = collect_values(full_rows, full_active, full_depth, combo.fallback_scale)
            for target in list(full_values):
                if target in values:
                    del full_values[target]
            merge_values(values, full_values)
    return values, rescued


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0}


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 2:
        return float("nan")
    xs = pd.Series(x, dtype=float)
    ys = pd.Series(y, dtype=float)
    return float(xs.corr(ys, method="pearson"))


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) < 2:
        return float("nan")
    xs = pd.Series(x, dtype=float)
    ys = pd.Series(y, dtype=float)
    return float(xs.corr(ys, method="spearman"))


def score_sample(sample: PreparedSample, combo: Combo) -> tuple[dict[str, object], dict[str, float]]:
    values, rescued = values_for_sample(sample, combo)
    pred = normalize(values)
    truth = sample.truth
    truth_set = set(truth)
    pred_set = set(pred)
    tp = pred_set & truth_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set
    precision = len(tp) / len(pred_set) if pred_set else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    keys = sorted(truth)
    y_true = [truth[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    sample_row = {
        "combo_id": combo.combo_id,
        "dataset": sample.dataset,
        "sample_id": sample.sample_id,
        "truth_taxa": len(truth_set),
        "pred_taxa": len(pred_set),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "missing_truth_mass": sum(truth[k] for k in fn),
        "fp_pred_mass": sum(pred.get(k, 0.0) for k in fp),
        "pred_sum_on_truth": sum(y_pred),
        "pred_sum_all": sum(pred.values()),
        "pearson": pearson(y_pred, y_true),
        "spearman": spearman(y_pred, y_true),
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) / max(len(keys), 1) * 100.0,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) * 100.0,
        "rescued_count": len(rescued),
    }
    return sample_row, pred


def combo_to_record(combo: Combo) -> dict[str, object]:
    rec = combo.__dict__.copy()
    rec["af_floor"] = af_floor(combo)
    return rec


def build_combos(max_combos: int | None = None, preset: str = "broad") -> list[Combo]:
    if preset == "focused":
        ani_values = [0.945, 0.95, 0.955, 0.96, 0.965, 0.97]
        xny_values = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
        af_defs = [
            ("fixed", 0.25),
            ("fixed", 0.28),
            ("fixed", 0.30),
            ("fixed", 0.32),
            ("fixed", 0.35),
            ("formula", 22.0),
            ("formula", 24.0),
            ("formula", 26.0),
            ("formula", 27.0),
            ("formula", 30.0),
        ]
        delta_mean_values = [2.0, 3.0, 4.0, 5.0]
        delta_vmr_values = [10.0, 20.0, 50.0, 100.0]
        delta_max_values = [0.015, 0.02, 0.025, 0.03, 0.04]
        depth_defs = [(10.0, 1.0), (10.0, 1.05), (20.0, 1.0), (20.0, 1.05), (30.0, 1.05)]
        fallback_defs = [
            (0, 0.0, 1.0),
            (100, 700.0, 0.4),
            (200, 500.0, 0.4),
            (200, 700.0, 0.4),
            (300, 700.0, 0.4),
            (500, 1000.0, 0.4),
        ]
        rescue_modes = ["strict", "off", "moderate"]
    else:
        ani_values = [0.945, 0.95, 0.955, 0.96, 0.965]
        xny_values = [5.0, 10.0, 15.0, 20.0, 30.0, 50.0]
        af_defs = [
            ("fixed", 0.30),
            ("fixed", 0.35),
            ("fixed", 0.40),
            ("fixed", 0.45),
            ("formula", 22.0),
            ("formula", 24.0),
            ("formula", 27.0),
            ("formula", 32.0),
        ]
        delta_mean_values = [2.0, 3.0, 4.0]
        delta_vmr_values = [20.0, 50.0, 100.0]
        delta_max_values = [0.02, 0.03, 0.04, 0.05]
        depth_defs = [(10.0, 1.0), (10.0, 1.05), (20.0, 1.0), (20.0, 1.05)]
        fallback_defs = [
            (0, 0.0, 1.0),
            (100, 700.0, 0.4),
            (200, 500.0, 0.4),
            (200, 700.0, 0.4),
            (300, 700.0, 0.4),
            (500, 1000.0, 0.4),
        ]
        rescue_modes = ["strict", "off", "moderate"]

    combos: list[Combo] = []
    for ani in ani_values:
        for xny in xny_values:
            for af_mode, af_value in af_defs:
                for delta_mean in delta_mean_values:
                    for delta_vmr in delta_vmr_values:
                        for delta_max in delta_max_values:
                            for median_cutoff, abundance_af_exp in depth_defs:
                                for fallback_cutoff, fallback_xny, fallback_scale in fallback_defs:
                                    for rescue_mode in rescue_modes:
                                        raw = (
                                            ani,
                                            xny,
                                            af_mode,
                                            af_value,
                                            delta_mean,
                                            delta_vmr,
                                            delta_max,
                                            median_cutoff,
                                            abundance_af_exp,
                                            fallback_cutoff,
                                            fallback_xny,
                                            fallback_scale,
                                            rescue_mode,
                                        )
                                        digest = hashlib.sha1(repr(raw).encode()).hexdigest()[:12]
                                        combos.append(
                                            Combo(
                                                combo_id=f"combo_{digest}",
                                                ani_min=ani,
                                                xny_min=xny,
                                                af_mode=af_mode,
                                                af_value=af_value,
                                                delta_mean_min=delta_mean,
                                                delta_vmr_min=delta_vmr,
                                                delta_max=delta_max,
                                                median_cutoff=median_cutoff,
                                                abundance_af_exp=abundance_af_exp,
                                                fallback_marker_cutoff=fallback_cutoff,
                                                fallback_xny_min=fallback_xny,
                                                fallback_scale=fallback_scale,
                                                rescue_mode=rescue_mode,
                                            )
                                        )
    if max_combos is not None and max_combos > 0 and len(combos) > max_combos:
        # Deterministic downsampling preserves all known-default combinations below.
        combos = sorted(combos, key=lambda c: c.combo_id)[:max_combos]

    known = [
        Combo("known_current_ctx_marker", 0.95, 15.0, "fixed", 0.40, 3.0, 50.0, 0.03, 20.0, 1.05, 0, 0.0, 1.0, "strict"),
        Combo("known_hybrid_xny700", 0.95, 15.0, "fixed", 0.40, 3.0, 50.0, 0.03, 20.0, 1.05, 200, 700.0, 0.4, "strict"),
        Combo("known_formula_len24", 0.95, 10.0, "formula", 24.0, 3.0, 50.0, 0.03, 20.0, 1.05, 200, 700.0, 0.4, "strict"),
    ]
    by_id = {c.combo_id: c for c in combos}
    for combo in known:
        by_id[combo.combo_id] = combo
    return list(by_id.values())


def summarize_combo(combo: Combo, sample_rows: list[dict[str, object]], mouse_ani: pd.DataFrame) -> dict[str, object]:
    df = pd.DataFrame(sample_rows)
    rec = combo_to_record(combo)
    for prefix, sub in [("all", df), ("mouse", df.loc[df["dataset"] == "mouse_gtdb"]), ("cami3", df.loc[df["dataset"] == "cami3_ncbi"])]:
        rec[f"{prefix}_F1"] = float(sub["F1"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_L1"] = float(sub["l1_pct_points"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_Pearson"] = float(sub["pearson"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_TP"] = float(sub["TP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FP"] = float(sub["FP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FN"] = float(sub["FN"].mean()) if len(sub) else float("nan")
    rec.update(mouse_ani_for_combo(mouse_ani, combo.combo_id))
    return rec


def mouse_ani_for_combo(mouse_ani: pd.DataFrame, combo_id: str) -> dict[str, object]:
    if mouse_ani.empty:
        return {"mouse0_ani_n": 0, "mouse0_ani_mae": float("nan"), "mouse0_ani_pearson": float("nan")}
    sub = mouse_ani.loc[mouse_ani["combo_id"] == combo_id].copy()
    if sub.empty:
        return {"mouse0_ani_n": 0, "mouse0_ani_mae": float("nan"), "mouse0_ani_pearson": float("nan")}
    err = sub["pred_ani"].astype(float) - sub["ANIm_ANI"].astype(float)
    return {
        "mouse0_ani_n": int(len(sub)),
        "mouse0_ani_mae": float(err.abs().mean()),
        "mouse0_ani_bias": float(err.mean()),
        "mouse0_ani_pearson": float(sub["pred_ani"].astype(float).corr(sub["ANIm_ANI"].astype(float), method="pearson")),
    }


def load_mouse_ani_truth() -> pd.DataFrame:
    if not MOUSE_ANI_PATH.exists():
        return pd.DataFrame()
    rows = pd.read_csv(MOUSE_ANI_PATH, sep="\t")
    rows = rows[["gtdb_species", "ANIm_ANI", "minco_naive_ANI_coden11", "minco_aafANI_coden11"]].copy()
    rows["gtdb_species"] = rows["gtdb_species"].astype(str)
    for col in ["ANIm_ANI", "minco_naive_ANI_coden11", "minco_aafANI_coden11"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    return rows.dropna(subset=["ANIm_ANI"])


def append_mouse_ani_rows(
    out_rows: list[dict[str, object]],
    combo_id: str,
    sample: PreparedSample,
    pred: dict[str, float],
    ani_truth: pd.DataFrame,
    estimator: str,
) -> None:
    if sample.dataset != "mouse_gtdb" or sample.sample_id != 0 or ani_truth.empty:
        return
    pred_species = pd.DataFrame({"gtdb_species": list(pred)})
    if pred_species.empty:
        return
    joined = pred_species.merge(ani_truth, on="gtdb_species", how="inner")
    pred_col = "minco_aafANI_coden11" if estimator == "aaf" else "minco_naive_ANI_coden11"
    for row in joined.itertuples(index=False):
        value = getattr(row, pred_col)
        if not math.isfinite(float(value)):
            continue
        out_rows.append(
            {
                "combo_id": combo_id,
                "gtdb_species": row.gtdb_species,
                "ANIm_ANI": float(row.ANIm_ANI),
                "pred_ani": float(value),
                "ani_estimator": estimator,
            }
        )


def score_sylph_baselines(samples: list[PreparedSample]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    by_accession, by_core = mouse_hybrid.truth.load_gtdb_metadata(mouse_hybrid.truth.GTDB_METADATA)
    for sample_id in mouse_hybrid.SAMPLES:
        truth = mouse_hybrid.ev.load_truth_profile(sample_id)
        sylph = mouse_hybrid.more.load_sylph_rows(mouse_hybrid.more.sylph_profile_path(sample_id), by_accession, by_core)
        values = mouse_hybrid.rescue.sylph_values(sylph)
        pseudo = PreparedSample("mouse_gtdb", sample_id, truth, pd.DataFrame(), pd.DataFrame())
        pred = normalize(values)
        pred_set = set(pred)
        truth_set = set(truth)
        tp = pred_set & truth_set
        fp = pred_set - truth_set
        fn = truth_set - pred_set
        precision = len(tp) / len(pred_set) if pred_set else 0.0
        recall = len(tp) / len(truth_set) if truth_set else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        keys = sorted(truth)
        rows.append(
            {
                "method": "sylph",
                "dataset": pseudo.dataset,
                "sample_id": sample_id,
                "F1": f1,
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "l1_pct_points": sum(abs(pred.get(k, 0.0) - truth[k]) for k in keys) * 100.0,
            }
        )

    by_accession, by_core = cami3_base.truth.load_gtdb_metadata(cami3_base.truth.GTDB_METADATA)
    for sample_id, paths in cami3_base.SAMPLES.items():
        truth_rows = cami3_base.load_truth(paths["truth"])
        truth = dict(zip(truth_rows["ncbi_species_taxid"].astype(str), truth_rows["truth_abundance_bacterial_norm"].astype(float)))
        sylph = cami3_base.load_sylph(paths["sylph"], by_accession, by_core)
        selected = cami3_base.select_predictions("sylph", sylph)
        collapsed = cami3_base.collapse_predictions(selected)
        pred_raw = dict(zip(collapsed["ncbi_species_taxid"].astype(str), collapsed["pred_abundance_raw"].astype(float)))
        pred = normalize(pred_raw)
        pred_set = set(pred)
        truth_set = set(truth)
        tp = pred_set & truth_set
        fp = pred_set - truth_set
        fn = truth_set - pred_set
        precision = len(tp) / len(pred_set) if pred_set else 0.0
        recall = len(tp) / len(truth_set) if truth_set else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        keys = sorted(truth)
        rows.append(
            {
                "method": "sylph",
                "dataset": "cami3_ncbi",
                "sample_id": sample_id,
                "F1": f1,
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "l1_pct_points": sum(abs(pred.get(k, 0.0) - truth[k]) for k in keys) * 100.0,
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    global RESULTS
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-combos", type=int, default=60000)
    parser.add_argument("--preset", choices=["broad", "focused"], default="broad")
    parser.add_argument("--ref-mode", choices=["ctxmarker_full", "ctxobj"], default="ctxmarker_full")
    parser.add_argument("--ani-estimator", choices=["naive", "aaf"], default="aaf")
    args = parser.parse_args(argv)

    RESULTS = RESULTS / args.ref_mode
    RESULTS.mkdir(parents=True, exist_ok=True)
    if args.ref_mode == "ctxobj":
        samples = load_mouse_ctxobj_samples() + load_cami3_ctxobj_samples()
    else:
        samples = load_mouse_samples() + load_cami3_samples()
    combos = build_combos(args.max_combos, args.preset)
    ani_truth = load_mouse_ani_truth()

    summary_rows: list[dict[str, object]] = []
    sample_metric_rows: list[dict[str, object]] = []
    saved_ani_rows: list[dict[str, object]] = []

    for i, combo in enumerate(combos, start=1):
        rows_for_combo: list[dict[str, object]] = []
        ani_rows_for_combo: list[dict[str, object]] = []
        for sample in samples:
            sample_row, pred = score_sample(sample, combo)
            rows_for_combo.append(sample_row)
            sample_metric_rows.append(sample_row)
            append_mouse_ani_rows(ani_rows_for_combo, combo.combo_id, sample, pred, ani_truth, args.ani_estimator)
        if combo.combo_id in {"known_current_ctx_marker", "known_hybrid_xny700"} or i <= 100:
            saved_ani_rows.extend(ani_rows_for_combo)
        summary_rows.append(summarize_combo(combo, rows_for_combo, pd.DataFrame(ani_rows_for_combo)))
        if i % 1000 == 0:
            print(f"scored {i}/{len(combos)} combos", file=sys.stderr)

    summary = pd.DataFrame(summary_rows)
    sort_cols = ["all_F1", "all_L1", "mouse0_ani_mae", "cami3_F1", "mouse_F1"]
    summary = summary.sort_values(sort_cols, ascending=[False, True, True, False, False])
    summary.to_csv(RESULTS / "combo_search_summary.tsv", sep="\t", index=False)
    summary.head(200).to_csv(RESULTS / "combo_search_top200.tsv", sep="\t", index=False)
    pd.DataFrame(sample_metric_rows).to_csv(RESULTS / "combo_search_sample_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(saved_ani_rows).to_csv(RESULTS / "combo_search_mouse0_ani_matched_subset.tsv", sep="\t", index=False)

    sylph_samples = score_sylph_baselines(samples)
    sylph_mean = (
        sylph_samples.groupby(["method", "dataset"], as_index=False)[["F1", "TP", "FP", "FN", "l1_pct_points"]]
        .mean()
    )
    sylph_overall = (
        sylph_samples.groupby("method", as_index=False)[["F1", "TP", "FP", "FN", "l1_pct_points"]]
        .mean()
        .assign(dataset="all")
    )
    pd.concat([sylph_mean, sylph_overall], ignore_index=True).to_csv(RESULTS / "sylph_baseline.tsv", sep="\t", index=False)

    print(summary.head(20).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
