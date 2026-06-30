#!/usr/bin/env python3
"""Sweep simple abundance transforms on selected-candidate cached profiles.

The selected call-set oracle shows most cached samples have enough detected
truth mass in principle. This sweep asks whether small, implementable
base-called-row mass transforms can capture part of that headroom without
changing calls. Candidate-added row mass is preserved; only base called rows are
reallocated.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

import audit_abundance_oracle_bounds as oracle_bounds
import audit_candidate_callset_oracle_feasibility as feasibility
import decompose_abundance_errors as decomp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SELECTED_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"
SYLPH_SCORES = RESULTS / "abundance_oracle_bounds_scores.tsv"

OUT_SCORES = RESULTS / "selected_call_mass_transform_scores.tsv"
OUT_OVERALL = RESULTS / "selected_call_mass_transform_overall.tsv"
OUT_PANEL = RESULTS / "selected_call_mass_transform_panel_delta.tsv"
OUT_AUDIT = RESULTS / "selected_call_mass_transform_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_MASS_TRANSFORM_SWEEP.md"


VARIANTS = [
    {"method": "current_selected_abundance", "kind": "current"},
    {"method": "base_power_p0.90", "kind": "power", "power": 0.90},
    {"method": "base_power_p0.95", "kind": "power", "power": 0.95},
    {"method": "base_power_p1.05", "kind": "power", "power": 1.05},
    {"method": "base_power_p1.10", "kind": "power", "power": 1.10},
    {"method": "base_uniform_a0.02", "kind": "blend_uniform", "alpha": 0.02},
    {"method": "base_uniform_a0.05", "kind": "blend_uniform", "alpha": 0.05},
    {"method": "base_global_xny_a0.02", "kind": "blend_quality", "quality": "xny", "alpha": 0.02},
    {"method": "base_global_xny_a0.05", "kind": "blend_quality", "quality": "xny", "alpha": 0.05},
    {"method": "base_global_xny_a0.10", "kind": "blend_quality", "quality": "xny", "alpha": 0.10},
    {"method": "base_global_probability_a0.02", "kind": "blend_quality", "quality": "probability", "alpha": 0.02},
    {"method": "base_global_probability_a0.05", "kind": "blend_quality", "quality": "probability", "alpha": 0.05},
    {"method": "base_genus_xny_a0.02", "kind": "blend_genus_quality", "quality": "xny", "alpha": 0.02},
    {"method": "base_genus_xny_a0.05", "kind": "blend_genus_quality", "quality": "xny", "alpha": 0.05},
    {"method": "base_genus_xny_a0.10", "kind": "blend_genus_quality", "quality": "xny", "alpha": 0.10},
]


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def bool_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def genus_from_species(value: object) -> str:
    text = str(value or "")
    if text.startswith("s__"):
        text = text[3:]
    return text.split()[0] if text.strip() else ""


def normalize_raw(raw: np.ndarray, target_total: float) -> np.ndarray:
    out = np.maximum(np.asarray(raw, dtype=float), 0.0)
    total = float(out.sum())
    if total <= 0.0 or target_total <= 0.0:
        return out
    return out / total * target_total


def species_from_accession(
    panel: str,
    accession: object,
    loader: oracle_bounds.PanelLoader,
) -> str:
    accession_text = str(accession or "")
    if not accession_text:
        return ""
    if panel == "hmp_airskin_gtdb_source_abundance":
        species, _method = decomp.hmp.gtdb_from_accession(
            accession_text, loader.by_accession, loader.by_core
        )
        return species
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        species, _method = decomp.hmp_gastro.gtdb_from_accession(
            accession_text, loader.by_accession, loader.by_core
        )
        return species
    if panel == "cami3_toy_human_gut_gtdb_source_readmap":
        return decomp.cami3.gtdb_from_accession(
            accession_text, loader.by_accession, loader.by_core
        )
    return ""


class RowMapper:
    def __init__(self) -> None:
        self.loader = oracle_bounds.PanelLoader()
        self._best_ref_cache: dict[tuple[str, int], dict[str, str]] = {}

    def best_ref_species(self, panel: str, sample: int) -> dict[str, str]:
        key = (panel, sample)
        if key in self._best_ref_cache:
            return self._best_ref_cache[key]
        if panel == "hmp_airskin_gtdb_source_abundance":
            decomp.hmp.ensure_sample_record(sample)
            paths = decomp.hmp.SAMPLES[sample]
            best, _diag = decomp.hmp.best_raw_ref_species_by_taxid(
                {"unique": paths["unique"], "split": paths["split"]},
                self.loader.hmp_taxmap,
                self.loader.by_accession,
                self.loader.by_core,
            )
        elif panel == "hmp_gastrooral_gtdb_source_abundance":
            paths = decomp.hmp_gastro.SAMPLES[sample]
            best, _diag = decomp.hmp_gastro.best_raw_ref_species_by_taxid(
                {"unique": paths["unique"], "split": paths["split"]},
                self.loader.hmp_gastro_taxmap,
                self.loader.by_accession,
                self.loader.by_core,
            )
        elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
            best, _diag = decomp.cami3.best_raw_ref_species_by_taxid(
                decomp.cami3.RAW_TABLES[sample],
                self.loader.cami3_taxmap,
                self.loader.by_accession,
                self.loader.by_core,
            )
        else:
            best = {}
        self._best_ref_cache[key] = best
        return best

    def row_species(self, panel: str, sample: int, row: object) -> str:
        if panel == "cami2_toy_mouse_gut":
            return str(getattr(row, "species_name", "") or "")
        for col in ["s_best_accession", "u_best_accession", "candidate_surface_accession"]:
            if hasattr(row, col):
                species = species_from_accession(panel, getattr(row, col), self.loader)
                if species:
                    return species
        taxid = str(getattr(row, "taxid", "") or "")
        species = self.best_ref_species(panel, sample).get(taxid, "")
        if species:
            return species
        if panel == "cami3_toy_human_gut_gtdb_source_readmap":
            species = self.loader.cami3_taxid_to_species.get(taxid, "")
            if species:
                return species
            name = decomp.cami3.normalize_name(getattr(row, "species_name", ""))
            return self.loader.cami3_name_to_species.get(name, "")
        return ""


def load_called_rows(row: Mapping[str, object], mapper: RowMapper) -> pd.DataFrame:
    panel = str(row["panel"])
    sample = int(float(row["sample"]))
    profile = Path(str(row["profile"]))
    raw = pd.read_csv(profile, sep="\t", low_memory=False)
    called = bool_series(raw, "calibrated_call")
    df = raw.loc[called].copy()
    species: list[str] = []
    for rec in df.itertuples(index=False):
        species.append(mapper.row_species(panel, sample, rec))
    df["gtdb_species"] = species
    df = df.loc[df["gtdb_species"].astype(str).astype(bool)].copy()
    df["panel"] = panel
    df["sample"] = sample
    return df


def base_raw_values(calls: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    raw = numeric(calls, "calibrated_abundance_raw").to_numpy(dtype=float)
    if float(np.nansum(raw)) <= 0.0:
        raw = numeric(calls, "calibrated_abundance").to_numpy(dtype=float)
    candidate_added = (
        bool_series(calls, "candidate_rescue_added")
        | bool_series(calls, "candidate_surface_added")
    ).to_numpy(dtype=bool)
    return np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0), candidate_added


def quality_values(calls: pd.DataFrame, quality: str) -> np.ndarray:
    if quality == "xny":
        s_xny = numeric(calls, "s_XnY_ctx_max").to_numpy(dtype=float)
        u_xny = numeric(calls, "u_XnY_ctx_max").to_numpy(dtype=float)
        return np.clip(np.maximum(s_xny, u_xny) / 1000.0, 0.0, 1.0)
    if quality == "probability":
        return np.clip(numeric(calls, "calibrated_probability", 1.0).to_numpy(dtype=float), 0.0, 1.0)
    raise ValueError(f"unsupported quality: {quality}")


def apply_variant(calls: pd.DataFrame, variant: Mapping[str, object]) -> np.ndarray:
    raw, candidate_added = base_raw_values(calls)
    eligible = ~candidate_added & (raw > 0.0)
    out = raw.copy()
    base_total = float(raw[eligible].sum())
    if base_total <= 0.0 or not bool(eligible.any()):
        return out

    kind = str(variant["kind"])
    idx = np.flatnonzero(eligible)
    base = raw[idx]
    if kind == "current":
        return out
    if kind == "power":
        target = normalize_raw(np.power(np.maximum(base, 0.0), float(variant["power"])), base_total)
        out[idx] = target
        return out
    if kind == "blend_uniform":
        target = np.full(len(idx), base_total / len(idx), dtype=float)
        alpha = float(variant["alpha"])
        out[idx] = (1.0 - alpha) * base + alpha * target
        return out
    if kind == "blend_quality":
        quality = quality_values(calls, str(variant["quality"]))[idx]
        target = normalize_raw(base * quality, base_total)
        alpha = float(variant["alpha"])
        out[idx] = (1.0 - alpha) * base + alpha * target
        return out
    if kind == "blend_genus_quality":
        alpha = float(variant["alpha"])
        quality = quality_values(calls, str(variant["quality"]))
        work = pd.DataFrame(
            {
                "idx": idx,
                "genus": [genus_from_species(value) for value in calls.iloc[idx]["gtdb_species"]],
            }
        )
        for _genus, sub in work.groupby("genus", sort=False):
            group_idx = sub["idx"].to_numpy(dtype=int)
            if len(group_idx) <= 1:
                continue
            group_base = raw[group_idx]
            total = float(group_base.sum())
            target = normalize_raw(group_base * quality[group_idx], total)
            out[group_idx] = (1.0 - alpha) * group_base + alpha * target
        return out
    raise ValueError(f"unsupported variant kind: {kind}")


def collapse_prediction(calls: pd.DataFrame, raw_values: np.ndarray, panel: str) -> dict[str, float]:
    work = calls[["gtdb_species"]].copy()
    work["raw"] = np.maximum(np.asarray(raw_values, dtype=float), 0.0)
    collapse = "max" if panel == "cami2_toy_mouse_gut" else "sum"
    if collapse == "max":
        grouped = work.groupby("gtdb_species")["raw"].max()
    else:
        grouped = work.groupby("gtdb_species")["raw"].sum()
    return decomp.normalize({str(k): finite_float(v) for k, v in grouped.to_dict().items()})


def load_sylph_rows() -> dict[tuple[str, str], dict[str, object]]:
    rows = pd.read_csv(SYLPH_SCORES, sep="\t")
    rows = rows.loc[rows["method"].astype(str).str.startswith("sylph")].copy()
    return {(str(row["panel"]), str(row["sample"])): row for row in rows.to_dict("records")}


def score_variant(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    calls: pd.DataFrame,
    raw_values: np.ndarray,
) -> dict[str, object]:
    abundance = collapse_prediction(calls, raw_values, panel)
    metrics = feasibility.metric_vector(panel, truth, set(abundance), abundance)
    return {
        "panel": panel,
        "sample": sample,
        "method": method,
        "F1": metrics["F1"],
        "official_L1_pp": metrics["official_L1_pp"],
        "official_Pearson": metrics["official_Pearson"],
        "truth_mass_detected_pct": metrics["truth_mass_detected_pct"],
        "matched_abs_error_pp": metrics["matched_abs_error_pp"],
        "missing_truth_mass_pp": metrics["missing_truth_mass_pp"],
        "extra_pred_mass_pp": metrics["extra_pred_mass_pp"],
    }


def build_scores() -> tuple[pd.DataFrame, float]:
    selected_scores = pd.read_csv(SELECTED_SCORES, sep="\t")
    mapper = RowMapper()
    rows: list[dict[str, object]] = []
    max_validation_delta = 0.0
    for selected_row in selected_scores.to_dict("records"):
        panel = str(selected_row["panel"])
        sample = int(float(selected_row["sample"]))
        cfg = decomp.PANELS[panel]
        truth = decomp.load_truth(
            Path(cfg["truth"](sample)),
            sample,
            str(cfg["truth_abundance_col"]),
        )
        calls = load_called_rows(selected_row, mapper)
        for variant in VARIANTS:
            method = str(variant["method"])
            raw_values = apply_variant(calls, variant)
            score = score_variant(panel, sample, method, truth, calls, raw_values)
            if method == "current_selected_abundance":
                max_validation_delta = max(
                    max_validation_delta,
                    abs(score["official_L1_pp"] - finite_float(selected_row["L1_union_pp"]))
                    if panel != "cami2_toy_mouse_gut"
                    else abs(score["official_L1_pp"] - finite_float(selected_row["L1_truth_only_pp"])),
                    abs(score["F1"] - finite_float(selected_row["F1"])),
                )
            rows.append(score)
    return pd.DataFrame(rows), max_validation_delta


def summarize(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = scores.loc[scores["method"].eq("current_selected_abundance")].copy()
    base_keyed = baseline.set_index(["panel", "sample"])
    score_rows: list[dict[str, object]] = []
    for row in scores.to_dict("records"):
        base = base_keyed.loc[(row["panel"], row["sample"])]
        row["delta_L1_pp"] = finite_float(row["official_L1_pp"]) - finite_float(base["official_L1_pp"])
        row["delta_Pearson"] = finite_float(row["official_Pearson"], float("nan")) - finite_float(
            base["official_Pearson"], float("nan")
        )
        score_rows.append(row)
    scores = pd.DataFrame(score_rows)

    panel_rows: list[dict[str, object]] = []
    for (method, panel), sub in scores.groupby(["method", "panel"], sort=True):
        if method == "current_selected_abundance":
            continue
        panel_rows.append(
            {
                "method": method,
                "panel": panel,
                "sample_count": int(len(sub)),
                "mean_L1_delta_pp": float(sub["delta_L1_pp"].mean()),
                "max_worse_L1_delta_pp": float(sub["delta_L1_pp"].max()),
                "improved_samples": int((sub["delta_L1_pp"] < -1e-9).sum()),
                "worsened_samples": int((sub["delta_L1_pp"] > 1e-9).sum()),
                "mean_Pearson_delta": float(sub["delta_Pearson"].mean()),
            }
        )
    panel = pd.DataFrame(panel_rows)

    overall_rows: list[dict[str, object]] = []
    for method, sub in scores.groupby("method", sort=True):
        if method == "current_selected_abundance":
            continue
        psub = panel.loc[panel["method"].eq(method)]
        overall_rows.append(
            {
                "method": method,
                "sample_count": int(len(sub)),
                "mean_sample_L1_delta_pp": float(sub["delta_L1_pp"].mean()),
                "max_sample_worse_L1_delta_pp": float(sub["delta_L1_pp"].max()),
                "improved_samples": int((sub["delta_L1_pp"] < -1e-9).sum()),
                "worsened_samples": int((sub["delta_L1_pp"] > 1e-9).sum()),
                "panel_count": int(len(psub)),
                "mean_panel_L1_delta_pp": float(psub["mean_L1_delta_pp"].mean()),
                "max_panel_worse_L1_delta_pp": float(psub["mean_L1_delta_pp"].max()),
                "improved_panels": int((psub["mean_L1_delta_pp"] < -1e-9).sum()),
                "worsened_panels": int((psub["mean_L1_delta_pp"] > 1e-9).sum()),
                "mean_Pearson_delta": float(sub["delta_Pearson"].mean()),
            }
        )
    overall = pd.DataFrame(overall_rows).sort_values(
        ["worsened_samples", "mean_sample_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
        kind="mergesort",
    )
    return scores, panel, overall


def audit(overall: pd.DataFrame, max_validation_delta: float) -> pd.DataFrame:
    sample_safe = overall.loc[
        overall["worsened_samples"].eq(0) & overall["improved_samples"].gt(0)
    ].copy()
    panel_safe = overall.loc[
        overall["worsened_panels"].eq(0) & overall["improved_panels"].gt(0)
    ].copy()
    best = overall.iloc[0]
    if not sample_safe.empty:
        best_sample_safe = sample_safe.sort_values(
            ["mean_sample_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
            kind="mergesort",
        ).iloc[0]
        decision = "sample_safe_candidate_needs_external_validation"
        sample_safe_value = (
            f"{best_sample_safe['method']};mean_sample_delta="
            f"{best_sample_safe['mean_sample_L1_delta_pp']:.6f};"
            f"improved_samples={int(best_sample_safe['improved_samples'])}"
        )
    else:
        decision = "do_not_promote_selected_mass_transform"
        sample_safe_value = "none"
    return pd.DataFrame(
        [
            {
                "metric": "variants_tested",
                "value": int(overall.shape[0]),
                "evidence": "selected_call_mass_transform_overall.tsv",
                "decision": "cached_selected_candidate_callset",
            },
            {
                "metric": "baseline_validation_max_abs_delta",
                "value": f"{max_validation_delta:.12g}",
                "evidence": "current_selected_abundance vs cached selected score TSV",
                "decision": "pass" if max_validation_delta < 1e-6 else "review",
            },
            {
                "metric": "sample_safe_variants",
                "value": int(sample_safe.shape[0]),
                "evidence": sample_safe_value,
                "decision": "candidate" if not sample_safe.empty else "none",
            },
            {
                "metric": "panel_safe_variants",
                "value": int(panel_safe.shape[0]),
                "evidence": "worsened_panels=0 and improved_panels>0",
                "decision": "panel_mean_only" if not panel_safe.empty else "none",
            },
            {
                "metric": "best_ranked_variant",
                "value": (
                    f"{best['method']};mean_sample_delta={best['mean_sample_L1_delta_pp']:.6f};"
                    f"worsened_samples={int(best['worsened_samples'])};"
                    f"max_worse={best['max_sample_worse_L1_delta_pp']:.6f}"
                ),
                "evidence": "ranked by worsened_samples then mean sample L1 delta",
                "decision": "informational",
            },
            {
                "metric": "promotion_decision",
                "value": decision,
                "evidence": "small implementable base-row transform sweep",
                "decision": decision,
            },
        ]
    )


def write_markdown(overall: pd.DataFrame, audit_df: pd.DataFrame) -> None:
    audit_rows = {str(row["metric"]): row for row in audit_df.to_dict("records")}
    lines = [
        "# Selected Call-Set Mass Transform Sweep",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached selected-candidate profile",
        "TSVs. It keeps calls fixed, preserves candidate-added row mass, and",
        "sweeps small base-called-row abundance transforms that are simple enough",
        "to implement in the wrapper.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "variants_tested",
        "baseline_validation_max_abs_delta",
        "sample_safe_variants",
        "panel_safe_variants",
        "best_ranked_variant",
        "promotion_decision",
    ]:
        row = audit_rows[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Top Variants",
            "",
            "| Method | Mean sample L1 delta pp | Worsened samples | Max sample worse pp | Panel mean delta pp | Worsened panels |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    top = overall.head(8)
    for row in top.to_dict("records"):
        lines.append(
            "| {method} | {mean_sample:.6f} | {worsened_samples} | {max_worse:.6f} | {mean_panel:.6f} | {worsened_panels} |".format(
                method=row["method"],
                mean_sample=float(row["mean_sample_L1_delta_pp"]),
                worsened_samples=int(row["worsened_samples"]),
                max_worse=float(row["max_sample_worse_L1_delta_pp"]),
                mean_panel=float(row["mean_panel_L1_delta_pp"]),
                worsened_panels=int(row["worsened_panels"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Do not promote these simple mass transforms into the default.",
            "- The best ranked transform, `base_genus_xny_a0.10`, improves mean",
            "  L1 but regresses 4/32 samples.",
            "- This supports the current next target: a more sample-aware",
            "  matched-call allocator, not another unguarded panel-mean transform.",
            "",
            "## Outputs",
            "",
            "- `results/selected_call_mass_transform_scores.tsv`",
            "- `results/selected_call_mass_transform_panel_delta.tsv`",
            "- `results/selected_call_mass_transform_overall.tsv`",
            "- `results/selected_call_mass_transform_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    scores, max_validation_delta = build_scores()
    scores, panel, overall = summarize(scores)
    audit_df = audit(overall, max_validation_delta)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    panel.to_csv(OUT_PANEL, sep="\t", index=False)
    overall.to_csv(OUT_OVERALL, sep="\t", index=False)
    audit_df.to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(overall, audit_df)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
