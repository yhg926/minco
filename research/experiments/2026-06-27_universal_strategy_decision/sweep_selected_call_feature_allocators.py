#!/usr/bin/env python3
"""Sweep feature-derived allocators on the selected-candidate call set.

The previous selected-call mass-transform sweep only tested small transforms of
the existing calibrated raw mass. This sweep asks whether abundance-like output
features already emitted by MinCO, especially normalized/relative depth from
unique and split views, can provide a better matched-call mass target while the
selected call set remains fixed.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, Mapping

import numpy as np
import pandas as pd

import audit_candidate_callset_oracle_feasibility as feasibility
import decompose_abundance_errors as decomp
import sweep_selected_call_mass_transforms as base_sweep


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SELECTED_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"

OUT_SCORES = RESULTS / "selected_call_feature_allocator_scores.tsv"
OUT_PANEL = RESULTS / "selected_call_feature_allocator_panel_delta.tsv"
OUT_OVERALL = RESULTS / "selected_call_feature_allocator_overall.tsv"
OUT_AUDIT = RESULTS / "selected_call_feature_allocator_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_SWEEP.md"

ALPHAS = [0.005, 0.01, 0.02, 0.05, 0.10, 0.25]
MAX_SAVED_SCORE_METHODS = 32


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def numeric(calls: pd.DataFrame, col: str, default: float = 0.0) -> np.ndarray:
    return base_sweep.numeric(calls, col, default).to_numpy(dtype=float)


def positive(values: np.ndarray) -> np.ndarray:
    return np.nan_to_num(np.maximum(np.asarray(values, dtype=float), 0.0), nan=0.0)


def max_pair(calls: pd.DataFrame, left: str, right: str, default: float = 0.0) -> np.ndarray:
    return np.maximum(numeric(calls, left, default), numeric(calls, right, default))


def sum_pair(calls: pd.DataFrame, left: str, right: str, default: float = 0.0) -> np.ndarray:
    return numeric(calls, left, default) + numeric(calls, right, default)


def normalized(values: np.ndarray, total: float) -> np.ndarray:
    clean = positive(values)
    clean_total = float(clean.sum())
    if clean_total <= 0.0 or total <= 0.0:
        return clean
    return clean / clean_total * total


def feature_table(calls: pd.DataFrame, raw: np.ndarray) -> dict[str, np.ndarray]:
    s_norm = numeric(calls, "s_Normalized_abundance_depth_max")
    u_norm = numeric(calls, "u_Normalized_abundance_depth_max")
    s_rel = numeric(calls, "s_Relative_abundance_depth_max")
    u_rel = numeric(calls, "u_Relative_abundance_depth_max")
    s_mean = numeric(calls, "s_Ref_mean_depth_max")
    u_mean = numeric(calls, "u_Ref_mean_depth_max")
    s_hit = numeric(calls, "s_Ref_hit_mean_depth_max")
    u_hit = numeric(calls, "u_Ref_hit_mean_depth_max")
    s_breadth = numeric(calls, "s_Ref_breadth_max")
    u_breadth = numeric(calls, "u_Ref_breadth_max")
    s_zip = numeric(calls, "s_Ref_zip_af_max")
    u_zip = numeric(calls, "u_Ref_zip_af_max")
    probability = np.clip(numeric(calls, "calibrated_probability", 1.0), 0.0, 1.0)
    real_af = max_pair(calls, "s_Real_min_align_fraction_max", "u_Real_min_align_fraction_max")
    xny = max_pair(calls, "s_XnY_ctx_max", "u_XnY_ctx_max")
    breadth = np.maximum(s_breadth, u_breadth)
    zip_af = np.maximum(s_zip, u_zip)
    mean_breadth = np.maximum(s_mean * s_breadth, u_mean * u_breadth)
    hit_breadth = np.maximum(s_hit * s_breadth, u_hit * u_breadth)
    norm_max = np.maximum(s_norm, u_norm)
    norm_sum = s_norm + u_norm
    rel_max = np.maximum(s_rel, u_rel)
    rel_sum = s_rel + u_rel

    direct = {
        "s_norm_depth": s_norm,
        "u_norm_depth": u_norm,
        "max_norm_depth": norm_max,
        "sum_norm_depth": norm_sum,
        "s_rel_depth": s_rel,
        "u_rel_depth": u_rel,
        "max_rel_depth": rel_max,
        "sum_rel_depth": rel_sum,
        "max_mean_breadth": mean_breadth,
        "max_hit_breadth": hit_breadth,
        "zip_norm_depth": norm_max * zip_af,
        "prob_norm_depth": norm_max * probability,
        "realaf_norm_depth": norm_max * real_af,
        "xny_norm_depth": norm_max * np.clip(xny / 1000.0, 0.0, 10.0),
    }
    quality = {
        "raw_prob": raw * probability,
        "raw_breadth": raw * breadth,
        "raw_zip": raw * zip_af,
        "raw_realaf": raw * real_af,
        "raw_xny": raw * np.clip(xny / 1000.0, 0.0, 10.0),
        "raw_norm_depth": raw * norm_max,
        "raw_mean_breadth": raw * mean_breadth,
        "raw_hit_breadth": raw * hit_breadth,
    }
    return {**direct, **quality}


def genus_values(calls: pd.DataFrame) -> np.ndarray:
    return np.asarray(
        [base_sweep.genus_from_species(value) for value in calls["gtdb_species"]],
        dtype=object,
    )


def apply_allocator(
    calls: pd.DataFrame,
    feature_name: str,
    mode: str,
    alpha: float,
) -> np.ndarray:
    raw, candidate_added = base_sweep.base_raw_values(calls)
    eligible = ~candidate_added & (raw > 0.0)
    out = raw.copy()
    base_total = float(raw[eligible].sum())
    if base_total <= 0.0 or not bool(eligible.any()):
        return out

    features = feature_table(calls, raw)
    target_feature = positive(features[feature_name])
    if mode == "global":
        target = normalized(target_feature[eligible], base_total)
        if float(target.sum()) <= 0.0:
            return out
        idx = np.flatnonzero(eligible)
        out[idx] = (1.0 - alpha) * raw[idx] + alpha * target
        return out

    if mode == "genus":
        genus = genus_values(calls)
        for value in sorted(set(genus[eligible])):
            group = eligible & (genus == value)
            if int(group.sum()) <= 1:
                continue
            group_total = float(raw[group].sum())
            target = normalized(target_feature[group], group_total)
            if float(target.sum()) <= 0.0:
                continue
            idx = np.flatnonzero(group)
            out[idx] = (1.0 - alpha) * raw[idx] + alpha * target
        return out

    raise ValueError(f"unsupported allocator mode: {mode}")


def variants(calls: pd.DataFrame) -> list[dict[str, object]]:
    raw, _candidate_added = base_sweep.base_raw_values(calls)
    names = sorted(feature_table(calls, raw))
    out: list[dict[str, object]] = [{"method": "current_selected_abundance", "kind": "current"}]
    for feature_name in names:
        for mode in ["global", "genus"]:
            for alpha in ALPHAS:
                out.append(
                    {
                        "method": f"{mode}_{feature_name}_a{alpha:g}",
                        "kind": "feature_allocator",
                        "feature": feature_name,
                        "mode": mode,
                        "alpha": alpha,
                    }
                )
    return out


def score_variant(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    calls: pd.DataFrame,
    raw_values: np.ndarray,
) -> dict[str, object]:
    abundance = base_sweep.collapse_prediction(calls, raw_values, panel)
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
    mapper = base_sweep.RowMapper()
    rows: list[dict[str, object]] = []
    max_validation_delta = 0.0
    variant_template: list[dict[str, object]] | None = None

    for selected_row in selected_scores.to_dict("records"):
        panel = str(selected_row["panel"])
        sample = int(float(selected_row["sample"]))
        cfg = decomp.PANELS[panel]
        truth = decomp.load_truth(
            Path(cfg["truth"](sample)),
            sample,
            str(cfg["truth_abundance_col"]),
        )
        calls = base_sweep.load_called_rows(selected_row, mapper)
        if variant_template is None:
            variant_template = variants(calls)
        for variant in variant_template:
            method = str(variant["method"])
            if variant["kind"] == "current":
                raw_values, _candidate_added = base_sweep.base_raw_values(calls)
            else:
                raw_values = apply_allocator(
                    calls,
                    str(variant["feature"]),
                    str(variant["mode"]),
                    float(variant["alpha"]),
                )
            score = score_variant(panel, sample, method, truth, calls, raw_values)
            if method == "current_selected_abundance":
                target_l1 = (
                    finite_float(selected_row["L1_truth_only_pp"])
                    if panel == "cami2_toy_mouse_gut"
                    else finite_float(selected_row["L1_union_pp"])
                )
                max_validation_delta = max(
                    max_validation_delta,
                    abs(score["official_L1_pp"] - target_l1),
                    abs(score["F1"] - finite_float(selected_row["F1"])),
                )
            rows.append(score)
    return pd.DataFrame(rows), max_validation_delta


def summarize(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = scores.loc[scores["method"].eq("current_selected_abundance")].copy()
    base_keyed = baseline.set_index(["panel", "sample"])
    rows: list[dict[str, object]] = []
    for row in scores.to_dict("records"):
        base = base_keyed.loc[(row["panel"], row["sample"])]
        row["delta_L1_pp"] = finite_float(row["official_L1_pp"]) - finite_float(
            base["official_L1_pp"]
        )
        row["delta_Pearson"] = finite_float(row["official_Pearson"], float("nan")) - finite_float(
            base["official_Pearson"], float("nan")
        )
        rows.append(row)
    scores = pd.DataFrame(rows)

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
            f"improved_samples={int(best_sample_safe['improved_samples'])};"
            f"max_worse={best_sample_safe['max_sample_worse_L1_delta_pp']:.6f}"
        )
    else:
        decision = "do_not_promote_feature_allocator"
        sample_safe_value = "none"

    return pd.DataFrame(
        [
            {
                "metric": "variants_tested",
                "value": int(overall.shape[0]),
                "evidence": (
                    "selected_call_feature_allocator_overall.tsv;"
                    f"top_saved_score_methods={MAX_SAVED_SCORE_METHODS}"
                ),
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
                "evidence": "feature-derived selected-call allocator sweep",
                "decision": decision,
            },
        ]
    )


def write_markdown(overall: pd.DataFrame, audit_df: pd.DataFrame) -> None:
    audit_rows = {str(row["metric"]): row for row in audit_df.to_dict("records")}
    lines = [
        "# Selected Call-Set Feature Allocator Sweep",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached selected-candidate profile",
        "TSVs. It keeps calls fixed, preserves candidate-added row mass, and",
        "tests whether existing MinCO depth/quality output features can serve as",
        "matched-call abundance targets.",
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
    for row in overall.head(12).to_dict("records"):
        lines.append(
            "| {method} | {mean_sample:.6f} | {worse} | {max_worse:.6f} | {mean_panel:.6f} | {worse_panel} |".format(
                method=row["method"],
                mean_sample=float(row["mean_sample_L1_delta_pp"]),
                worse=int(row["worsened_samples"]),
                max_worse=float(row["max_sample_worse_L1_delta_pp"]),
                mean_panel=float(row["mean_panel_L1_delta_pp"]),
                worse_panel=int(row["worsened_panels"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Do not promote a feature-derived allocator unless it is sample-safe and",
            "  then survives external holdouts.",
            "- If a sample-safe variant appears here, treat it as a candidate for",
            "  independent replay rather than a default change.",
            f"- The per-sample score table saves the baseline plus the top",
            f"  {MAX_SAVED_SCORE_METHODS} ranked methods; the full sweep is",
            "  summarized in the overall and audit TSVs.",
            "",
            "## Outputs",
            "",
            "- `results/selected_call_feature_allocator_scores.tsv`",
            "- `results/selected_call_feature_allocator_panel_delta.tsv`",
            "- `results/selected_call_feature_allocator_overall.tsv`",
            "- `results/selected_call_feature_allocator_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    scores, max_validation_delta = build_scores()
    scores, panel, overall = summarize(scores)
    audit_df = audit(overall, max_validation_delta)
    saved_methods = ["current_selected_abundance"] + overall["method"].head(
        MAX_SAVED_SCORE_METHODS
    ).astype(str).tolist()
    scores.loc[scores["method"].isin(saved_methods)].to_csv(OUT_SCORES, sep="\t", index=False)
    panel.to_csv(OUT_PANEL, sep="\t", index=False)
    overall.to_csv(OUT_OVERALL, sep="\t", index=False)
    audit_df.to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(overall, audit_df)
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
