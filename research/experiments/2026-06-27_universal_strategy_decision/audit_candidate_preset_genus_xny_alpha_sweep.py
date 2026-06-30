#!/usr/bin/env python3
"""Sweep conservative selected-default genus-XnY abundance blend strengths."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import audit_candidate_preset_genus_xny_blend as blend


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_candidate_preset_genus_xny_alpha_sweep_20260629")

VARIANTS = [
    {"method": "genus_xny_a0.02", "alpha": 0.02, "cap": None},
    {"method": "genus_xny_a0.05", "alpha": 0.05, "cap": None},
    {"method": "genus_xny_a0.10", "alpha": 0.10, "cap": None},
    {"method": "genus_xny_a0.15", "alpha": 0.15, "cap": None},
    {"method": "genus_xny_a0.25_cap0.005", "alpha": 0.25, "cap": 0.005},
    {"method": "genus_xny_a0.25", "alpha": 0.25, "cap": None},
]


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(0.0, index=df.index)
    return pd.to_numeric(df[column], errors="coerce").fillna(0.0)


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def apply_variant_to_profile(
    panel: str,
    sample: int,
    src: Path,
    dst: Path,
    method: str,
    alpha: float,
    cap_sample_base_frac: float | None,
) -> dict[str, object]:
    df = pd.read_csv(src, sep="\t", low_memory=False)
    called = bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
    rescue_added = bool_series(df.get("candidate_rescue_added", pd.Series(False, index=df.index)))
    surface_added = bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    candidate_added = rescue_added | surface_added
    raw = numeric_series(df, "calibrated_abundance_raw")
    adjusted = raw.copy()

    eligible = called & ~candidate_added & (raw > 0.0)
    base_total = float(raw.loc[eligible].sum())
    quality = numeric_series(df, "s_XnY_ctx_max").clip(lower=0.0, upper=1000.0) / 1000.0
    work = pd.DataFrame(
        {
            "idx": df.index[eligible],
            "genus": [blend.genus_from_name(value) for value in df.loc[eligible, "species_name"]],
            "base": raw.loc[eligible].to_numpy(dtype=float),
            "quality": quality.loc[eligible].to_numpy(dtype=float),
        }
    )

    adjusted_rows = 0
    capped_genus_count = 0
    min_effective_alpha = alpha
    for _genus, sub in work.groupby("genus", sort=False):
        if len(sub) <= 1:
            continue
        idx = sub["idx"].to_numpy()
        base_values = sub["base"].to_numpy(dtype=float)
        weighted = base_values * sub["quality"].to_numpy(dtype=float)
        total = float(base_values.sum())
        weighted_sum = float(weighted.sum())
        if total <= 0.0 or weighted_sum <= 0.0:
            continue
        target_delta = weighted / weighted_sum * total - base_values
        effective_alpha = alpha
        if cap_sample_base_frac is not None and base_total > 0.0:
            max_abs_delta = float(abs(target_delta).max())
            cap_abs_delta = cap_sample_base_frac * base_total
            if max_abs_delta > 0.0 and alpha * max_abs_delta > cap_abs_delta:
                effective_alpha = cap_abs_delta / max_abs_delta
                capped_genus_count += 1
        min_effective_alpha = min(min_effective_alpha, effective_alpha)
        blended = base_values + effective_alpha * target_delta
        adjusted.loc[idx] = blended
        adjusted_rows += int((abs(blended - base_values) > 1e-12).sum())

    df["calibrated_abundance_raw"] = adjusted
    called_mass = float(adjusted.loc[called].sum())
    df["calibrated_abundance"] = 0.0
    if called_mass > 0.0:
        df.loc[called, "calibrated_abundance"] = adjusted.loc[called] / called_mass
    df["abundance_rule"] = (
        "posthoc_base_called_raw_then_within_genus_xny_alpha_sweep"
    )
    df["abundance_genus_xny_blend_alpha"] = alpha
    df["abundance_genus_xny_quality"] = "clip(s_XnY_ctx_max/1000,0,1)"
    df["abundance_genus_xny_cap_sample_base_frac"] = (
        "" if cap_sample_base_frac is None else cap_sample_base_frac
    )

    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dst, sep="\t", index=False)

    return {
        "method": method,
        "panel": panel,
        "sample": sample,
        "source_profile": str(src),
        "adjusted_profile": str(dst),
        "alpha": alpha,
        "cap_sample_base_frac": "" if cap_sample_base_frac is None else cap_sample_base_frac,
        "called_rows": int(called.sum()),
        "candidate_added_rows": int(candidate_added.sum()),
        "eligible_base_rows": int(eligible.sum()),
        "adjusted_base_rows": adjusted_rows,
        "capped_genus_count": capped_genus_count,
        "min_effective_alpha": min_effective_alpha,
        "base_raw_mass_before": float(raw.loc[called & ~candidate_added].sum()),
        "base_raw_mass_after": float(adjusted.loc[called & ~candidate_added].sum()),
        "candidate_raw_mass_before": float(raw.loc[called & candidate_added].sum()),
        "candidate_raw_mass_after": float(adjusted.loc[called & candidate_added].sum()),
        "called_raw_mass_before": float(raw.loc[called].sum()),
        "called_raw_mass_after": called_mass,
    }


def build_variant_profiles(
    variant: dict[str, object],
    paths: dict[tuple[str, int], Path],
) -> tuple[dict[tuple[str, int], Path], pd.DataFrame]:
    method = str(variant["method"])
    alpha = float(variant["alpha"])
    cap = variant["cap"]
    cap_value = None if cap is None else float(cap)
    adjusted_paths: dict[tuple[str, int], Path] = {}
    rows: list[dict[str, object]] = []
    for (panel, sample), src in sorted(paths.items()):
        dst = OUT_DIR / method / f"{panel}.sample{sample}.{method}.tsv"
        rows.append(apply_variant_to_profile(panel, sample, src, dst, method, alpha, cap_value))
        adjusted_paths[(panel, sample)] = dst
    return adjusted_paths, pd.DataFrame(rows)


def annotate_method(df: pd.DataFrame, method: str) -> pd.DataFrame:
    out = df.copy()
    out["method"] = method
    return out


def summarize_by_variant(panel_deltas: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method, sub in panel_deltas.groupby("method", sort=True):
        rows.append(
            {
                "method": method,
                "panel_count": int(sub.shape[0]),
                "mean_panel_L1_delta_pp": float(sub["mean_L1_delta_pp"].mean()),
                "max_panel_worse_L1_delta_pp": float(sub["mean_L1_delta_pp"].max()),
                "improved_panels": int((sub["mean_L1_delta_pp"] < -1e-9).sum()),
                "worsened_panels": int((sub["mean_L1_delta_pp"] > 1e-9).sum()),
                "sample_count": int(sub["sample_count"].sum()),
                "improved_samples": int(sub["improved_samples"].sum()),
                "worsened_samples": int(sub["worsened_samples"].sum()),
                "max_sample_worse_L1_delta_pp": float(sub["max_worse_L1_delta_pp"].max()),
                "mean_panel_Pearson_delta": float(sub["mean_Pearson_delta"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["worsened_samples", "mean_panel_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
        kind="mergesort",
    )


def audit(overall: pd.DataFrame) -> pd.DataFrame:
    sample_safe = overall.loc[
        overall["worsened_samples"].eq(0) & overall["improved_samples"].gt(0)
    ].copy()
    panel_safe = overall.loc[
        overall["worsened_panels"].eq(0) & overall["improved_panels"].gt(0)
    ].copy()
    best = overall.sort_values(
        ["worsened_samples", "mean_panel_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
        kind="mergesort",
    ).iloc[0]
    if not sample_safe.empty:
        best_sample_safe = sample_safe.sort_values(
            ["mean_panel_L1_delta_pp", "max_sample_worse_L1_delta_pp"],
            kind="mergesort",
        ).iloc[0]
        decision = "sample_safe_alpha_candidate_needs_wrapper_validation"
        sample_safe_value = (
            f"{best_sample_safe['method']};mean_panel_delta="
            f"{best_sample_safe['mean_panel_L1_delta_pp']:.6f};"
            f"max_worse={best_sample_safe['max_sample_worse_L1_delta_pp']:.6f};"
            f"improved_samples={int(best_sample_safe['improved_samples'])}"
        )
    else:
        decision = "do_not_promote_alpha_sweep"
        sample_safe_value = "none"
    return pd.DataFrame(
        [
            {
                "metric": "variants_tested",
                "value": int(overall.shape[0]),
                "evidence": "candidate_preset_genus_xny_alpha_sweep_overall.tsv",
                "decision": "cached_selected_default_posthoc_sweep",
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
                    f"{best['method']};mean_panel_delta={best['mean_panel_L1_delta_pp']:.6f};"
                    f"worsened_samples={int(best['worsened_samples'])};"
                    f"max_worse={best['max_sample_worse_L1_delta_pp']:.6f}"
                ),
                "evidence": "ranked by worsened_samples then mean_panel_L1_delta_pp",
                "decision": "informational",
            },
            {
                "metric": "promotion_decision",
                "value": decision,
                "evidence": "alpha/cap sweep over selected-default profiles",
                "decision": decision,
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    baseline_paths = blend.path_map_from_baseline()

    all_scores: list[pd.DataFrame] = []
    all_summaries: list[pd.DataFrame] = []
    all_comparisons: list[pd.DataFrame] = []
    all_panel_deltas: list[pd.DataFrame] = []
    all_metadata: list[pd.DataFrame] = []
    for variant in VARIANTS:
        method = str(variant["method"])
        adjusted_paths, metadata = build_variant_profiles(variant, baseline_paths)
        blend.METHOD = method
        scores = annotate_method(blend.score_adjusted(adjusted_paths), method)
        summary = annotate_method(blend.summarize(scores), method)
        comparison = annotate_method(blend.compare_to_baseline(scores), method)
        panel_deltas = annotate_method(blend.panel_delta_summary(comparison), method)
        all_scores.append(scores)
        all_summaries.append(summary)
        all_comparisons.append(comparison)
        all_panel_deltas.append(panel_deltas)
        all_metadata.append(metadata)

    scores_df = pd.concat(all_scores, ignore_index=True, sort=False)
    summary_df = pd.concat(all_summaries, ignore_index=True, sort=False)
    comparison_df = pd.concat(all_comparisons, ignore_index=True, sort=False)
    panel_df = pd.concat(all_panel_deltas, ignore_index=True, sort=False)
    metadata_df = pd.concat(all_metadata, ignore_index=True, sort=False)
    overall_df = summarize_by_variant(panel_df)
    audit_df = audit(overall_df)

    metadata_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_metadata.tsv", sep="\t", index=False)
    scores_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_summary.tsv", sep="\t", index=False)
    comparison_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_vs_baseline.tsv", sep="\t", index=False)
    panel_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_panel_delta.tsv", sep="\t", index=False)
    overall_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_overall.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / "candidate_preset_genus_xny_alpha_sweep_audit.tsv", sep="\t", index=False)

    print(overall_df.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
