#!/usr/bin/env python3
"""Audit genus-XnY abundance blending on selected candidate-preset outputs.

This is a cached post-processing audit: it reads the 32 selected-default
candidate-preset profiles from validate_candidate_preset_replay.py, applies
the same base-row within-genus XnY abundance blend used by the wrapper option,
and scores the adjusted profiles with the existing GTDB panel scorers.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import score_adaptive_call_filter_wrapper_validation as adaptive
import validate_integrated_candidate_surface_hmp_abundance_policy as hmp_abundance


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_candidate_preset_genus_xny_blend_20260629")

BASELINE_SCORES = RESULTS / "candidate_preset_replay_scores.tsv"
METHOD = "default_candidate_preset_genus_xny_blend_a0.25"
BASELINE_METHOD = "default_candidate_preset_normalized_depth_alpha2"
ALPHA = 0.25


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(0.0, index=df.index)
    return pd.to_numeric(df[column], errors="coerce").fillna(0.0)


def genus_from_name(name: object) -> str:
    text = str(name or "").strip()
    if text.startswith("s__"):
        text = text[3:]
    parts = text.split()
    if parts:
        return parts[0]
    return text


def path_map_from_baseline() -> dict[tuple[str, int], Path]:
    if not BASELINE_SCORES.exists():
        raise SystemExit(f"missing baseline scores: {BASELINE_SCORES}")
    scores = pd.read_csv(BASELINE_SCORES, sep="\t")
    paths: dict[tuple[str, int], Path] = {}
    for row in scores.itertuples(index=False):
        panel = str(row.panel)
        sample = int(row.sample)
        profile = Path(str(row.profile))
        if profile.exists():
            paths[(panel, sample)] = profile
    if len(paths) != 32:
        missing = 32 - len(paths)
        raise SystemExit(f"expected 32 candidate-preset profiles, found {len(paths)}; missing={missing}")
    return paths


def apply_blend_to_profile(panel: str, sample: int, src: Path, dst: Path) -> dict[str, object]:
    df = pd.read_csv(src, sep="\t", low_memory=False)
    called = bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
    rescue_added = bool_series(df.get("candidate_rescue_added", pd.Series(False, index=df.index)))
    surface_added = bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    candidate_added = rescue_added | surface_added
    raw = numeric_series(df, "calibrated_abundance_raw")
    adjusted = raw.copy()

    eligible = called & ~candidate_added & (raw > 0.0)
    quality = numeric_series(df, "s_XnY_ctx_max").clip(lower=0.0, upper=1000.0) / 1000.0
    work = pd.DataFrame(
        {
            "idx": df.index[eligible],
            "genus": [genus_from_name(value) for value in df.loc[eligible, "species_name"]],
            "base": raw.loc[eligible].to_numpy(dtype=float),
            "quality": quality.loc[eligible].to_numpy(dtype=float),
        }
    )
    adjusted_rows = 0
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
        reallocated = weighted / weighted_sum * total
        blended = (1.0 - ALPHA) * base_values + ALPHA * reallocated
        adjusted.loc[idx] = blended
        adjusted_rows += int((abs(blended - base_values) > 1e-12).sum())

    df["calibrated_abundance_raw"] = adjusted
    called_mass = float(adjusted.loc[called].sum())
    df["calibrated_abundance"] = 0.0
    if called_mass > 0.0:
        df.loc[called, "calibrated_abundance"] = adjusted.loc[called] / called_mass
    df["abundance_rule"] = "posthoc_base_called_raw_then_within_genus_xny_blend"
    df["abundance_genus_xny_blend_alpha"] = ALPHA
    df["abundance_genus_xny_quality"] = "clip(s_XnY_ctx_max/1000,0,1)"

    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dst, sep="\t", index=False)

    return {
        "panel": panel,
        "sample": sample,
        "source_profile": str(src),
        "adjusted_profile": str(dst),
        "called_rows": int(called.sum()),
        "candidate_added_rows": int(candidate_added.sum()),
        "eligible_base_rows": int(eligible.sum()),
        "adjusted_base_rows": adjusted_rows,
        "base_raw_mass_before": float(raw.loc[called & ~candidate_added].sum()),
        "base_raw_mass_after": float(adjusted.loc[called & ~candidate_added].sum()),
        "candidate_raw_mass_before": float(raw.loc[called & candidate_added].sum()),
        "candidate_raw_mass_after": float(adjusted.loc[called & candidate_added].sum()),
        "called_raw_mass_before": float(raw.loc[called].sum()),
        "called_raw_mass_after": called_mass,
    }


def build_adjusted_profiles(paths: dict[tuple[str, int], Path]) -> tuple[dict[tuple[str, int], Path], pd.DataFrame]:
    adjusted_paths: dict[tuple[str, int], Path] = {}
    rows: list[dict[str, object]] = []
    for (panel, sample), src in sorted(paths.items()):
        dst = OUT_DIR / f"{panel}.sample{sample}.genus_xny_a{ALPHA:g}.tsv"
        rows.append(apply_blend_to_profile(panel, sample, src, dst))
        adjusted_paths[(panel, sample)] = dst
    return adjusted_paths, pd.DataFrame(rows)


def score_adjusted(paths: dict[tuple[str, int], Path]) -> pd.DataFrame:
    adaptive.METHOD = METHOD
    hmp_abundance.METHOD = METHOD
    emitted_panels = {"cami2_toy_mouse_gut", "cami3_toy_human_gut_gtdb_source_readmap"}
    emitted = {
        key: path
        for key, path in paths.items()
        if key[0] in emitted_panels
    }
    hmp_paths = {key: path for key, path in paths.items() if key[0] not in emitted_panels}
    emitted_scores, _emitted_meta = adaptive.score_profiles(emitted)
    mapper = raw_audit.RawMapper()
    hmp_rows: list[dict[str, object]] = []
    for (panel, sample), path in hmp_paths.items():
        row = hmp_abundance.score_profile(panel, int(sample), path, mapper)
        row["method"] = METHOD
        row["panel"] = panel
        row["sample"] = str(sample)
        row["profile"] = str(path)
        hmp_rows.append(row)
    return pd.concat([emitted_scores, pd.DataFrame(hmp_rows)], ignore_index=True, sort=False)


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, sub in scores.groupby("panel", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        l1_col = adaptive.official_l1_col(panel)
        pearson_col = adaptive.official_pearson_col(panel)
        rows.append(
            {
                "panel": panel,
                "method": METHOD,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "official_L1_pp": sub[l1_col].mean(),
                "official_Pearson": sub[pearson_col].mean(),
            }
        )
    return pd.DataFrame(rows)


def compare_to_baseline(scores: pd.DataFrame) -> pd.DataFrame:
    baseline = pd.read_csv(BASELINE_SCORES, sep="\t")
    rows: list[dict[str, object]] = []
    for row in scores.itertuples(index=False):
        panel = str(row.panel)
        sample = int(row.sample)
        ref = baseline.loc[
            baseline["panel"].astype(str).eq(panel)
            & baseline["sample"].astype(int).eq(sample)
            & baseline["method"].astype(str).eq(BASELINE_METHOD)
        ]
        if ref.empty:
            rows.append({"panel": panel, "sample": sample, "status": "missing_baseline"})
            continue
        ref_row = ref.iloc[0]
        l1_col = adaptive.official_l1_col(panel)
        pearson_col = adaptive.official_pearson_col(panel)
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "status": "compared",
                "TP_delta": int(row.TP) - int(ref_row["TP"]),
                "FP_delta": int(row.FP) - int(ref_row["FP"]),
                "FN_delta": int(row.FN) - int(ref_row["FN"]),
                "F1_delta": finite(row.F1) - finite(ref_row["F1"]),
                "official_L1_delta_pp": finite(getattr(row, l1_col)) - finite(ref_row[l1_col]),
                "official_Pearson_delta": finite(getattr(row, pearson_col)) - finite(ref_row[pearson_col]),
            }
        )
    return pd.DataFrame(rows)


def panel_delta_summary(comparison: pd.DataFrame) -> pd.DataFrame:
    compared = comparison.loc[comparison["status"].astype(str).eq("compared")].copy()
    rows: list[dict[str, object]] = []
    for panel, sub in compared.groupby("panel", sort=True):
        rows.append(
            {
                "panel": panel,
                "method": METHOD,
                "sample_count": int(len(sub)),
                "mean_L1_delta_pp": sub["official_L1_delta_pp"].mean(),
                "max_worse_L1_delta_pp": sub["official_L1_delta_pp"].max(),
                "improved_samples": int((sub["official_L1_delta_pp"] < -1e-9).sum()),
                "worsened_samples": int((sub["official_L1_delta_pp"] > 1e-9).sum()),
                "mean_Pearson_delta": sub["official_Pearson_delta"].mean(),
                "max_count_delta": int(sub[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max()),
                "max_F1_delta": float(sub["F1_delta"].abs().max()),
            }
        )
    return pd.DataFrame(rows)


def audit(comparison: pd.DataFrame, panel_deltas: pd.DataFrame) -> pd.DataFrame:
    compared = comparison.loc[comparison["status"].astype(str).eq("compared")].copy()
    max_count_delta = int(compared[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max())
    max_f1_delta = float(compared["F1_delta"].abs().max())
    mean_l1_delta = float(compared["official_L1_delta_pp"].mean())
    max_worse_l1 = float(compared["official_L1_delta_pp"].max())
    worsened_samples = int((compared["official_L1_delta_pp"] > 1e-9).sum())
    improved_samples = int((compared["official_L1_delta_pp"] < -1e-9).sum())
    worsened_panels = int((panel_deltas["mean_L1_delta_pp"] > 1e-9).sum())
    decision = (
        "promotable_abundance_component_pending_wrapper_rerun"
        if max_count_delta == 0
        and max_f1_delta < 1e-12
        and max_worse_l1 <= 0.0
        and worsened_panels == 0
        else "do_not_promote_abundance_blend"
    )
    return pd.DataFrame(
        [
            {
                "metric": "profiles_scored",
                "value": int(compared.shape[0]),
                "evidence": str(OUT_DIR),
                "decision": "cached_candidate_preset_posthoc_audit",
            },
            {
                "metric": "max_call_count_delta_vs_candidate_preset",
                "value": max_count_delta,
                "evidence": str(BASELINE_SCORES),
                "decision": "pass" if max_count_delta == 0 and max_f1_delta < 1e-12 else "fail",
            },
            {
                "metric": "mean_L1_delta_vs_candidate_preset_pp",
                "value": mean_l1_delta,
                "evidence": "candidate_preset_genus_xny_blend_vs_baseline.tsv",
                "decision": "improves" if mean_l1_delta < 0.0 else "worse_or_equal",
            },
            {
                "metric": "max_worse_L1_delta_vs_candidate_preset_pp",
                "value": max_worse_l1,
                "evidence": "candidate_preset_genus_xny_blend_vs_baseline.tsv",
                "decision": "sample_safe" if max_worse_l1 <= 0.0 else "has_sample_regression",
            },
            {
                "metric": "sample_L1_direction",
                "value": f"improved={improved_samples};worsened={worsened_samples}",
                "evidence": "candidate_preset_genus_xny_blend_vs_baseline.tsv",
                "decision": "sample_safe" if worsened_samples == 0 else "mixed",
            },
            {
                "metric": "panel_L1_direction",
                "value": f"worsened_panels={worsened_panels}",
                "evidence": "candidate_preset_genus_xny_blend_panel_delta.tsv",
                "decision": "panel_safe" if worsened_panels == 0 else "mixed",
            },
            {
                "metric": "promotion_decision",
                "value": decision,
                "evidence": "posthoc wrapper-equivalent abundance-only audit; no raw profiling rerun",
                "decision": decision,
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    baseline_paths = path_map_from_baseline()
    adjusted_paths, metadata = build_adjusted_profiles(baseline_paths)
    scores = score_adjusted(adjusted_paths)
    summary = summarize(scores)
    comparison = compare_to_baseline(scores)
    panel_deltas = panel_delta_summary(comparison)
    audit_df = audit(comparison, panel_deltas)

    metadata.to_csv(RESULTS / "candidate_preset_genus_xny_blend_metadata.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "candidate_preset_genus_xny_blend_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "candidate_preset_genus_xny_blend_summary.tsv", sep="\t", index=False)
    comparison.to_csv(
        RESULTS / "candidate_preset_genus_xny_blend_vs_baseline.tsv",
        sep="\t",
        index=False,
    )
    panel_deltas.to_csv(
        RESULTS / "candidate_preset_genus_xny_blend_panel_delta.tsv",
        sep="\t",
        index=False,
    )
    audit_df.to_csv(RESULTS / "candidate_preset_genus_xny_blend_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nPANEL DELTAS")
    print(panel_deltas.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
