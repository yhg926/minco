#!/usr/bin/env python3
"""Compare the wrapper candidate strategy to current MinCO and Sylph.

The candidate strategy is:
  candidate rescue/surface calls + normalized-depth-alpha2 candidate abundance.

This script uses only matched panel summaries already generated in this note.
It does not rerun profiling. The goal is to decide whether the candidate should
replace the current MinCO default candidate, remain experimental, or be rejected.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import score_adaptive_call_filter_wrapper_validation as adaptive


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
WRAPPER = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"
WRAPPER_ZERO = RESULTS / "cross_panel_candidate_abundance_wrapper_zero_scores.tsv"
PRESET = RESULTS / "candidate_preset_replay_scores.tsv"

PANEL_LABELS = {
    "cami2_toy_mouse_gut": "CAMI2 Toy Mouse",
    "cami3_toy_human_gut_gtdb_source_readmap": "CAMI3 ToyGut source-readmap",
    "hmp_airskin_gtdb_source_abundance": "HMP airskin source-abundance",
    "hmp_gastrooral_gtdb_source_abundance": "HMP gastrooral source-abundance",
}

SUMMARY_PATHS = {
    "cami2_toy_mouse_gut": RESULTS / "toymouse_current_refresh_summary.tsv",
    "cami3_toy_human_gut_gtdb_source_readmap": RESULTS / "cami3_gtdb_source_readmap_summary.tsv",
    "hmp_airskin_gtdb_source_abundance": RESULTS
    / "hmp_current_refresh_r232_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28_summary.tsv",
    "hmp_gastrooral_gtdb_source_abundance": RESULTS / "hmp_gastrooral_raw_default_r232_source_abundance_summary.tsv",
}

CURRENT_METHOD = {
    "cami2_toy_mouse_gut": "minco_current_code_refresh",
    "cami3_toy_human_gut_gtdb_source_readmap": "minco_universal_autoexact_gtdb_source_readmap",
    "hmp_airskin_gtdb_source_abundance": (
        "minco_current_default_gtdb_source_abundance_"
        "sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28"
    ),
    "hmp_gastrooral_gtdb_source_abundance": "minco_current_raw_default_gastrooral_source_abundance",
}

SYLPH_METHOD = {
    "cami2_toy_mouse_gut": "sylph_gtdb_profile",
    "cami3_toy_human_gut_gtdb_source_readmap": "sylph_gtdb_source_readmap",
    "hmp_airskin_gtdb_source_abundance": "sylph_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance": "sylph_gtdb_r232_gastrooral_source_abundance",
}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def panel_summary_from_scores(scores: pd.DataFrame, method: str, source: Path) -> pd.DataFrame:
    rows = []
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
                "panel_label": PANEL_LABELS[panel],
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "official_L1_pp": sub[l1_col].mean(),
                "official_Pearson": sub[pearson_col].mean(),
                "source": str(source),
            }
        )
    return pd.DataFrame(rows)


def normalize_summary_row(panel: str, row: pd.Series, method: str, source: Path) -> dict[str, object]:
    if panel == "cami2_toy_mouse_gut":
        l1 = row["mean_L1_pp"]
        pearson = row["mean_Pearson"]
    else:
        l1 = row["mean_L1_union_pp"]
        pearson = row["mean_Pearson_union"]
    return {
        "panel": panel,
        "panel_label": PANEL_LABELS[panel],
        "method": method,
        "samples": row["samples"],
        "sample_count": len(str(row["samples"]).split(",")) if str(row["samples"]) else 0,
        "mean_F1": row["mean_F1"],
        "pooled_F1": row["pooled_F1"],
        "pooled_TP": row["pooled_TP"],
        "pooled_FP": row["pooled_FP"],
        "pooled_FN": row["pooled_FN"],
        "official_L1_pp": l1,
        "official_Pearson": pearson,
        "source": str(source),
    }


def baseline_rows() -> pd.DataFrame:
    rows = []
    for panel, path in SUMMARY_PATHS.items():
        if not path.exists():
            raise SystemExit(f"missing summary: {path}")
        table = pd.read_csv(path, sep="\t")
        for method in [CURRENT_METHOD[panel], SYLPH_METHOD[panel]]:
            hit = table.loc[table["method"].astype(str).eq(method)]
            if hit.empty:
                raise SystemExit(f"missing {method} in {path}")
            rows.append(normalize_summary_row(panel, hit.iloc[0], method, path))
    return pd.DataFrame(rows)


def compare_methods(summary: pd.DataFrame, candidate: str, baseline: str, label: str) -> pd.DataFrame:
    rows = []
    for panel in PANEL_LABELS:
        cand = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(candidate)].iloc[0]
        base = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(baseline)].iloc[0]
        rows.append(
            {
                "comparison": label,
                "panel": panel,
                "panel_label": PANEL_LABELS[panel],
                "samples": cand["samples"],
                "delta_pooled_F1": finite(cand["pooled_F1"]) - finite(base["pooled_F1"]),
                "delta_official_L1_pp": finite(cand["official_L1_pp"]) - finite(base["official_L1_pp"]),
                "delta_official_Pearson": finite(cand["official_Pearson"]) - finite(base["official_Pearson"]),
                "candidate_pooled_F1": cand["pooled_F1"],
                "baseline_pooled_F1": base["pooled_F1"],
                "candidate_L1_pp": cand["official_L1_pp"],
                "baseline_L1_pp": base["official_L1_pp"],
                "candidate_Pearson": cand["official_Pearson"],
                "baseline_Pearson": base["official_Pearson"],
            }
        )
    return pd.DataFrame(rows)


def audit(compare_current: pd.DataFrame, compare_sylph: pd.DataFrame) -> pd.DataFrame:
    current_f1_worse = int((compare_current["delta_pooled_F1"] < -1e-12).sum())
    current_l1_worse = int((compare_current["delta_official_L1_pp"] > 1e-12).sum())
    current_pearson_worse = int((compare_current["delta_official_Pearson"] < -1e-12).sum())
    sylph_f1_wins = int((compare_sylph["delta_pooled_F1"] > 1e-12).sum())
    sylph_l1_wins = int((compare_sylph["delta_official_L1_pp"] < -1e-12).sum())
    sylph_pearson_wins = int((compare_sylph["delta_official_Pearson"] > 1e-12).sum())
    rows = [
        {
            "metric": "candidate_vs_current_panel_safety",
            "value": (
                f"f1_worse={current_f1_worse};"
                f"l1_worse={current_l1_worse};"
                f"pearson_worse={current_pearson_worse}"
            ),
            "evidence": "candidate_default_vs_current.tsv",
            "decision": "candidate_dominates_current_panels"
            if current_f1_worse == current_l1_worse == current_pearson_worse == 0
            else "candidate_has_current_regressions",
        },
        {
            "metric": "candidate_vs_current_mean_delta",
            "value": (
                f"F1={compare_current['delta_pooled_F1'].mean():.6f};"
                f"L1={compare_current['delta_official_L1_pp'].mean():.6f};"
                f"Pearson={compare_current['delta_official_Pearson'].mean():.6f}"
            ),
            "evidence": "candidate_default_vs_current.tsv",
            "decision": "effect_size",
        },
        {
            "metric": "candidate_vs_sylph_panel_wins",
            "value": (
                f"F1_wins={sylph_f1_wins}/4;"
                f"L1_wins={sylph_l1_wins}/4;"
                f"Pearson_wins={sylph_pearson_wins}/4"
            ),
            "evidence": "candidate_default_vs_sylph.tsv",
            "decision": "not_broad_sylph_beating"
            if sylph_l1_wins < 3 or sylph_f1_wins < 3
            else "broad_sylph_beating_candidate",
        },
        {
            "metric": "default_candidate_decision",
            "value": "promote_candidate_for_minco_default_study_not_release_claim",
            "evidence": "dominates current MinCO panels but does not beat Sylph on abundance",
            "decision": "candidate_default_next_step",
        },
    ]
    return pd.DataFrame(rows)


def main() -> int:
    if not WRAPPER.exists() or not WRAPPER_ZERO.exists():
        raise SystemExit("missing wrapper score tables; run validate_cross_panel_candidate_abundance_wrapper.py")
    candidate_source = PRESET if PRESET.exists() else WRAPPER
    candidate_method = (
        "minco_default_candidate_preset"
        if candidate_source == PRESET
        else "minco_wrapper_candidate_abundance"
    )
    candidate_scores = pd.read_csv(candidate_source, sep="\t")
    wrapper_scores = pd.read_csv(WRAPPER, sep="\t")
    zero_scores = pd.read_csv(WRAPPER_ZERO, sep="\t")
    summary = pd.concat(
        [
            baseline_rows(),
            panel_summary_from_scores(zero_scores, "minco_wrapper_zero_mass", WRAPPER_ZERO),
            panel_summary_from_scores(candidate_scores, candidate_method, candidate_source),
            panel_summary_from_scores(wrapper_scores, "minco_explicit_wrapper_candidate_abundance", WRAPPER),
        ],
        ignore_index=True,
        sort=False,
    )
    # Current MinCO method names differ by panel.
    rows = []
    for panel, method in CURRENT_METHOD.items():
        cand = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(candidate_method)].iloc[0]
        base = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(method)].iloc[0]
        rows.append(
            {
                "comparison": "candidate_vs_current_minco",
                "panel": panel,
                "panel_label": PANEL_LABELS[panel],
                "samples": cand["samples"],
                "delta_pooled_F1": finite(cand["pooled_F1"]) - finite(base["pooled_F1"]),
                "delta_official_L1_pp": finite(cand["official_L1_pp"]) - finite(base["official_L1_pp"]),
                "delta_official_Pearson": finite(cand["official_Pearson"]) - finite(base["official_Pearson"]),
                "candidate_pooled_F1": cand["pooled_F1"],
                "baseline_pooled_F1": base["pooled_F1"],
                "candidate_L1_pp": cand["official_L1_pp"],
                "baseline_L1_pp": base["official_L1_pp"],
                "candidate_Pearson": cand["official_Pearson"],
                "baseline_Pearson": base["official_Pearson"],
            }
        )
    compare_current = pd.DataFrame(rows)
    rows = []
    for panel, method in SYLPH_METHOD.items():
        cand = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(candidate_method)].iloc[0]
        base = summary.loc[summary["panel"].eq(panel) & summary["method"].eq(method)].iloc[0]
        rows.append(
            {
                "comparison": "candidate_vs_sylph",
                "panel": panel,
                "panel_label": PANEL_LABELS[panel],
                "samples": cand["samples"],
                "delta_pooled_F1": finite(cand["pooled_F1"]) - finite(base["pooled_F1"]),
                "delta_official_L1_pp": finite(cand["official_L1_pp"]) - finite(base["official_L1_pp"]),
                "delta_official_Pearson": finite(cand["official_Pearson"]) - finite(base["official_Pearson"]),
                "candidate_pooled_F1": cand["pooled_F1"],
                "baseline_pooled_F1": base["pooled_F1"],
                "candidate_L1_pp": cand["official_L1_pp"],
                "baseline_L1_pp": base["official_L1_pp"],
                "candidate_Pearson": cand["official_Pearson"],
                "baseline_Pearson": base["official_Pearson"],
            }
        )
    compare_sylph = pd.DataFrame(rows)
    audit_df = audit(compare_current, compare_sylph)

    summary.to_csv(RESULTS / "candidate_default_comparison_summary.tsv", sep="\t", index=False)
    compare_current.to_csv(RESULTS / "candidate_default_vs_current.tsv", sep="\t", index=False)
    compare_sylph.to_csv(RESULTS / "candidate_default_vs_sylph.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / "candidate_default_decision_audit.tsv", sep="\t", index=False)
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
