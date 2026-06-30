#!/usr/bin/env python3
"""Score wrapper-emitted candidate-rescue profiles.

This validates the off-by-default `--candidate-rescue-switch
emitted-ani90-xny100-br01-af70` implementation in
`scripts/minco_profile_calibrated.py`. The switch replays the emitted-profile
part of the cross-panel zero-mass candidate rescue diagnostic. HMP raw-cache
candidates from the offline sweep are not available to the wrapper, so exact
agreement with the offline cross-panel rule is not required on those samples.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import score_adaptive_call_filter_wrapper_validation as adaptive


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
VALIDATION_DIR = Path("/tmp/minco_candidate_rescue_validation_20260629")
METHOD = "wrapper_candidate_rescue_emitted_ani90_xny100_br01_af70"
OFFLINE_METHOD = "cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass"
CURRENT_METHOD = "current_default"


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def path_map() -> dict[tuple[str, int], Path]:
    return {
        ("cami2_toy_mouse_gut", 5): VALIDATION_DIR / "toymouse_sample5.tsv",
        ("cami2_toy_mouse_gut", 6): VALIDATION_DIR / "toymouse_sample6.tsv",
        ("cami2_toy_mouse_gut", 7): VALIDATION_DIR / "toymouse_sample7.tsv",
        ("hmp_airskin_gtdb_source_abundance", 6): VALIDATION_DIR / "hmp_airskin_sample6.tsv",
        ("hmp_airskin_gtdb_source_abundance", 11): VALIDATION_DIR / "hmp_airskin_sample11.tsv",
        ("hmp_gastrooral_gtdb_source_abundance", 0): VALIDATION_DIR / "hmp_gastrooral_sample0.tsv",
        ("hmp_gastrooral_gtdb_source_abundance", 6): VALIDATION_DIR / "hmp_gastrooral_sample6.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 0): VALIDATION_DIR / "cami3_sample0.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 1): VALIDATION_DIR / "cami3_sample1.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 2): VALIDATION_DIR / "cami3_sample2.tsv",
    }


def profile_candidate_metadata(path: Path) -> dict[str, object]:
    cols = {
        "calibrated_call",
        "calibrated_abundance",
        "calibrated_abundance_raw",
        "candidate_rescue_added",
        "candidate_rescue_switch",
        "candidate_rescue_applied",
        "candidate_rescue_added_n",
        "candidate_rescue_input_call_n",
        "candidate_rescue_output_call_n",
        "candidate_rescue_zero_mass",
        "candidate_rescue_rule",
        "candidate_rescue_ani_min",
        "candidate_rescue_xny_min",
        "candidate_rescue_breadth_min",
        "candidate_rescue_real_af_min",
        "candidate_rescue_ani_source_rule",
    }
    raw = pd.read_csv(path, sep="\t", usecols=lambda col: col in cols, low_memory=False)
    if raw.empty:
        return {
            "profile": str(path),
            "candidate_rescue_switch": "",
            "candidate_rescue_applied": False,
            "candidate_rescue_added_n": 0,
            "candidate_rescue_zero_mass": False,
            "candidate_rescue_zero_mass_violations": 0,
        }
    first = raw.iloc[0]
    called = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    if "candidate_rescue_added" in raw.columns:
        added = raw["candidate_rescue_added"].astype(str).str.lower().isin({"true", "1", "yes"})
    else:
        added = pd.Series([False] * len(raw), index=raw.index)
    abundance_raw = pd.to_numeric(raw.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    abundance = pd.to_numeric(raw.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    zero_mass_violations = int(((abundance_raw.abs() > 1e-15) | (abundance.abs() > 1e-15))[added].sum())
    return {
        "profile": str(path),
        "candidate_rescue_switch": str(first.get("candidate_rescue_switch", "")),
        "candidate_rescue_applied": bool_value(first.get("candidate_rescue_applied", False)),
        "candidate_rescue_added_n": int(finite(first.get("candidate_rescue_added_n", added.sum()))),
        "candidate_rescue_observed_added_n": int(added.sum()),
        "candidate_rescue_input_call_n": int(finite(first.get("candidate_rescue_input_call_n", 0))),
        "candidate_rescue_output_call_n": int(finite(first.get("candidate_rescue_output_call_n", called.sum()))),
        "candidate_rescue_observed_call_n": int(called.sum()),
        "candidate_rescue_zero_mass": bool_value(first.get("candidate_rescue_zero_mass", False)),
        "candidate_rescue_zero_mass_violations": zero_mass_violations,
        "candidate_rescue_rule": str(first.get("candidate_rescue_rule", "")),
        "candidate_rescue_ani_min": finite(first.get("candidate_rescue_ani_min", 0.0)),
        "candidate_rescue_xny_min": finite(first.get("candidate_rescue_xny_min", 0.0)),
        "candidate_rescue_breadth_min": finite(first.get("candidate_rescue_breadth_min", 0.0)),
        "candidate_rescue_real_af_min": finite(first.get("candidate_rescue_real_af_min", 0.0)),
        "candidate_rescue_ani_source_rule": str(first.get("candidate_rescue_ani_source_rule", "")),
    }


def candidate_metadata(paths: dict[tuple[str, int], Path]) -> pd.DataFrame:
    rows = []
    for (panel, sample), path in paths.items():
        row = profile_candidate_metadata(path)
        row.update({"panel": panel, "sample": str(sample)})
        rows.append(row)
    return pd.DataFrame(rows)


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
                "samples": ",".join(map(str, sorted(sub["sample"].unique(), key=str))),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "official_L1_pp": sub[l1_col].mean(),
                "official_Pearson": sub[pearson_col].mean(),
                "rescued_rows": int(sub["candidate_rescue_observed_added_n"].sum()),
                "zero_mass_violations": int(sub["candidate_rescue_zero_mass_violations"].sum()),
            }
        )
    return pd.DataFrame(rows)


def compare_to_offline(scores: pd.DataFrame) -> pd.DataFrame:
    offline_path = Path("/tmp/minco_cross_panel_candidate_rescue_scores.tsv")
    if not offline_path.exists():
        return pd.DataFrame(
            [
                {
                    "panel": "",
                    "sample": "",
                    "comparison_status": "missing_offline_scores",
                    "offline_scores": str(offline_path),
                }
            ]
        )
    offline = pd.read_csv(offline_path, sep="\t")
    rows: list[dict[str, object]] = []
    for row in scores.itertuples(index=False):
        panel = str(getattr(row, "panel"))
        sample = str(getattr(row, "sample"))
        expected = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(str).eq(sample)
            & offline["method"].astype(str).eq(OFFLINE_METHOD)
        ]
        current = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(str).eq(sample)
            & offline["method"].astype(str).eq(CURRENT_METHOD)
        ]
        if expected.empty or current.empty:
            rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "comparison_status": "missing_expected_or_current",
                    "offline_scores": str(offline_path),
                }
            )
            continue
        expected_row = expected.iloc[0]
        current_row = current.iloc[0]
        l1_col = adaptive.official_l1_col(panel)
        pearson_col = adaptive.official_pearson_col(panel)
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "comparison_status": "compared",
                "offline_scores": str(offline_path),
                "candidate_rescue_added_n": int(getattr(row, "candidate_rescue_observed_added_n")),
                "wrapper_TP": int(getattr(row, "TP")),
                "wrapper_FP": int(getattr(row, "FP")),
                "wrapper_FN": int(getattr(row, "FN")),
                "offline_TP": int(expected_row["TP"]),
                "offline_FP": int(expected_row["FP"]),
                "offline_FN": int(expected_row["FN"]),
                "TP_delta_vs_offline": int(getattr(row, "TP")) - int(expected_row["TP"]),
                "FP_delta_vs_offline": int(getattr(row, "FP")) - int(expected_row["FP"]),
                "FN_delta_vs_offline": int(getattr(row, "FN")) - int(expected_row["FN"]),
                "F1_delta_vs_offline": finite(getattr(row, "F1")) - finite(expected_row["F1"]),
                "official_L1_delta_vs_offline_pp": finite(getattr(row, l1_col)) - finite(expected_row[l1_col]),
                "official_Pearson_delta_vs_offline": finite(getattr(row, pearson_col)) - finite(expected_row[pearson_col]),
                "F1_delta_vs_current": finite(getattr(row, "F1")) - finite(current_row["F1"]),
                "official_L1_delta_vs_current_pp": finite(getattr(row, l1_col)) - finite(current_row[l1_col]),
                "official_Pearson_delta_vs_current": finite(getattr(row, pearson_col)) - finite(current_row[pearson_col]),
                "offline_F1_delta_vs_current": finite(expected_row["F1"]) - finite(current_row["F1"]),
                "offline_L1_delta_vs_current_pp": finite(expected_row[l1_col]) - finite(current_row[l1_col]),
                "offline_Pearson_delta_vs_current": finite(expected_row[pearson_col]) - finite(current_row[pearson_col]),
            }
        )
    return pd.DataFrame(rows)


def audit(scores: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    compared = validation.loc[validation["comparison_status"].astype(str).eq("compared")].copy()
    rows = [
        {
            "metric": "wrapper_profiles",
            "value": int(scores.shape[0]),
            "evidence": str(VALIDATION_DIR),
            "decision": "cached_table_wrapper_replay",
        },
        {
            "metric": "rescued_rows",
            "value": int(scores["candidate_rescue_observed_added_n"].sum()),
            "evidence": "candidate_rescue_added column",
            "decision": "candidate_switch_behavior",
        },
        {
            "metric": "zero_mass_violations",
            "value": int(scores["candidate_rescue_zero_mass_violations"].sum()),
            "evidence": "calibrated_abundance_raw and calibrated_abundance for candidate_rescue_added rows",
            "decision": "pass" if int(scores["candidate_rescue_zero_mass_violations"].sum()) == 0 else "fail",
        },
    ]
    if compared.empty:
        rows.append(
            {
                "metric": "offline_comparison",
                "value": "missing",
                "evidence": "/tmp/minco_cross_panel_candidate_rescue_scores.tsv",
                "decision": "diagnostic_incomplete",
            }
        )
    else:
        max_count_delta = int(
            compared[["TP_delta_vs_offline", "FP_delta_vs_offline", "FN_delta_vs_offline"]]
            .abs()
            .max()
            .max()
        )
        max_f1_delta = float(compared["F1_delta_vs_offline"].abs().max())
        rows.extend(
            [
                {
                    "metric": "max_abs_call_delta_vs_offline_cross_panel_rule",
                    "value": f"counts={max_count_delta};F1={max_f1_delta:.6g}",
                    "evidence": "candidate_rescue_wrapper_validation_vs_offline.tsv",
                    "decision": "expected_partial_replay_review" if max_count_delta else "wrapper_matches_offline_subset",
                },
                {
                    "metric": "mean_wrapper_delta_vs_current",
                    "value": (
                        f"F1={compared['F1_delta_vs_current'].mean():.6f};"
                        f"L1_pp={compared['official_L1_delta_vs_current_pp'].mean():.6f};"
                        f"Pearson={compared['official_Pearson_delta_vs_current'].mean():.6f}"
                    ),
                    "evidence": "candidate_rescue_wrapper_validation_vs_offline.tsv",
                    "decision": "wrapper_candidate_effect_on_validation_subset",
                },
                {
                    "metric": "mean_offline_delta_vs_current",
                    "value": (
                        f"F1={compared['offline_F1_delta_vs_current'].mean():.6f};"
                        f"L1_pp={compared['offline_L1_delta_vs_current_pp'].mean():.6f};"
                        f"Pearson={compared['offline_Pearson_delta_vs_current'].mean():.6f}"
                    ),
                    "evidence": "candidate_rescue_wrapper_validation_vs_offline.tsv",
                    "decision": "offline_candidate_effect_on_validation_subset",
                },
            ]
        )
    rows.append(
        {
            "metric": "promotion_decision",
            "value": "historical_zero_mass_validation_not_final_preset",
            "evidence": "emitted-profile replay lacks HMP raw-cache side channel and abundance policy; final candidate preset is validated separately",
            "decision": "superseded_by_candidate_preset_replay",
        }
    )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    paths = path_map()
    missing = [path for path in paths.values() if not path.exists()]
    if missing:
        raise SystemExit("missing wrapper validation profiles:\n" + "\n".join(map(str, missing)))

    adaptive.METHOD = METHOD
    scores, _adaptive_metadata = adaptive.score_profiles(paths)
    metadata = candidate_metadata(paths)
    scores = scores.merge(metadata, on=["panel", "sample"], how="left")
    panel_summary = summarize(scores)
    validation = compare_to_offline(scores)
    audit_df = audit(scores, validation)

    scores.to_csv(RESULTS / "candidate_rescue_wrapper_validation_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(
        RESULTS / "candidate_rescue_wrapper_validation_panel_summary.tsv",
        sep="\t",
        index=False,
    )
    metadata.to_csv(
        RESULTS / "candidate_rescue_wrapper_validation_metadata.tsv",
        sep="\t",
        index=False,
    )
    validation.to_csv(
        RESULTS / "candidate_rescue_wrapper_validation_vs_offline.tsv",
        sep="\t",
        index=False,
    )
    audit_df.to_csv(
        RESULTS / "candidate_rescue_wrapper_validation_audit.tsv",
        sep="\t",
        index=False,
    )

    print(panel_summary.to_string(index=False))
    print("\nVALIDATION")
    print(validation.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
