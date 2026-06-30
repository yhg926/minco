#!/usr/bin/env python3
"""Validate accession-level raw-side rescue rows on all cached HMP panels.

This extends the focused gastrooral replay to all HMP source-abundance samples
used by the cross-panel rescue audit. It appends zero-mass accession-level
raw-side candidate rows to current MinCO profiles and checks whether the output
matches the offline raw-cache component of the selected cross-panel rule.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

import apply_raw_side_candidate_rescue_to_profile as apply_rescue
import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp
import sweep_hmp_raw_candidate_rescue as raw_rescue


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_raw_side_candidate_surface_hmp_20260629")
PANELS = [
    "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance",
]
METHOD = "raw_side_candidate_surface_hmp"
OFFLINE_METHOD = "cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass"


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def score_profile(panel: str, sample: int, profile: Path, mapper: raw_audit.RawMapper) -> dict[str, object]:
    cfg = decomp.PANELS[panel]
    truth = decomp.load_truth(
        Path(cfg["truth"](sample)),
        sample,
        str(cfg["truth_abundance_col"]),
    )
    pred = raw_rescue.load_current_pred_gtdb(
        panel,
        profile,
        str(cfg["minco_collapse"]),
        mapper,
    )
    return raw_rescue.score(
        panel,
        sample,
        METHOD,
        truth,
        set(pred),
        pred,
        "profile with appended zero-mass accession-level raw-side rescue rows",
    )


def profile_zero_mass_audit(profile: Path) -> dict[str, object]:
    df = pd.read_csv(
        profile,
        sep="\t",
        usecols=lambda col: col
        in {
            "raw_side_rescue_added",
            "calibrated_abundance",
            "calibrated_abundance_raw",
            "raw_side_rescue_gtdb_species",
        },
        low_memory=False,
    )
    added = bool_series(df["raw_side_rescue_added"])
    abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    return {
        "added_rows": int(added.sum()),
        "zero_mass_violations": int(((abundance.abs() > 1e-15) | (raw.abs() > 1e-15))[added].sum()),
        "added_species": ",".join(
            df.loc[added, "raw_side_rescue_gtdb_species"].astype(str).sort_values().tolist()
        ),
    }


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
        sample = int(getattr(row, "sample"))
        expected = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(int).eq(sample)
            & offline["method"].astype(str).eq(OFFLINE_METHOD)
        ]
        current = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(int).eq(sample)
            & offline["method"].astype(str).eq("current_default")
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
        exp = expected.iloc[0]
        cur = current.iloc[0]
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "comparison_status": "compared",
                "offline_scores": str(offline_path),
                "added_rows": int(getattr(row, "added_rows")),
                "TP_delta_vs_offline": int(getattr(row, "TP")) - int(exp["TP"]),
                "FP_delta_vs_offline": int(getattr(row, "FP")) - int(exp["FP"]),
                "FN_delta_vs_offline": int(getattr(row, "FN")) - int(exp["FN"]),
                "F1_delta_vs_offline": finite(getattr(row, "F1")) - finite(exp["F1"]),
                "L1_delta_vs_offline_pp": finite(getattr(row, "L1_union_pp")) - finite(exp["L1_union_pp"]),
                "Pearson_delta_vs_offline": finite(getattr(row, "Pearson_union")) - finite(exp["Pearson_union"]),
                "F1_delta_vs_current": finite(getattr(row, "F1")) - finite(cur["F1"]),
                "L1_delta_vs_current_pp": finite(getattr(row, "L1_union_pp")) - finite(cur["L1_union_pp"]),
                "Pearson_delta_vs_current": finite(getattr(row, "Pearson_union")) - finite(cur["Pearson_union"]),
                "offline_F1_delta_vs_current": finite(exp["F1"]) - finite(cur["F1"]),
            }
        )
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
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "added_rows": int(sub["added_rows"].sum()),
                "zero_mass_violations": int(sub["zero_mass_violations"].sum()),
            }
        )
    return pd.DataFrame(rows)


def overall(summary: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    compared = validation.loc[validation["comparison_status"].astype(str).eq("compared")].copy()
    if compared.empty:
        max_count_delta = ""
        max_f1_delta = ""
        match = False
        mean_delta_f1 = float("nan")
    else:
        max_count_delta_value = int(
            compared[["TP_delta_vs_offline", "FP_delta_vs_offline", "FN_delta_vs_offline"]]
            .abs()
            .max()
            .max()
        )
        max_f1_delta_value = float(compared["F1_delta_vs_offline"].abs().max())
        match = max_count_delta_value == 0 and max_f1_delta_value < 1e-12
        max_count_delta = str(max_count_delta_value)
        max_f1_delta = f"{max_f1_delta_value:.3g}"
        mean_delta_f1 = float(compared["F1_delta_vs_current"].mean())
    return pd.DataFrame(
        [
            {
                "metric": "profiles_replayed",
                "value": int(summary["sample_count"].sum()),
                "evidence": str(OUT_DIR),
                "decision": "all_cached_hmp_replay",
            },
            {
                "metric": "added_rows",
                "value": int(summary["added_rows"].sum()),
                "evidence": "raw_side_rescue_added rows",
                "decision": "candidate_surface_behavior",
            },
            {
                "metric": "zero_mass_violations",
                "value": int(summary["zero_mass_violations"].sum()),
                "evidence": "calibrated_abundance and calibrated_abundance_raw on added rows",
                "decision": "pass" if int(summary["zero_mass_violations"].sum()) == 0 else "fail",
            },
            {
                "metric": "max_delta_vs_offline_raw_cache_rule",
                "value": f"counts={max_count_delta};F1={max_f1_delta}",
                "evidence": "raw_side_candidate_surface_hmp_vs_offline.tsv",
                "decision": "matches_offline_rule" if match else "review_mismatch",
            },
            {
                "metric": "mean_delta_F1_vs_current",
                "value": f"{mean_delta_f1:.6f}",
                "evidence": "raw_side_candidate_surface_hmp_vs_offline.tsv",
                "decision": "candidate_surface_effect",
            },
            {
                "metric": "promotion_decision",
                "value": "experimental_not_default",
                "evidence": "postprocessed cached HMP profiles only; needs integrated raw-side candidate surface and nonzero abundance policy",
                "decision": "do_not_change_current_default",
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mapper = raw_audit.RawMapper()
    apply_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    for panel in PANELS:
        cfg = decomp.PANELS[panel]
        for sample_raw in cfg["samples"]:
            sample = int(sample_raw)
            profile = Path(cfg["minco"](sample))
            raw_best = raw_rescue.cache_path(panel, sample)
            out = OUT_DIR / f"{panel}.sample{sample}.raw_side_rescue.tsv"
            apply_rows.append(
                apply_rescue.apply_raw_side_rescue(profile, raw_best, panel, sample, out, mapper)
            )
            score = score_profile(panel, sample, out, mapper)
            score.update(profile_zero_mass_audit(out))
            score["profile"] = str(out)
            score_rows.append(score)
    apply_df = pd.DataFrame(apply_rows)
    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    validation = compare_to_offline(scores)
    audit = overall(summary, validation)

    apply_df.to_csv(RESULTS / "raw_side_candidate_surface_hmp_apply.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "raw_side_candidate_surface_hmp_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "raw_side_candidate_surface_hmp_summary.tsv", sep="\t", index=False)
    validation.to_csv(RESULTS / "raw_side_candidate_surface_hmp_vs_offline.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "raw_side_candidate_surface_hmp_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nVALIDATION")
    print(validation.to_string(index=False))
    print("\nAUDIT")
    print(audit.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
