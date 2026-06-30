#!/usr/bin/env python3
"""Validate accession-level raw-side rescue profile rows on HMP gastrooral.

This focused replay targets the two HMP gastrooral samples where emitted-profile
candidate rescue missed raw-cache true positives. It writes small postprocessed
profiles under /tmp, scores them in the same GTDB source-abundance namespace,
and compares the call set to the offline cross-panel raw-cache rescue rule.
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
OUT_DIR = Path("/tmp/minco_raw_side_candidate_surface_hmp_gastrooral_20260629")
PANEL = "hmp_gastrooral_gtdb_source_abundance"
SAMPLES = [0, 6]
METHOD = "raw_side_candidate_surface_gastrooral"
OFFLINE_METHOD = "cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass"


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


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
    added = df["raw_side_rescue_added"].astype(str).str.lower().isin({"true", "1", "yes"})
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
                    "panel": PANEL,
                    "sample": "",
                    "comparison_status": "missing_offline_scores",
                    "offline_scores": str(offline_path),
                }
            ]
        )
    offline = pd.read_csv(offline_path, sep="\t")
    rows: list[dict[str, object]] = []
    for row in scores.itertuples(index=False):
        sample = int(getattr(row, "sample"))
        expected = offline.loc[
            offline["panel"].astype(str).eq(PANEL)
            & offline["sample"].astype(int).eq(sample)
            & offline["method"].astype(str).eq(OFFLINE_METHOD)
        ]
        current = offline.loc[
            offline["panel"].astype(str).eq(PANEL)
            & offline["sample"].astype(int).eq(sample)
            & offline["method"].astype(str).eq("current_default")
        ]
        if expected.empty or current.empty:
            rows.append(
                {
                    "panel": PANEL,
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
                "panel": PANEL,
                "sample": sample,
                "comparison_status": "compared",
                "offline_scores": str(offline_path),
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
    tp = int(scores["TP"].sum())
    fp = int(scores["FP"].sum())
    fn = int(scores["FN"].sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return pd.DataFrame(
        [
            {
                "panel": PANEL,
                "method": METHOD,
                "samples": ",".join(map(str, SAMPLES)),
                "mean_F1": scores["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": scores["L1_union_pp"].mean(),
                "mean_Pearson_union": scores["Pearson_union"].mean(),
                "added_rows": int(scores["added_rows"].sum()),
                "zero_mass_violations": int(scores["zero_mass_violations"].sum()),
            }
        ]
    )


def audit(summary: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    compared = validation.loc[validation["comparison_status"].astype(str).eq("compared")]
    if compared.empty:
        match = False
        max_count_delta = ""
        max_f1_delta = ""
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
    return pd.DataFrame(
        [
            {
                "metric": "profiles_replayed",
                "value": len(SAMPLES),
                "evidence": str(OUT_DIR),
                "decision": "focused_gap_replay",
            },
            {
                "metric": "added_rows",
                "value": int(summary["added_rows"].iloc[0]),
                "evidence": "raw_side_rescue_added rows",
                "decision": "candidate_surface_behavior",
            },
            {
                "metric": "zero_mass_violations",
                "value": int(summary["zero_mass_violations"].iloc[0]),
                "evidence": "calibrated_abundance and calibrated_abundance_raw on added rows",
                "decision": "pass" if int(summary["zero_mass_violations"].iloc[0]) == 0 else "fail",
            },
            {
                "metric": "max_delta_vs_offline_raw_cache_rule",
                "value": f"counts={max_count_delta};F1={max_f1_delta}",
                "evidence": "raw_side_candidate_surface_gastrooral_vs_offline.tsv",
                "decision": "matches_offline_rule" if match else "review_mismatch",
            },
            {
                "metric": "promotion_decision",
                "value": "experimental_not_default",
                "evidence": "focused HMP gastrooral replay only; needs integrated raw-side candidate surface",
                "decision": "do_not_change_current_default",
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mapper = raw_audit.RawMapper()
    apply_rows = []
    score_rows = []
    for sample in SAMPLES:
        cfg = decomp.PANELS[PANEL]
        profile = Path(cfg["minco"](sample))
        raw_best = raw_rescue.cache_path(PANEL, sample)
        out = OUT_DIR / f"hmp_gastrooral_sample{sample}.raw_side_rescue.tsv"
        apply_rows.append(
            apply_rescue.apply_raw_side_rescue(profile, raw_best, PANEL, sample, out, mapper)
        )
        score = score_profile(PANEL, sample, out, mapper)
        score.update(profile_zero_mass_audit(out))
        score["profile"] = str(out)
        score_rows.append(score)
    scores = pd.DataFrame(score_rows)
    apply_df = pd.DataFrame(apply_rows)
    summary = summarize(scores)
    validation = compare_to_offline(scores)
    audit_df = audit(summary, validation)

    apply_df.to_csv(RESULTS / "raw_side_candidate_surface_gastrooral_apply.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "raw_side_candidate_surface_gastrooral_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "raw_side_candidate_surface_gastrooral_summary.tsv", sep="\t", index=False)
    validation.to_csv(
        RESULTS / "raw_side_candidate_surface_gastrooral_vs_offline.tsv",
        sep="\t",
        index=False,
    )
    audit_df.to_csv(RESULTS / "raw_side_candidate_surface_gastrooral_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nVALIDATION")
    print(validation.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
