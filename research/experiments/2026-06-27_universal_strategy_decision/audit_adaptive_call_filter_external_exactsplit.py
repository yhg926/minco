#!/usr/bin/env python3
"""Stress-test the adaptive call-filter candidate on external exact-split panels.

The four-panel adaptive call-filter audit found an in-sample/LOPO candidate:
apply `min_xny_ge_25` only when the sample-level median max-XnY support is at
least 253. This script applies that same output-only rule to cached current
exact-split plant-associated and strainmadness profiles. These panels are
diagnostic, not release-grade, but they are outside the four GTDB panels used
to select the rule.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

PLANT_EXP = REPO_ROOT / "research/experiments/2026-06-25_cami2_plant_samples3_5_minco_vs_sylph"
sys.path.insert(0, str(PLANT_EXP))
import score_plant_samples3_5 as plant_score  # noqa: E402


EXACT_DIR = Path("/tmp/minco_exactsplit_universal_20260626")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
PLANT_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read")
STRAIN_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_strain_20260621/short_read")

SWITCH = "lopo-min-xny25"
METHOD_CURRENT = "current_exactsplit_calls"
METHOD_FILTERED = "adaptive_call_filter_lopo_min_xny25"
MAX_XNY_MEDIAN_THRESHOLD = 253.0
MIN_XNY_THRESHOLD = 25.0


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def min_positive(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    out = np.maximum(left, right)
    both = (left > 0.0) & (right > 0.0)
    out[both] = np.minimum(left[both], right[both])
    return out


def profile_path(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return EXACT_DIR / f"plant_holdout{sample}_universal_exactsplit_strategy.tsv"
    if dataset == "strainmadness":
        return EXACT_DIR / f"strain{sample}_universal_exactsplit_strategy.tsv"
    raise KeyError(dataset)


def truth_path(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return PLANT_TRUTH_DIR / f"taxonomic_profile_{sample}.txt"
    if dataset == "strainmadness":
        return STRAIN_TRUTH_DIR / f"taxonomic_profile_{sample}.txt"
    raise KeyError(dataset)


def score_file(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return EXACT_DIR / f"plant_holdout{sample}_exactsplit_score.tsv"
    if dataset == "strainmadness":
        return EXACT_DIR / f"strain{sample}_exactsplit_score.tsv"
    raise KeyError(dataset)


def sample_key(dataset: str, sample: int) -> str:
    return f"plant_holdout{sample}" if dataset == "plant_holdout" else f"strain{sample}"


def load_profile(path: Path, scope_by_taxid: dict[str, str]) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t", low_memory=False)
    if "calibrated_call" not in rows.columns:
        raise ValueError(f"{path} lacks calibrated_call")
    rows["taxid"] = rows["taxid"].astype(str)
    rows = rows.loc[
        rows["taxid"].map(lambda taxid: scope_by_taxid.get(str(taxid), "other") == "bacteria")
    ].copy()
    if rows.empty:
        raise ValueError(f"{path} has no bacteria-scope rows")
    return rows


def current_mask(rows: pd.DataFrame) -> np.ndarray:
    return rows["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()


def adaptive_mask(rows: pd.DataFrame, call: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    s_xny = numeric(rows, "s_XnY_ctx_max").to_numpy(dtype=float)
    u_xny = numeric(rows, "u_XnY_ctx_max").to_numpy(dtype=float)
    max_xny = np.maximum(s_xny, u_xny)
    min_xny = min_positive(s_xny, u_xny)
    called_n = int(call.sum())
    if called_n == 0:
        return call.copy(), {
            "adaptive_call_filter_switch": SWITCH,
            "adaptive_call_filter_applied": False,
            "adaptive_call_filter_removed_n": 0,
            "adaptive_call_filter_input_call_n": 0,
            "adaptive_call_filter_output_call_n": 0,
            "adaptive_call_filter_max_xny_median": 0.0,
        }
    median_max_xny = float(np.median(max_xny[call]))
    applies = median_max_xny >= MAX_XNY_MEDIAN_THRESHOLD
    filtered = call & (min_xny >= MIN_XNY_THRESHOLD) if applies else call.copy()
    return filtered, {
        "adaptive_call_filter_switch": SWITCH,
        "adaptive_call_filter_applied": applies,
        "adaptive_call_filter_removed_n": int(called_n - filtered.sum()),
        "adaptive_call_filter_input_call_n": called_n,
        "adaptive_call_filter_output_call_n": int(filtered.sum()),
        "adaptive_call_filter_max_xny_median": median_max_xny,
    }


def score_mask(
    dataset: str,
    sample: int,
    rows: pd.DataFrame,
    mask: np.ndarray,
    method: str,
) -> dict[str, object]:
    calls = rows.loc[mask].copy()
    gold = plant_score.bacteria_gold(truth_path(dataset, sample))
    pred_taxids = set(calls["taxid"].astype(str))
    row = plant_score.score_method(
        sample,
        method,
        gold,
        pred_taxids,
        calls,
        "calibrated_abundance",
    )
    row["dataset"] = dataset
    row["sample_key"] = sample_key(dataset, sample)
    row["profile_path"] = str(profile_path(dataset, sample))
    return row


def reference_current_score(dataset: str, sample: int) -> dict[str, str]:
    path = score_file(dataset, sample)
    rows = pd.read_csv(path, sep="\t")
    matches = rows.loc[rows["method"].astype(str).eq("minco_universal_exactsplit")]
    if matches.empty:
        raise ValueError(f"{path} lacks minco_universal_exactsplit")
    return matches.iloc[0].to_dict()


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    current = scores.loc[
        scores["method"].eq(METHOD_CURRENT),
        ["dataset", "sample", "F1", "bacteria_scope_l1", "tp_abundance_pearson"],
    ].rename(
        columns={
            "F1": "current_F1",
            "bacteria_scope_l1": "current_l1",
            "tp_abundance_pearson": "current_pearson",
        }
    )
    scores = scores.merge(current, on=["dataset", "sample"], how="left")
    scores["delta_current_F1"] = scores["F1"] - scores["current_F1"]
    scores["delta_current_l1"] = scores["bacteria_scope_l1"] - scores["current_l1"]
    scores["delta_current_pearson"] = scores["tp_abundance_pearson"] - scores["current_pearson"]

    rows: list[dict[str, object]] = []
    for (dataset, method), sub in scores.groupby(["dataset", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "dataset": dataset,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "mean_F1": float(sub["F1"].mean()),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_delta_current_F1": float(sub["delta_current_F1"].mean()),
                "min_delta_current_F1": float(sub["delta_current_F1"].min()),
                "worsened_samples_F1": int((sub["delta_current_F1"] < -1e-12).sum()),
                "improved_samples_F1": int((sub["delta_current_F1"] > 1e-12).sum()),
                "mean_l1": float(sub["bacteria_scope_l1"].mean()),
                "mean_delta_current_l1": float(sub["delta_current_l1"].mean()),
                "max_worse_l1": float(sub["delta_current_l1"].max()),
                "mean_pearson": float(sub["tp_abundance_pearson"].mean()),
                "mean_delta_current_pearson": float(sub["delta_current_pearson"].mean()),
                "switched_samples": int(sub["adaptive_call_filter_applied"].astype(bool).sum())
                if "adaptive_call_filter_applied" in sub.columns
                else 0,
                "removed_rows": int(sub["adaptive_call_filter_removed_n"].sum())
                if "adaptive_call_filter_removed_n" in sub.columns
                else 0,
            }
        )
    return pd.DataFrame(rows)


def validation_rows(scores: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    current = scores.loc[scores["method"].eq(METHOD_CURRENT)].copy()
    for item in current.itertuples(index=False):
        ref = reference_current_score(str(getattr(item, "dataset")), int(getattr(item, "sample")))
        rows.append(
            {
                "dataset": getattr(item, "dataset"),
                "sample": getattr(item, "sample"),
                "TP_delta": int(getattr(item, "TP")) - int(float(ref["TP"])),
                "FP_delta": int(getattr(item, "FP")) - int(float(ref["FP"])),
                "FN_delta": int(getattr(item, "FN")) - int(float(ref["FN"])),
                "F1_delta": float(getattr(item, "F1")) - float(ref["F1"]),
                "L1_delta": float(getattr(item, "bacteria_scope_l1")) - float(ref["L1"]),
                "Pearson_delta": float(getattr(item, "tp_abundance_pearson")) - float(ref["Pearson"]),
                "source": str(score_file(str(getattr(item, "dataset")), int(getattr(item, "sample")))),
            }
        )
    return rows


def finite_max_abs(values: Iterable[object]) -> float:
    clean = []
    for value in values:
        try:
            out = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(out):
            clean.append(abs(out))
    return max(clean) if clean else 0.0


def audit(summary: pd.DataFrame, validation: pd.DataFrame) -> list[dict[str, object]]:
    filtered = summary.loc[summary["method"].eq(METHOD_FILTERED)].copy()
    if filtered.empty:
        return []
    mean_delta_f1 = float(filtered["mean_delta_current_F1"].mean())
    worsened_datasets = int((filtered["mean_delta_current_F1"] < -1e-12).sum())
    improved_datasets = int((filtered["mean_delta_current_F1"] > 1e-12).sum())
    switched = int(filtered["switched_samples"].sum())
    removed = int(filtered["removed_rows"].sum())
    max_count_delta = int(
        validation[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max()
    )
    max_f1_delta = finite_max_abs(validation["F1_delta"])
    max_l1_delta = finite_max_abs(validation["L1_delta"])
    max_pearson_delta = finite_max_abs(validation["Pearson_delta"])
    promoted = worsened_datasets == 0 and improved_datasets > 0
    return [
        {
            "metric": "external_panels",
            "value": "plant_holdout3-5,strainmadness0-2",
            "evidence": "cached exact-split current profile tables under /tmp/minco_exactsplit_universal_20260626",
            "decision": "diagnostic_nonrelease_stress_test",
        },
        {
            "metric": "baseline_validation",
            "value": (
                f"counts={max_count_delta};F1={max_f1_delta:.3g};"
                f"L1={max_l1_delta:.3g};Pearson={max_pearson_delta:.3g}"
            ),
            "evidence": "adaptive_call_filter_external_exactsplit_validation.tsv",
            "decision": "matches_cached_exactsplit_scores"
            if max_count_delta == 0 and max_f1_delta < 1e-9
            else "baseline_mismatch_review",
        },
        {
            "metric": "candidate_effect",
            "value": (
                f"mean_delta_F1={mean_delta_f1:.6f};"
                f"improved_datasets={improved_datasets};"
                f"worsened_datasets={worsened_datasets};"
                f"switched_samples={switched};removed_rows={removed}"
            ),
            "evidence": "adaptive_call_filter_external_exactsplit_summary.tsv",
            "decision": "external_supports_candidate" if promoted else "external_stress_rejects_default",
        },
        {
            "metric": "promotion_decision",
            "value": (
                "candidate_requires_release_grade_external_validation"
                if promoted
                else "reject_adaptive_call_filter_as_default"
            ),
            "evidence": "adaptive_call_filter_external_exactsplit_audit.tsv",
            "decision": "current_call_gate_remains_default",
        },
    ]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    samples = [
        ("plant_holdout", 3),
        ("plant_holdout", 4),
        ("plant_holdout", 5),
        ("strainmadness", 0),
        ("strainmadness", 1),
        ("strainmadness", 2),
    ]
    required = [TAXMAP]
    for dataset, sample in samples:
        required.extend([profile_path(dataset, sample), truth_path(dataset, sample), score_file(dataset, sample)])
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required external exact-split inputs:\n" + "\n".join(missing))

    taxmap = plant_score.parse_species_taxmap(TAXMAP)
    scope_by_taxid = plant_score.taxid_scope_map(taxmap)
    score_rows: list[dict[str, object]] = []
    for dataset, sample in samples:
        rows = load_profile(profile_path(dataset, sample), scope_by_taxid)
        call = current_mask(rows)
        filtered, details = adaptive_mask(rows, call)
        current_row = score_mask(dataset, sample, rows, call, METHOD_CURRENT)
        filtered_row = score_mask(dataset, sample, rows, filtered, METHOD_FILTERED)
        for row in [current_row, filtered_row]:
            row.update(details)
            if row["method"] == METHOD_CURRENT:
                row["adaptive_call_filter_applied"] = False
                row["adaptive_call_filter_removed_n"] = 0
                row["adaptive_call_filter_output_call_n"] = details["adaptive_call_filter_input_call_n"]
            score_rows.append(row)

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    validation = pd.DataFrame(validation_rows(scores))
    audit_df = pd.DataFrame(audit(summary, validation))

    scores.to_csv(RESULTS / "adaptive_call_filter_external_exactsplit_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "adaptive_call_filter_external_exactsplit_summary.tsv", sep="\t", index=False)
    validation.to_csv(
        RESULTS / "adaptive_call_filter_external_exactsplit_validation.tsv",
        sep="\t",
        index=False,
    )
    write_tsv(
        RESULTS / "adaptive_call_filter_external_exactsplit_audit.tsv",
        audit_df.to_dict("records"),
        ["metric", "value", "evidence", "decision"],
    )

    print(audit_df.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
