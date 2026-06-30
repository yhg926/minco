#!/usr/bin/env python3
"""External exact-split stress test for the guarded feature allocator.

This is a cached-profile replay, not a raw-read rerun. It applies the
implemented `guarded-genus-hit-breadth-a002` abundance switch to the current
exact-split plant-associated and strainmadness profile tables. The call set is
fixed; only per-called-row abundance mass is changed.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

sys.path.insert(0, str(REPO_ROOT))
from scripts import minco_profile_calibrated as wrapper  # noqa: E402

PLANT_EXP = REPO_ROOT / "research/experiments/2026-06-25_cami2_plant_samples3_5_minco_vs_sylph"
sys.path.insert(0, str(PLANT_EXP))
import score_plant_samples3_5 as plant_score  # noqa: E402


EXACT_DIR = Path("/tmp/minco_exactsplit_universal_20260626")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
PLANT_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read")
STRAIN_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_strain_20260621/short_read")

METHOD_CURRENT = "current_exactsplit_abundance"
METHOD_SWITCH = "guarded_feature_allocator_exactsplit"
SWITCH = wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002

OUT_SCORES = RESULTS / "feature_allocator_external_exactsplit_scores.tsv"
OUT_SUMMARY = RESULTS / "feature_allocator_external_exactsplit_summary.tsv"
OUT_VALIDATION = RESULTS / "feature_allocator_external_exactsplit_validation.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_external_exactsplit_audit.tsv"
OUT_MD = NOTE_DIR / "ABUNDANCE_FEATURE_ALLOCATOR_EXTERNAL_EXACTSPLIT.md"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def truth_path(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return PLANT_TRUTH_DIR / f"taxonomic_profile_{sample}.txt"
    if dataset == "strainmadness":
        return STRAIN_TRUTH_DIR / f"taxonomic_profile_{sample}.txt"
    raise KeyError(dataset)


def profile_path(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return EXACT_DIR / f"plant_holdout{sample}_universal_exactsplit_strategy.tsv"
    if dataset == "strainmadness":
        return EXACT_DIR / f"strain{sample}_universal_exactsplit_strategy.tsv"
    raise KeyError(dataset)


def score_file(dataset: str, sample: int) -> Path:
    if dataset == "plant_holdout":
        return EXACT_DIR / f"plant_holdout{sample}_exactsplit_score.tsv"
    if dataset == "strainmadness":
        return EXACT_DIR / f"strain{sample}_exactsplit_score.tsv"
    raise KeyError(dataset)


def sample_key(dataset: str, sample: int) -> str:
    return f"plant_holdout{sample}" if dataset == "plant_holdout" else f"strain{sample}"


def load_profile(path: Path, scope_by_taxid: Mapping[str, str]) -> pd.DataFrame:
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


def current_raw(rows: pd.DataFrame) -> np.ndarray:
    col = "calibrated_abundance_raw" if "calibrated_abundance_raw" in rows.columns else "calibrated_abundance"
    raw = pd.to_numeric(rows[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return np.where(raw.to_numpy(dtype=float) > 0.0, raw.to_numpy(dtype=float), 0.0)


def score_calls(
    dataset: str,
    sample: int,
    rows: pd.DataFrame,
    mask: np.ndarray,
    method: str,
    abundance_col: str,
) -> dict[str, object]:
    calls = rows.loc[mask].copy()
    gold = plant_score.bacteria_gold(truth_path(dataset, sample))
    pred_taxids = set(calls["taxid"].astype(str))
    row = plant_score.score_method(sample, method, gold, pred_taxids, calls, abundance_col)
    row["dataset"] = dataset
    row["sample_key"] = sample_key(dataset, sample)
    row["profile_path"] = str(profile_path(dataset, sample))
    return row


def reference_current_score(dataset: str, sample: int) -> dict[str, object]:
    rows = pd.read_csv(score_file(dataset, sample), sep="\t")
    matches = rows.loc[rows["method"].astype(str).eq("minco_universal_exactsplit")]
    if matches.empty:
        raise ValueError(f"{score_file(dataset, sample)} lacks minco_universal_exactsplit")
    return matches.iloc[0].to_dict()


def evaluate_sample(
    dataset: str,
    sample: int,
    scope_by_taxid: Mapping[str, str],
    taxmap: Mapping[str, Mapping[str, str]],
) -> list[dict[str, object]]:
    rows = load_profile(profile_path(dataset, sample), scope_by_taxid)
    call = current_mask(rows)
    raw = current_raw(rows)
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        rows,
        call,
        raw,
        SWITCH,
        taxmap,
    )

    rows = rows.copy()
    rows["feature_allocator_raw"] = adjusted
    current = score_calls(dataset, sample, rows, call, METHOD_CURRENT, "calibrated_abundance")
    switched = score_calls(dataset, sample, rows, call, METHOD_SWITCH, "feature_allocator_raw")
    for row in [current, switched]:
        row.update(details)
    current["abundance_feature_allocator_applied"] = False
    current["abundance_feature_allocator_adjusted_rows_n"] = 0
    return [current, switched]


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
                "improved_samples_l1": int((sub["delta_current_l1"] < -1e-12).sum()),
                "worsened_samples_l1": int((sub["delta_current_l1"] > 1e-12).sum()),
                "mean_pearson": float(sub["tp_abundance_pearson"].mean()),
                "mean_delta_current_pearson": float(sub["delta_current_pearson"].mean()),
                "switched_samples": int(
                    sub["abundance_feature_allocator_applied"].astype(str).str.lower().isin({"true", "1"}).sum()
                ),
                "adjusted_rows": int(pd.to_numeric(
                    sub["abundance_feature_allocator_adjusted_rows_n"], errors="coerce"
                ).fillna(0.0).sum()),
                "mean_base_mass_multi_genus_frac": float(pd.to_numeric(
                    sub["abundance_feature_allocator_base_mass_multi_genus_frac"], errors="coerce"
                ).fillna(0.0).mean()),
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
                "F1_delta": finite_float(getattr(item, "F1")) - finite_float(ref["F1"]),
                "L1_delta": finite_float(getattr(item, "bacteria_scope_l1")) - finite_float(ref["L1"]),
                "Pearson_delta": finite_float(getattr(item, "tp_abundance_pearson"), float("nan"))
                - finite_float(ref["Pearson"], float("nan")),
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


def audit_rows(summary: pd.DataFrame, validation: pd.DataFrame) -> list[dict[str, object]]:
    switched = summary.loc[summary["method"].eq(METHOD_SWITCH)].copy()
    if switched.empty:
        return []
    mean_delta_l1 = float(switched["mean_delta_current_l1"].mean())
    worsened_datasets = int((switched["mean_delta_current_l1"] > 1e-12).sum())
    improved_datasets = int((switched["mean_delta_current_l1"] < -1e-12).sum())
    worsened_samples = int(switched["worsened_samples_l1"].sum())
    improved_samples = int(switched["improved_samples_l1"].sum())
    max_sample_worse = float(switched["max_worse_l1"].max())
    switched_samples = int(switched["switched_samples"].sum())
    adjusted_rows = int(switched["adjusted_rows"].sum())
    max_count_delta = int(validation[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max())
    max_f1_delta = finite_max_abs(validation["F1_delta"])
    max_l1_delta = finite_max_abs(validation["L1_delta"])
    max_pearson_delta = finite_max_abs(validation["Pearson_delta"])
    promoted = (
        worsened_datasets == 0
        and worsened_samples == 0
        and improved_samples > 0
        and max_count_delta == 0
        and max_f1_delta < 1e-9
    )
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
            "evidence": str(OUT_VALIDATION.relative_to(NOTE_DIR)),
            "decision": "matches_cached_exactsplit_scores"
            if max_count_delta == 0 and max_f1_delta < 1e-9
            else "baseline_mismatch_review",
        },
        {
            "metric": "candidate_effect",
            "value": (
                f"mean_delta_l1={mean_delta_l1:.6f};"
                f"improved_datasets={improved_datasets};worsened_datasets={worsened_datasets};"
                f"improved_samples={improved_samples};worsened_samples={worsened_samples};"
                f"max_sample_worse={max_sample_worse:.6f};"
                f"switched_samples={switched_samples};adjusted_rows={adjusted_rows}"
            ),
            "evidence": str(OUT_SUMMARY.relative_to(NOTE_DIR)),
            "decision": "external_supports_candidate" if promoted else "external_stress_rejects_default",
        },
        {
            "metric": "promotion_decision",
            "value": (
                "candidate_requires_release_grade_external_validation"
                if promoted
                else "keep_feature_allocator_experimental_off_by_default"
            ),
            "evidence": str(OUT_AUDIT.relative_to(NOTE_DIR)),
            "decision": "current_default_unchanged",
        },
    ]


def write_markdown(summary: pd.DataFrame, audit: pd.DataFrame) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit.to_dict("records")}
    lines = [
        "# Feature Allocator External Exact-Split Stress Test",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-profile replay applies the implemented experimental",
        "`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002`",
        "to plant-associated samples 3-5 and strainmadness samples 0-2.",
        "The call set is fixed; only abundance mass is adjusted.",
        "",
        "These panels are diagnostic, not release-grade GTDB holdouts, because",
        "the scorers use local bacteria-scope CAMI truth namespaces.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in ["baseline_validation", "candidate_effect", "promotion_decision"]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Dataset | Method | Samples | Mean F1 | Mean L1 | Delta L1 | Worsened L1 samples | Switched | Adjusted rows |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.to_dict("records"):
        lines.append(
            "| {dataset} | {method} | {samples} | {f1:.6f} | {l1:.6f} | {delta:.6f} | {worse} | {switched} | {adjusted} |".format(
                dataset=row["dataset"],
                method=row["method"],
                samples=row["samples"],
                f1=float(row["mean_F1"]),
                l1=float(row["mean_l1"]),
                delta=float(row["mean_delta_current_l1"]),
                worse=int(row["worsened_samples_l1"]),
                switched=int(row["switched_samples"]),
                adjusted=int(row["adjusted_rows"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Keep the current default unchanged.",
            "- Leave the guarded feature allocator experimental/off by default.",
            "- It remains useful as a candidate for a future same-namespace",
            "  release-grade holdout because this replay tests only diagnostic",
            "  external exact-split panels.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SCORES.relative_to(NOTE_DIR)}`",
            f"- `{OUT_SUMMARY.relative_to(NOTE_DIR)}`",
            f"- `{OUT_VALIDATION.relative_to(NOTE_DIR)}`",
            f"- `{OUT_AUDIT.relative_to(NOTE_DIR)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


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
    rows: list[dict[str, object]] = []
    for dataset, sample in samples:
        rows.extend(evaluate_sample(dataset, sample, scope_by_taxid, taxmap))

    scores = pd.DataFrame(rows)
    summary = summarize(scores)
    validation = pd.DataFrame(validation_rows(scores))
    audit = pd.DataFrame(audit_rows(summary, validation))

    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    validation.to_csv(OUT_VALIDATION, sep="\t", index=False)
    audit.to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(summary, audit)

    print(audit.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
