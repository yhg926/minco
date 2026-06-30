#!/usr/bin/env python3
"""External exact-split stress test for supervised abundance candidates.

The four-panel supervised abundance audit found small RF/blend candidates. This
script trains the predefined RF model on the four-panel audit rows, then applies
several blend strengths to cached plant-associated and strainmadness
exact-split MinCO profiles. The call set is kept fixed; only per-called-row
abundance mass is changed. These panels are diagnostic rather than
release-grade, but they are outside the supervised abundance training panels.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import audit_supervised_abundance_calibrator as base


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

PLANT_EXP = REPO_ROOT / "research/experiments/2026-06-25_cami2_plant_samples3_5_minco_vs_sylph"
sys.path.insert(0, str(PLANT_EXP))
import score_plant_samples3_5 as plant_score  # noqa: E402


TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
EXACT_DIR = Path("/tmp/minco_exactsplit_universal_20260626")
PLANT_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read")
STRAIN_TRUTH_DIR = Path("/mnt/new3T/minco_cami2_strain_20260621/short_read")

RF_BLENDS = [0.25, 0.5, 0.75]
RF_SPEC = base.ModelSpec("rf_log_leaf3", "rf", 3.0)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def numeric_frame(rows: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=rows.index)
    for col in features:
        if col in rows.columns:
            out[col] = pd.to_numeric(rows[col], errors="coerce")
        else:
            out[col] = np.nan
    return out


def load_called_rows(path: Path, scope_by_taxid: dict[str, str]) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    if "calibrated_call" not in rows.columns:
        raise ValueError(f"{path} lacks calibrated_call")
    calls = rows.loc[rows["calibrated_call"].astype(str).str.lower().eq("true")].copy()
    calls["taxid"] = calls["taxid"].astype(str)
    calls = calls.loc[calls["taxid"].map(lambda taxid: scope_by_taxid.get(str(taxid), "other") == "bacteria")].copy()
    if calls.empty:
        raise ValueError(f"{path} has no called bacteria rows")
    return calls


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


def score_rows(dataset: str, sample: int, calls: pd.DataFrame, method: str, abundance_col: str) -> dict[str, object]:
    gold = plant_score.bacteria_gold(truth_path(dataset, sample))
    pred_taxids = set(calls["taxid"].astype(str))
    row = plant_score.score_method(sample, method, gold, pred_taxids, calls, abundance_col)
    row["dataset"] = dataset
    row["sample_key"] = f"plant_holdout{sample}" if dataset == "plant_holdout" else f"strain{sample}"
    return row


def blend_label(blend: float) -> str:
    return f"{blend:g}"


def train_rf_model() -> tuple[object, list[str]]:
    calls, _sample_meta = base.load_call_tables()
    features = base.feature_columns(calls)
    x_train = calls[features].apply(pd.to_numeric, errors="coerce")
    y_train = np.log1p(calls["target_row_abundance"].to_numpy(dtype=float) * base.TARGET_SCALE)
    model = base.make_model(RF_SPEC)
    model.fit(x_train, y_train)
    return model, features


def evaluate_panel(model, features: list[str], taxmap: dict[str, object]) -> tuple[pd.DataFrame, pd.DataFrame]:
    scope_by_taxid = plant_score.taxid_scope_map(taxmap)
    samples = [
        ("plant_holdout", 3),
        ("plant_holdout", 4),
        ("plant_holdout", 5),
        ("strainmadness", 0),
        ("strainmadness", 1),
        ("strainmadness", 2),
    ]
    records: list[dict[str, object]] = []
    for dataset, sample in samples:
        path = profile_path(dataset, sample)
        calls = load_called_rows(path, scope_by_taxid)
        current_raw = base.current_raw(calls).to_numpy(dtype=float)
        pred_log = model.predict(numeric_frame(calls, features))
        pred_raw = np.maximum(np.expm1(pred_log) / base.TARGET_SCALE, 0.0)
        for blend in RF_BLENDS:
            calls[f"rf_blend{blend_label(blend)}_raw"] = blend * pred_raw + (1.0 - blend) * current_raw
        calls["rf_full_raw"] = pred_raw
        calls["current_raw_for_external"] = current_raw

        records.append(score_rows(dataset, sample, calls, "current_calibrated_abundance", "current_raw_for_external"))
        for blend in RF_BLENDS:
            label = blend_label(blend)
            records.append(score_rows(dataset, sample, calls, f"rf_log_leaf3_blend{label}", f"rf_blend{label}_raw"))
        records.append(score_rows(dataset, sample, calls, "rf_log_leaf3_full", "rf_full_raw"))

    scores = pd.DataFrame(records)
    current = scores.loc[
        scores["method"].eq("current_calibrated_abundance"),
        ["dataset", "sample", "bacteria_scope_l1", "tp_abundance_pearson"],
    ].rename(
        columns={
            "bacteria_scope_l1": "current_bacteria_scope_l1",
            "tp_abundance_pearson": "current_tp_abundance_pearson",
        }
    )
    scores = scores.merge(current, on=["dataset", "sample"], how="left")
    scores["delta_current_l1"] = scores["bacteria_scope_l1"] - scores["current_bacteria_scope_l1"]
    scores["delta_current_tp_pearson"] = scores["tp_abundance_pearson"] - scores["current_tp_abundance_pearson"]

    summary_rows = []
    for (dataset, method), sub in scores.groupby(["dataset", "method"], sort=True):
        summary_rows.append(
            {
                "dataset": dataset,
                "method": method,
                "samples": ",".join(sorted(sub["sample"].astype(str).unique())),
                "mean_F1": float(sub["F1"].mean()),
                "mean_l1": float(sub["bacteria_scope_l1"].mean()),
                "mean_delta_current_l1": float(sub["delta_current_l1"].mean()),
                "max_sample_worse_l1": float(sub["delta_current_l1"].max()),
                "improved_samples": int((sub["delta_current_l1"] < -1e-9).sum()),
                "worsened_samples": int((sub["delta_current_l1"] > 1e-9).sum()),
                "mean_tp_abundance_pearson": float(sub["tp_abundance_pearson"].mean()),
                "mean_delta_current_tp_pearson": float(sub["delta_current_tp_pearson"].mean()),
            }
        )
    return scores, pd.DataFrame(summary_rows)


def audit(summary: pd.DataFrame) -> list[dict[str, object]]:
    candidates = summary.loc[summary["method"].str.startswith(f"{RF_SPEC.name}_blend")].copy()
    method_summary = []
    for method, sub in candidates.groupby("method", sort=True):
        method_summary.append(
            {
                "method": method,
                "mean_delta": float(sub["mean_delta_current_l1"].mean()),
                "max_worse": float(sub["mean_delta_current_l1"].max()),
                "improved": int((sub["mean_delta_current_l1"] < -1e-9).sum()),
                "worsened": int((sub["mean_delta_current_l1"] > 1e-9).sum()),
            }
        )
    strict_safe = [row for row in method_summary if row["worsened"] == 0 and row["improved"] > 0]
    best = min(method_summary, key=lambda row: (row["mean_delta"], row["max_worse"])) if method_summary else {}
    predefined = next((row for row in method_summary if row["method"] == "rf_log_leaf3_blend0.75"), {})
    decision = "external_stress_supports_candidate" if strict_safe else "external_stress_rejects_candidate_as_default"
    return [
        {
            "metric": "external_panels",
            "value": "plant_holdout3-5,strainmadness0-2",
            "evidence": "cached exact-split current profile tables under /tmp/minco_exactsplit_universal_20260626",
            "decision": "diagnostic_nonrelease_stress_test",
        },
        {
            "metric": "candidate_model",
            "value": RF_SPEC.name,
            "evidence": "trained on four-panel supervised abundance audit rows only",
            "decision": "predefined_no_external_tuning",
        },
        {
            "metric": "tested_blends",
            "value": ",".join(blend_label(blend) for blend in RF_BLENDS),
            "evidence": "supervised_abundance_external_exactsplit_summary.tsv",
            "decision": "fixed_before_external_scoring",
        },
        {
            "metric": "strict_external_safe_blends",
            "value": len(strict_safe),
            "evidence": "supervised_abundance_external_exactsplit_summary.tsv",
            "decision": "must_be_positive_for_default_candidate",
        },
        {
            "metric": "best_external_mean_delta",
            "value": (
                f"{best.get('method', 'NA')};mean_delta_l1={best.get('mean_delta', float('nan')):.6f};"
                f"worsened_datasets={best.get('worsened', 0)};max_dataset_worse={best.get('max_worse', float('nan')):.6f}"
                if best
                else "NA"
            ),
            "evidence": "supervised_abundance_external_exactsplit_summary.tsv",
            "decision": "diagnostic_only_nonrelease_external_panel",
        },
        {
            "metric": "candidate_external_delta",
            "value": (
                f"method={predefined.get('method', 'rf_log_leaf3_blend0.75')};"
                f"mean_delta_l1={predefined.get('mean_delta', float('nan')):.6f};"
                f"improved_datasets={predefined.get('improved', 0)};"
                f"worsened_datasets={predefined.get('worsened', 0)};"
                f"max_dataset_worse={predefined.get('max_worse', float('nan')):.6f}"
            ),
            "evidence": "supervised_abundance_external_exactsplit_summary.tsv",
            "decision": decision,
        },
        {
            "metric": "promotion_decision",
            "value": (
                "candidate_requires_release_grade_external_validation"
                if strict_safe
                else "reject_supervised_abundance_external_default"
            ),
            "evidence": "supervised_abundance_external_exactsplit_audit.tsv",
            "decision": "do_not_promote_from_nonrelease_external_panels",
        },
    ]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    required = [TAXMAP]
    for dataset, sample in [
        ("plant_holdout", 3),
        ("plant_holdout", 4),
        ("plant_holdout", 5),
        ("strainmadness", 0),
        ("strainmadness", 1),
        ("strainmadness", 2),
    ]:
        required.extend([profile_path(dataset, sample), truth_path(dataset, sample)])
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required external exact-split inputs:\n" + "\n".join(missing))

    taxmap = plant_score.parse_species_taxmap(TAXMAP)
    model, features = train_rf_model()
    scores, summary = evaluate_panel(model, features, taxmap)
    audit_rows = audit(summary)

    scores.to_csv(RESULTS / "supervised_abundance_external_exactsplit_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "supervised_abundance_external_exactsplit_summary.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "supervised_abundance_external_exactsplit_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )

    print(pd.DataFrame(audit_rows).to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
