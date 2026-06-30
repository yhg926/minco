#!/usr/bin/env python3
"""Leave-one-panel-out audit for supervised abundance calibration.

This diagnostic trains row-level abundance models from MinCO output features.
It keeps the species call set fixed and predicts only per-called-row abundance
mass. Panel labels, sample IDs, species names, and reference IDs are excluded
from model features; truth is used only for training and evaluation.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import score_cami3_gtdb_taxid_transfer as taxid_score
import sweep_cross_panel_abundance_variants as sweep


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

TARGET_SCALE = 1_000_000.0
CURRENT_METHOD = "current_calibrated_abundance"


@dataclass(frozen=True)
class ModelSpec:
    name: str
    kind: str
    alpha: float = 1.0


MODEL_SPECS = [
    ModelSpec("ridge_log_alpha1", "ridge", 1.0),
    ModelSpec("ridge_log_alpha10", "ridge", 10.0),
    ModelSpec("hgb_log_l2_0.1", "hgb", 0.1),
    ModelSpec("hgb_log_l2_1", "hgb", 1.0),
    ModelSpec("rf_log_leaf3", "rf", 3.0),
]
BLENDS = [1.0, 0.75, 0.5, 0.25]


TEXT_OR_ID_COLUMNS = {
    "taxid",
    "species_name",
    "profile_species_name",
    "gtdb_species",
    "u_best_accession",
    "u_best_ref",
    "u_best_ref_annotation",
    "s_best_accession",
    "s_best_ref",
    "s_best_ref_annotation",
    "profile_strategy",
    "train_pool",
    "scope",
    "adaptive_mode",
    "tail_rescue_abundance_rule",
    "low_extra_split_rescue_rule",
    "raw_unique_fallback_rule",
    "auto_exact_split_low_extra_mode",
    "auto_exact_split_guard_rule",
    "auto_exact_split_path",
    "auto_exact_split_source",
    "auto_exact_split_unavailable_reason",
    "abundance_rule",
    "abundance_genus_xny_quality",
    "calibrated_call",
}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def current_raw(calls: pd.DataFrame) -> pd.Series:
    raw = pd.to_numeric(calls.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    if float(raw.sum()) <= 0.0:
        raw = pd.to_numeric(calls.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    return raw.clip(lower=0.0)


def sample_current_method(panel: str) -> str:
    if panel == "cami2_toy_mouse_gut":
        return "minco_current_code_refresh"
    if panel == "hmp_airskin_gtdb_source_abundance":
        return "minco_current_default_gtdb_source_abundance_exact6_sample28"
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        return "minco_current_universal_gastrooral_source_abundance"
    if panel == "cami3_toy_human_gut_gtdb_source_readmap":
        return "minco_universal_autoexact_gtdb_source_readmap"
    raise KeyError(panel)


def load_truth_df(panel: str, sample: int) -> pd.DataFrame:
    if panel == "cami2_toy_mouse_gut":
        return sweep.toy_truth(sample)
    if panel == "hmp_airskin_gtdb_source_abundance":
        return sweep.hmp_truth(sample)
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        return sweep.hmp_gastro_truth(sample)
    if panel == "cami3_toy_human_gut_gtdb_source_readmap":
        return sweep.cami3_truth(sample)
    raise KeyError(panel)


def load_call_tables() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    toy_mod = sweep.decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    hmp_taxmap = sweep.hmp.parse_species_taxmap(sweep.hmp.TAXMAP)
    hmp_gastro_taxmap = sweep.hmp_gastro.parse_species_taxmap(sweep.hmp_gastro.TAXMAP)
    _wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _map_diag = sweep.cami3.build_transfer_maps()
    cami3_taxmap = sweep.cami3.parse_species_taxmap(sweep.cami3.TAXMAP)

    sample_records = [
        ("cami2_toy_mouse_gut", 5),
        ("cami2_toy_mouse_gut", 6),
        ("cami2_toy_mouse_gut", 7),
        ("hmp_airskin_gtdb_source_abundance", 6),
        ("hmp_airskin_gtdb_source_abundance", 11),
        ("hmp_airskin_gtdb_source_abundance", 28),
        ("hmp_gastrooral_gtdb_source_abundance", 0),
        ("hmp_gastrooral_gtdb_source_abundance", 6),
        ("cami3_toy_human_gut_gtdb_source_readmap", 0),
        ("cami3_toy_human_gut_gtdb_source_readmap", 1),
        ("cami3_toy_human_gut_gtdb_source_readmap", 2),
    ]

    frames = []
    sample_meta = []
    for panel, sample in sample_records:
        if panel == "cami2_toy_mouse_gut":
            calls, collapse, path = sweep.load_toy_calls(sample)
        elif panel == "hmp_airskin_gtdb_source_abundance":
            calls, collapse, path = sweep.load_hmp_calls(sample, hmp_taxmap, by_accession, by_core)
        elif panel == "hmp_gastrooral_gtdb_source_abundance":
            calls, collapse, path = sweep.load_hmp_gastro_calls(
                sample,
                hmp_gastro_taxmap,
                by_accession,
                by_core,
            )
        elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
            calls, collapse, path = sweep.load_cami3_calls(
                sample,
                cami3_taxid_to_species,
                cami3_name_to_species,
                cami3_taxmap,
                by_accession,
                by_core,
            )
        else:
            raise AssertionError(panel)

        calls = calls.copy()
        calls["panel"] = panel
        calls["sample"] = str(sample)
        calls["collapse_rule"] = collapse
        calls["profile"] = str(path)
        calls["row_id"] = [f"{panel}|{sample}|{idx}" for idx in range(len(calls))]
        raw = current_raw(calls)
        calls["current_raw_for_training"] = raw
        truth = load_truth_df(panel, sample)
        truth_map = {
            str(row.gtdb_species): finite(row.truth_abundance)
            for row in truth.itertuples(index=False)
        }
        species_raw_sum = calls.groupby("gtdb_species")["current_raw_for_training"].transform("sum")
        species_count = calls.groupby("gtdb_species")["current_raw_for_training"].transform("size")
        calls["truth_abundance"] = calls["gtdb_species"].map(truth_map).fillna(0.0)
        proportional = np.where(
            species_raw_sum.to_numpy(dtype=float) > 0.0,
            calls["current_raw_for_training"].to_numpy(dtype=float) / species_raw_sum.to_numpy(dtype=float),
            1.0 / species_count.to_numpy(dtype=float),
        )
        calls["target_row_abundance"] = calls["truth_abundance"].to_numpy(dtype=float) * proportional
        frames.append(calls)
        sample_meta.append(
            {
                "panel": panel,
                "sample": str(sample),
                "collapse_rule": collapse,
                "profile": str(path),
                "called_rows": int(len(calls)),
                "called_species": int(calls["gtdb_species"].nunique()),
                "truth_species": int(truth["gtdb_species"].nunique()),
            }
        )
    return pd.concat(frames, ignore_index=True, sort=False), sample_meta


def feature_columns(calls: pd.DataFrame) -> list[str]:
    columns = []
    for col in calls.columns:
        if col in TEXT_OR_ID_COLUMNS or col.startswith("Unnamed"):
            continue
        if col in {
            "panel",
            "sample",
            "collapse_rule",
            "profile",
            "row_id",
            "truth_abundance",
            "target_row_abundance",
        }:
            continue
        converted = pd.to_numeric(calls[col], errors="coerce")
        if converted.notna().sum() == 0 or converted.nunique(dropna=True) <= 1:
            continue
        columns.append(col)
    return columns


def make_model(spec: ModelSpec):
    if spec.kind == "ridge":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), Ridge(alpha=spec.alpha))
    if spec.kind == "hgb":
        return make_pipeline(
            SimpleImputer(strategy="median"),
            HistGradientBoostingRegressor(
                loss="squared_error",
                learning_rate=0.05,
                max_iter=200,
                max_leaf_nodes=8,
                l2_regularization=spec.alpha,
                random_state=17,
            ),
        )
    if spec.kind == "rf":
        return make_pipeline(
            SimpleImputer(strategy="median"),
            RandomForestRegressor(
                n_estimators=300,
                min_samples_leaf=int(spec.alpha),
                max_features=0.65,
                n_jobs=1,
                random_state=17,
            ),
        )
    raise ValueError(spec)


def collapse_prediction(sample_calls: pd.DataFrame, raw_values: np.ndarray) -> pd.DataFrame:
    collapse = str(sample_calls["collapse_rule"].iloc[0])
    pred = sweep.collapse_prediction(sample_calls, raw_values, collapse)
    return pred


def score_sample(panel: str, sample: str, method: str, calls: pd.DataFrame, raw_values: np.ndarray) -> dict[str, object]:
    pred = collapse_prediction(calls, raw_values)
    truth = load_truth_df(panel, int(sample))
    row = taxid_score.score_prediction(int(sample), method, pred, truth, {})
    row["panel"] = panel
    row["sample"] = sample
    row["official_L1_pp"] = (
        row["L1_truth_only_pp"] if panel == "cami2_toy_mouse_gut" else row["L1_union_pp"]
    )
    row["official_Pearson"] = (
        row["Pearson_truth_only"] if panel == "cami2_toy_mouse_gut" else row["Pearson_union"]
    )
    return row


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    current = scores.loc[
        scores["method"].eq(CURRENT_METHOD),
        ["panel", "sample", "official_L1_pp", "official_Pearson"],
    ].rename(
        columns={
            "official_L1_pp": "current_official_L1_pp",
            "official_Pearson": "current_official_Pearson",
        }
    )
    work = scores.merge(current, on=["panel", "sample"], how="left")
    work["delta_current_L1_pp"] = work["official_L1_pp"] - work["current_official_L1_pp"]
    work["delta_current_Pearson"] = work["official_Pearson"] - work["current_official_Pearson"]

    rows = []
    for (method, panel), sub in work.groupby(["method", "panel"], sort=True):
        rows.append(
            {
                "method": method,
                "panel": panel,
                "samples": ",".join(sorted(sub["sample"].astype(str).unique())),
                "mean_F1": float(sub["F1"].mean()),
                "mean_official_L1_pp": float(sub["official_L1_pp"].mean()),
                "mean_delta_current_L1_pp": float(sub["delta_current_L1_pp"].mean()),
                "max_sample_worse_L1_pp": float(sub["delta_current_L1_pp"].max()),
                "improved_samples": int((sub["delta_current_L1_pp"] < -1e-9).sum()),
                "worsened_samples": int((sub["delta_current_L1_pp"] > 1e-9).sum()),
                "mean_official_Pearson": float(sub["official_Pearson"].mean()),
                "mean_delta_current_Pearson": float(sub["delta_current_Pearson"].mean()),
            }
        )
    panel_summary = pd.DataFrame(rows)
    overall_rows = []
    for method, sub in panel_summary.groupby("method", sort=True):
        overall_rows.append(
            {
                "method": method,
                "panel_count": int(sub["panel"].nunique()),
                "samples": ",".join(sorted(work.loc[work["method"].eq(method), "sample"].astype(str).unique())),
                "mean_official_L1_pp": float(sub["mean_official_L1_pp"].mean()),
                "mean_delta_current_L1_pp": float(sub["mean_delta_current_L1_pp"].mean()),
                "max_panel_worse_L1_pp": float(sub["mean_delta_current_L1_pp"].max()),
                "improved_panels": int((sub["mean_delta_current_L1_pp"] < -1e-9).sum()),
                "worsened_panels": int((sub["mean_delta_current_L1_pp"] > 1e-9).sum()),
                "mean_official_Pearson": float(sub["mean_official_Pearson"].mean()),
                "mean_delta_current_Pearson": float(sub["mean_delta_current_Pearson"].mean()),
            }
        )
    return work, panel_summary, pd.DataFrame(overall_rows).sort_values(
        ["mean_delta_current_L1_pp", "max_panel_worse_L1_pp", "method"],
        kind="mergesort",
    )


def audit(overall: pd.DataFrame) -> list[dict[str, object]]:
    noncurrent = overall.loc[~overall["method"].eq(CURRENT_METHOD)].copy()
    strict = noncurrent.loc[
        noncurrent["worsened_panels"].eq(0) & noncurrent["improved_panels"].gt(0)
    ].copy()
    best = noncurrent.iloc[0]
    if strict.empty:
        decision_value = "reject_supervised_abundance_calibrator"
        decision = "current_calibrated_abundance_remains_default"
    else:
        best_strict = strict.sort_values(
            ["mean_delta_current_L1_pp", "max_panel_worse_L1_pp"],
            kind="mergesort",
        ).iloc[0]
        decision_value = (
            "candidate_requires_independent_holdout"
            if float(best_strict["mean_delta_current_L1_pp"]) < 0.0
            else "reject_supervised_abundance_calibrator"
        )
        decision = "do_not_promote_without_independent_holdout"
    return [
        {
            "metric": "model_specs",
            "value": len(MODEL_SPECS),
            "evidence": ",".join(spec.name for spec in MODEL_SPECS),
            "decision": "leave_one_panel_out_training",
        },
        {
            "metric": "blend_values",
            "value": len(BLENDS),
            "evidence": ",".join(map(str, BLENDS)),
            "decision": "blend_model_raw_with_current_raw",
        },
        {
            "metric": "strict_panel_safe_models",
            "value": int(strict.shape[0]),
            "evidence": "supervised_abundance_calibrator_overall.tsv",
            "decision": "must_be_positive_for_candidate_status",
        },
        {
            "metric": "best_model",
            "value": (
                f"{best['method']};mean_delta={best['mean_delta_current_L1_pp']:.6f};"
                f"worsened_panels={int(best['worsened_panels'])};"
                f"max_worse={best['max_panel_worse_L1_pp']:.6f}"
            ),
            "evidence": "supervised_abundance_calibrator_overall.tsv",
            "decision": "held_out_panel_result",
        },
        {
            "metric": "promotion_decision",
            "value": decision_value,
            "evidence": "supervised_abundance_calibrator_overall.tsv",
            "decision": decision,
        },
    ]


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    calls, sample_meta = load_call_tables()
    features = feature_columns(calls)
    score_rows = []

    for (panel, sample), sample_calls in calls.groupby(["panel", "sample"], sort=True):
        score_rows.append(
            score_sample(
                str(panel),
                str(sample),
                CURRENT_METHOD,
                sample_calls,
                current_raw(sample_calls).to_numpy(dtype=float),
            )
        )

    X_all = calls[features].apply(pd.to_numeric, errors="coerce")
    y_all = np.log1p(calls["target_row_abundance"].to_numpy(dtype=float) * TARGET_SCALE)
    panels = sorted(calls["panel"].astype(str).unique())

    for holdout_panel in panels:
        train_mask = ~calls["panel"].astype(str).eq(holdout_panel)
        test_mask = calls["panel"].astype(str).eq(holdout_panel)
        X_train = X_all.loc[train_mask]
        y_train = y_all[train_mask.to_numpy()]
        X_test = X_all.loc[test_mask]
        test_calls = calls.loc[test_mask].copy()
        current_test_raw = current_raw(test_calls).to_numpy(dtype=float)
        for spec in MODEL_SPECS:
            model = make_model(spec)
            model.fit(X_train, y_train)
            pred_log = model.predict(X_test)
            pred_raw = np.maximum(np.expm1(pred_log) / TARGET_SCALE, 0.0)
            for blend in BLENDS:
                raw = blend * pred_raw + (1.0 - blend) * current_test_raw
                method = f"lopo_{spec.name}_blend{blend:g}"
                for (panel, sample), sample_calls in test_calls.groupby(["panel", "sample"], sort=True):
                    idx = sample_calls.index
                    local_raw = raw[test_calls.index.get_indexer(idx)]
                    score_rows.append(score_sample(str(panel), str(sample), method, sample_calls, local_raw))

    scores = pd.DataFrame(score_rows)
    scores_with_delta, panel_summary, overall = summarize(scores)
    audit_rows = audit(overall)

    scores_with_delta.to_csv(RESULTS / "supervised_abundance_calibrator_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(RESULTS / "supervised_abundance_calibrator_panel_summary.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / "supervised_abundance_calibrator_overall.tsv", sep="\t", index=False)
    pd.DataFrame(sample_meta).to_csv(RESULTS / "supervised_abundance_calibrator_samples.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "supervised_abundance_calibrator_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )

    print(pd.DataFrame(audit_rows).to_string(index=False))
    print("\nTOP OVERALL")
    print(overall.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
