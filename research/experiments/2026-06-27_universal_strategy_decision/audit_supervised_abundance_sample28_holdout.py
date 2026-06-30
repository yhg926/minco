#!/usr/bin/env python3
"""Post-hoc sample28 holdout for the supervised abundance calibrator.

This trains abundance models with HMP airskin sample28 excluded, then scores
sample28 as a single held-out sample. It is stricter than including sample28 in
the LOPO panel audit, but it is still a same-panel holdout rather than a new
independent dataset, so it cannot promote a default by itself.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import audit_supervised_abundance_calibrator as base


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"

HOLDOUT_PANEL = "hmp_airskin_gtdb_source_abundance"
HOLDOUT_SAMPLE = "28"
PREDEFINED_SPEC = "rf_log_leaf3"
PREDEFINED_BLEND = 0.75


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_hmp_sample28_sylph() -> dict[str, str]:
    path = RESULTS / "hmp_airskin28_r232_source_abundance_summary.tsv"
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        if row.get("method") == "sylph_gtdb_source_abundance":
            return row
    return {}


def add_current_deltas(scores: pd.DataFrame) -> pd.DataFrame:
    current = scores.loc[scores["method"].eq(base.CURRENT_METHOD)].iloc[0]
    out = scores.copy()
    out["delta_current_L1_pp"] = out["official_L1_pp"].astype(float) - float(current["official_L1_pp"])
    out["delta_current_Pearson"] = out["official_Pearson"].astype(float) - float(current["official_Pearson"])
    return out.sort_values(["delta_current_L1_pp", "method"], kind="mergesort")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    calls, _sample_meta = base.load_call_tables()
    features = base.feature_columns(calls)
    holdout_mask = (
        calls["panel"].astype(str).eq(HOLDOUT_PANEL)
        & calls["sample"].astype(str).eq(HOLDOUT_SAMPLE)
    )
    if int(holdout_mask.sum()) == 0:
        raise SystemExit(f"missing holdout {HOLDOUT_PANEL} sample {HOLDOUT_SAMPLE}")

    train_mask = ~holdout_mask
    X_all = calls[features].apply(pd.to_numeric, errors="coerce")
    y_all = np.log1p(calls["target_row_abundance"].to_numpy(dtype=float) * base.TARGET_SCALE)
    X_train = X_all.loc[train_mask]
    y_train = y_all[train_mask.to_numpy()]
    X_test = X_all.loc[holdout_mask]
    test_calls = calls.loc[holdout_mask].copy()
    current_test_raw = base.current_raw(test_calls).to_numpy(dtype=float)

    score_rows = [
        base.score_sample(
            HOLDOUT_PANEL,
            HOLDOUT_SAMPLE,
            base.CURRENT_METHOD,
            test_calls,
            current_test_raw,
        )
    ]

    for spec in base.MODEL_SPECS:
        model = base.make_model(spec)
        model.fit(X_train, y_train)
        pred_log = model.predict(X_test)
        pred_raw = np.maximum(np.expm1(pred_log) / base.TARGET_SCALE, 0.0)
        for blend in base.BLENDS:
            raw = blend * pred_raw + (1.0 - blend) * current_test_raw
            method = f"sample28_holdout_{spec.name}_blend{blend:g}"
            score_rows.append(base.score_sample(HOLDOUT_PANEL, HOLDOUT_SAMPLE, method, test_calls, raw))

    scores = add_current_deltas(pd.DataFrame(score_rows))
    scores.to_csv(RESULTS / "supervised_abundance_sample28_holdout_scores.tsv", sep="\t", index=False)

    current = scores.loc[scores["method"].eq(base.CURRENT_METHOD)].iloc[0]
    predefined_method = f"sample28_holdout_{PREDEFINED_SPEC}_blend{PREDEFINED_BLEND:g}"
    predefined = scores.loc[scores["method"].eq(predefined_method)].iloc[0]
    noncurrent = scores.loc[~scores["method"].eq(base.CURRENT_METHOD)]
    best = noncurrent.iloc[0]
    sylph = read_hmp_sample28_sylph()
    sylph_l1 = sylph.get("mean_L1_union_pp", "")
    sylph_f1 = sylph.get("mean_F1", "")
    predefined_decision = (
        "survives_same_panel_holdout"
        if float(predefined["delta_current_L1_pp"]) < 0.0
        else "does_not_improve_same_panel_holdout"
    )
    audit_rows = [
        {
            "metric": "holdout",
            "value": f"{HOLDOUT_PANEL}:{HOLDOUT_SAMPLE}",
            "evidence": "sample28 excluded from training rows",
            "decision": "same_panel_holdout_not_release_promotion",
        },
        {
            "metric": "train_rows",
            "value": int(train_mask.sum()),
            "evidence": f"holdout_rows={int(holdout_mask.sum())};features={len(features)}",
            "decision": "no_holdout_rows_in_training",
        },
        {
            "metric": "current_default",
            "value": (
                f"F1={float(current['F1']):.6f};"
                f"L1={float(current['official_L1_pp']):.6f};"
                f"Pearson={float(current['official_Pearson']):.6f}"
            ),
            "evidence": "supervised_abundance_sample28_holdout_scores.tsv",
            "decision": "baseline",
        },
        {
            "metric": "predefined_rf_blend075",
            "value": (
                f"F1={float(predefined['F1']):.6f};"
                f"L1={float(predefined['official_L1_pp']):.6f};"
                f"delta_L1={float(predefined['delta_current_L1_pp']):.6f};"
                f"delta_Pearson={float(predefined['delta_current_Pearson']):.6f}"
            ),
            "evidence": "supervised_abundance_sample28_holdout_scores.tsv",
            "decision": predefined_decision,
        },
        {
            "metric": "best_sample28_candidate",
            "value": (
                f"{best['method']};"
                f"L1={float(best['official_L1_pp']):.6f};"
                f"delta_L1={float(best['delta_current_L1_pp']):.6f};"
                f"delta_Pearson={float(best['delta_current_Pearson']):.6f}"
            ),
            "evidence": "supervised_abundance_sample28_holdout_scores.tsv",
            "decision": "diagnostic_only_selection_used_holdout",
        },
        {
            "metric": "sylph_reference",
            "value": f"F1={sylph_f1};L1={sylph_l1}",
            "evidence": "hmp_airskin28_r232_source_abundance_summary.tsv",
            "decision": "external_baseline_unchanged",
        },
        {
            "metric": "promotion_decision",
            "value": "candidate_requires_new_independent_dataset",
            "evidence": "supervised_abundance_sample28_holdout_audit.tsv",
            "decision": "do_not_promote_from_single_same_panel_holdout",
        },
    ]
    write_tsv(
        RESULTS / "supervised_abundance_sample28_holdout_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )

    print(pd.DataFrame(audit_rows).to_string(index=False))
    print("\nTOP SAMPLE28 HOLDOUT")
    print(scores.head(10).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
