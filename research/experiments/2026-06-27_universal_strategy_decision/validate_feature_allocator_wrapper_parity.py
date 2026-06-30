#!/usr/bin/env python3
"""Validate the guarded feature allocator implementation on cached profiles.

This is a wrapper-function parity check, not a raw-read rerun. It loads the
selected-candidate profile TSVs, applies the implemented
`guarded-genus-hit-breadth-a002` abundance switch to the emitted called rows,
and scores the result with the same GTDB helpers used by the offline sweeps.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
sys.path.insert(0, str(ROOT))

import audit_candidate_callset_oracle_feasibility as feasibility
import decompose_abundance_errors as decomp
import sweep_selected_call_mass_transforms as base_sweep
from scripts import minco_profile_calibrated as wrapper


RESULTS = EXP / "results"

SELECTED_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"

OUT_SCORES = RESULTS / "feature_allocator_wrapper_parity_scores.tsv"
OUT_SUMMARY = RESULTS / "feature_allocator_wrapper_parity_summary.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_wrapper_parity_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_WRAPPER_PARITY.md"

METHOD_BASE = "current_selected_abundance"
METHOD_SWITCH = "wrapper_guarded_genus_hit_breadth_a002"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def official_l1(row: Mapping[str, object], panel: str) -> float:
    if panel == "cami2_toy_mouse_gut":
        return finite_float(row.get("L1_truth_only_pp"))
    return finite_float(row.get("L1_union_pp"))


def official_pearson(row: Mapping[str, object], panel: str) -> float:
    if panel == "cami2_toy_mouse_gut":
        return finite_float(row.get("Pearson_truth_only"), float("nan"))
    return finite_float(row.get("Pearson_union"), float("nan"))


def score_raw(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    calls: pd.DataFrame,
    raw_values: np.ndarray,
    details: Mapping[str, object] | None = None,
) -> dict[str, object]:
    abundance = base_sweep.collapse_prediction(calls, raw_values, panel)
    metrics = feasibility.metric_vector(panel, truth, set(abundance), abundance)
    out = {
        "panel": panel,
        "sample": sample,
        "method": method,
        "F1": metrics["F1"],
        "official_L1_pp": metrics["official_L1_pp"],
        "official_Pearson": metrics["official_Pearson"],
        "truth_mass_detected_pct": metrics["truth_mass_detected_pct"],
        "matched_abs_error_pp": metrics["matched_abs_error_pp"],
        "missing_truth_mass_pp": metrics["missing_truth_mass_pp"],
        "extra_pred_mass_pp": metrics["extra_pred_mass_pp"],
    }
    if details:
        for key, value in details.items():
            if key.startswith("abundance_feature_allocator_"):
                out[key] = value
    return out


def allocator_taxmap_from_profile(calls: pd.DataFrame) -> Mapping[str, Mapping[str, str]] | None:
    if "candidate_surface_taxmap" not in calls.columns:
        return None
    values = [
        str(value)
        for value in calls["candidate_surface_taxmap"].dropna().astype(str).unique()
        if str(value) and str(value).lower() != "nan"
    ]
    for value in values:
        path = Path(value)
        if path.is_file():
            return wrapper.parse_species_taxmap(path)
    return None


def build_scores() -> tuple[pd.DataFrame, float]:
    selected = pd.read_csv(SELECTED_SCORES, sep="\t")
    mapper = base_sweep.RowMapper()
    rows: list[dict[str, object]] = []
    max_validation_delta = 0.0
    for selected_row in selected.to_dict("records"):
        panel = str(selected_row["panel"])
        sample = int(float(selected_row["sample"]))
        cfg = decomp.PANELS[panel]
        truth = decomp.load_truth(
            Path(cfg["truth"](sample)),
            sample,
            str(cfg["truth_abundance_col"]),
        )
        calls = base_sweep.load_called_rows(selected_row, mapper)
        allocator_taxmap = allocator_taxmap_from_profile(calls)
        raw, _candidate_added = base_sweep.base_raw_values(calls)
        base_score = score_raw(panel, sample, METHOD_BASE, truth, calls, raw)
        target_l1 = official_l1(selected_row, panel)
        max_validation_delta = max(
            max_validation_delta,
            abs(base_score["official_L1_pp"] - target_l1),
            abs(base_score["F1"] - finite_float(selected_row.get("F1"))),
        )
        adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
            calls,
            np.ones(len(calls), dtype=bool),
            raw,
            wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002,
            allocator_taxmap,
        )
        switch_score = score_raw(panel, sample, METHOD_SWITCH, truth, calls, adjusted, details)
        rows.append(base_score)
        rows.append(switch_score)
    return pd.DataFrame(rows), max_validation_delta


def summarize(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = scores.loc[scores["method"].eq(METHOD_BASE)].set_index(["panel", "sample"])
    rows: list[dict[str, object]] = []
    for row in scores.to_dict("records"):
        if row["method"] == METHOD_BASE:
            row["delta_L1_pp"] = 0.0
            row["delta_Pearson"] = 0.0
        else:
            base_row = base.loc[(row["panel"], row["sample"])]
            row["delta_L1_pp"] = finite_float(row["official_L1_pp"]) - finite_float(
                base_row["official_L1_pp"]
            )
            row["delta_Pearson"] = finite_float(row["official_Pearson"], float("nan")) - finite_float(
                base_row["official_Pearson"], float("nan")
            )
        rows.append(row)
    scores = pd.DataFrame(rows)
    switch = scores.loc[scores["method"].eq(METHOD_SWITCH)].copy()
    panel_rows: list[dict[str, object]] = []
    for panel, sub in switch.groupby("panel", sort=True):
        panel_rows.append(
            {
                "panel": panel,
                "sample_count": int(len(sub)),
                "switched_samples": int(
                    sub["abundance_feature_allocator_applied"].astype(str).str.lower().isin({"true", "1"}).sum()
                ),
                "mean_L1_delta_pp": float(sub["delta_L1_pp"].mean()),
                "max_worse_L1_delta_pp": float(sub["delta_L1_pp"].max()),
                "improved_samples": int((sub["delta_L1_pp"] < -1e-9).sum()),
                "worsened_samples": int((sub["delta_L1_pp"] > 1e-9).sum()),
                "mean_Pearson_delta": float(sub["delta_Pearson"].mean()),
            }
        )
    panel = pd.DataFrame(panel_rows)
    overall = {
        "panel": "all",
        "sample_count": int(len(switch)),
        "switched_samples": int(
            switch["abundance_feature_allocator_applied"].astype(str).str.lower().isin({"true", "1"}).sum()
        ),
        "mean_L1_delta_pp": float(switch["delta_L1_pp"].mean()),
        "max_worse_L1_delta_pp": float(switch["delta_L1_pp"].max()),
        "improved_samples": int((switch["delta_L1_pp"] < -1e-9).sum()),
        "worsened_samples": int((switch["delta_L1_pp"] > 1e-9).sum()),
        "mean_Pearson_delta": float(switch["delta_Pearson"].mean()),
    }
    panel = pd.concat([panel, pd.DataFrame([overall])], ignore_index=True)
    return scores, panel


def audit_rows(summary: pd.DataFrame, max_validation_delta: float) -> list[dict[str, object]]:
    overall = summary.loc[summary["panel"].eq("all")].iloc[0]
    decision = (
        "candidate_needs_independent_holdout"
        if int(overall["worsened_samples"]) == 0 and int(overall["improved_samples"]) > 0
        else "do_not_promote_wrapper_feature_allocator"
    )
    return [
        {
            "metric": "baseline_validation_max_abs_delta",
            "value": f"{max_validation_delta:.12g}",
            "evidence": "current_selected_abundance vs cached selected score TSV",
            "decision": "pass" if max_validation_delta < 1e-6 else "review",
        },
        {
            "metric": "switched_samples",
            "value": int(overall["switched_samples"]),
            "evidence": str(OUT_SCORES.relative_to(EXP)),
            "decision": "wrapper_guard_applied",
        },
        {
            "metric": "mean_L1_delta_pp",
            "value": f"{float(overall['mean_L1_delta_pp']):.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "abundance_delta",
        },
        {
            "metric": "worsened_samples",
            "value": int(overall["worsened_samples"]),
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "sample_safety",
        },
        {
            "metric": "max_worse_L1_delta_pp",
            "value": f"{float(overall['max_worse_L1_delta_pp']):.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "sample_safety",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": "implemented switch on cached selected profiles",
            "decision": decision,
        },
    ]


def write_markdown(summary: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    audit_map = {str(row["metric"]): row for row in audit}
    lines = [
        "# Feature Allocator Wrapper Parity",
        "",
        "Date: 2026-06-29",
        "",
        "This note validates the implemented experimental switch",
        "`--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002`",
        "on cached selected-candidate profiles. It is not a raw-read rerun and",
        "does not change the default.",
        "The validator passes accession/taxmap labels to the allocator for",
        "genus grouping when the cached profile records a candidate-surface",
        "taxmap sidecar.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "baseline_validation_max_abs_delta",
        "switched_samples",
        "mean_L1_delta_pp",
        "worsened_samples",
        "max_worse_L1_delta_pp",
        "promotion_decision",
    ]:
        row = audit_map[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Panel Summary",
            "",
            "| Panel | Samples | Switched | Mean L1 delta pp | Max worse pp | Improved | Worsened |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.to_dict("records"):
        lines.append(
            "| {panel} | {sample_count} | {switched_samples} | {mean:.6f} | {max_worse:.6f} | {improved} | {worse} |".format(
                panel=row["panel"],
                sample_count=int(row["sample_count"]),
                switched_samples=int(row["switched_samples"]),
                mean=float(row["mean_L1_delta_pp"]),
                max_worse=float(row["max_worse_L1_delta_pp"]),
                improved=int(row["improved_samples"]),
                worse=int(row["worsened_samples"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Keep the selected default unchanged.",
            "- Treat this implemented switch as an independent-holdout candidate.",
            "- The next validation step is a same-namespace holdout replay before",
            "  considering any default promotion.",
            "",
            "## Outputs",
            "",
            "- `results/feature_allocator_wrapper_parity_scores.tsv`",
            "- `results/feature_allocator_wrapper_parity_summary.tsv`",
            "- `results/feature_allocator_wrapper_parity_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    scores, max_validation_delta = build_scores()
    scores, summary = summarize(scores)
    audit = audit_rows(summary, max_validation_delta)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    pd.DataFrame(audit).to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(summary, audit)
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
