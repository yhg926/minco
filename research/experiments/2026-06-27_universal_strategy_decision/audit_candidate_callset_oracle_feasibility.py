#!/usr/bin/env python3
"""Audit whether selected-candidate calls need allocation or call recovery first.

This diagnostic keeps the selected candidate preset call set fixed and assigns
truth-aware abundance to detected true species. It is an oracle bound, not an
implementable strategy. The point is to decide whether the next default work
should focus on abundance allocation over already called rows or on recovering
missing calls before abundance can match the reference comparator.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_abundance_oracle_bounds as oracle_bounds
import decompose_abundance_errors as decomp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

SELECTED_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"
SYLPH_SCORES = RESULTS / "abundance_oracle_bounds_scores.tsv"

OUT_SAMPLE = RESULTS / "candidate_callset_oracle_feasibility.tsv"
OUT_SUMMARY = RESULTS / "candidate_callset_oracle_feasibility_summary.tsv"
OUT_AUDIT = RESULTS / "candidate_callset_oracle_feasibility_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_NEXT_TARGET.md"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def official_l1_kind(panel: str) -> str:
    return "truth_only" if panel == "cami2_toy_mouse_gut" else "union"


def official_l1(row: Mapping[str, object], panel: str) -> float:
    if official_l1_kind(panel) == "truth_only":
        return finite_float(row.get("L1_truth_only_pp"))
    return finite_float(row.get("L1_union_pp"))


def official_pearson(row: Mapping[str, object], panel: str) -> float:
    if official_l1_kind(panel) == "truth_only":
        return finite_float(row.get("Pearson_truth_only"), float("nan"))
    return finite_float(row.get("Pearson_union"), float("nan"))


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) < 2:
        return float("nan")
    left_series = pd.Series(left, dtype=float)
    right_series = pd.Series(right, dtype=float)
    if left_series.nunique(dropna=True) < 2 or right_series.nunique(dropna=True) < 2:
        return float("nan")
    return float(left_series.corr(right_series, method="pearson"))


def metric_vector(
    panel: str,
    truth: Mapping[str, float],
    call_species: set[str],
    abundance: Mapping[str, float],
) -> dict[str, float]:
    truth_norm = decomp.normalize(truth)
    pred_species = set(call_species)
    truth_species = set(truth_norm)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species

    union_species = sorted(truth_species | pred_species)
    truth_order = sorted(truth_species)
    union_true = [truth_norm.get(species, 0.0) for species in union_species]
    union_pred = [abundance.get(species, 0.0) for species in union_species]
    truth_true = [truth_norm.get(species, 0.0) for species in truth_order]
    truth_pred = [abundance.get(species, 0.0) for species in truth_order]

    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    union_l1 = sum(abs(p - t) for p, t in zip(union_pred, union_true)) * 100.0
    truth_l1 = sum(abs(p - t) for p, t in zip(truth_pred, truth_true)) * 100.0
    matched_abs = sum(abs(abundance.get(species, 0.0) - truth_norm.get(species, 0.0)) for species in tp) * 100.0
    missing = sum(truth_norm.get(species, 0.0) for species in fn) * 100.0
    extra = sum(abundance.get(species, 0.0) for species in fp) * 100.0
    detected_truth = sum(truth_norm.get(species, 0.0) for species in tp) * 100.0

    return {
        "TP": float(len(tp)),
        "FP": float(len(fp)),
        "FN": float(len(fn)),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "union_L1_pp": union_l1,
        "truth_only_L1_pp": truth_l1,
        "official_L1_pp": truth_l1 if official_l1_kind(panel) == "truth_only" else union_l1,
        "Pearson_union": pearson(union_pred, union_true),
        "Pearson_truth_only": pearson(truth_pred, truth_true),
        "official_Pearson": pearson(truth_pred, truth_true)
        if official_l1_kind(panel) == "truth_only"
        else pearson(union_pred, union_true),
        "matched_abs_error_pp": matched_abs,
        "missing_truth_mass_pp": missing,
        "extra_pred_mass_pp": extra,
        "truth_mass_detected_pct": detected_truth,
    }


def selected_oracle_abundance(
    truth: Mapping[str, float],
    call_species: set[str],
) -> tuple[dict[str, float], float]:
    truth_norm = decomp.normalize(truth)
    detected = set(truth_norm) & set(call_species)
    detected_mass = sum(truth_norm.get(species, 0.0) for species in detected)
    if detected_mass <= 0.0:
        return {}, 0.0
    return (
        {species: truth_norm.get(species, 0.0) / detected_mass for species in detected},
        detected_mass * 100.0,
    )


def load_sylph_rows() -> dict[tuple[str, str], dict[str, object]]:
    rows = pd.read_csv(SYLPH_SCORES, sep="\t")
    rows = rows.loc[rows["method"].astype(str).str.startswith("sylph")].copy()
    out: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows.to_dict("records"):
        out[(str(row["panel"]), str(row["sample"]))] = row
    return out


def classify(
    selected_minus_sylph_l1: float,
    oracle_minus_sylph_l1: float,
    selected_minus_sylph_f1: float,
    truth_mass_detected_pct: float,
) -> str:
    if selected_minus_sylph_l1 <= 0.0:
        return "selected_l1_already_not_worse"
    if oracle_minus_sylph_l1 <= 0.0:
        return "allocation_only_can_close_l1_gap"
    if selected_minus_sylph_f1 < 0.0 or truth_mass_detected_pct < 98.0:
        return "call_recovery_required_before_l1_can_match"
    return "allocation_plus_residual_call_or_fp_control_required"


def build_rows() -> tuple[list[dict[str, object]], float]:
    selected = pd.read_csv(SELECTED_SCORES, sep="\t")
    sylph = load_sylph_rows()
    loader = oracle_bounds.PanelLoader()
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
        selected_abundance, _rule = loader.minco_pred(
            panel,
            sample,
            Path(str(selected_row["profile"])),
            str(cfg["minco_collapse"]),
        )
        call_species = set(selected_abundance)
        baseline_metrics = metric_vector(panel, truth, call_species, selected_abundance)
        selected_l1 = official_l1(selected_row, panel)
        selected_pearson = official_pearson(selected_row, panel)
        max_validation_delta = max(
            max_validation_delta,
            abs(selected_l1 - baseline_metrics["official_L1_pp"]),
            abs(finite_float(selected_row["F1"]) - baseline_metrics["F1"]),
        )

        oracle_abundance, detected_mass = selected_oracle_abundance(truth, call_species)
        oracle_metrics = metric_vector(panel, truth, call_species, oracle_abundance)
        sylph_row = sylph.get((panel, str(sample)))
        if sylph_row is None:
            raise ValueError(f"missing Sylph score row for {panel} sample {sample}")
        sylph_l1 = finite_float(sylph_row.get("official_L1_pp"))
        sylph_f1 = finite_float(sylph_row.get("F1"))
        sylph_pearson = finite_float(sylph_row.get("official_Pearson"), float("nan"))
        selected_f1 = finite_float(selected_row["F1"])
        oracle_l1 = oracle_metrics["official_L1_pp"]
        selected_minus_sylph_l1 = selected_l1 - sylph_l1
        oracle_minus_sylph_l1 = oracle_l1 - sylph_l1
        selected_minus_sylph_f1 = selected_f1 - sylph_f1

        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "selected_F1": selected_f1,
                "sylph_F1": sylph_f1,
                "selected_minus_sylph_F1": selected_minus_sylph_f1,
                "selected_L1_pp": selected_l1,
                "sylph_L1_pp": sylph_l1,
                "selected_minus_sylph_L1_pp": selected_minus_sylph_l1,
                "selected_Pearson": selected_pearson,
                "sylph_Pearson": sylph_pearson,
                "selected_minus_sylph_Pearson": selected_pearson - sylph_pearson,
                "selected_oracle_L1_pp": oracle_l1,
                "selected_oracle_minus_sylph_L1_pp": oracle_minus_sylph_l1,
                "allocation_headroom_L1_pp": selected_l1 - oracle_l1,
                "truth_mass_detected_pct": detected_mass,
                "selected_TP": int(round(baseline_metrics["TP"])),
                "selected_FP": int(round(baseline_metrics["FP"])),
                "selected_FN": int(round(baseline_metrics["FN"])),
                "selected_missing_truth_mass_pp": baseline_metrics["missing_truth_mass_pp"],
                "selected_extra_pred_mass_pp": baseline_metrics["extra_pred_mass_pp"],
                "selected_matched_abs_error_pp": baseline_metrics["matched_abs_error_pp"],
                "selected_oracle_missing_truth_mass_pp": oracle_metrics["missing_truth_mass_pp"],
                "feasibility_class": classify(
                    selected_minus_sylph_l1,
                    oracle_minus_sylph_l1,
                    selected_minus_sylph_f1,
                    detected_mass,
                ),
            }
        )

    return rows, max_validation_delta


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    df = pd.DataFrame(rows)
    out: list[dict[str, object]] = []
    for panel, sub in df.groupby("panel", sort=True):
        counts = sub["feasibility_class"].value_counts().to_dict()
        out.append(
            {
                "panel": panel,
                "sample_count": int(len(sub)),
                "mean_selected_minus_sylph_L1_pp": float(sub["selected_minus_sylph_L1_pp"].mean()),
                "mean_selected_oracle_minus_sylph_L1_pp": float(
                    sub["selected_oracle_minus_sylph_L1_pp"].mean()
                ),
                "mean_allocation_headroom_L1_pp": float(sub["allocation_headroom_L1_pp"].mean()),
                "mean_truth_mass_detected_pct": float(sub["truth_mass_detected_pct"].mean()),
                "selected_l1_not_worse_samples": int(
                    counts.get("selected_l1_already_not_worse", 0)
                ),
                "allocation_only_can_close_samples": int(
                    counts.get("allocation_only_can_close_l1_gap", 0)
                ),
                "call_recovery_required_samples": int(
                    counts.get("call_recovery_required_before_l1_can_match", 0)
                ),
                "residual_required_samples": int(
                    counts.get("allocation_plus_residual_call_or_fp_control_required", 0)
                ),
            }
        )
    out.append(
        {
            "panel": "all",
            "sample_count": int(len(df)),
            "mean_selected_minus_sylph_L1_pp": float(df["selected_minus_sylph_L1_pp"].mean()),
            "mean_selected_oracle_minus_sylph_L1_pp": float(
                df["selected_oracle_minus_sylph_L1_pp"].mean()
            ),
            "mean_allocation_headroom_L1_pp": float(df["allocation_headroom_L1_pp"].mean()),
            "mean_truth_mass_detected_pct": float(df["truth_mass_detected_pct"].mean()),
            "selected_l1_not_worse_samples": int(
                (df["feasibility_class"] == "selected_l1_already_not_worse").sum()
            ),
            "allocation_only_can_close_samples": int(
                (df["feasibility_class"] == "allocation_only_can_close_l1_gap").sum()
            ),
            "call_recovery_required_samples": int(
                (df["feasibility_class"] == "call_recovery_required_before_l1_can_match").sum()
            ),
            "residual_required_samples": int(
                (
                    df["feasibility_class"]
                    == "allocation_plus_residual_call_or_fp_control_required"
                ).sum()
            ),
        }
    )
    return out


def write_markdown(summary_rows: list[dict[str, object]], audit_rows: list[dict[str, str]]) -> None:
    lines = [
        "# Abundance Next-Target Feasibility",
        "",
        "Date: 2026-06-29",
        "",
        "This note is generated from cached selected-candidate profile TSVs and",
        "cached truth/comparator scores. It keeps the selected candidate call set",
        "fixed, then computes a truth-aware abundance oracle over detected true",
        "species. The oracle is not an implementable strategy; it is a feasibility",
        "bound for deciding the next algorithmic target.",
        "",
        "## Summary",
        "",
        "| Panel | Samples | Selected-Sylph L1 pp | Oracle-Sylph L1 pp | Allocation headroom pp | Detected truth % | Already not worse | Allocation-only | Call recovery first | Residual work |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            "| {panel} | {sample_count} | {sel:.6f} | {oracle:.6f} | {headroom:.6f} | {truth:.6f} | {not_worse} | {alloc} | {call} | {resid} |".format(
                panel=row["panel"],
                sample_count=row["sample_count"],
                sel=float(row["mean_selected_minus_sylph_L1_pp"]),
                oracle=float(row["mean_selected_oracle_minus_sylph_L1_pp"]),
                headroom=float(row["mean_allocation_headroom_L1_pp"]),
                truth=float(row["mean_truth_mass_detected_pct"]),
                not_worse=row["selected_l1_not_worse_samples"],
                alloc=row["allocation_only_can_close_samples"],
                call=row["call_recovery_required_samples"],
                resid=row["residual_required_samples"],
            )
        )

    audit = {row["metric"]: row for row in audit_rows}
    lines.extend(
        [
            "",
            "## Decision",
            "",
            f"- Baseline validation max delta: `{audit['baseline_validation_max_abs_delta']['value']}`.",
            f"- Allocation-only can close the selected-candidate L1 gap on `{audit['allocation_only_can_close_samples']['value']}` of `{audit['evaluated_samples']['value']}` cached samples.",
            f"- Call recovery is required first on `{audit['call_recovery_required_samples']['value']}` cached samples.",
            f"- Residual work remains after allocation on `{audit['residual_required_samples']['value']}` cached samples.",
            "- Therefore the next default-strategy work should not be another",
            "  panel-mean-only abundance rescaling. It should combine a sample-safe",
            "  matched-call allocator with high-confidence call recovery for samples",
            "  whose detected truth mass is not enough.",
            "",
            "## Outputs",
            "",
            "- `results/candidate_callset_oracle_feasibility.tsv`",
            "- `results/candidate_callset_oracle_feasibility_summary.tsv`",
            "- `results/candidate_callset_oracle_feasibility_audit.tsv`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows, max_validation_delta = build_rows()
    summary_rows = summarize(rows)
    df = pd.DataFrame(rows)
    audit_rows = [
        {
            "metric": "evaluated_samples",
            "value": str(len(rows)),
            "evidence": str(SELECTED_SCORES.relative_to(EXP)),
            "decision": "selected_candidate_callset",
        },
        {
            "metric": "baseline_validation_max_abs_delta",
            "value": f"{max_validation_delta:.12g}",
            "evidence": "selected profile rescoring versus cached selected-candidate score TSV",
            "decision": "pass" if max_validation_delta < 1e-6 else "review",
        },
        {
            "metric": "allocation_only_can_close_samples",
            "value": str(int((df["feasibility_class"] == "allocation_only_can_close_l1_gap").sum())),
            "evidence": str(OUT_SAMPLE.relative_to(EXP)),
            "decision": "allocation_target",
        },
        {
            "metric": "call_recovery_required_samples",
            "value": str(int((df["feasibility_class"] == "call_recovery_required_before_l1_can_match").sum())),
            "evidence": str(OUT_SAMPLE.relative_to(EXP)),
            "decision": "call_recovery_target",
        },
        {
            "metric": "residual_required_samples",
            "value": str(
                int(
                    (
                        df["feasibility_class"]
                        == "allocation_plus_residual_call_or_fp_control_required"
                    ).sum()
                )
            ),
            "evidence": str(OUT_SAMPLE.relative_to(EXP)),
            "decision": "mixed_target",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "prioritize_allocator_plus_call_recovery",
        },
    ]

    pd.DataFrame(rows).to_csv(OUT_SAMPLE, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT_SUMMARY, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(summary_rows, audit_rows)


if __name__ == "__main__":
    main()
