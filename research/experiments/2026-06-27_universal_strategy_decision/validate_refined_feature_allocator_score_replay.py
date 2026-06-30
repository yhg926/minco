#!/usr/bin/env python3
"""Score replay for the refined guarded feature allocator switch.

This is a cached-profile replay. It keeps call sets fixed, applies the
implemented `guarded-genus-hit-breadth-a002-xny230` abundance switch, and
scores the adjusted abundance values through the same selected-profile and
external diagnostic helpers used by earlier audits.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
sys.path.insert(0, str(ROOT))

import decompose_abundance_errors as decomp  # noqa: E402
import sweep_selected_call_mass_transforms as base_sweep  # noqa: E402
from scripts import minco_profile_calibrated as wrapper  # noqa: E402
import validate_feature_allocator_external_exactsplit as external  # noqa: E402
import validate_feature_allocator_wrapper_parity as selected_parity  # noqa: E402


SELECTED_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"
COMBINED_GUARD_AUDIT = RESULTS / "feature_allocator_combined_guard_audit.tsv"

OUT_SCORES = RESULTS / "feature_allocator_refined_score_replay_scores.tsv"
OUT_SUMMARY = RESULTS / "feature_allocator_refined_score_replay_summary.tsv"
OUT_VALIDATION = RESULTS / "feature_allocator_refined_score_replay_validation.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_refined_score_replay_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_REFINED_SCORE_REPLAY.md"

SWITCH = wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230

METHOD_SELECTED_BASE = "current_selected_abundance"
METHOD_SELECTED_REFINED = "wrapper_guarded_genus_hit_breadth_a002_xny230"
METHOD_EXTERNAL_BASE = "current_exactsplit_abundance"
METHOD_EXTERNAL_REFINED = "refined_feature_allocator_exactsplit"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def finite_max_abs(values: Iterable[object]) -> float:
    clean = []
    for value in values:
        out = finite_float(value, float("nan"))
        if math.isfinite(out):
            clean.append(abs(out))
    return max(clean) if clean else 0.0


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def parse_metric(text: object, key: str, default: float = 0.0) -> float:
    match = re.search(rf"(?:^|;){re.escape(key)}=([-+0-9.eE]+)", str(text))
    return finite_float(match.group(1), default) if match else default


def add_common(
    row: dict[str, object],
    evidence_group: str,
    panel: str,
    sample: int | str,
    source_profile: object,
    details: Mapping[str, object] | None = None,
) -> dict[str, object]:
    out = dict(row)
    out["evidence_group"] = evidence_group
    out["panel"] = panel
    out["sample"] = sample
    out["source_profile"] = str(source_profile)
    if details:
        for key, value in details.items():
            if key.startswith("abundance_feature_allocator_"):
                out[key] = value
    return out


def selected_score_rows() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    selected = pd.read_csv(SELECTED_SCORES, sep="\t")
    mapper = base_sweep.RowMapper()
    rows: list[dict[str, object]] = []
    validation: list[dict[str, object]] = []

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
        raw, _candidate_added = base_sweep.base_raw_values(calls)
        base = selected_parity.score_raw(panel, sample, METHOD_SELECTED_BASE, truth, calls, raw)
        target_l1 = selected_parity.official_l1(selected_row, panel)
        validation.append(
            {
                "evidence_group": "selected_cached_wrapper",
                "panel": panel,
                "sample": sample,
                "count_delta": 0,
                "F1_delta": finite_float(base["F1"]) - finite_float(selected_row.get("F1")),
                "L1_delta": finite_float(base["official_L1_pp"]) - target_l1,
                "Pearson_delta": finite_float(base["official_Pearson"], float("nan"))
                - selected_parity.official_pearson(selected_row, panel),
                "source": str(SELECTED_SCORES.relative_to(EXP)),
            }
        )
        taxmap = selected_parity.allocator_taxmap_from_profile(calls)
        adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
            calls,
            np.ones(len(calls), dtype=bool),
            raw,
            SWITCH,
            taxmap,
        )
        refined = selected_parity.score_raw(
            panel,
            sample,
            METHOD_SELECTED_REFINED,
            truth,
            calls,
            adjusted,
            details,
        )
        source_profile = selected_row.get("profile", "")
        for score in [base, refined]:
            score["L1"] = score.get("official_L1_pp", 0.0)
            score["Pearson"] = score.get("official_Pearson", float("nan"))
            rows.append(
                add_common(
                    score,
                    "selected_cached_wrapper",
                    panel,
                    sample,
                    source_profile,
                    details if score["method"] == METHOD_SELECTED_REFINED else None,
                )
            )
    return rows, validation


def external_score_rows() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    samples = [
        ("plant_holdout", 3),
        ("plant_holdout", 4),
        ("plant_holdout", 5),
        ("strainmadness", 0),
        ("strainmadness", 1),
        ("strainmadness", 2),
    ]
    required = [external.TAXMAP]
    for dataset, sample in samples:
        required.extend(
            [
                external.profile_path(dataset, sample),
                external.truth_path(dataset, sample),
                external.score_file(dataset, sample),
            ]
        )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required external score replay inputs:\n" + "\n".join(missing))

    taxmap = external.plant_score.parse_species_taxmap(external.TAXMAP)
    scope_by_taxid = external.plant_score.taxid_scope_map(taxmap)
    rows: list[dict[str, object]] = []
    validation: list[dict[str, object]] = []

    for dataset, sample in samples:
        profile_path = external.profile_path(dataset, sample)
        profile = external.load_profile(profile_path, scope_by_taxid)
        call = external.current_mask(profile)
        raw = external.current_raw(profile)
        adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
            profile,
            call,
            raw,
            SWITCH,
            taxmap,
        )
        work = profile.copy()
        work["feature_allocator_refined_raw"] = adjusted
        base = external.score_calls(
            dataset,
            sample,
            work,
            call,
            METHOD_EXTERNAL_BASE,
            "calibrated_abundance",
        )
        refined = external.score_calls(
            dataset,
            sample,
            work,
            call,
            METHOD_EXTERNAL_REFINED,
            "feature_allocator_refined_raw",
        )
        ref = external.reference_current_score(dataset, sample)
        validation.append(
            {
                "evidence_group": "external_exactsplit_diagnostic",
                "panel": dataset,
                "sample": sample,
                "count_delta": max(
                    abs(int(base["TP"]) - int(float(ref["TP"]))),
                    abs(int(base["FP"]) - int(float(ref["FP"]))),
                    abs(int(base["FN"]) - int(float(ref["FN"]))),
                ),
                "F1_delta": finite_float(base["F1"]) - finite_float(ref["F1"]),
                "L1_delta": finite_float(base["bacteria_scope_l1"]) - finite_float(ref["L1"]),
                "Pearson_delta": finite_float(base["tp_abundance_pearson"], float("nan"))
                - finite_float(ref["Pearson"], float("nan")),
                "source": str(external.score_file(dataset, sample)),
            }
        )
        for score in [base, refined]:
            score["L1"] = score.get("bacteria_scope_l1", 0.0)
            score["Pearson"] = score.get("tp_abundance_pearson", float("nan"))
            rows.append(
                add_common(
                    score,
                    "external_exactsplit_diagnostic",
                    dataset,
                    sample,
                    profile_path,
                    details if score["method"] == METHOD_EXTERNAL_REFINED else None,
                )
            )
    return rows, validation


def add_deltas(scores: pd.DataFrame) -> pd.DataFrame:
    scores = scores.copy()
    base_methods = {METHOD_SELECTED_BASE, METHOD_EXTERNAL_BASE}
    base = scores.loc[scores["method"].isin(base_methods)].set_index(
        ["evidence_group", "panel", "sample"]
    )
    out_rows: list[dict[str, object]] = []
    for row in scores.to_dict("records"):
        key = (row["evidence_group"], row["panel"], row["sample"])
        if row["method"] in base_methods:
            row["delta_F1"] = 0.0
            row["delta_L1"] = 0.0
            row["delta_Pearson"] = 0.0
        else:
            base_row = base.loc[key]
            row["delta_F1"] = finite_float(row.get("F1")) - finite_float(base_row.get("F1"))
            row["delta_L1"] = finite_float(row.get("L1")) - finite_float(base_row.get("L1"))
            row["delta_Pearson"] = finite_float(row.get("Pearson"), float("nan")) - finite_float(
                base_row.get("Pearson"), float("nan")
            )
        out_rows.append(row)
    return pd.DataFrame(out_rows)


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    refined = scores.loc[scores["method"].isin({METHOD_SELECTED_REFINED, METHOD_EXTERNAL_REFINED})].copy()
    rows: list[dict[str, object]] = []
    groups = [(name, sub) for name, sub in refined.groupby("evidence_group", sort=True)]
    groups.append(("all", refined))
    for group_name, sub in groups:
        rows.append(
            {
                "evidence_group": group_name,
                "profile_count": int(len(sub)),
                "switched_profiles": int(
                    sub.get("abundance_feature_allocator_applied", pd.Series([], dtype=object))
                    .astype(str)
                    .str.lower()
                    .isin({"true", "1"})
                    .sum()
                ),
                "mean_F1_delta": float(sub["delta_F1"].mean()),
                "min_F1_delta": float(sub["delta_F1"].min()),
                "mean_L1_delta": float(sub["delta_L1"].mean()),
                "max_worse_L1_delta": float(sub["delta_L1"].max()),
                "improved_L1_profiles": int((sub["delta_L1"] < -1e-9).sum()),
                "worsened_L1_profiles": int((sub["delta_L1"] > 1e-9).sum()),
                "mean_Pearson_delta": float(sub["delta_Pearson"].mean()),
            }
        )
    return pd.DataFrame(rows)


def build_audit(summary: pd.DataFrame, validation: pd.DataFrame) -> list[dict[str, object]]:
    combined_guard = pd.read_csv(COMBINED_GUARD_AUDIT, sep="\t")
    guard_by_metric = {str(row["metric"]): row for row in combined_guard.to_dict("records")}
    guard_text = guard_by_metric.get("best_strict_rule", {}).get("value", "")
    guard_mean = parse_metric(guard_text, "mean_L1_delta_pp")
    guard_selected_mean = parse_metric(guard_text, "selected_mean_L1_delta_pp")
    guard_external_mean = parse_metric(guard_text, "external_mean_L1_delta_pp")

    overall = summary.loc[summary["evidence_group"].eq("all")].iloc[0]
    selected = summary.loc[summary["evidence_group"].eq("selected_cached_wrapper")].iloc[0]
    external_row = summary.loc[summary["evidence_group"].eq("external_exactsplit_diagnostic")].iloc[0]
    count_max = int(pd.to_numeric(validation["count_delta"], errors="coerce").fillna(0).max())
    f1_max = finite_max_abs(validation["F1_delta"])
    l1_max = finite_max_abs(validation["L1_delta"])
    pearson_max = finite_max_abs(validation["Pearson_delta"])
    replay_guard_delta = max(
        abs(float(overall["mean_L1_delta"]) - guard_mean),
        abs(float(selected["mean_L1_delta"]) - guard_selected_mean),
        abs(float(external_row["mean_L1_delta"]) - guard_external_mean),
    )
    guard_match_tolerance = 1e-6
    pass_replay = (
        count_max == 0
        and f1_max < 1e-9
        and int(overall["worsened_L1_profiles"]) == 0
        and int(overall["improved_L1_profiles"]) > 0
        and replay_guard_delta < guard_match_tolerance
    )
    decision = (
        "refined_guard_score_replay_pass_candidate_needs_independent_holdout"
        if pass_replay
        else "refined_guard_score_replay_review"
    )
    return [
        {
            "metric": "baseline_validation",
            "value": (
                f"counts={count_max};F1={f1_max:.3g};L1={l1_max:.3g};"
                f"Pearson={pearson_max:.3g}"
            ),
            "evidence": str(OUT_VALIDATION.relative_to(EXP)),
            "decision": "pass" if count_max == 0 and f1_max < 1e-9 else "review",
        },
        {
            "metric": "score_replay_effect",
            "value": (
                f"samples={int(overall['profile_count'])};"
                f"switched={int(overall['switched_profiles'])};"
                f"improved={int(overall['improved_L1_profiles'])};"
                f"worsened={int(overall['worsened_L1_profiles'])};"
                f"mean_L1_delta_pp={float(overall['mean_L1_delta']):.6f};"
                f"selected_mean_L1_delta_pp={float(selected['mean_L1_delta']):.6f};"
                f"external_mean_L1_delta_pp={float(external_row['mean_L1_delta']):.6f};"
                f"max_worse_L1_delta_pp={float(overall['max_worse_L1_delta']):.6f}"
            ),
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "fixed_call_abundance_replay",
        },
        {
            "metric": "score_replay_matches_guard_estimate",
            "value": f"max_abs_mean_delta={replay_guard_delta:.12g}",
            "evidence": f"{OUT_SUMMARY.relative_to(EXP)};{COMBINED_GUARD_AUDIT.relative_to(EXP)}",
            "decision": "pass" if replay_guard_delta < guard_match_tolerance else "review",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "current_default_unchanged",
        },
    ]


def write_markdown(summary: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Refined Feature Allocator Score Replay",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-profile replay validates the implemented opt-in switch",
        f"`--abundance-feature-allocator-switch {SWITCH}`. The call set is fixed;",
        "only abundance mass is adjusted, and the default remains unchanged.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "baseline_validation",
        "score_replay_effect",
        "score_replay_matches_guard_estimate",
        "promotion_decision",
    ]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Evidence group | Profiles | Switched | Mean F1 delta | Mean L1 delta | Max worse L1 | Improved L1 | Worsened L1 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.to_dict("records"):
        lines.append(
            "| {group} | {n} | {switched} | {f1:.6f} | {l1:.6f} | {max_worse:.6f} | {improved} | {worse} |".format(
                group=row["evidence_group"],
                n=int(row["profile_count"]),
                switched=int(row["switched_profiles"]),
                f1=float(row["mean_F1_delta"]),
                l1=float(row["mean_L1_delta"]),
                max_worse=float(row["max_worse_L1_delta"]),
                improved=int(row["improved_L1_profiles"]),
                worse=int(row["worsened_L1_profiles"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Keep the current default unchanged.",
            "- Treat the refined allocator as the strongest opt-in abundance candidate.",
            "- Require independent holdout validation before any default promotion.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SCORES.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_VALIDATION.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    selected_rows, selected_validation = selected_score_rows()
    external_rows, external_validation = external_score_rows()
    scores = add_deltas(pd.DataFrame(selected_rows + external_rows))
    validation = pd.DataFrame(selected_validation + external_validation)
    summary = summarize(scores)
    audit = build_audit(summary, validation)

    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    validation.to_csv(OUT_VALIDATION, sep="\t", index=False)
    pd.DataFrame(audit).to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(summary, audit)

    print(pd.DataFrame(audit).to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
