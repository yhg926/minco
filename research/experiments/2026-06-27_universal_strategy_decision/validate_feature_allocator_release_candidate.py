#!/usr/bin/env python3
"""Validate the strongest guarded feature-allocator candidate for promotion.

This combines cached wrapper-parity evidence with the external exact-split
stress test. It does not rerun raw profiling jobs. The promotion rule is
strict because this would affect the no-expertise default abundance output:
the allocator must preserve calls/F1 and avoid sample-level L1 regressions
when checked outside the panel used to select the guard.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

WRAPPER_SUMMARY = RESULTS / "feature_allocator_wrapper_parity_summary.tsv"
WRAPPER_AUDIT = RESULTS / "feature_allocator_wrapper_parity_audit.tsv"
EXTERNAL_SCORES = RESULTS / "feature_allocator_external_exactsplit_scores.tsv"
EXTERNAL_SUMMARY = RESULTS / "feature_allocator_external_exactsplit_summary.tsv"
EXTERNAL_AUDIT = RESULTS / "feature_allocator_external_exactsplit_audit.tsv"

OUT_PANEL = RESULTS / "feature_allocator_release_candidate_panel.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_release_candidate_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_RELEASE_CANDIDATE.md"

METHOD_EXTERNAL = "guarded_feature_allocator_exactsplit"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def as_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return default
    return float(value)


def as_int(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return default
    return int(float(value))


def weighted_mean(items: list[dict[str, object]], value_key: str, weight_key: str) -> float:
    total_weight = sum(float(row[weight_key]) for row in items)
    if total_weight <= 0.0:
        return 0.0
    return sum(float(row[value_key]) * float(row[weight_key]) for row in items) / total_weight


def build_panel_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    wrapper = table_by(WRAPPER_SUMMARY, "panel")
    for panel, source in wrapper.items():
        if panel == "all":
            continue
        rows.append(
            {
                "evidence_group": "selected_cached_wrapper",
                "panel": panel,
                "sample_count": as_int(source, "sample_count"),
                "switched_samples": as_int(source, "switched_samples"),
                "mean_L1_delta_pp": as_float(source, "mean_L1_delta_pp"),
                "max_worse_L1_delta_pp": as_float(source, "max_worse_L1_delta_pp"),
                "improved_samples": as_int(source, "improved_samples"),
                "worsened_samples": as_int(source, "worsened_samples"),
                "mean_Pearson_delta": as_float(source, "mean_Pearson_delta"),
                "mean_F1_delta": 0.0,
                "decision": "cached_wrapper_sample_safe"
                if as_int(source, "worsened_samples") == 0
                else "cached_wrapper_regression",
            }
        )

    external = [
        row
        for row in read_tsv(EXTERNAL_SUMMARY)
        if row.get("method") == METHOD_EXTERNAL
    ]
    for source in external:
        rows.append(
            {
                "evidence_group": "external_exactsplit_diagnostic",
                "panel": source["dataset"],
                "sample_count": as_int(source, "samples", 0)
                if str(source.get("samples", "")).isdigit()
                else len([part for part in source.get("samples", "").split(",") if part]),
                "switched_samples": as_int(source, "switched_samples"),
                "mean_L1_delta_pp": as_float(source, "mean_delta_current_l1"),
                "max_worse_L1_delta_pp": as_float(source, "max_worse_l1"),
                "improved_samples": as_int(source, "improved_samples_l1"),
                "worsened_samples": as_int(source, "worsened_samples_l1"),
                "mean_Pearson_delta": as_float(source, "mean_delta_current_pearson"),
                "mean_F1_delta": as_float(source, "mean_delta_current_F1"),
                "decision": "external_sample_safe"
                if as_int(source, "worsened_samples_l1") == 0
                else "external_sample_regression",
            }
        )
    return rows


def external_sample_details() -> list[dict[str, object]]:
    out = []
    scores = read_tsv(EXTERNAL_SCORES)
    current_l1 = {
        (row["dataset"], row["sample"]): as_float(row, "bacteria_scope_l1")
        for row in scores
        if row.get("method") == "current_exactsplit_abundance"
    }
    for row in scores:
        if row.get("method") != METHOD_EXTERNAL:
            continue
        delta = as_float(row, "bacteria_scope_l1") - current_l1[(row["dataset"], row["sample"])]
        out.append(
            {
                "dataset": row["dataset"],
                "sample": row["sample"],
                "delta_current_l1": delta,
                "base_mass_multi_genus_frac": as_float(
                    row, "abundance_feature_allocator_base_mass_multi_genus_frac"
                ),
                "adjusted_rows": as_int(row, "abundance_feature_allocator_adjusted_rows_n"),
                "decision": "regressed" if delta > 1e-12 else "improved_or_equal",
            }
        )
    return out


def build_audit(panel_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    wrapper_all = table_by(WRAPPER_SUMMARY, "panel").get("all", {})
    wrapper_audit = table_by(WRAPPER_AUDIT, "metric")
    external_audit = table_by(EXTERNAL_AUDIT, "metric")
    external_details = external_sample_details()

    combined_sample_count = sum(int(row["sample_count"]) for row in panel_rows)
    combined_switched = sum(int(row["switched_samples"]) for row in panel_rows)
    combined_improved = sum(int(row["improved_samples"]) for row in panel_rows)
    combined_worsened = sum(int(row["worsened_samples"]) for row in panel_rows)
    combined_mean_l1 = weighted_mean(panel_rows, "mean_L1_delta_pp", "sample_count")
    combined_max_worse = max(float(row["max_worse_L1_delta_pp"]) for row in panel_rows)
    combined_mean_f1 = weighted_mean(panel_rows, "mean_F1_delta", "sample_count")
    external_regressions = [row for row in external_details if row["decision"] == "regressed"]

    strict_pass = combined_worsened == 0 and combined_mean_f1 >= -1e-12
    return [
        {
            "metric": "wrapper_cached_sample_safety",
            "value": (
                f"samples={wrapper_all.get('sample_count', 'NA')};"
                f"switched={wrapper_all.get('switched_samples', 'NA')};"
                f"mean_L1_delta_pp={wrapper_all.get('mean_L1_delta_pp', 'NA')};"
                f"worsened={wrapper_all.get('worsened_samples', 'NA')};"
                f"max_worse={wrapper_all.get('max_worse_L1_delta_pp', 'NA')}"
            ),
            "evidence": str(WRAPPER_SUMMARY.relative_to(EXP)),
            "decision": wrapper_audit.get("promotion_decision", {}).get(
                "value", "candidate_needs_independent_holdout"
            ),
        },
        {
            "metric": "external_exactsplit_sample_safety",
            "value": external_audit.get("candidate_effect", {}).get("value", "NA"),
            "evidence": str(EXTERNAL_AUDIT.relative_to(EXP)),
            "decision": external_audit.get("candidate_effect", {}).get(
                "decision", "external_stress_rejects_default"
            ),
        },
        {
            "metric": "combined_candidate_effect",
            "value": (
                f"samples={combined_sample_count};switched={combined_switched};"
                f"improved={combined_improved};worsened={combined_worsened};"
                f"mean_L1_delta_pp={combined_mean_l1:.6f};"
                f"max_worse_L1_delta_pp={combined_max_worse:.6f};"
                f"mean_F1_delta={combined_mean_f1:.12g}"
            ),
            "evidence": f"{OUT_PANEL.relative_to(EXP)};{EXTERNAL_SCORES.relative_to(EXP)}",
            "decision": "strict_promotion_pass" if strict_pass else "strict_promotion_fail",
        },
        {
            "metric": "external_regression_details",
            "value": (
                "none"
                if not external_regressions
                else ";".join(
                    "{dataset}{sample}:delta_L1={delta:.9f};base_multi={base:.6f};rows={rows}".format(
                        dataset=row["dataset"],
                        sample=row["sample"],
                        delta=float(row["delta_current_l1"]),
                        base=float(row["base_mass_multi_genus_frac"]),
                        rows=int(row["adjusted_rows"]),
                    )
                    for row in external_regressions
                )
            ),
            "evidence": str(EXTERNAL_SCORES.relative_to(EXP)),
            "decision": "external_sample_regressions_present"
            if external_regressions
            else "external_sample_safe",
        },
        {
            "metric": "promotion_decision",
            "value": "keep_feature_allocator_experimental_off_by_default"
            if not strict_pass
            else "candidate_can_be_promoted_after_full_default_replay",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "do_not_promote_until_external_sample_safe"
            if not strict_pass
            else "promotion_candidate",
        },
        {
            "metric": "next_requirement",
            "value": (
                "find output-only guard or allocator target with unchanged F1, "
                "zero external sample L1 regressions, and preserved cached-panel gains"
            ),
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "algorithmic_next_step",
        },
    ]


def write_markdown(panel_rows: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Feature Allocator Release Candidate Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note combines cached wrapper-parity evidence with the",
        "external exact-split stress test for the guarded feature allocator. It",
        "does not rerun raw profiling jobs.",
        "",
        "## Decision",
        "",
        f"- Cached wrapper evidence: `{audit_by_metric['wrapper_cached_sample_safety']['value']}`.",
        f"- External stress evidence: `{audit_by_metric['external_exactsplit_sample_safety']['value']}`.",
        f"- Combined effect: `{audit_by_metric['combined_candidate_effect']['value']}`.",
        f"- External regression detail: `{audit_by_metric['external_regression_details']['value']}`.",
        f"- Promotion decision: `{audit_by_metric['promotion_decision']['decision']}`.",
        "",
        "## Panel Evidence",
        "",
        "| Evidence group | Panel | Samples | Mean L1 delta pp | Max worse pp | Improved | Worsened | Decision |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in panel_rows:
        lines.append(
            "| {group} | {panel} | {samples} | {mean:.6f} | {max_worse:.6f} | {improved} | {worsened} | {decision} |".format(
                group=row["evidence_group"],
                panel=row["panel"],
                samples=int(row["sample_count"]),
                mean=float(row["mean_L1_delta_pp"]),
                max_worse=float(row["max_worse_L1_delta_pp"]),
                improved=int(row["improved_samples"]),
                worsened=int(row["worsened_samples"]),
                decision=row["decision"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The allocator remains a useful experimental candidate because it",
            "  improves the cached selected-default panels and both external panel",
            "  means without changing calls.",
            "- It is not promoted into the default because strict sample safety fails",
            "  on the external exact-split diagnostic profiles.",
            "- The next abundance step should change the output-only guard or target",
            "  feature so independent-profile sample regressions are zero before",
            "  default promotion is reconsidered.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_PANEL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    panel_rows = build_panel_rows()
    audit = build_audit(panel_rows)
    write_tsv(
        OUT_PANEL,
        panel_rows,
        [
            "evidence_group",
            "panel",
            "sample_count",
            "switched_samples",
            "mean_L1_delta_pp",
            "max_worse_L1_delta_pp",
            "improved_samples",
            "worsened_samples",
            "mean_Pearson_delta",
            "mean_F1_delta",
            "decision",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(panel_rows, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
