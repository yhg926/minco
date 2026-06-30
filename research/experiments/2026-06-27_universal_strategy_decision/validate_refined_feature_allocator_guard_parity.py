#!/usr/bin/env python3
"""Validate wrapper guard parity for the refined feature allocator switch."""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
sys.path.insert(0, str(ROOT))

from scripts import minco_profile_calibrated as wrapper  # noqa: E402
import audit_feature_allocator_combined_guard_refinement as combined_guard  # noqa: E402
import validate_feature_allocator_wrapper_parity as selected_parity  # noqa: E402


FEATURES_TSV = RESULTS / "feature_allocator_combined_guard_features.tsv"
SELECTED_SCORES = RESULTS / "feature_allocator_wrapper_parity_scores.tsv"
EXTERNAL_SCORES = RESULTS / "feature_allocator_external_exactsplit_scores.tsv"
EXTERNAL_TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")

OUT_APPLY = RESULTS / "feature_allocator_refined_guard_parity_apply.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_refined_guard_parity_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_REFINED_GUARD_PARITY.md"


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    return wrapper.finite(value, default)


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def load_profile(path: Path) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    rows = pd.read_csv(path, sep="\t", low_memory=False)
    call = wrapper.boolean_column(rows, "calibrated_call").to_numpy(dtype=bool)
    if not call.any():
        call = np.ones(len(rows), dtype=bool)
    raw_col = "calibrated_abundance_raw" if "calibrated_abundance_raw" in rows.columns else "calibrated_abundance"
    raw = pd.to_numeric(rows[raw_col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    raw_values = np.where(raw.to_numpy(dtype=float) > 0.0, raw.to_numpy(dtype=float), 0.0)
    return rows, call, raw_values


def old_switch_status() -> dict[tuple[str, str], bool]:
    out: dict[tuple[str, str], bool] = {}
    selected = pd.read_csv(SELECTED_SCORES, sep="\t")
    selected = selected.loc[selected["method"].eq("wrapper_guarded_genus_hit_breadth_a002")]
    for row in selected.to_dict("records"):
        out[(str(row["panel"]), str(row["sample"]))] = is_true(
            row.get("abundance_feature_allocator_applied", False)
        )
    external = pd.read_csv(EXTERNAL_SCORES, sep="\t")
    external = external.loc[external["method"].eq("guarded_feature_allocator_exactsplit")]
    for row in external.to_dict("records"):
        out[(str(row["dataset"]), str(int(float(row["sample"]))))] = is_true(
            row.get("abundance_feature_allocator_applied", False)
        )
    return out


def external_taxmap() -> dict[str, dict[str, str]]:
    return wrapper.parse_species_taxmap(EXTERNAL_TAXMAP) if EXTERNAL_TAXMAP.is_file() else {}


def build_apply_rows() -> list[dict[str, object]]:
    if not FEATURES_TSV.exists():
        combined_guard.main()
    features = pd.read_csv(FEATURES_TSV, sep="\t")
    old_applied = old_switch_status()
    ext_taxmap = external_taxmap()
    rows: list[dict[str, object]] = []
    for feature in features.to_dict("records"):
        panel = str(feature["panel"])
        sample = str(feature["sample"])
        path = Path(str(feature["profile"]))
        profile, call, raw = load_profile(path)
        taxmap = None
        if str(feature.get("evidence_group", "")) == "selected_cached_wrapper":
            taxmap = selected_parity.allocator_taxmap_from_profile(profile.loc[call].copy())
        elif ext_taxmap:
            taxmap = ext_taxmap
        _adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
            profile,
            call,
            raw,
            wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230,
            taxmap,
        )
        rule_selected = finite(feature.get("s_xny_median")) >= wrapper.ABUNDANCE_FEATURE_ALLOCATOR_S_XNY_MEDIAN_MIN
        expected_applied = bool(rule_selected and old_applied.get((panel, sample), False))
        actual_applied = bool(details.get("abundance_feature_allocator_applied", False))
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "evidence_group": feature.get("evidence_group", ""),
                "profile": str(path),
                "rule_selected": rule_selected,
                "old_switch_applied": old_applied.get((panel, sample), False),
                "expected_refined_applied": expected_applied,
                "wrapper_refined_applied": actual_applied,
                "apply_match": expected_applied == actual_applied,
                "s_xny_median_feature": finite(feature.get("s_xny_median")),
                "s_xny_median_wrapper": finite(details.get("abundance_feature_allocator_s_xny_median")),
                "base_mass_multi_genus_frac_wrapper": finite(
                    details.get("abundance_feature_allocator_base_mass_multi_genus_frac")
                ),
                "adjusted_rows": int(details.get("abundance_feature_allocator_adjusted_rows_n", 0)),
            }
        )
    return rows


def build_audit(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    total = len(rows)
    mismatches = [row for row in rows if not bool(row["apply_match"])]
    applied = [row for row in rows if bool(row["wrapper_refined_applied"])]
    selected_applied = [
        row
        for row in applied
        if str(row.get("evidence_group")) == "selected_cached_wrapper"
    ]
    external_applied = [
        row
        for row in applied
        if str(row.get("evidence_group")) == "external_exactsplit_diagnostic"
    ]
    max_xny_delta = max(
        abs(float(row["s_xny_median_feature"]) - float(row["s_xny_median_wrapper"]))
        for row in rows
    )
    return [
        {
            "metric": "profile_count",
            "value": total,
            "evidence": str(OUT_APPLY.relative_to(EXP)),
            "decision": "cached_profiles_only",
        },
        {
            "metric": "wrapper_apply_matches_refined_rule",
            "value": f"mismatches={len(mismatches)};max_s_xny_median_delta={max_xny_delta:.12g}",
            "evidence": str(OUT_APPLY.relative_to(EXP)),
            "decision": "pass" if not mismatches and max_xny_delta < 1e-9 else "review",
        },
        {
            "metric": "refined_switch_applied_profiles",
            "value": (
                f"total={len(applied)};selected_cached={len(selected_applied)};"
                f"external={len(external_applied)}"
            ),
            "evidence": str(OUT_APPLY.relative_to(EXP)),
            "decision": "wrapper_refined_guard_active",
        },
        {
            "metric": "promotion_decision",
            "value": (
                "refined_guard_wrapper_validated_candidate_not_default"
                if not mismatches and max_xny_delta < 1e-9
                else "review_refined_guard_wrapper"
            ),
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": (
                "guard_parity_only_not_default"
                if not mismatches and max_xny_delta < 1e-9
                else "review"
            ),
        },
    ]


def write_markdown(audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Refined Feature Allocator Guard Parity",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note validates the implemented refined allocator switch",
        "against the cached combined guard-refinement rule. It does not rerun raw",
        "profiling jobs.",
        "",
        "## Decision",
        "",
        f"- Profiles: `{audit_by_metric['profile_count']['value']}`.",
        f"- Rule parity: `{audit_by_metric['wrapper_apply_matches_refined_rule']['value']}` "
        f"({audit_by_metric['wrapper_apply_matches_refined_rule']['decision']}).",
        f"- Applied profiles: `{audit_by_metric['refined_switch_applied_profiles']['value']}`.",
        f"- Promotion decision: `{audit_by_metric['promotion_decision']['decision']}`.",
        "",
        "## Outputs",
        "",
        f"- `{OUT_APPLY.relative_to(EXP)}`",
        f"- `{OUT_AUDIT.relative_to(EXP)}`",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows = build_apply_rows()
    write_tsv(
        OUT_APPLY,
        rows,
        [
            "panel",
            "sample",
            "evidence_group",
            "profile",
            "rule_selected",
            "old_switch_applied",
            "expected_refined_applied",
            "wrapper_refined_applied",
            "apply_match",
            "s_xny_median_feature",
            "s_xny_median_wrapper",
            "base_mass_multi_genus_frac_wrapper",
            "adjusted_rows",
        ],
    )
    audit = build_audit(rows)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
