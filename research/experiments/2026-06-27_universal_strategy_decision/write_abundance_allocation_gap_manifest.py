#!/usr/bin/env python3
"""Write a compact abundance-gap manifest from cached decision TSVs."""

from __future__ import annotations

import csv
from pathlib import Path


EXP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXP_DIR / "results"

DELTA_TSV = RESULTS_DIR / "abundance_error_decomposition_delta.tsv"
VARIANT_TSV = RESULTS_DIR / "cross_panel_abundance_variant_overall.tsv"
CANDIDATE_VS_SYLPH_TSV = RESULTS_DIR / "candidate_default_vs_sylph.tsv"
CANDIDATE_VS_CURRENT_TSV = RESULTS_DIR / "candidate_default_vs_current.tsv"
BLEND_AUDIT_TSV = RESULTS_DIR / "candidate_preset_genus_xny_blend_audit.tsv"
BLEND_PANEL_TSV = RESULTS_DIR / "candidate_preset_genus_xny_blend_panel_delta.tsv"
BLEND_GUARD_AUDIT_TSV = RESULTS_DIR / "candidate_preset_genus_xny_blend_guard_audit.tsv"
BLEND_ALPHA_AUDIT_TSV = RESULTS_DIR / "candidate_preset_genus_xny_alpha_sweep_audit.tsv"
BLEND_ALPHA_OVERALL_TSV = RESULTS_DIR / "candidate_preset_genus_xny_alpha_sweep_overall.tsv"

OUT_TSV = RESULTS_DIR / "abundance_allocation_gap_manifest.tsv"
OUT_MD = EXP_DIR / "ABUNDANCE_ALLOCATION_GAP.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fnum(value: str) -> float:
    return float(value) if value not in {"", "NA", "nan"} else float("nan")


def fmt(value: str | float, digits: int = 6) -> str:
    if isinstance(value, str):
        value = fnum(value)
    return f"{value:.{digits}f}"


def best_cached_allocator(rows: list[dict[str, str]]) -> dict[str, str]:
    candidates = [
        row
        for row in rows
        if row.get("decision") == "promotable_by_L1_if_validated_from_raw"
    ]
    if not candidates:
        return {}
    return min(candidates, key=lambda row: fnum(row["mean_delta_current_L1_pp"]))


def audit_metric(rows: list[dict[str, str]], metric: str) -> dict[str, str]:
    for row in rows:
        if row.get("metric") == metric:
            return row
    return {}


def write_tsv(rows: list[dict[str, str]]) -> None:
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "record_type",
        "panel",
        "main_gap_component",
        "selected_default_minus_sylph_L1_pp",
        "legacy_minco_minus_sylph_L1_pp",
        "legacy_matched_abs_error_gap_pp",
        "legacy_missing_truth_mass_gap_pp",
        "legacy_extra_pred_mass_gap_pp",
        "selected_default_minus_sylph_F1",
        "selected_default_minus_sylph_Pearson",
        "cached_allocator",
        "cached_allocator_mean_delta_current_L1_pp",
        "cached_allocator_max_worse_current_L1_pp",
        "cached_allocator_improved_panel_count",
        "decision",
    ]
    with OUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(
    panel_rows: list[dict[str, str]],
    allocator_row: dict[str, str],
    blend_row: dict[str, str],
    blend_panel_rows: list[dict[str, str]],
    guard_row: dict[str, str],
    alpha_row: dict[str, str],
    alpha_overall_rows: list[dict[str, str]],
    default_row: dict[str, str],
) -> None:
    lines = [
        "# Abundance Allocation Gap Manifest",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached decision TSVs. It records the",
        "remaining abundance gap for the selected MinCO default without rerunning",
        "raw profiling jobs.",
        "",
        "## Source Tables",
        "",
        f"- `{DELTA_TSV.relative_to(EXP_DIR)}`",
        f"- `{VARIANT_TSV.relative_to(EXP_DIR)}`",
        f"- `{CANDIDATE_VS_SYLPH_TSV.relative_to(EXP_DIR)}`",
        f"- `{CANDIDATE_VS_CURRENT_TSV.relative_to(EXP_DIR)}`",
        f"- `{BLEND_AUDIT_TSV.relative_to(EXP_DIR)}`",
        f"- `{BLEND_PANEL_TSV.relative_to(EXP_DIR)}`",
        f"- `{BLEND_GUARD_AUDIT_TSV.relative_to(EXP_DIR)}`",
        f"- `{BLEND_ALPHA_AUDIT_TSV.relative_to(EXP_DIR)}`",
        f"- `{BLEND_ALPHA_OVERALL_TSV.relative_to(EXP_DIR)}`",
        "",
        "## Panel Gap",
        "",
        "| Panel | Main gap | Selected default - Sylph L1 pp | Legacy matched-gap pp | Legacy missing-gap pp | Legacy extra-gap pp |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in panel_rows:
        lines.append(
            "| {panel} | {main_gap_component} | {selected_l1} | {matched} | {missing} | {extra} |".format(
                panel=row["panel"],
                main_gap_component=row["main_gap_component"],
                selected_l1=row["selected_default_minus_sylph_L1_pp"],
                matched=row["legacy_matched_abs_error_gap_pp"],
                missing=row["legacy_missing_truth_mass_gap_pp"],
                extra=row["legacy_extra_pred_mass_gap_pp"],
            )
        )

    lines.extend(
        [
            "",
            "## Cached Allocation Candidate",
            "",
            "| Method | Mean L1 delta vs current pp | Max worse pp | Improved panels | Status |",
            "|---|---:|---:|---:|---|",
        ]
    )
    if allocator_row:
        lines.append(
            "| {method} | {mean_delta} | {max_worse} | {improved} | {status} |".format(
                method=allocator_row["cached_allocator"],
                mean_delta=allocator_row["cached_allocator_mean_delta_current_L1_pp"],
                max_worse=allocator_row["cached_allocator_max_worse_current_L1_pp"],
                improved=allocator_row["cached_allocator_improved_panel_count"],
                status=allocator_row["decision"],
            )
        )
    else:
        lines.append("| NA | NA | NA | NA | no cached candidate |")

    lines.extend(
        [
            "",
            "## Selected-Default Blend Audit",
            "",
            "| Method | Mean L1 delta vs selected default pp | Max worse pp | Sample direction | Decision |",
            "|---|---:|---:|---|---|",
        ]
    )
    if blend_row:
        lines.append(
            "| {method} | {mean_delta} | {max_worse} | {direction} | {status} |".format(
                method=blend_row["cached_allocator"],
                mean_delta=blend_row["cached_allocator_mean_delta_current_L1_pp"],
                max_worse=blend_row["cached_allocator_max_worse_current_L1_pp"],
                direction=blend_row["cached_allocator_improved_panel_count"],
                status=blend_row["decision"],
            )
        )
    else:
        lines.append("| NA | NA | NA | NA | no selected-default blend audit |")

    if blend_panel_rows:
        lines.extend(
            [
                "",
                "| Panel | Mean L1 delta pp | Max worse pp | Improved samples | Worsened samples |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for row in blend_panel_rows:
            lines.append(
                "| {panel} | {mean_delta} | {max_worse} | {improved} | {worsened} |".format(
                    panel=row["panel"],
                    mean_delta=fmt(row["mean_L1_delta_pp"]),
                    max_worse=fmt(row["max_worse_L1_delta_pp"]),
                    improved=row["improved_samples"],
                    worsened=row["worsened_samples"],
                )
            )

    lines.extend(
        [
            "",
            "## Adaptive Blend Guard Audit",
            "",
            "| Tested rules | Sample-safe in-panel | LOPO result | Decision |",
            "|---:|---:|---|---|",
        ]
    )
    if guard_row:
        lines.append(
            "| {tested} | {sample_safe} | {lopo} | {decision} |".format(
                tested=guard_row.get("tested_guard_rules", "NA"),
                sample_safe=guard_row.get("sample_safe_guards", "NA"),
                lopo=guard_row.get("lopo_guard_result", "NA"),
                decision=guard_row.get("decision", "NA"),
            )
        )
    else:
        lines.append("| NA | NA | NA | no guard audit |")

    lines.extend(
        [
            "",
            "## Conservative Alpha Sweep",
            "",
            "| Variants | Sample-safe variants | Best ranked variant | Decision |",
            "|---:|---:|---|---|",
        ]
    )
    if alpha_row:
        lines.append(
            "| {tested} | {sample_safe} | {best} | {decision} |".format(
                tested=alpha_row.get("variants_tested", "NA"),
                sample_safe=alpha_row.get("sample_safe_variants", "NA"),
                best=alpha_row.get("best_ranked_variant", "NA"),
                decision=alpha_row.get("decision", "NA"),
            )
        )
    else:
        lines.append("| NA | NA | NA | no alpha sweep |")
    if alpha_overall_rows:
        lines.extend(
            [
                "",
                "| Method | Mean panel L1 delta pp | Worsened samples | Max sample worse pp |",
                "|---|---:|---:|---:|",
            ]
        )
        for row in alpha_overall_rows:
            lines.append(
                "| {method} | {mean_delta} | {worsened} | {max_worse} |".format(
                    method=row["method"],
                    mean_delta=fmt(row["mean_panel_L1_delta_pp"]),
                    worsened=row["worsened_samples"],
                    max_worse=fmt(row["max_sample_worse_L1_delta_pp"]),
                )
            )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            f"- Selected default status: `{default_row['decision']}`.",
            "- The candidate preset improves MinCO versus its previous default on the",
            "  matched panels, but it still loses abundance L1 to Sylph on all four",
            "  cached comparison panels.",
            "- The fixed-call allocator sweep has a small safe cached signal, but the",
            "  selected-default posthoc audit has individual-sample regressions, so",
            "  it is not promoted.",
            "- Output-derived blend guards can be sample-safe in-panel, but the",
            "  leave-one-panel-out audit still has held-out sample regressions, so",
            "  no adaptive blend guard is promoted.",
            "- Lower-alpha and capped genus-XnY sweeps improve all panel means, but",
            "  no nonzero tested variant is sample-safe, so this family remains",
            "  off by default.",
            "- Next abundance work should target matched-call mass allocation first.",
            "  The selected-call oracle feasibility audit in",
            "  `ABUNDANCE_NEXT_TARGET.md` narrows call recovery to one cached",
            "  insufficient-call sample under the selected candidate preset.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    delta_rows = read_tsv(DELTA_TSV)
    variant_rows = read_tsv(VARIANT_TSV)
    vs_sylph_rows = read_tsv(CANDIDATE_VS_SYLPH_TSV)
    vs_current_rows = read_tsv(CANDIDATE_VS_CURRENT_TSV)
    blend_audit_rows = read_tsv(BLEND_AUDIT_TSV) if BLEND_AUDIT_TSV.exists() else []
    blend_panel_rows = read_tsv(BLEND_PANEL_TSV) if BLEND_PANEL_TSV.exists() else []
    blend_guard_rows = read_tsv(BLEND_GUARD_AUDIT_TSV) if BLEND_GUARD_AUDIT_TSV.exists() else []
    blend_alpha_rows = read_tsv(BLEND_ALPHA_AUDIT_TSV) if BLEND_ALPHA_AUDIT_TSV.exists() else []
    blend_alpha_overall_rows = (
        read_tsv(BLEND_ALPHA_OVERALL_TSV) if BLEND_ALPHA_OVERALL_TSV.exists() else []
    )

    vs_sylph_by_panel = {row["panel"]: row for row in vs_sylph_rows}

    panel_rows: list[dict[str, str]] = []
    for row in delta_rows:
        panel = row["panel"]
        selected = vs_sylph_by_panel.get(panel, {})
        panel_rows.append(
            {
                "record_type": "panel_gap",
                "panel": panel,
                "main_gap_component": row["main_gap_component"],
                "selected_default_minus_sylph_L1_pp": fmt(
                    selected.get("delta_official_L1_pp", "nan")
                ),
                "legacy_minco_minus_sylph_L1_pp": fmt(
                    row["minco_minus_sylph_union_L1_pp"]
                ),
                "legacy_matched_abs_error_gap_pp": fmt(
                    row["minco_minus_sylph_matched_abs_error_pp"]
                ),
                "legacy_missing_truth_mass_gap_pp": fmt(
                    row["minco_minus_sylph_missing_truth_mass_pp"]
                ),
                "legacy_extra_pred_mass_gap_pp": fmt(
                    row["minco_minus_sylph_extra_pred_mass_pp"]
                ),
                "selected_default_minus_sylph_F1": fmt(
                    selected.get("delta_pooled_F1", "nan")
                ),
                "selected_default_minus_sylph_Pearson": fmt(
                    selected.get("delta_official_Pearson", "nan")
                ),
                "cached_allocator": "",
                "cached_allocator_mean_delta_current_L1_pp": "",
                "cached_allocator_max_worse_current_L1_pp": "",
                "cached_allocator_improved_panel_count": "",
                "decision": "abundance_gap_not_closed",
            }
        )

    best_allocator = best_cached_allocator(variant_rows)
    allocator_manifest = {
        "record_type": "cached_allocator",
        "panel": "all_cached_panels",
        "main_gap_component": "",
        "selected_default_minus_sylph_L1_pp": "",
        "legacy_minco_minus_sylph_L1_pp": "",
        "legacy_matched_abs_error_gap_pp": "",
        "legacy_missing_truth_mass_gap_pp": "",
        "legacy_extra_pred_mass_gap_pp": "",
        "selected_default_minus_sylph_F1": "",
        "selected_default_minus_sylph_Pearson": "",
        "cached_allocator": best_allocator.get("method", ""),
        "cached_allocator_mean_delta_current_L1_pp": fmt(
            best_allocator.get("mean_delta_current_L1_pp", "nan")
        )
        if best_allocator
        else "",
        "cached_allocator_max_worse_current_L1_pp": fmt(
            best_allocator.get("max_worse_current_L1_pp", "nan")
        )
        if best_allocator
        else "",
        "cached_allocator_improved_panel_count": best_allocator.get(
            "improved_panel_count", ""
        ),
        "decision": "cached_signal_not_promoted_without_raw_validation"
        if best_allocator
        else "no_cached_allocator_candidate",
    }

    promotion = audit_metric(blend_audit_rows, "promotion_decision")
    mean_l1 = audit_metric(blend_audit_rows, "mean_L1_delta_vs_candidate_preset_pp")
    max_worse = audit_metric(blend_audit_rows, "max_worse_L1_delta_vs_candidate_preset_pp")
    direction = audit_metric(blend_audit_rows, "sample_L1_direction")
    blend_manifest = {
        "record_type": "selected_default_blend_audit",
        "panel": "all_cached_profile_samples",
        "main_gap_component": "",
        "selected_default_minus_sylph_L1_pp": "",
        "legacy_minco_minus_sylph_L1_pp": "",
        "legacy_matched_abs_error_gap_pp": "",
        "legacy_missing_truth_mass_gap_pp": "",
        "legacy_extra_pred_mass_gap_pp": "",
        "selected_default_minus_sylph_F1": "",
        "selected_default_minus_sylph_Pearson": "",
        "cached_allocator": "default_candidate_preset_genus_xny_blend_a0.25"
        if blend_audit_rows
        else "",
        "cached_allocator_mean_delta_current_L1_pp": fmt(mean_l1.get("value", "nan"))
        if mean_l1
        else "",
        "cached_allocator_max_worse_current_L1_pp": fmt(max_worse.get("value", "nan"))
        if max_worse
        else "",
        "cached_allocator_improved_panel_count": direction.get("value", ""),
        "decision": promotion.get("decision", "selected_default_blend_audit_missing"),
    }
    guard_metrics = {row.get("metric", ""): row for row in blend_guard_rows}
    guard_manifest = {
        "record_type": "selected_default_blend_guard_audit",
        "panel": "all_cached_profile_samples",
        "main_gap_component": "",
        "selected_default_minus_sylph_L1_pp": "",
        "legacy_minco_minus_sylph_L1_pp": "",
        "legacy_matched_abs_error_gap_pp": "",
        "legacy_missing_truth_mass_gap_pp": "",
        "legacy_extra_pred_mass_gap_pp": "",
        "selected_default_minus_sylph_F1": "",
        "selected_default_minus_sylph_Pearson": "",
        "cached_allocator": "adaptive_default_candidate_preset_genus_xny_blend_guard"
        if blend_guard_rows
        else "",
        "cached_allocator_mean_delta_current_L1_pp": "",
        "cached_allocator_max_worse_current_L1_pp": "",
        "cached_allocator_improved_panel_count": guard_metrics.get("lopo_guard_result", {}).get(
            "value", ""
        ),
        "decision": guard_metrics.get("promotion_decision", {}).get(
            "decision", "selected_default_blend_guard_audit_missing"
        ),
    }
    guard_summary = {
        "tested_guard_rules": guard_metrics.get("tested_guard_rules", {}).get("value", ""),
        "sample_safe_guards": guard_metrics.get("sample_safe_guards", {}).get("value", ""),
        "lopo_guard_result": guard_metrics.get("lopo_guard_result", {}).get("value", ""),
        "decision": guard_manifest["decision"],
    }
    alpha_metrics = {row.get("metric", ""): row for row in blend_alpha_rows}
    alpha_manifest = {
        "record_type": "selected_default_blend_alpha_sweep",
        "panel": "all_cached_profile_samples",
        "main_gap_component": "",
        "selected_default_minus_sylph_L1_pp": "",
        "legacy_minco_minus_sylph_L1_pp": "",
        "legacy_matched_abs_error_gap_pp": "",
        "legacy_missing_truth_mass_gap_pp": "",
        "legacy_extra_pred_mass_gap_pp": "",
        "selected_default_minus_sylph_F1": "",
        "selected_default_minus_sylph_Pearson": "",
        "cached_allocator": "selected_default_genus_xny_alpha_cap_sweep"
        if blend_alpha_rows
        else "",
        "cached_allocator_mean_delta_current_L1_pp": "",
        "cached_allocator_max_worse_current_L1_pp": "",
        "cached_allocator_improved_panel_count": alpha_metrics.get("best_ranked_variant", {}).get(
            "value", ""
        ),
        "decision": alpha_metrics.get("promotion_decision", {}).get(
            "decision", "selected_default_blend_alpha_sweep_missing"
        ),
    }
    alpha_summary = {
        "variants_tested": alpha_metrics.get("variants_tested", {}).get("value", ""),
        "sample_safe_variants": alpha_metrics.get("sample_safe_variants", {}).get("value", ""),
        "best_ranked_variant": alpha_metrics.get("best_ranked_variant", {}).get("value", ""),
        "decision": alpha_manifest["decision"],
    }

    mean_delta_current_l1 = sum(
        fnum(row["delta_official_L1_pp"]) for row in vs_current_rows
    ) / len(vs_current_rows)
    mean_delta_sylph_l1 = sum(
        fnum(row["delta_official_L1_pp"]) for row in vs_sylph_rows
    ) / len(vs_sylph_rows)
    default_manifest = {
        "record_type": "default_decision",
        "panel": "matched_four_panel_summary",
        "main_gap_component": "matched_allocation_primary_except_cami3",
        "selected_default_minus_sylph_L1_pp": fmt(mean_delta_sylph_l1),
        "legacy_minco_minus_sylph_L1_pp": "",
        "legacy_matched_abs_error_gap_pp": "",
        "legacy_missing_truth_mass_gap_pp": "",
        "legacy_extra_pred_mass_gap_pp": "",
        "selected_default_minus_sylph_F1": "",
        "selected_default_minus_sylph_Pearson": "",
        "cached_allocator": "",
        "cached_allocator_mean_delta_current_L1_pp": fmt(mean_delta_current_l1),
        "cached_allocator_max_worse_current_L1_pp": "",
        "cached_allocator_improved_panel_count": "",
        "decision": "selected_default_is_best_minco_candidate_not_abundance_release_claim",
    }

    manifest_rows = panel_rows + [
        allocator_manifest,
        blend_manifest,
        guard_manifest,
        alpha_manifest,
        default_manifest,
    ]
    write_tsv(manifest_rows)
    write_markdown(
        panel_rows,
        allocator_manifest,
        blend_manifest,
        blend_panel_rows,
        guard_summary,
        alpha_summary,
        blend_alpha_overall_rows,
        default_manifest,
    )


if __name__ == "__main__":
    main()
