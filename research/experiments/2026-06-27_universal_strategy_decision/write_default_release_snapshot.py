#!/usr/bin/env python3
"""Write a consolidated snapshot for the current MinCO default strategy.

The snapshot is intentionally conservative. It summarizes cached evidence for
the selected candidate preset, compares it with the previous MinCO default and
Sylph on the matched GTDB panels, records diagnostic exact-split evidence, and
keeps experimental abundance switches out of the default unless their audits
are sample-safe.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

VS_CURRENT = RESULTS / "candidate_default_vs_current.tsv"
VS_SYLPH = RESULTS / "candidate_default_vs_sylph.tsv"
DEFAULT_MANIFEST = RESULTS / "current_default_strategy_manifest.tsv"
EXACT_DIAG = RESULTS / "cached_exactsplit_pairwise_summary.tsv"
ALLOCATOR_EXTERNAL = RESULTS / "feature_allocator_external_exactsplit_audit.tsv"
RUNTIME = RESULTS / "runtime_memory_minco_vs_sylph.tsv"
CAMI3_POLICY = RESULTS / "cami3_binomial_fallback_truth_policy_audit.tsv"
CAMI3_CANDIDATE_SUMMARY = RESULTS / "cami3_binomial_fallback_candidate_default_summary.tsv"
CAMI3_CANDIDATE_DELTA = RESULTS / "cami3_binomial_fallback_candidate_default_delta.tsv"

OUT_PANEL = RESULTS / "default_release_snapshot_panel_metrics.tsv"
OUT_DIAG = RESULTS / "default_release_snapshot_diagnostic_metrics.tsv"
OUT_RUNTIME = RESULTS / "default_release_snapshot_runtime.tsv"
OUT_AUDIT = RESULTS / "default_release_snapshot_audit.tsv"
OUT_MD = EXP / "DEFAULT_RELEASE_SNAPSHOT.md"


PANEL_GRADE = {
    "cami2_toy_mouse_gut": ("release", "clean GTDB species comparison"),
    "cami3_toy_human_gut_gtdb_source_readmap": (
        "diagnostic",
        "source-readmap GTDB transfer with unmapped read rows",
    ),
    "hmp_airskin_gtdb_source_abundance": (
        "release",
        "same-release GTDB source-abundance comparison",
    ),
    "hmp_gastrooral_gtdb_source_abundance": (
        "release",
        "same-release GTDB source-abundance comparison",
    ),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def metric_table(path: Path) -> dict[str, dict[str, str]]:
    return {row["metric"] if "metric" in row else row["item"]: row for row in read_tsv(path)}


def grade_for_panel(panel: str) -> tuple[str, str]:
    if panel == "cami3_toy_human_gut_gtdb_source_readmap" and cami3_policy_accepted():
        return "release", "accepted exact-binomial source-readmap GTDB truth"
    return PANEL_GRADE.get(panel, ("diagnostic", "grade not explicitly mapped"))


def cami3_policy_accepted() -> bool:
    if not CAMI3_POLICY.exists():
        return False
    policy = metric_table(CAMI3_POLICY)
    return (
        policy.get("truth_policy_decision", {}).get("value")
        == "accept_exact_binomial_fallback_for_release_truth"
    )


def apply_cami3_binomial_override(row: dict[str, object]) -> dict[str, object]:
    if row["panel"] != "cami3_toy_human_gut_gtdb_source_readmap":
        return row
    if not (cami3_policy_accepted() and CAMI3_CANDIDATE_SUMMARY.exists() and CAMI3_CANDIDATE_DELTA.exists()):
        return row
    summary = {item["method"]: item for item in read_tsv(CAMI3_CANDIDATE_SUMMARY)}
    minco = summary.get("minco_candidate_default_binomial_fallback", {})
    sylph = summary.get("sylph_binomial_fallback", {})
    delta_rows = [
        item
        for item in read_tsv(CAMI3_CANDIDATE_DELTA)
        if item["method"] == "minco_candidate_default_binomial_fallback"
    ]
    mean_delta_f1 = sum(finite(item["delta_F1"]) for item in delta_rows) / len(delta_rows)
    mean_delta_l1 = sum(finite(item["delta_L1_union_pp"]) for item in delta_rows) / len(delta_rows)
    mean_delta_pearson = sum(finite(item["delta_Pearson_union"]) for item in delta_rows) / len(delta_rows)
    candidate_f1 = finite(minco.get("pooled_F1"))
    sylph_f1 = finite(sylph.get("pooled_F1"))
    candidate_l1 = finite(minco.get("mean_L1_union_pp"))
    sylph_l1 = finite(sylph.get("mean_L1_union_pp"))
    candidate_pearson = finite(minco.get("mean_Pearson_union"))
    sylph_pearson = finite(sylph.get("mean_Pearson_union"))
    row.update(
        {
            "grade": "release",
            "grade_caveat": "accepted exact-binomial source-readmap GTDB truth",
            "candidate_pooled_F1": candidate_f1,
            "sylph_pooled_F1": sylph_f1,
            "delta_F1_vs_sylph": candidate_f1 - sylph_f1,
            "candidate_L1_pp": candidate_l1,
            "sylph_L1_pp": sylph_l1,
            "delta_L1_pp_vs_sylph": candidate_l1 - sylph_l1,
            "candidate_Pearson": candidate_pearson,
            "sylph_Pearson": sylph_pearson,
            "delta_Pearson_vs_sylph": candidate_pearson - sylph_pearson,
            "wins_F1_vs_sylph": int(candidate_f1 > sylph_f1),
            "wins_L1_vs_sylph": int(candidate_l1 < sylph_l1),
            "wins_Pearson_vs_sylph": int(candidate_pearson > sylph_pearson),
            "delta_F1_vs_previous_minco": mean_delta_f1,
            "delta_L1_pp_vs_previous_minco": mean_delta_l1,
            "delta_Pearson_vs_previous_minco": mean_delta_pearson,
        }
    )
    return row


def build_panel_rows() -> list[dict[str, object]]:
    current_by_panel = {row["panel"]: row for row in read_tsv(VS_CURRENT)}
    rows: list[dict[str, object]] = []
    for row in read_tsv(VS_SYLPH):
        panel = row["panel"]
        grade, caveat = grade_for_panel(panel)
        current = current_by_panel.get(panel, {})
        f1_delta = finite(row["delta_pooled_F1"])
        l1_delta = finite(row["delta_official_L1_pp"])
        pearson_delta = finite(row["delta_official_Pearson"])
        rows.append(
            apply_cami3_binomial_override(
                {
                "panel": panel,
                "panel_label": row["panel_label"],
                "samples": row["samples"],
                "grade": grade,
                "grade_caveat": caveat,
                "candidate_pooled_F1": finite(row["candidate_pooled_F1"]),
                "sylph_pooled_F1": finite(row["baseline_pooled_F1"]),
                "delta_F1_vs_sylph": f1_delta,
                "candidate_L1_pp": finite(row["candidate_L1_pp"]),
                "sylph_L1_pp": finite(row["baseline_L1_pp"]),
                "delta_L1_pp_vs_sylph": l1_delta,
                "candidate_Pearson": finite(row["candidate_Pearson"]),
                "sylph_Pearson": finite(row["baseline_Pearson"]),
                "delta_Pearson_vs_sylph": pearson_delta,
                "wins_F1_vs_sylph": int(f1_delta > 0.0),
                "wins_L1_vs_sylph": int(l1_delta < 0.0),
                "wins_Pearson_vs_sylph": int(pearson_delta > 0.0),
                "delta_F1_vs_previous_minco": finite(current.get("delta_pooled_F1", "")),
                "delta_L1_pp_vs_previous_minco": finite(
                    current.get("delta_official_L1_pp", "")
                ),
                "delta_Pearson_vs_previous_minco": finite(
                    current.get("delta_official_Pearson", "")
                ),
            }
            )
        )
    return rows


def build_diag_rows() -> list[dict[str, object]]:
    rows = []
    for row in read_tsv(EXACT_DIAG):
        rows.append(
            {
                "dataset": row["dataset"],
                "samples": row["samples"],
                "comparable_samples": int(float(row["comparable_samples"])),
                "mean_delta_F1_minco_minus_sylph": finite(
                    row["mean_delta_F1_minco_minus_sylph"]
                ),
                "mean_delta_L1_minco_minus_sylph": finite(
                    row["mean_delta_L1_minco_minus_sylph"]
                ),
                "mean_delta_Pearson_minco_minus_sylph": finite(
                    row["mean_delta_Pearson_minco_minus_sylph"]
                ),
                "F1_wins": int(float(row["F1_wins"])),
                "F1_losses": int(float(row["F1_losses"])),
                "L1_wins": int(float(row["L1_wins"])),
                "L1_losses": int(float(row["L1_losses"])),
                "Pearson_wins": int(float(row["Pearson_wins"])),
                "Pearson_losses": int(float(row["Pearson_losses"])),
                "grade": "diagnostic_nonrelease",
            }
        )
    return rows


def preferred_minco_runtime(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    candidates: dict[tuple[str, str], list[tuple[int, dict[str, str]]]] = defaultdict(list)
    for row in rows:
        method = row["method"]
        if not method.startswith("MinCO"):
            continue
        if "current concurrent default wrapper" in method:
            rank = 0
        elif "current strategy components" in method:
            rank = 1
        elif "default wrapper/autoexact" in method:
            rank = 2
        else:
            continue
        candidates[(row["dataset"], row["sample"])].append((rank, row))
    return {key: sorted(values, key=lambda item: item[0])[0][1] for key, values in candidates.items()}


def build_runtime_rows() -> tuple[list[dict[str, object]], dict[str, float]]:
    rows = read_tsv(RUNTIME)
    minco = preferred_minco_runtime(rows)
    sylph = {
        (row["dataset"], row["sample"]): row
        for row in rows
        if row["method"] == "Sylph sketch+profile"
    }
    out: list[dict[str, object]] = []
    for key in sorted(set(minco) & set(sylph)):
        m = minco[key]
        s = sylph[key]
        m_sec = finite(m["seconds"])
        s_sec = finite(s["seconds"])
        m_rss = finite(m["peak_rss_gib"])
        s_rss = finite(s["peak_rss_gib"])
        out.append(
            {
                "dataset": key[0],
                "sample": key[1],
                "minco_method": m["method"],
                "sylph_method": s["method"],
                "minco_seconds": m_sec,
                "sylph_seconds": s_sec,
                "seconds_ratio_minco_over_sylph": m_sec / s_sec if s_sec else 0.0,
                "minco_peak_rss_gib": m_rss,
                "sylph_peak_rss_gib": s_rss,
                "rss_ratio_minco_over_sylph": m_rss / s_rss if s_rss else 0.0,
                "source_caveat": m.get("source_caveat", ""),
            }
        )
    summary = {
        "comparable_runtime_samples": float(len(out)),
        "mean_minco_seconds": sum(row["minco_seconds"] for row in out) / len(out) if out else 0.0,
        "mean_sylph_seconds": sum(row["sylph_seconds"] for row in out) / len(out) if out else 0.0,
        "mean_minco_peak_rss_gib": sum(row["minco_peak_rss_gib"] for row in out) / len(out)
        if out
        else 0.0,
        "mean_sylph_peak_rss_gib": sum(row["sylph_peak_rss_gib"] for row in out) / len(out)
        if out
        else 0.0,
    }
    summary["mean_seconds_ratio_minco_over_sylph"] = (
        summary["mean_minco_seconds"] / summary["mean_sylph_seconds"]
        if summary["mean_sylph_seconds"]
        else 0.0
    )
    summary["mean_rss_ratio_minco_over_sylph"] = (
        summary["mean_minco_peak_rss_gib"] / summary["mean_sylph_peak_rss_gib"]
        if summary["mean_sylph_peak_rss_gib"]
        else 0.0
    )
    return out, summary


def build_audit(
    panel_rows: list[dict[str, object]],
    diag_rows: list[dict[str, object]],
    runtime_summary: dict[str, float],
) -> list[dict[str, object]]:
    release = [row for row in panel_rows if row["grade"] == "release"]
    diagnostic = [row for row in panel_rows if row["grade"] != "release"]
    f1_wins = sum(int(row["wins_F1_vs_sylph"]) for row in panel_rows)
    l1_wins = sum(int(row["wins_L1_vs_sylph"]) for row in panel_rows)
    pearson_wins = sum(int(row["wins_Pearson_vs_sylph"]) for row in panel_rows)
    current_f1_worse = sum(1 for row in panel_rows if float(row["delta_F1_vs_previous_minco"]) < -1e-12)
    current_l1_worse = sum(1 for row in panel_rows if float(row["delta_L1_pp_vs_previous_minco"]) > 1e-12)
    current_pearson_worse = sum(
        1 for row in panel_rows if float(row["delta_Pearson_vs_previous_minco"]) < -1e-12
    )
    manifest = metric_table(DEFAULT_MANIFEST)
    allocator = metric_table(ALLOCATOR_EXTERNAL)
    return [
        {
            "metric": "selected_default",
            "value": manifest.get("selected_profile_preset", {}).get("value", "NA"),
            "evidence": "current_default_strategy_manifest.tsv",
            "decision": "user_default_entrypoint"
            if manifest.get("selected_profile_preset", {}).get("value") == "candidate"
            else "review_default_entrypoint",
        },
        {
            "metric": "matched_panel_scope",
            "value": (
                f"panels={len(panel_rows)};release={len(release)};"
                f"diagnostic={len(diagnostic)}"
            ),
            "evidence": "candidate_default_vs_sylph.tsv",
            "decision": "cached_gtdb_default_evidence",
        },
        {
            "metric": "candidate_vs_previous_minco",
            "value": (
                f"F1_regressions={current_f1_worse};L1_regressions={current_l1_worse};"
                f"Pearson_regressions={current_pearson_worse}"
            ),
            "evidence": "candidate_default_vs_current.tsv",
            "decision": "candidate_dominates_previous_cached_default"
            if current_f1_worse == current_l1_worse == current_pearson_worse == 0
            else "review_candidate_regressions",
        },
        {
            "metric": "candidate_vs_sylph",
            "value": (
                f"F1_wins={f1_wins}/{len(panel_rows)};"
                f"L1_wins={l1_wins}/{len(panel_rows)};"
                f"Pearson_wins={pearson_wins}/{len(panel_rows)}"
            ),
            "evidence": "candidate_default_vs_sylph.tsv",
            "decision": "not_broad_sylph_beating",
        },
        {
            "metric": "external_exactsplit_diagnostic",
            "value": (
                f"datasets={len(diag_rows)};"
                f"F1_win_datasets={sum(1 for row in diag_rows if row['F1_wins'] > row['F1_losses'])};"
                f"L1_win_datasets={sum(1 for row in diag_rows if row['L1_wins'] > row['L1_losses'])}"
            ),
            "evidence": "cached_exactsplit_pairwise_summary.tsv",
            "decision": "diagnostic_only_nonrelease",
        },
        {
            "metric": "guarded_allocator_external_status",
            "value": allocator.get("candidate_effect", {}).get("value", "NA"),
            "evidence": "feature_allocator_external_exactsplit_audit.tsv",
            "decision": allocator.get("promotion_decision", {}).get(
                "value", "keep_feature_allocator_experimental_off_by_default"
            ),
        },
        {
            "metric": "runtime_tradeoff",
            "value": (
                f"samples={runtime_summary['comparable_runtime_samples']:.0f};"
                f"seconds_ratio={runtime_summary['mean_seconds_ratio_minco_over_sylph']:.3f};"
                f"rss_ratio={runtime_summary['mean_rss_ratio_minco_over_sylph']:.3f}"
            ),
            "evidence": "runtime_memory_minco_vs_sylph.tsv",
            "decision": "minco_slower_lower_memory_in_cached_timed_runs",
        },
        {
            "metric": "release_decision",
            "value": (
                "selected_candidate_default;experimental_allocators_off;"
                "do_not_claim_broad_sylph_beating"
            ),
            "evidence": "default_release_snapshot_audit.tsv",
            "decision": "stable_default_candidate_not_final_universal_claim",
        },
    ]


def write_markdown(
    panel_rows: list[dict[str, object]],
    diag_rows: list[dict[str, object]],
    runtime_rows: list[dict[str, object]],
    runtime_summary: dict[str, float],
    audit_rows: list[dict[str, object]],
) -> None:
    audit = {str(row["metric"]): row for row in audit_rows}
    lines = [
        "# Default Release Snapshot",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note consolidates cached evidence for the current MinCO",
        "default. It separates the selected user default from broader release",
        "claims. No raw profiling jobs are run by this snapshot.",
        "",
        "## Decision",
        "",
        f"- Selected default: `{audit['selected_default']['value']}` preset via `scripts/minco_profile_default.py`.",
        "- Keep experimental abundance allocators off by default.",
        "- Current evidence supports the candidate preset as the best cached MinCO",
        "  default versus the previous MinCO default, but not a broad Sylph-beating",
        "  abundance claim.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit_rows:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Matched Panels",
            "",
            "| Panel | Grade | Samples | Delta F1 vs Sylph | Delta L1 pp vs Sylph | Delta Pearson vs Sylph |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for row in panel_rows:
        lines.append(
            "| {label} | {grade} | {samples} | {f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                label=row["panel_label"],
                grade=row["grade"],
                samples=row["samples"],
                f1=float(row["delta_F1_vs_sylph"]),
                l1=float(row["delta_L1_pp_vs_sylph"]),
                pearson=float(row["delta_Pearson_vs_sylph"]),
            )
        )

    lines.extend(
        [
            "",
            "## Diagnostic Exact-Split Panels",
            "",
            "| Dataset | Samples | Mean delta F1 | Mean delta L1 | Mean delta Pearson |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in diag_rows:
        lines.append(
            "| {dataset} | {samples} | {f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                dataset=row["dataset"],
                samples=row["samples"],
                f1=float(row["mean_delta_F1_minco_minus_sylph"]),
                l1=float(row["mean_delta_L1_minco_minus_sylph"]),
                pearson=float(row["mean_delta_Pearson_minco_minus_sylph"]),
            )
        )

    lines.extend(
        [
            "",
            "## Runtime",
            "",
            "| Samples | Mean MinCO seconds | Mean Sylph seconds | Seconds ratio | Mean MinCO RSS GiB | Mean Sylph RSS GiB | RSS ratio |",
            "|---:|---:|---:|---:|---:|---:|---:|",
            "| {n:.0f} | {ms:.2f} | {ss:.2f} | {sr:.3f} | {mr:.3f} | {srss:.3f} | {rr:.3f} |".format(
                n=runtime_summary["comparable_runtime_samples"],
                ms=runtime_summary["mean_minco_seconds"],
                ss=runtime_summary["mean_sylph_seconds"],
                sr=runtime_summary["mean_seconds_ratio_minco_over_sylph"],
                mr=runtime_summary["mean_minco_peak_rss_gib"],
                srss=runtime_summary["mean_sylph_peak_rss_gib"],
                rr=runtime_summary["mean_rss_ratio_minco_over_sylph"],
            ),
            "",
            "## Outputs",
            "",
            f"- `{OUT_PANEL.relative_to(EXP)}`",
            f"- `{OUT_DIAG.relative_to(EXP)}`",
            f"- `{OUT_RUNTIME.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    panel_rows = build_panel_rows()
    diag_rows = build_diag_rows()
    runtime_rows, runtime_summary = build_runtime_rows()
    audit = build_audit(panel_rows, diag_rows, runtime_summary)

    write_tsv(
        OUT_PANEL,
        panel_rows,
        [
            "panel",
            "panel_label",
            "samples",
            "grade",
            "grade_caveat",
            "candidate_pooled_F1",
            "sylph_pooled_F1",
            "delta_F1_vs_sylph",
            "candidate_L1_pp",
            "sylph_L1_pp",
            "delta_L1_pp_vs_sylph",
            "candidate_Pearson",
            "sylph_Pearson",
            "delta_Pearson_vs_sylph",
            "wins_F1_vs_sylph",
            "wins_L1_vs_sylph",
            "wins_Pearson_vs_sylph",
            "delta_F1_vs_previous_minco",
            "delta_L1_pp_vs_previous_minco",
            "delta_Pearson_vs_previous_minco",
        ],
    )
    write_tsv(
        OUT_DIAG,
        diag_rows,
        [
            "dataset",
            "samples",
            "comparable_samples",
            "mean_delta_F1_minco_minus_sylph",
            "mean_delta_L1_minco_minus_sylph",
            "mean_delta_Pearson_minco_minus_sylph",
            "F1_wins",
            "F1_losses",
            "L1_wins",
            "L1_losses",
            "Pearson_wins",
            "Pearson_losses",
            "grade",
        ],
    )
    write_tsv(
        OUT_RUNTIME,
        runtime_rows,
        [
            "dataset",
            "sample",
            "minco_method",
            "sylph_method",
            "minco_seconds",
            "sylph_seconds",
            "seconds_ratio_minco_over_sylph",
            "minco_peak_rss_gib",
            "sylph_peak_rss_gib",
            "rss_ratio_minco_over_sylph",
            "source_caveat",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(panel_rows, diag_rows, runtime_rows, runtime_summary, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
