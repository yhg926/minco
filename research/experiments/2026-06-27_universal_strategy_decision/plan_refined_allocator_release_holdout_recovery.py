#!/usr/bin/env python3
"""Plan release-grade holdout recovery for the refined allocator.

The refined guarded allocator has strong cached score-replay support, but it is
not default because the current cache has no unused release-grade holdout
profile set. This script turns that blocker into ranked recovery actions using
only local file-presence audits and cached summaries.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

INVENTORY = RESULTS / "feature_allocator_refined_independent_holdout_inventory.tsv"
INVENTORY_AUDIT = RESULTS / "feature_allocator_refined_independent_holdout_inventory_audit.tsv"
CAMI3_EXTENSION_DETAIL = RESULTS / "cami3_source_readmap_extension_cache_detail.tsv"
CAMI3_EXTENSION_AUDIT = RESULTS / "cami3_source_readmap_extension_cache_audit.tsv"
CAMI3_RECOVERY_INPUTS_AUDIT = RESULTS / "cami3_source_readmap_recovery_inputs_audit.tsv"
CAMI3_AFTER_RECOVERY_AUDIT = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
MARINE_AUDIT = RESULTS / "feature_allocator_refined_marine_diagnostic_audit.tsv"
SOURCE_PROFILE_AUDIT = RESULTS / "cami3_source_profile_binomial_extension_audit.tsv"

OUT_PLAN = RESULTS / "feature_allocator_refined_holdout_recovery_plan.tsv"
OUT_RUNTIME = RESULTS / "feature_allocator_refined_holdout_recovery_runtime.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_refined_holdout_recovery_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_REFINED_HOLDOUT_RECOVERY_PLAN.md"


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


def table_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def count_true(rows: list[dict[str, str]], key: str) -> int:
    return sum(str(row.get(key, "")).strip().lower() == "true" for row in rows)


def parse_count(value: str) -> tuple[int, int]:
    left, right = str(value).split("/", 1)
    return int(left), int(right)


def bool_from_text(text: str, prefix: str) -> bool:
    return f"{prefix}=missing" not in text and f"{prefix}=" in text


def parse_time_log(path: Path) -> dict[str, str]:
    out = {"wall_clock": "", "max_rss_kb": "", "exit_status": ""}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            out["wall_clock"] = line.rsplit("):", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            out["max_rss_kb"] = line.rsplit(":", 1)[-1].strip()
        elif "Exit status" in line:
            out["exit_status"] = line.rsplit(":", 1)[-1].strip()
    return out


def cami3_runtime_rows(extension_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    out = []
    for row in extension_rows:
        sample = int(row["sample"])
        profile = Path(row["selected_default_profile_path_or_checked"].split(";", 1)[0])
        time_log = profile.with_suffix(".time.log")
        if not time_log.is_file():
            # The replay script writes time logs beside the output using this
            # explicit naming convention.
            time_log = profile.parent / (
                f"cami3_toy_human_gut_gtdb_source_readmap.sample{sample}."
                "candidate_preset.time.log"
            )
        timing = parse_time_log(time_log)
        out.append(
            {
                "route_id": "cami3_source_readmap_extension_3_5",
                "sample": sample,
                "profile": str(profile),
                "time_log": str(time_log),
                "wall_clock": timing["wall_clock"],
                "max_rss_kb": timing["max_rss_kb"],
                "exit_status": timing["exit_status"],
            }
        )
    return out


def build_plan() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    inventory = read_tsv(INVENTORY)
    inventory_audit = table_by(INVENTORY_AUDIT, "metric")
    extension_detail = read_tsv(CAMI3_EXTENSION_DETAIL)
    extension_audit = table_by(CAMI3_EXTENSION_AUDIT, "metric")
    recovery_inputs_audit = (
        table_by(CAMI3_RECOVERY_INPUTS_AUDIT, "metric")
        if CAMI3_RECOVERY_INPUTS_AUDIT.exists()
        else {}
    )
    after_recovery_audit = (
        table_by(CAMI3_AFTER_RECOVERY_AUDIT, "metric")
        if CAMI3_AFTER_RECOVERY_AUDIT.exists()
        else {}
    )
    marine_audit = table_by(MARINE_AUDIT, "metric")
    source_profile_audit = table_by(SOURCE_PROFILE_AUDIT, "metric") if SOURCE_PROFILE_AUDIT.exists() else {}

    cami3_extension = [
        row for row in extension_detail if row["scope"] == "extension_candidate_subset"
    ]
    cami3_samples = ",".join(row["sample"] for row in cami3_extension)
    cami3_source = extension_audit["extension_samples_with_source_readmap_truth"]["value"]
    cami3_selected = extension_audit["extension_samples_with_selected_default_profiles"]["value"]
    cami3_raw = extension_audit["extension_samples_with_raw_tables"]["value"]
    cami3_baseline = extension_audit["extension_samples_with_sylph_profiles"]["value"]
    cami3_blockers = extension_audit["extension_blockers"]["value"]
    cami3_recovery_artifacts = recovery_inputs_audit.get("rescoring_artifacts", {}).get(
        "value", "not_run"
    )
    cami3_recovery_inputs = recovery_inputs_audit.get("recovery_blockers", {}).get(
        "value", "not_run"
    )
    cami3_size_estimate = recovery_inputs_audit.get("compressed_archive_size_estimate", {}).get(
        "value", "not_run"
    )
    cami3_after_recovery_status = after_recovery_audit.get("post_recovery_input_status", {}).get(
        "value", "not_run"
    )
    cami3_after_recovery_promotion = after_recovery_audit.get("promotion_decision", {}).get(
        "value", "not_run"
    )
    cami3_after_recovery_decision = after_recovery_audit.get("promotion_decision", {}).get(
        "decision", "not_run"
    )
    cami3_after_recovery_comparison = after_recovery_audit.get(
        "candidate_minco_vs_sylph", {}
    ).get("value", "not_run")
    cami3_after_recovery_effect = after_recovery_audit.get("refined_allocator_effect", {}).get(
        "value", "not_run"
    )
    cami3_after_recovery_scored = cami3_after_recovery_promotion == "post_recovery_extension_scored"
    cami3_after_recovery_blocker = after_recovery_audit.get("blocking_reason", {}).get(
        "value", "none" if cami3_after_recovery_scored else "not_run"
    )
    cami3_runtime = cami3_runtime_rows(cami3_extension)

    hmp_omitted = [
        row for row in inventory if row["candidate_set"].startswith("hmp_airskin_omitted_sample")
    ]
    hmp_samples = ",".join(row["samples"] for row in hmp_omitted)
    hmp_truth = sum(row["truth_status"] == "truth_file_present" for row in hmp_omitted)
    hmp_minco = sum(bool_from_text(row["evidence"], "minco") for row in hmp_omitted)
    hmp_sylph = sum(bool_from_text(row["evidence"], "sylph") for row in hmp_omitted)
    hmp_reads = sum(bool_from_text(row["evidence"], "read_input") for row in hmp_omitted)

    marine_effect = marine_audit["refined_allocator_marine_effect"]["value"]
    marine_scope = marine_audit["diagnostic_scope"]["value"]
    inventory_decision = inventory_audit["promotion_decision"]["value"]
    source_profile_scope = source_profile_audit.get("source_profile_truth_scope", {}).get(
        "value", "not_run"
    )
    source_profile_effect = source_profile_audit.get("refined_allocator_effect", {}).get(
        "value", "not_run"
    )
    source_profile_decision = source_profile_audit.get("promotion_decision", {}).get(
        "decision", "not_run"
    )

    plan = [
        {
            "priority_rank": 10,
            "route_id": "cami3_source_readmap_extension_3_5",
            "samples": cami3_samples,
            "sample_count": len(cami3_extension),
            "route_type": "recover_missing_truth_cache",
            "current_status": (
                f"post_recovery_scored={cami3_after_recovery_scored};"
                f"candidate_vs_sylph={cami3_after_recovery_comparison};"
                f"refined_effect={cami3_after_recovery_effect}"
                if cami3_after_recovery_scored
                else (
                    f"source_truth={cami3_source};selected_default={cami3_selected};"
                    f"raw_tables={cami3_raw};baseline={cami3_baseline};"
                    f"recovery_inputs={cami3_recovery_inputs}"
                )
            ),
            "remaining_blocker": (
                "extension_result_does_not_support_default_promotion"
                if cami3_after_recovery_scored
                else cami3_blockers
            ),
            "local_work_done": (
                "readmap_only_recovery_and_scoring_complete"
                if cami3_after_recovery_scored
                else "selected_default_replay_profiles_present"
            ),
            "next_action": (
                "Do not promote the refined allocator from this route; the independent "
                "extension result has no refined-allocator gain and favors the external "
                "baseline on F1 and abundance. Move to the next independent holdout route."
                if cami3_after_recovery_scored
                else (
                    "Recover or generate source-readmap truth cache for samples3-5, "
                    "apply the accepted exact-binomial GTDB fallback policy, and rescore "
                    "the already generated candidate-preset profiles in the same namespace. "
                    f"Current input audit: {cami3_recovery_artifacts}; {cami3_size_estimate}. "
                    "After readmap recovery, run "
                    "score_cami3_source_readmap_extension_after_recovery.py to score "
                    "selected-default MinCO, the refined allocator, and Sylph."
                )
            ),
            "release_use_rule": (
                "Use for refined-allocator promotion only if profile-scope mapped truth "
                "coverage is release-grade and the allocator preserves F1 with no "
                "sample-level L1 regressions on the independent samples."
            ),
            "decision": (
                "completed_negative_holdout_result"
                if cami3_after_recovery_scored
                else "top_feasible_recovery_path"
            ),
            "evidence": (
                f"{CAMI3_EXTENSION_DETAIL.relative_to(EXP)};"
                f"{CAMI3_EXTENSION_AUDIT.relative_to(EXP)};"
                f"{CAMI3_RECOVERY_INPUTS_AUDIT.relative_to(EXP)};"
                f"{CAMI3_AFTER_RECOVERY_AUDIT.relative_to(EXP)}"
            ),
            "post_recovery_runner_status": (
                f"{cami3_after_recovery_status};blocker={cami3_after_recovery_blocker};"
                f"promotion={cami3_after_recovery_promotion};decision={cami3_after_recovery_decision}"
            ),
        },
        {
            "priority_rank": 20,
            "route_id": "hmp_airskin_omitted_profile_pair_recovery",
            "samples": hmp_samples,
            "sample_count": len(hmp_omitted),
            "route_type": "recover_profile_pairs",
            "current_status": (
                f"truth={hmp_truth}/{len(hmp_omitted)};minco={hmp_minco}/{len(hmp_omitted)};"
                f"sylph={hmp_sylph}/{len(hmp_omitted)};read_input={hmp_reads}/{len(hmp_omitted)}"
            ),
            "remaining_blocker": "missing_minco_or_sylph_profile_pairs",
            "local_work_done": "truth_files_present",
            "next_action": (
                "Recover read inputs, then run same-release selected-default MinCO and "
                "baseline profiles for the omitted sample IDs."
            ),
            "release_use_rule": (
                "Use only after both profile types exist locally in the same GTDB source-"
                "abundance namespace."
            ),
            "decision": "blocked_by_missing_local_reads_and_profiles",
            "evidence": str(INVENTORY.relative_to(EXP)),
        },
        {
            "priority_rank": 25,
            "route_id": "cami3_source_profile_binomial_extension",
            "samples": "3,4,5",
            "sample_count": 3,
            "route_type": "tested_local_truth_substitute",
            "current_status": f"{source_profile_scope};{source_profile_effect}",
            "remaining_blocker": "insufficient_source_profile_gtdb_truth_mapping",
            "local_work_done": "source_profile_binomial_extension_scored",
            "next_action": (
                "Do not use the taxonomic-profile source rows as a release substitute; "
                "continue with per-read source-readmap recovery."
            ),
            "release_use_rule": (
                "Rejected for release promotion unless a separate policy and stronger "
                "GTDB mapping can raise in-scope truth coverage to the release threshold."
            ),
            "decision": source_profile_decision,
            "evidence": str(SOURCE_PROFILE_AUDIT.relative_to(EXP)),
        },
        {
            "priority_rank": 30,
            "route_id": "marine_diagnostic_truth_upgrade",
            "samples": "marine0,marine2",
            "sample_count": 2,
            "route_type": "upgrade_truth_mapping",
            "current_status": f"{marine_scope};{marine_effect}",
            "remaining_blocker": "partial_truth_mapping_nonrelease",
            "local_work_done": "independent_diagnostic_replay_supports_allocator",
            "next_action": (
                "Improve the GTDB species truth transfer before using marine as a "
                "release-grade abundance holdout."
            ),
            "release_use_rule": (
                "Keep diagnostic until mapped Bacteria/Archaea truth mass reaches the "
                "release threshold and selected-default profiles are in the same namespace."
            ),
            "decision": "diagnostic_only",
            "evidence": str(MARINE_AUDIT.relative_to(EXP)),
        },
        {
            "priority_rank": 40,
            "route_id": "new_clean_gtdb_holdout",
            "samples": "",
            "sample_count": 0,
            "route_type": "new_release_panel",
            "current_status": inventory_decision,
            "remaining_blocker": "no_unused_release_grade_profile_set",
            "local_work_done": "none",
            "next_action": (
                "Add a new clean GTDB truth panel only if CAMI3 or HMP recovery cannot "
                "produce independent release-grade profile pairs."
            ),
            "release_use_rule": "Same promotion rule as other independent release holdouts.",
            "decision": "fallback_after_local_recovery_routes",
            "evidence": str(INVENTORY_AUDIT.relative_to(EXP)),
        },
    ]
    return plan, cami3_runtime


def build_audit(plan: list[dict[str, object]], runtime: list[dict[str, object]]) -> list[dict[str, object]]:
    top = plan[0]
    source_profile = next(
        (row for row in plan if row["route_id"] == "cami3_source_profile_binomial_extension"),
        {},
    )
    max_rss = [
        int(row["max_rss_kb"])
        for row in runtime
        if str(row.get("max_rss_kb", "")).isdigit()
    ]
    exits = [str(row.get("exit_status", "")) for row in runtime if row.get("exit_status", "") != ""]
    return [
        {
            "metric": "top_recovery_route",
            "value": top["route_id"],
            "evidence": str(OUT_PLAN.relative_to(EXP)),
            "decision": top["decision"],
        },
        {
            "metric": "completed_local_step",
            "value": str(top["current_status"]),
            "evidence": (
                str(CAMI3_AFTER_RECOVERY_AUDIT.relative_to(EXP))
                if top["decision"] == "completed_negative_holdout_result"
                else str(CAMI3_EXTENSION_AUDIT.relative_to(EXP))
            ),
            "decision": (
                "readmap_recovery_and_scoring_complete"
                if top["decision"] == "completed_negative_holdout_result"
                else "selected_default_replay_complete_for_cami3_extension"
            ),
        },
        {
            "metric": "remaining_blocker",
            "value": str(top["remaining_blocker"]),
            "evidence": (
                str(CAMI3_AFTER_RECOVERY_AUDIT.relative_to(EXP))
                if top["decision"] == "completed_negative_holdout_result"
                else str(CAMI3_EXTENSION_AUDIT.relative_to(EXP))
            ),
            "decision": (
                "route_completed_but_negative_for_default_promotion"
                if top["decision"] == "completed_negative_holdout_result"
                else "recover_source_readmap_truth_cache_before_release_scoring"
            ),
        },
        {
            "metric": "source_readmap_recovery_inputs",
            "value": str(top["current_status"]),
            "evidence": str(CAMI3_RECOVERY_INPUTS_AUDIT.relative_to(EXP)),
            "decision": (
                "readmaps_restored_and_scored"
                if top["decision"] == "completed_negative_holdout_result"
                else "readmaps_missing_but_rescoring_artifacts_ready"
            ),
        },
        {
            "metric": "post_recovery_runner",
            "value": str(top.get("post_recovery_runner_status", "")),
            "evidence": str(CAMI3_AFTER_RECOVERY_AUDIT.relative_to(EXP)),
            "decision": (
                "post_recovery_extension_scored_negative"
                if top["decision"] == "completed_negative_holdout_result"
                else "runner_ready_blocked_only_by_readmaps"
            ),
        },
        {
            "metric": "tested_local_substitute",
            "value": str(source_profile.get("current_status", "")),
            "evidence": str(SOURCE_PROFILE_AUDIT.relative_to(EXP)),
            "decision": str(source_profile.get("decision", "not_run")),
        },
        {
            "metric": "candidate_replay_runtime",
            "value": (
                f"samples={len(runtime)};exit_statuses={','.join(exits)};"
                f"max_rss_kb={max(max_rss) if max_rss else ''}"
            ),
            "evidence": str(OUT_RUNTIME.relative_to(EXP)),
            "decision": "local_table_mode_replay_completed",
        },
        {
            "metric": "default_decision",
            "value": "candidate_preset_default_with_refined_allocator_opt_in",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "do_not_promote_refined_allocator_without_independent_release_holdout",
        },
    ]


def write_markdown(plan: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Refined Allocator Release Holdout Recovery Plan",
        "",
        "Date: 2026-06-29",
        "",
        "This generated note ranks local recovery routes for the refined guarded",
        "abundance allocator. It does not promote the allocator or change the",
        "selected default.",
        "",
        "## Decision",
        "",
        f"- Top route: `{audit_by_metric['top_recovery_route']['value']}` "
        f"({audit_by_metric['top_recovery_route']['decision']}).",
        f"- Completed local step: `{audit_by_metric['completed_local_step']['value']}`.",
        f"- Remaining blocker: `{audit_by_metric['remaining_blocker']['value']}`.",
        f"- Source-readmap recovery inputs: `{audit_by_metric['source_readmap_recovery_inputs']['value']}`.",
        f"- Post-recovery scorer: `{audit_by_metric['post_recovery_runner']['value']}`.",
        f"- Tested local substitute: `{audit_by_metric['tested_local_substitute']['value']}` "
        f"({audit_by_metric['tested_local_substitute']['decision']}).",
        f"- Replay runtime: `{audit_by_metric['candidate_replay_runtime']['value']}`.",
        f"- Default decision: `{audit_by_metric['default_decision']['decision']}`.",
        "",
        "## Ranked Routes",
        "",
        "| Rank | Route | Samples | Status | Next Action |",
        "|---:|---|---|---|---|",
    ]
    for row in plan:
        lines.append(
            "| {rank} | `{route}` | {samples} | {status} | {action} |".format(
                rank=row["priority_rank"],
                route=row["route_id"],
                samples=row["samples"],
                status=str(row["current_status"]).replace("|", "/"),
                action=str(row["next_action"]).replace("|", "/"),
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_PLAN.relative_to(EXP)}`",
            f"- `{OUT_RUNTIME.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    plan, runtime = build_plan()
    audit = build_audit(plan, runtime)
    write_tsv(
        OUT_PLAN,
        plan,
        [
            "priority_rank",
            "route_id",
            "samples",
            "sample_count",
            "route_type",
            "current_status",
            "remaining_blocker",
            "local_work_done",
            "next_action",
            "release_use_rule",
            "decision",
            "evidence",
            "post_recovery_runner_status",
        ],
    )
    write_tsv(
        OUT_RUNTIME,
        runtime,
        ["route_id", "sample", "profile", "time_log", "wall_clock", "max_rss_kb", "exit_status"],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(plan, audit)
    print("\n".join(f"{row['priority_rank']}\t{row['route_id']}\t{row['decision']}" for row in plan))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
