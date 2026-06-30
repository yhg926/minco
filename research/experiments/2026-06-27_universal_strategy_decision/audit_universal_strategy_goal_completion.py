#!/usr/bin/env python3
"""Summarize whether the active universal-strategy goal is complete.

This audit is intentionally lightweight: it consumes existing decision TSVs and
does not rerun profiling. Its purpose is to keep the current goal state
machine-readable after long sessions and restarts.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

OBJECTIVE = RESULTS / "objective_audit.tsv"
RELEASE_GATE = RESULTS / "universal_strategy_release_gate.tsv"
ABUNDANCE = RESULTS / "abundance_release_blocker_audit.tsv"
HOLDOUT = RESULTS / "holdout_gap_action_plan_audit.tsv"
CAMI3_RECOVERY = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
MANIFEST = RESULTS / "current_default_strategy_manifest.tsv"
NEXT_ROUTES = RESULTS / "universal_strategy_next_evidence_routes_audit.tsv"
REMAINING_ROUTES = RESULTS / "remaining_holdout_route_options_audit.tsv"

OUT_TSV = RESULTS / "universal_strategy_goal_completion_audit.tsv"
OUT_MD = EXP / "UNIVERSAL_STRATEGY_GOAL_COMPLETION_AUDIT.md"


FIELDS = [
    "check_id",
    "requirement",
    "status",
    "observed",
    "evidence",
    "decision",
    "next_action",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def keyed(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def optional_keyed(path: Path, key: str) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    return keyed(path, key)


def write_tsv(path: Path, rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def row(
    check_id: str,
    requirement: str,
    status: str,
    observed: str,
    evidence: str,
    decision: str,
    next_action: str,
) -> dict[str, str]:
    return {
        "check_id": check_id,
        "requirement": requirement,
        "status": status,
        "observed": observed,
        "evidence": evidence,
        "decision": decision,
        "next_action": next_action,
    }


def build_rows() -> list[dict[str, str]]:
    objective = keyed(OBJECTIVE, "requirement")
    gate = keyed(RELEASE_GATE, "gate_id")
    abundance = keyed(ABUNDANCE, "metric")
    holdout = keyed(HOLDOUT, "metric")
    cami3 = keyed(CAMI3_RECOVERY, "metric")
    manifest = keyed(MANIFEST, "item")
    routes = keyed(NEXT_ROUTES, "metric")
    remaining = optional_keyed(REMAINING_ROUTES, "metric")
    remaining_status = remaining.get("promotion_decision", {}).get("decision", "NA")
    remaining_top = remaining.get("top_remaining_route", {}).get("value", "NA")
    remaining_local = remaining.get("local_viable_uncompleted_routes", {}).get(
        "value", "NA"
    )
    cami3_scored = (
        cami3.get("promotion_decision", {}).get("value") == "post_recovery_extension_scored"
    )
    top_route = routes["top_next_evidence_route"]["value"]
    completed_routes = {
        item.strip()
        for item in routes["completed_evidence_routes"]["value"].split(",")
        if item.strip() and item.strip() != "none"
    }
    hmp_completed = "hmp_airskin_omitted_samples_2_8_12_26_27" in completed_routes
    cami3_completed = "cami3_source_readmap_samples3_5_completed" in completed_routes
    if top_route == "marine_truth_upgrade":
        route_blocker = routes.get("top_route_blocker_detail", {}).get("value", "NA")
        route_decision = routes.get("top_next_evidence_route", {}).get("decision", "NA")
        if route_decision == "truth_upgrade_profiles_needed":
            route_next_action = (
                "Use the marine truth-upgrade route next; the HMP omitted-sample "
                "profile-pair route and CAMI3 readmap-recovery route are already "
                "scored negative for promotion. Marine truth mapping now clears the "
                "threshold, so the current route state is selected-default profile "
                f"scoring needed. Detail: {route_blocker}."
            )
            release_gap_next_action = (
                "Run selected-default same-namespace profile scoring for the marine "
                "route. HMP omitted and CAMI3 samples3-5 are scored negative, and "
                "marine truth mapping has cleared the threshold under the strict "
                "setup source-name or local assembly-summary rule. "
                f"Remaining-route audit: top={remaining_top}; "
                f"local_viable_uncompleted={remaining_local}; decision={remaining_status}."
            )
            cami3_next_action = (
                "Do not repeat the completed CAMI3 samples3-5 recovery unless "
                "inputs or strategy change; continue with marine selected-default "
                "profile scoring."
            )
        else:
            route_next_action = (
                "Use the marine truth-upgrade route next; the HMP omitted-sample "
                "profile-pair route and CAMI3 readmap-recovery route are already "
                "scored negative for promotion. Current route state: "
                f"{route_decision}; blocker/detail: {route_blocker}."
            )
            release_gap_next_action = (
                "Recover or generate additional same-namespace GTDB truth/profile "
                "pairs. HMP omitted and CAMI3 samples3-5 are scored negative, so "
                "the current nearest route is marine truth upgrade. The next step "
                "is stronger source-specific mapping or another clean holdout "
                "rather than local threshold tuning. "
                f"Remaining-route audit: top={remaining_top}; "
                f"local_viable_uncompleted={remaining_local}; decision={remaining_status}."
            )
            cami3_next_action = (
                "Do not repeat the completed CAMI3 samples3-5 recovery unless "
                "inputs or strategy change; continue with the marine truth-upgrade "
                "route and its stricter source-mapping requirement."
            )
    elif top_route == "refined_allocator_candidate":
        route_blocker = routes.get("top_route_blocker_detail", {}).get("value", "NA")
        route_next_action = (
            "Marine, HMP omitted, and CAMI3 samples3-5 routes are scored negative "
            "for default promotion. The remaining route is the refined allocator "
            "candidate, which still needs an independent release-grade holdout "
            f"before promotion. Current blocker/detail: {route_blocker}."
        )
        release_gap_next_action = (
            "Do not claim completion from current local routes. The next useful "
            "work is either a strategy change that addresses the negative marine/HMP/CAMI3 "
            "evidence or a new clean same-namespace holdout, then rerun the release gate. "
            f"Remaining-route audit: top={remaining_top}; "
            f"local_viable_uncompleted={remaining_local}; decision={remaining_status}."
        )
        cami3_next_action = (
            "Do not repeat the completed CAMI3 samples3-5 recovery unless inputs "
            "or strategy change; continue with a new strategy or new holdout."
        )
    else:
        route_blocker = routes.get("top_route_blocker_detail", {}).get("value", "NA")
        route_next_action = (
            "Use the HMP omitted-sample profile-pair recovery path next. Current "
            f"route blocker/counts: {route_blocker}. The CAMI3 readmap-recovery "
            "route is already scored and negative for allocator promotion."
        )
        release_gap_next_action = (
            "Recover or generate additional same-namespace GTDB truth/profile "
            "pairs; the CAMI3 samples3-5 route is scored, so the current "
            "nearest route is the HMP omitted-sample profile-pair recovery path "
            "after scratch headroom is restored."
        )
        cami3_next_action = (
            "Do not repeat the completed CAMI3 samples3-5 recovery unless "
            "inputs or strategy change; continue with the HMP omitted-sample "
            "profile-pair recovery path next."
            if cami3_scored
            else "Restore only the missing CAMI3 samples3-5 readmap files, then rerun score_cami3_source_readmap_extension_after_recovery.py."
        )
    next_route_ok = (
        (cami3_completed and not hmp_completed and top_route == "hmp_airskin_omitted_samples_2_8_12_26_27")
        or (cami3_completed and hmp_completed and top_route == "marine_truth_upgrade")
        or (cami3_completed and hmp_completed and top_route == "refined_allocator_candidate")
    )

    rows = [
        row(
            "default_entrypoint",
            "User can run one current default without choosing a dataset mode.",
            "pass"
            if gate["default_profile_preset"]["status"] == "pass"
            and gate["default_launcher_strategy"]["status"] == "pass"
            else "fail",
            "preset="
            + gate["default_profile_preset"]["observed"]
            + ";strategy="
            + gate["default_launcher_strategy"]["observed"],
            "results/universal_strategy_release_gate.tsv;scripts/minco_profile_default.py",
            "selected_candidate_default",
            "Keep scripts/minco_profile_default.py --profile-preset candidate as the documented no-manual entrypoint.",
        ),
        row(
            "best_minco_f1_default",
            "F1-priority default is supported versus previous MinCO cached panels.",
            "pass" if gate["candidate_vs_previous_minco"]["status"] == "pass" else "fail",
            gate["candidate_vs_previous_minco"]["observed"],
            gate["candidate_vs_previous_minco"]["evidence"],
            gate["candidate_vs_previous_minco"]["decision"],
            "Do not retune the species boundary from a small panel; use more clean holdouts for release evidence.",
        ),
        row(
            "claim_boundary",
            "Current evidence does not overclaim broad external-baseline superiority.",
            "pass" if gate["sylph_claim_boundary"]["status"] == "pass" else "fail",
            gate["sylph_claim_boundary"]["observed"],
            gate["sylph_claim_boundary"]["evidence"],
            gate["sylph_claim_boundary"]["decision"],
            "Keep wording as best current MinCO default, not broad external-baseline beating.",
        ),
        row(
            "release_holdout_gap",
            "Clean multi-dataset release holdouts are complete enough for a final universal claim.",
            "expected_gap" if gate["release_grade_bundle"]["status"] == "expected_gap" else gate["release_grade_bundle"]["status"],
            gate["release_grade_bundle"]["observed"]
            + ";"
            + holdout["panel_counts_by_grade"]["value"],
            "results/universal_strategy_release_gate.tsv;results/holdout_gap_action_plan_audit.tsv;results/remaining_holdout_route_options_audit.tsv",
            gate["release_grade_bundle"]["decision"],
            release_gap_next_action,
        ),
        row(
            "abundance_gap",
            "Abundance L1/Pearson is strong enough to promote a broad abundance claim.",
            "expected_gap"
            if gate["abundance_release_claim"]["status"] == "expected_gap"
            else gate["abundance_release_claim"]["status"],
            abundance["abundance_release_decision"]["value"],
            "results/abundance_release_blocker_audit.tsv;ABUNDANCE_RELEASE_BLOCKER.md",
            abundance["abundance_release_decision"]["decision"],
            "Keep the refined allocator opt-in; validate on independent release-grade holdouts before promotion.",
        ),
        row(
            "ani_reporting_boundary",
            "ANI/reporting is documented as diagnostic and not overpromoted.",
            "pass" if gate["ani_reporting"]["status"] == "pass" else "fail",
            manifest["ani_reporting"]["value"],
            gate["ani_reporting"]["evidence"],
            gate["ani_reporting"]["decision"],
            "Continue reporting diagnostic ANI separately from the default call gate.",
        ),
        row(
            "cami3_extension_path",
            "CAMI3 post-recovery scorer is ready or already scored for the independent extension.",
            "pass" if gate["cami3_extension_runner_contract"]["status"] == "pass" else "fail",
            gate["cami3_extension_runner_contract"]["observed"],
            gate["cami3_extension_runner_contract"]["evidence"],
            gate["cami3_extension_runner_contract"]["decision"],
            cami3_next_action,
        ),
        row(
            "next_evidence_route",
            "The next route toward closing expected gaps is explicit and non-circular.",
            "pass"
            if next_route_ok
            else "fail",
            "top="
            + routes["top_next_evidence_route"]["value"]
            + ";completed="
            + routes["completed_evidence_routes"]["value"]
            + ";ready_without_external_restore="
            + routes["routes_ready_without_external_restore"]["value"],
            "results/universal_strategy_next_evidence_routes_audit.tsv;UNIVERSAL_STRATEGY_NEXT_EVIDENCE_ROUTES.md",
            routes["top_next_evidence_route"]["decision"],
            route_next_action,
        ),
        row(
            "objective_state",
            "The active research goal has a single final decision state.",
            "expected_gap",
            gate["summary"]["observed"]
            + ";objective="
            + objective["document stable universal strategy"]["status"],
            "results/universal_strategy_release_gate.tsv;results/objective_audit.tsv",
            "goal_not_complete_expected_gaps_remain",
            "Do local evidence work next; no final default-promotion claim until holdout and abundance gaps close.",
        ),
    ]
    return rows


def write_md(rows: list[dict[str, str]]) -> None:
    fail_count = sum(row["status"] == "fail" for row in rows)
    gap_count = sum(row["status"] == "expected_gap" for row in rows)
    pass_count = sum(row["status"] == "pass" for row in rows)
    if fail_count:
        completion = "blocked_by_failed_checks"
    elif gap_count:
        completion = "not_complete_expected_gaps_remain"
    else:
        completion = "complete"

    lines = [
        "# Universal Strategy Goal Completion Audit",
        "",
        "Generated by `audit_universal_strategy_goal_completion.py` from cached decision TSVs.",
        "",
        f"Completion decision: **{completion}**.",
        "",
        f"Summary: pass={pass_count}; expected_gap={gap_count}; fail={fail_count}.",
        "",
        "| Check | Status | Observed | Decision | Next Action |",
        "| --- | --- | --- | --- | --- |",
    ]
    for item in rows:
        lines.append(
            "| {check_id} | `{status}` | {observed} | {decision} | {next_action} |".format(
                **item
            )
        )
    lines.extend(
        [
            "",
            "The machine-readable audit is `results/universal_strategy_goal_completion_audit.tsv`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    write_tsv(OUT_TSV, rows)
    write_md(rows)
    print(f"wrote {OUT_TSV}")
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
