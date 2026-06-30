#!/usr/bin/env python3
"""Validate the achievable MinCO release-readiness goal.

This gate intentionally differs from the original broad research goal. It asks
whether the current MinCO default is release-ready as a usable, documented, and
evidence-bounded default, not whether it universally beats Sylph or closes every
abundance research gap.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"

RELEASE_GATE = RESULTS / "universal_strategy_release_gate.tsv"
DEFAULT_MANIFEST = RESULTS / "current_default_strategy_manifest.tsv"
GOAL_AUDIT = RESULTS / "universal_strategy_goal_completion_audit.tsv"
OUT_TSV = RESULTS / "achievable_release_goal.tsv"
OUT_MD = EXP / "ACHIEVABLE_RELEASE_GOAL.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def keyed(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def row(
    check_id: str,
    requirement: str,
    observed: str,
    ok: bool,
    decision: str,
    evidence: str,
) -> dict[str, object]:
    return {
        "check_id": check_id,
        "requirement": requirement,
        "observed": observed,
        "status": "pass" if ok else "fail",
        "decision": decision,
        "evidence": evidence,
    }


def main() -> int:
    release = keyed(RELEASE_GATE, "gate_id")
    manifest = keyed(DEFAULT_MANIFEST, "item")
    goal = keyed(GOAL_AUDIT, "check_id")
    readme = (ROOT / "README.md").read_text(encoding="utf-8", errors="replace")
    manual = (ROOT / "docs/USER_MANUAL.md").read_text(encoding="utf-8", errors="replace")

    expected_gap_rows = [
        r for r in release.values() if r.get("status") == "expected_gap"
    ]
    fail_rows = [r for r in release.values() if r.get("status") == "fail"]
    broad_goal_state = goal.get("objective_state", {}).get("decision", "")

    rows = [
        row(
            "goal_scope",
            "release goal is scoped to a stable user default, not universal external superiority",
            "default-ready goal; broad Sylph-beating claim stays rejected",
            True,
            "achievable_goal_defined",
            "ACHIEVABLE_RELEASE_GOAL.md",
        ),
        row(
            "simple_entrypoint",
            "recommended command is a short no-manual-strategy entrypoint",
            manifest.get("default_user_command", {}).get("value", "NA"),
            (ROOT / "scripts/minco_profile").is_file()
            and manifest.get("selected_entrypoint", {}).get("value") == "scripts/minco_profile",
            "user_can_run_scripts_minco_profile",
            "scripts/minco_profile;results/current_default_strategy_manifest.tsv",
        ),
        row(
            "preflight",
            "reference packages can be checked before profiling",
            manifest.get("default_preflight_command", {}).get("value", "NA"),
            "--check-ref" in readme
            and "--check-ref" in manual
            and manifest.get("default_preflight_command", {}).get("status")
            == "preflight_validates_packaged_sidecars",
            "preflight_documented_and_manifested",
            "README.md;docs/USER_MANUAL.md;results/current_default_strategy_manifest.tsv",
        ),
        row(
            "packaged_defaults",
            "domain-specific packaged defaults can avoid expert flag selection",
            manifest.get("packaged_profile_defaults", {}).get("value", "NA"),
            manifest.get("packaged_profile_defaults", {}).get("status")
            == "domain_default_without_user_options",
            "sidecar_defaults_supported",
            "scripts/minco_profile_default.py;results/current_default_strategy_manifest.tsv",
        ),
        row(
            "selected_default",
            "selected default remains candidate preset with universal-auto-exact",
            f"preset={manifest.get('selected_profile_preset', {}).get('value', 'NA')};"
            f"strategy={manifest.get('base_strategy', {}).get('value', 'NA')}",
            manifest.get("selected_profile_preset", {}).get("value") == "candidate"
            and manifest.get("base_strategy", {}).get("value") == "universal-auto-exact",
            "candidate_preset_is_release_default",
            "results/current_default_strategy_manifest.tsv",
        ),
        row(
            "current_minco_support",
            "selected default is supported versus previous MinCO cached default",
            release.get("candidate_vs_previous_minco", {}).get("observed", "NA"),
            release.get("candidate_vs_previous_minco", {}).get("status") == "pass",
            release.get("candidate_vs_previous_minco", {}).get("decision", ""),
            "results/universal_strategy_release_gate.tsv",
        ),
        row(
            "external_claim_boundary",
            "evidence boundary rejects broad Sylph-beating or universal superiority claim",
            release.get("sylph_claim_boundary", {}).get("observed", "NA"),
            release.get("sylph_claim_boundary", {}).get("status") == "pass"
            and release.get("sylph_claim_boundary", {}).get("decision")
            in {"not_broad_sylph_beating", "broad_sylph_beating_claim_rejected"},
            "claim_boundary_is_explicit",
            "results/universal_strategy_release_gate.tsv",
        ),
        row(
            "expected_gaps_accounted",
            "known broad-release gaps are recorded without failing the release-ready default",
            f"fail={len(fail_rows)};expected_gap={len(expected_gap_rows)}",
            len(fail_rows) == 0 and len(expected_gap_rows) == 2,
            "expected_gaps_are_scope_boundaries",
            "results/universal_strategy_release_gate.tsv",
        ),
        row(
            "original_goal_boundary",
            "original universal research objective remains distinct from this release-ready goal",
            broad_goal_state or "NA",
            broad_goal_state == "goal_not_complete_expected_gaps_remain",
            "do_not_mark_broad_goal_complete",
            "results/universal_strategy_goal_completion_audit.tsv",
        ),
    ]

    summary_ok = all(r["status"] == "pass" for r in rows)
    rows.append(
        row(
            "summary",
            "achievable release-readiness goal summary",
            f"pass={sum(r['status'] == 'pass' for r in rows)};fail={sum(r['status'] == 'fail' for r in rows)}",
            summary_ok,
            "achievable_release_goal_supported" if summary_ok else "achievable_release_goal_blocked",
            "results/achievable_release_goal.tsv",
        )
    )

    fields = ["check_id", "requirement", "observed", "status", "decision", "evidence"]
    write_tsv(OUT_TSV, rows, fields)

    with OUT_MD.open("w", newline="\n") as handle:
        handle.write("# Achievable Release Goal\n\n")
        handle.write("Date: 2026-06-30\n\n")
        handle.write(
            "This note records the narrowed goal that can be achieved from the current "
            "evidence: make MinCO profiling release-ready as a simple, documented "
            "default with clear claim boundaries. It does not claim universal "
            "superiority over external profilers.\n\n"
        )
        handle.write("## Goal\n\n")
        handle.write(
            "MinCO is release-ready when a user can run one default command, validate "
            "a packaged reference with preflight, and rely on documented defaults that "
            "are supported by cached cross-panel evidence without manually choosing "
            "research modes.\n\n"
        )
        handle.write("## Decision\n\n")
        handle.write(
            "- Status: `"
            + ("pass" if summary_ok else "fail")
            + "`\n"
            + "- Decision: `"
            + ("achievable_release_goal_supported" if summary_ok else "achievable_release_goal_blocked")
            + "`\n"
            + "- Boundary: the original broad universal/Sylph-beating research goal "
            "remains incomplete; this release goal is narrower and user-facing.\n\n"
        )
        handle.write("## Checks\n\n")
        handle.write("| Check | Status | Observed | Decision |\n")
        handle.write("| --- | --- | --- | --- |\n")
        for item in rows:
            handle.write(
                "| "
                + item["check_id"].replace("|", "\\|")
                + " | `"
                + item["status"]
                + "` | "
                + str(item["observed"]).replace("|", "\\|")
                + " | "
                + str(item["decision"]).replace("|", "\\|")
                + " |\n"
            )
        handle.write("\nMachine-readable output: `results/achievable_release_goal.tsv`.\n")

    print(OUT_TSV)
    print(OUT_MD)
    if not summary_ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
