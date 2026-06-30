#!/usr/bin/env python3
"""Audit MinCO read/context ambiguity outputs for abundance research.

This script is intentionally lightweight: it checks representative existing
outputs and source-code surfaces, then writes a TSV stating which outputs are
sufficient for true context-level ambiguity assignment or EM.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
OUT = EXP / "results/read_context_ambiguity_output_audit.tsv"

TRACK = Path(
    "/tmp/cami2_hmp_unseen_transfer_20260626/run/"
    "minco_sample11_split_track_s1000.tsv"
)
TRACK_SUMMARY = Path(
    "/tmp/cami2_hmp_unseen_transfer_20260626/run/"
    "minco_sample11_split_track_s1000.summary.tsv"
)
PROFILE = Path(
    "/tmp/cami2_hmp_unseen_transfer_20260626/run/"
    "minco_sample11_split_track_s1000.profile.tsv"
)


def header(path: Path) -> list[str]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        first = fh.readline().rstrip("\n")
    return first.split("\t") if first else []


def line_count(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("rb") as fh:
        return sum(1 for _ in fh)


def contains(path: Path, needle: str) -> bool:
    if not path.exists():
        return False
    return needle in path.read_text(errors="replace")


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    track_header = header(TRACK)
    profile_header = header(PROFILE)

    ani_wrapper = ROOT / "minco_core/src/command_ani_wrapper.c"
    profile_wrapper = ROOT / "minco_core/src/command_profile_wrapper.c"
    ani_core = ROOT / "minco_core/src/command_ani.c"
    calibrated = ROOT / "scripts/minco_profile_calibrated.py"

    edge_cli = (
        contains(ani_wrapper, "readwise-edge-out")
        and contains(profile_wrapper, "edge-out")
        and contains(ani_core, "candidate_refs")
    )
    calibrated_uses_edges = contains(calibrated, "readwise-edge-out") or contains(
        calibrated, "edge-out"
    )

    rows = [
        {
            "surface": "minco profile --track",
            "artifact": str(TRACK),
            "available": yesno(TRACK.exists()),
            "rows_observed": "" if line_count(TRACK) is None else str(line_count(TRACK) - 1),
            "key_columns": ",".join(
                c
                for c in [
                    "read_id",
                    "selected_ctx",
                    "target_ref_count",
                    "target_refs",
                    "ctx_offsets",
                    "assignment_status",
                ]
                if c in track_header
            ),
            "has_full_candidate_refs_per_context": "no",
            "has_selected_refs_per_context": "no",
            "has_read_level_offsets": yesno("ctx_offsets" in track_header),
            "sufficient_for_context_em": "no",
            "default_calibrated_uses_it": "no",
            "decision": "useful for read/taxa provenance, not enough for context-level EM",
        },
        {
            "surface": "minco profile --track-summary",
            "artifact": str(TRACK_SUMMARY),
            "available": yesno(TRACK_SUMMARY.exists()),
            "rows_observed": "" if line_count(TRACK_SUMMARY) is None else str(line_count(TRACK_SUMMARY)),
            "key_columns": "metric,value",
            "has_full_candidate_refs_per_context": "no",
            "has_selected_refs_per_context": "no",
            "has_read_level_offsets": "no",
            "sufficient_for_context_em": "no",
            "default_calibrated_uses_it": "no",
            "decision": "summary statistics only",
        },
        {
            "surface": "readwise profile table",
            "artifact": str(PROFILE),
            "available": yesno(PROFILE.exists()),
            "rows_observed": "" if line_count(PROFILE) is None else str(line_count(PROFILE) - 1),
            "key_columns": ",".join(
                c
                for c in [
                    "Ref",
                    "XnY_ctx",
                    "Ref_breadth",
                    "Ref_mean_depth",
                    "Ref_zip_af",
                    "Reliable_Ref_zip_af",
                ]
                if c in profile_header
            ),
            "has_full_candidate_refs_per_context": "no",
            "has_selected_refs_per_context": "no",
            "has_read_level_offsets": "no",
            "sufficient_for_context_em": "no",
            "default_calibrated_uses_it": "yes",
            "decision": "aggregated per-reference features; loses ambiguity groups",
        },
        {
            "surface": "context ambiguity edge sidecar",
            "artifact": "minco ani --readwise-edge-out / minco profile --edge-out",
            "available": yesno(edge_cli),
            "rows_observed": "runtime-dependent",
            "key_columns": (
                "edge_id,read_id,unit_id,qctx,gid,diff,best_diff,"
                "candidate_refs,selected_refs,selected_by_mode,cov_inc"
            ),
            "has_full_candidate_refs_per_context": "yes",
            "has_selected_refs_per_context": "yes",
            "has_read_level_offsets": "no",
            "sufficient_for_context_em": "yes",
            "default_calibrated_uses_it": yesno(calibrated_uses_edges),
            "decision": (
                "sufficient diagnostic substrate for EM; not yet consumed by "
                "the calibrated default"
            ),
        },
        {
            "surface": "scripts/minco_profile_calibrated.py",
            "artifact": str(calibrated),
            "available": yesno(calibrated.exists()),
            "rows_observed": "NA",
            "key_columns": "unique-table,split-table,RF/HGB features",
            "has_full_candidate_refs_per_context": "no",
            "has_selected_refs_per_context": "no",
            "has_read_level_offsets": "no",
            "sufficient_for_context_em": "no",
            "default_calibrated_uses_it": "yes",
            "decision": (
                "current default still uses repeated aggregate passes; next "
                "speed/abundance work needs a one-scan feature/edge producer"
            ),
        },
    ]

    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
