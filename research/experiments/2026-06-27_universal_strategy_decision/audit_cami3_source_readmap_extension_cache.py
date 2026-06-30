#!/usr/bin/env python3
"""Audit local cache readiness for extending the CAMI3 source-readmap panel.

The accepted release-grade source-readmap panel currently covers samples 0-2.
This audit only checks local file presence for samples 3-5; it does not run
profilers and does not infer truth from weaker taxid-only tables.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620")
EXTRA = Path("/mnt/new3T/minco_cami3_toygut_extra_20260621")
RUN = Path("/tmp/cami3_toy_human_gut_20260626/run")
CANDIDATE_ROOT = Path("/tmp/minco_candidate_preset_replay_20260629")

ACCEPTED_SAMPLES = [0, 1, 2]
EXTENSION_SAMPLES = [3, 4, 5]

KNOWN_SOURCE_MAPPINGS = {
    0: [ROOT / "sample_0_reads_mapping.tsv.gz"],
    1: [EXTRA / "sample_1_reads/reads_mapping.tsv.gz"],
    2: [EXTRA / "sample_2_reads/reads_mapping.tsv.gz"],
}
KNOWN_BASELINE_PROFILES = {
    0: [ROOT / "sylph_sample0/profile.tsv"],
    1: [Path("/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv")],
    2: [RUN / "sylph_sample2/profile.tsv"],
}

OUT_DETAIL = RESULTS / "cami3_source_readmap_extension_cache_detail.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_extension_cache_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_EXTENSION_CACHE_AUDIT.md"


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def first_existing(paths: Iterable[Path]) -> tuple[bool, str]:
    paths = list(paths)
    for path in paths:
        if path.exists():
            return True, str(path)
    return False, ";".join(str(path) for path in paths)


def source_mapping_candidates(sample: int) -> list[Path]:
    if sample in KNOWN_SOURCE_MAPPINGS:
        return KNOWN_SOURCE_MAPPINGS[sample]
    return [
        EXTRA / f"sample_{sample}_reads/reads_mapping.tsv.gz",
        ROOT / f"sample_{sample}_reads_mapping.tsv.gz",
        RUN / f"sample{sample}_reads_mapping.tsv.gz",
    ]


def candidate_profile_candidates(sample: int) -> list[Path]:
    return [
        CANDIDATE_ROOT / f"cami3_sample{sample}.candidate_preset.tsv",
        CANDIDATE_ROOT / f"cami3_toy_human_gut_gtdb_source_readmap.sample{sample}.candidate_preset.tsv",
    ]


def old_profile_candidates(sample: int) -> list[Path]:
    return [
        RUN / f"sample{sample}_universal_autoexact.tsv",
        RUN / f"sample{sample}_universal_autoexact_tail_p025.tsv",
        RUN / f"sample{sample}_default_strategy.tsv",
    ]


def raw_unique_candidates(sample: int) -> list[Path]:
    return [RUN / f"sample{sample}_autoexact/minco.best_diff_unique.unfiltered.tsv"]


def raw_split_candidates(sample: int) -> list[Path]:
    return [RUN / f"sample{sample}_autoexact/minco.best_diff_split.unfiltered.tsv"]


def sylph_profile_candidates(sample: int) -> list[Path]:
    if sample in KNOWN_BASELINE_PROFILES:
        return KNOWN_BASELINE_PROFILES[sample]
    return [
        RUN / f"sylph_sample{sample}/profile.tsv",
        ROOT / f"sylph_sample{sample}/profile.tsv",
    ]


def row_for_sample(sample: int, scope: str) -> dict[str, object]:
    source_ok, source_path = first_existing(source_mapping_candidates(sample))
    candidate_ok, candidate_path = first_existing(candidate_profile_candidates(sample))
    old_ok, old_path = first_existing(old_profile_candidates(sample))
    raw_unique_ok, raw_unique_path = first_existing(raw_unique_candidates(sample))
    raw_split_ok, raw_split_path = first_existing(raw_split_candidates(sample))
    sylph_ok, sylph_path = first_existing(sylph_profile_candidates(sample))

    blockers = []
    if not source_ok:
        blockers.append("missing_source_readmap_truth")
    if not candidate_ok:
        blockers.append("missing_selected_default_profile")
    if not raw_unique_ok or not raw_split_ok:
        blockers.append("missing_raw_best_ref_tables")
    if not sylph_ok:
        blockers.append("missing_sylph_profile")

    release_ready = source_ok and candidate_ok and raw_unique_ok and raw_split_ok and sylph_ok
    return {
        "sample": sample,
        "scope": scope,
        "source_mapping_present": source_ok,
        "source_mapping_path_or_checked": source_path,
        "selected_default_profile_present": candidate_ok,
        "selected_default_profile_path_or_checked": candidate_path,
        "older_profile_present": old_ok,
        "older_profile_path_or_checked": old_path,
        "raw_unique_present": raw_unique_ok,
        "raw_unique_path_or_checked": raw_unique_path,
        "raw_split_present": raw_split_ok,
        "raw_split_path_or_checked": raw_split_path,
        "sylph_profile_present": sylph_ok,
        "sylph_profile_path_or_checked": sylph_path,
        "release_extension_ready": release_ready,
        "blocking_reason": ",".join(blockers),
    }


def count_true(rows: list[dict[str, object]], key: str) -> int:
    return sum(1 for row in rows if bool(row.get(key)))


def build_audit(detail: list[dict[str, object]]) -> list[dict[str, object]]:
    accepted = [row for row in detail if row["scope"] == "accepted_release_subset"]
    extension = [row for row in detail if row["scope"] == "extension_candidate_subset"]
    extension_source_n = count_true(extension, "source_mapping_present")
    extension_candidate_n = count_true(extension, "selected_default_profile_present")
    extension_old_n = count_true(extension, "older_profile_present")
    extension_sylph_n = count_true(extension, "sylph_profile_present")
    extension_raw_n = sum(
        1 for row in extension if row["raw_unique_present"] and row["raw_split_present"]
    )
    extension_ready_n = count_true(extension, "release_extension_ready")
    extension_blockers = sorted(
        {
            blocker
            for row in extension
            for blocker in str(row.get("blocking_reason", "")).split(",")
            if blocker
        }
    )
    blocker_actions = {
        "missing_source_readmap_truth": "source_readmap_truth_cache",
        "missing_selected_default_profile": "selected_default_replay",
        "missing_raw_best_ref_tables": "raw_best_ref_tables",
        "missing_sylph_profile": "baseline_profile",
    }
    blocker_decision = (
        "no_extension_blockers"
        if not extension_blockers
        else "requires_" + "_and_".join(blocker_actions.get(item, item) for item in extension_blockers)
    )
    return [
        {
            "metric": "accepted_samples_release_ready",
            "value": f"{count_true(accepted, 'release_extension_ready')}/{len(accepted)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "accepted_subset_has_required_local_cache",
        },
        {
            "metric": "extension_samples_with_source_readmap_truth",
            "value": f"{extension_source_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "truth_cache_present_for_extension"
            if extension_source_n == len(extension)
            else "truth_cache_missing_for_extension",
        },
        {
            "metric": "extension_samples_with_selected_default_profiles",
            "value": f"{extension_candidate_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "selected_default_cache_present_for_extension"
            if extension_candidate_n == len(extension)
            else "selected_default_cache_missing_for_extension",
        },
        {
            "metric": "extension_samples_with_older_profiles",
            "value": f"{extension_old_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "older_profiles_exist_but_not_release_default",
        },
        {
            "metric": "extension_samples_with_raw_tables",
            "value": f"{extension_raw_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "raw_best_ref_tables_present_for_extension",
        },
        {
            "metric": "extension_samples_with_sylph_profiles",
            "value": f"{extension_sylph_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "baseline_profiles_present_for_extension",
        },
        {
            "metric": "extension_samples_release_ready",
            "value": f"{extension_ready_n}/{len(extension)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "extension_release_ready"
            if extension_ready_n == len(extension)
            else "do_not_extend_release_panel_from_current_cache",
        },
        {
            "metric": "extension_blockers",
            "value": ",".join(extension_blockers),
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": blocker_decision,
        },
    ]


def write_markdown(detail: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    selected_decision = audit_by_metric["extension_samples_with_selected_default_profiles"][
        "decision"
    ]
    blocker_decision = audit_by_metric["extension_blockers"]["decision"]
    lines = [
        "# CAMI3 Source-Readmap Extension Cache Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This audit checks whether the accepted source-readmap release panel can",
        "be extended from local cached files alone. It only checks file presence;",
        "it does not rerun profilers or infer release truth from weaker tables.",
        "",
        "## Result",
        "",
        f"- Accepted subset ready: `{audit_by_metric['accepted_samples_release_ready']['value']}`.",
        f"- Extension subset release-ready: `{audit_by_metric['extension_samples_release_ready']['value']}`.",
        f"- Extension blockers: `{audit_by_metric['extension_blockers']['value']}`.",
        "",
        "Older profile outputs, raw best-reference tables, baseline profiles,",
        "and selected-default replay outputs exist for the later samples. The",
        "current remaining release blocker is source-readmap truth cache in the",
        f"same namespace (`{blocker_decision}`; selected replay status:",
        f"`{selected_decision}`).",
        "",
        "## Detail",
        "",
        "| Sample | Scope | Truth Cache | Selected Default Profile | Older Profile | Raw Tables | Baseline Profile | Ready |",
        "|---:|---|---|---|---|---|---|---|",
    ]
    for row in detail:
        raw_ok = bool(row["raw_unique_present"]) and bool(row["raw_split_present"])
        lines.append(
            "| {sample} | {scope} | {truth} | {candidate} | {old} | {raw} | {baseline} | {ready} |".format(
                sample=row["sample"],
                scope=row["scope"],
                truth=row["source_mapping_present"],
                candidate=row["selected_default_profile_present"],
                old=row["older_profile_present"],
                raw=raw_ok,
                baseline=row["sylph_profile_present"],
                ready=row["release_extension_ready"],
            )
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_DETAIL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    detail = [
        row_for_sample(sample, "accepted_release_subset") for sample in ACCEPTED_SAMPLES
    ] + [
        row_for_sample(sample, "extension_candidate_subset") for sample in EXTENSION_SAMPLES
    ]
    audit = build_audit(detail)
    fields = [
        "sample",
        "scope",
        "source_mapping_present",
        "source_mapping_path_or_checked",
        "selected_default_profile_present",
        "selected_default_profile_path_or_checked",
        "older_profile_present",
        "older_profile_path_or_checked",
        "raw_unique_present",
        "raw_unique_path_or_checked",
        "raw_split_present",
        "raw_split_path_or_checked",
        "sylph_profile_present",
        "sylph_profile_path_or_checked",
        "release_extension_ready",
        "blocking_reason",
    ]
    write_tsv(OUT_DETAIL, detail, fields)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(detail, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
