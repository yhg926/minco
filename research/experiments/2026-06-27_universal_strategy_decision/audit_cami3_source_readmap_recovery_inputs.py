#!/usr/bin/env python3
"""Audit exact CAMI3 samples3-5 source-readmap recovery inputs.

This is a local file-presence and recovery-manifest audit only. It does not
download or extract data. The goal is to make the remaining CAMI3 extension
blocker concrete: which archives or extracted files are missing, where they
should land, and which already-generated MinCO/Sylph profiling artifacts are
available for rescoring after truth recovery.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

DOWNLOAD_LIST = Path("/tmp/cami3_toy_human_gut_20260626/downloads/CAMI3_toy_dataset_download.list")
PRIMARY_ROOT = Path("/mnt/new3T/minco_cami3_toygut_extra_20260621")
LEGACY_ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620")
TMP_ROOT = Path("/tmp/cami3_toy_human_gut_20260626")
RUN = TMP_ROOT / "run"
CANDIDATE_ROOT = Path("/tmp/minco_candidate_preset_replay_20260629")
KNOWN_ARCHIVES = [
    LEGACY_ROOT / "sample_0_reads.tar.gz",
    PRIMARY_ROOT / "sample_1_reads.tar.gz",
    PRIMARY_ROOT / "sample_2_reads.tar.gz",
]

SAMPLES = [3, 4, 5]

OUT_DETAIL = RESULTS / "cami3_source_readmap_recovery_inputs.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_recovery_inputs_audit.tsv"
OUT_COMMANDS = RESULTS / "cami3_source_readmap_readmap_only_recovery_commands.sh"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_RECOVERY_INPUTS.md"


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


def download_urls() -> dict[int, str]:
    out: dict[int, str] = {}
    if not DOWNLOAD_LIST.is_file():
        return out
    for line in DOWNLOAD_LIST.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        for sample in SAMPLES:
            token = f"sample_{sample}_reads.tar.gz"
            if token in line:
                out[sample] = line
    return out


def archive_candidates(sample: int) -> list[Path]:
    return [
        PRIMARY_ROOT / f"sample_{sample}_reads.tar.gz",
        LEGACY_ROOT / f"sample_{sample}_reads.tar.gz",
        TMP_ROOT / "downloads" / f"sample_{sample}_reads.tar.gz",
    ]


def extracted_dir_candidates(sample: int) -> list[Path]:
    return [
        PRIMARY_ROOT / f"sample_{sample}_reads",
        LEGACY_ROOT / f"sample_{sample}_reads",
        TMP_ROOT / f"sample_{sample}_reads",
        TMP_ROOT / "downloads" / f"sample_{sample}_reads",
    ]


def mapping_candidates(sample: int) -> list[Path]:
    paths = []
    for root in extracted_dir_candidates(sample):
        paths.append(root / "reads_mapping.tsv.gz")
    paths.extend(
        [
            LEGACY_ROOT / f"sample_{sample}_reads_mapping.tsv.gz",
            RUN / f"sample{sample}_reads_mapping.tsv.gz",
        ]
    )
    return paths


def reads_candidates(sample: int) -> list[Path]:
    return [root / "anonymous_reads.fq.gz" for root in extracted_dir_candidates(sample)]


def raw_unique(sample: int) -> Path:
    return RUN / f"sample{sample}_autoexact/minco.best_diff_unique.unfiltered.tsv"


def raw_split(sample: int) -> Path:
    return RUN / f"sample{sample}_autoexact/minco.best_diff_split.unfiltered.tsv"


def candidate_profile(sample: int) -> Path:
    return CANDIDATE_ROOT / (
        f"cami3_toy_human_gut_gtdb_source_readmap.sample{sample}.candidate_preset.tsv"
    )


def sylph_profile(sample: int) -> Path:
    return RUN / f"sylph_sample{sample}/profile.tsv"


def mean_existing_archive_size_bytes() -> int:
    sizes = [path.stat().st_size for path in KNOWN_ARCHIVES if path.is_file()]
    return int(sum(sizes) / len(sizes)) if sizes else 0


def readmap_only_recovery_commands(urls: dict[int, str]) -> list[str]:
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Stream-extract only reads_mapping.tsv.gz for CAMI3 samples3-5.",
        "# This avoids storing the full sample read archives or extracting read FASTQ.",
    ]
    for sample in SAMPLES:
        url = urls.get(sample, "")
        target_dir = PRIMARY_ROOT / f"sample_{sample}_reads"
        member = f"sample_{sample}_reads/reads_mapping.tsv.gz"
        lines.extend(
            [
                "",
                f"mkdir -p {target_dir}",
                (
                    f"curl -fL {url} | "
                    f"tar -xz -C {PRIMARY_ROOT} {member}"
                ),
            ]
        )
    lines.append("")
    return lines


def build_detail() -> list[dict[str, object]]:
    urls = download_urls()
    rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        archive_ok, archive_checked = first_existing(archive_candidates(sample))
        extracted_ok, extracted_checked = first_existing(extracted_dir_candidates(sample))
        mapping_ok, mapping_checked = first_existing(mapping_candidates(sample))
        reads_ok, reads_checked = first_existing(reads_candidates(sample))
        blockers = []
        if not mapping_ok:
            blockers.append("missing_reads_mapping")
            if not archive_ok and not urls.get(sample):
                blockers.append("missing_archive_source")
        rows.append(
            {
                "sample": sample,
                "archive_url": urls.get(sample, ""),
                "archive_present": archive_ok,
                "archive_path_or_checked": archive_checked,
                "extracted_dir_present": extracted_ok,
                "extracted_dir_path_or_checked": extracted_checked,
                "reads_mapping_present": mapping_ok,
                "reads_mapping_path_or_checked": mapping_checked,
                "anonymous_reads_present": reads_ok,
                "anonymous_reads_path_or_checked": reads_checked,
                "anonymous_reads_required_for_truth_rescore": False,
                "selected_default_profile_present": candidate_profile(sample).is_file(),
                "selected_default_profile": candidate_profile(sample),
                "raw_unique_present": raw_unique(sample).is_file(),
                "raw_unique": raw_unique(sample),
                "raw_split_present": raw_split(sample).is_file(),
                "raw_split": raw_split(sample),
                "sylph_profile_present": sylph_profile(sample).is_file(),
                "sylph_profile": sylph_profile(sample),
                "release_truth_recovery_ready": mapping_ok,
                "blocking_reason": ",".join(blockers),
            }
        )
    return rows


def build_audit(detail: list[dict[str, object]]) -> list[dict[str, object]]:
    archive_n = sum(bool(row["archive_present"]) for row in detail)
    mapping_n = sum(bool(row["reads_mapping_present"]) for row in detail)
    reads_n = sum(bool(row["anonymous_reads_present"]) for row in detail)
    profiles_n = sum(bool(row["selected_default_profile_present"]) for row in detail)
    raw_n = sum(bool(row["raw_unique_present"]) and bool(row["raw_split_present"]) for row in detail)
    sylph_n = sum(bool(row["sylph_profile_present"]) for row in detail)
    urls_n = sum(bool(row["archive_url"]) for row in detail)
    mean_size = mean_existing_archive_size_bytes()
    total_estimate = mean_size * len(SAMPLES)
    blockers = sorted(
        {
            blocker
            for row in detail
            for blocker in str(row.get("blocking_reason", "")).split(",")
            if blocker
        }
    )
    ready_for_truth_rescore = (
        mapping_n == len(detail)
        and profiles_n == len(detail)
        and raw_n == len(detail)
        and sylph_n == len(detail)
    )
    return [
        {
            "metric": "remote_manifest_urls",
            "value": f"{urls_n}/{len(detail)}",
            "evidence": str(DOWNLOAD_LIST),
            "decision": "download_urls_available" if urls_n == len(detail) else "download_urls_missing",
        },
        {
            "metric": "local_archives",
            "value": f"{archive_n}/{len(detail)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "archives_missing_locally" if archive_n < len(detail) else "archives_present",
        },
        {
            "metric": "local_readmaps",
            "value": f"{mapping_n}/{len(detail)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "readmaps_missing_locally" if mapping_n < len(detail) else "readmaps_present",
        },
        {
            "metric": "local_anonymous_reads",
            "value": f"{reads_n}/{len(detail)}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "not_required_for_truth_rescore"
            if reads_n < len(detail)
            else "reads_present_but_not_required_for_truth_rescore",
        },
        {
            "metric": "rescoring_artifacts",
            "value": f"candidate_profiles={profiles_n}/3;raw_tables={raw_n}/3;sylph_profiles={sylph_n}/3",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "rescoring_artifacts_ready"
            if profiles_n == raw_n == sylph_n == len(detail)
            else "rescoring_artifacts_incomplete",
        },
        {
            "metric": "compressed_archive_size_estimate",
            "value": f"mean_existing_archive_bytes={mean_size};samples3_5_estimated_bytes={total_estimate}",
            "evidence": ";".join(str(path) for path in KNOWN_ARCHIVES),
            "decision": "size_estimate_from_samples0_2",
        },
        {
            "metric": "readmap_only_recovery_plan",
            "value": f"commands={len(detail)}/{len(detail)}",
            "evidence": str(OUT_COMMANDS.relative_to(EXP)),
            "decision": "stream_extract_readmaps_without_storing_archives",
        },
        {
            "metric": "recovery_blockers",
            "value": ",".join(blockers),
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "no_recovery_blockers"
            if not blockers
            else "requires_readmap_restore_or_stream_extract",
        },
        {
            "metric": "promotion_decision",
            "value": "source_readmap_truth_recovery_inputs_ready"
            if ready_for_truth_rescore
            else "source_readmap_truth_recovery_inputs_missing_locally",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "post_recovery_scoring_ready"
            if ready_for_truth_rescore
            else "keep_refined_allocator_opt_in",
        },
    ]


def write_markdown(detail: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# CAMI3 Source-Readmap Recovery Inputs",
        "",
        "Date: 2026-06-29",
        "",
        "This audit records the exact local and remote inputs needed to recover",
        "per-read source mapping for CAMI3 samples3-5. It does not download or",
        "extract data.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Detail",
            "",
            "| Sample | URL present | Archive | Readmap | Reads optional | Candidate profile | Raw tables | Baseline profile |",
            "|---:|---|---|---|---|---|---|---|",
        ]
    )
    for row in detail:
        raw_ok = bool(row["raw_unique_present"]) and bool(row["raw_split_present"])
        lines.append(
            "| {sample} | {url} | {archive} | {mapping} | {reads} | {candidate} | {raw} | {sylph} |".format(
                sample=row["sample"],
                url=bool(row["archive_url"]),
                archive=row["archive_present"],
                mapping=row["reads_mapping_present"],
                reads=not bool(row["anonymous_reads_required_for_truth_rescore"]),
                candidate=row["selected_default_profile_present"],
                raw=raw_ok,
                sylph=row["sylph_profile_present"],
            )
        )
    ready_for_truth_rescore = (
        audit_by_metric["promotion_decision"]["decision"] == "post_recovery_scoring_ready"
    )
    lines.extend(["", "## Decision", ""])
    if ready_for_truth_rescore:
        lines.extend(
            [
                "- Samples3-5 now have extracted `reads_mapping.tsv.gz` files,",
                "  selected-default profiles, raw tables, and baseline profiles ready",
                "  for rescoring.",
                "- Anonymous reads are not needed for this truth-rescore path because",
                "  profiling artifacts already exist.",
            ]
        )
    else:
        lines.extend(
            [
                "- Samples3-5 have selected-default profiles, raw tables, and baseline",
                "  profiles ready for rescoring.",
                "- The missing release-truth input is the extracted `reads_mapping.tsv.gz`",
                "  for each sample. Anonymous reads are not needed for this truth-rescore",
                "  path because profiling artifacts already exist.",
            ]
        )
    lines.extend(
        [
            "- The generated command file streams each archive URL through `tar` and",
            "  extracts only `reads_mapping.tsv.gz`, avoiding storage of the full",
            "  sample read archives.",
            "- The source-profile fallback has already been tested and is not a release",
            "  substitute.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_DETAIL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            f"- `{OUT_COMMANDS.relative_to(EXP)}`",
            "",
            f"Promotion decision is `{audit_by_metric['promotion_decision']['decision']}`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    detail = build_detail()
    audit = build_audit(detail)
    OUT_COMMANDS.write_text("\n".join(readmap_only_recovery_commands(download_urls())), encoding="utf-8")
    detail_fields = [
        "sample",
        "archive_url",
        "archive_present",
        "archive_path_or_checked",
        "extracted_dir_present",
        "extracted_dir_path_or_checked",
        "reads_mapping_present",
        "reads_mapping_path_or_checked",
        "anonymous_reads_present",
        "anonymous_reads_path_or_checked",
        "anonymous_reads_required_for_truth_rescore",
        "selected_default_profile_present",
        "selected_default_profile",
        "raw_unique_present",
        "raw_unique",
        "raw_split_present",
        "raw_split",
        "sylph_profile_present",
        "sylph_profile",
        "release_truth_recovery_ready",
        "blocking_reason",
    ]
    write_tsv(OUT_DETAIL, detail, detail_fields)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(detail, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
