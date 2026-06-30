#!/usr/bin/env python3
"""Audit whether local marine files contain source-ID crosswalk metadata.

The marine readmaps identify reads by CAMI genome IDs and source sequence IDs.
This audit checks local files for the missing source genome/contig-to-assembly
mapping needed to turn those readmaps into release-grade GTDB truth. The setup
archive metadata, if extracted locally, is treated as present but not sufficient
by itself; the separate setup truth-upgrade audit decides whether it clears the
release threshold.
"""

from __future__ import annotations

import csv
from pathlib import Path
import subprocess


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
MARINE_WORK = Path("/tmp/cami2_marine_samples3_5_20260625")
SOURCE_READMAP_SUMMARY = RESULTS / "marine_source_readmap_feasibility_summary.tsv"
OUT_TSV = RESULTS / "marine_source_mapping_local_inventory.tsv"
OUT_AUDIT = RESULTS / "marine_source_mapping_local_inventory_audit.tsv"
OUT_MD = EXP / "MARINE_SOURCE_MAPPING_LOCAL_INVENTORY.md"

FIELDS = ["item", "status", "value", "evidence", "decision"]
AUDIT_FIELDS = ["metric", "value", "evidence", "decision"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def count_named_files(root: Path, names: set[str]) -> dict[str, int]:
    counts = {name: 0 for name in names}
    if not root.exists():
        return counts
    for path in root.rglob("*"):
        if path.is_file() and path.name in counts:
            counts[path.name] += 1
    return counts


def tar_members(path: Path) -> list[str]:
    if not path.exists():
        return []
    proc = subprocess.run(
        ["tar", "-tzf", str(path)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        return [f"ERROR:{proc.stderr.strip()}"]
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def build_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    local_counts = count_named_files(
        MARINE_WORK,
        {"genome_to_id.tsv", "metadata.tsv", "gsa_pooled_mapping.tsv.gz"},
    )
    rows.append(
        {
            "item": "local_marine_workdir",
            "status": str(MARINE_WORK.exists()).lower(),
            "value": str(MARINE_WORK),
            "evidence": str(MARINE_WORK),
            "decision": "workdir_available" if MARINE_WORK.exists() else "workdir_missing",
        }
    )
    for name, count in sorted(local_counts.items()):
        is_setup_metadata = name in {"genome_to_id.tsv", "metadata.tsv"}
        rows.append(
            {
                "item": f"local_marine_{name}",
                "status": "present" if count else "missing",
                "value": count,
                "evidence": str(MARINE_WORK),
                "decision": "local_setup_metadata_present"
                if count and is_setup_metadata
                else "not_a_blocker"
                if count
                else "source_crosswalk_not_found",
            }
        )

    if SOURCE_READMAP_SUMMARY.exists():
        for summary in read_tsv(SOURCE_READMAP_SUMMARY):
            archive = Path(summary["archive"])
            members = tar_members(archive)
            relevant = [
                member
                for member in members
                if any(
                    token in member.lower()
                    for token in ("map", "metadata", "genome", "profile", "gold", "tax")
                )
            ]
            rows.append(
                {
                    "item": f"sample{summary['sample']}_archive_relevant_members",
                    "status": "present" if relevant else "missing",
                    "value": ";".join(relevant[:20]),
                    "evidence": str(archive),
                    "decision": (
                        "readmap_only_no_source_crosswalk"
                        if relevant == [f"simulation_short_read/2018.08.15_09.49.32_sample_{summary['sample']}/reads/reads_mapping.tsv.gz"]
                        else "review_relevant_archive_members"
                    ),
                }
            )

    sibling_paths = [
        Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read/genome_to_id.tsv"),
        Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read/metadata.tsv"),
        Path("/mnt/new3T/minco_cami2_strain_20260621/short_read/genome_to_id.tsv"),
        Path("/mnt/new3T/minco_cami2_strain_20260621/short_read/metadata.tsv"),
    ]
    rows.append(
        {
            "item": "sibling_panel_crosswalk_files",
            "status": "present" if any(path.exists() for path in sibling_paths) else "missing",
            "value": ";".join(str(path) for path in sibling_paths if path.exists()),
            "evidence": "/mnt/new3T/minco_cami2_plant_20260621;/mnt/new3T/minco_cami2_strain_20260621",
            "decision": "sibling_crosswalks_exist_but_do_not_resolve_marine",
        }
    )
    return rows


def build_audit(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    local_missing = [
        row["item"]
        for row in rows
        if str(row["item"]).startswith("local_marine_")
        and row["decision"] == "source_crosswalk_not_found"
    ]
    setup_present = {
        str(row["item"])
        for row in rows
        if str(row["item"]).startswith("local_marine_")
        and row["decision"] == "local_setup_metadata_present"
    }
    readmap_only = [
        row["item"]
        for row in rows
        if row["decision"] == "readmap_only_no_source_crosswalk"
    ]
    return [
        {
            "metric": "local_marine_setup_metadata_status",
            "value": ",".join(sorted(setup_present)) or "none",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "setup_metadata_present"
            if {"local_marine_genome_to_id.tsv", "local_marine_metadata.tsv"} <= setup_present
            else "setup_metadata_incomplete",
        },
        {
            "metric": "local_marine_crosswalk_files_missing",
            "value": ",".join(local_missing) or "none",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "pooled_mapping_missing_but_setup_metadata_present"
            if {"local_marine_genome_to_id.tsv", "local_marine_metadata.tsv"} <= setup_present
            and local_missing
            else "source_crosswalk_not_found"
            if local_missing
            else "review_available_crosswalk",
        },
        {
            "metric": "archive_relevant_member_status",
            "value": f"readmap_only_archives={len(readmap_only)}",
            "evidence": str(OUT_TSV.relative_to(EXP)),
            "decision": "archives_do_not_embed_source_crosswalk"
            if readmap_only
            else "review_archive_members",
        },
        {
            "metric": "promotion_decision",
            "value": "do_not_promote_marine_from_inventory_alone"
            if {"local_marine_genome_to_id.tsv", "local_marine_metadata.tsv"} <= setup_present
            else "do_not_promote_marine_from_local_inventory",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "setup_metadata_present_requires_truth_upgrade_audit"
            if {"local_marine_genome_to_id.tsv", "local_marine_metadata.tsv"} <= setup_present
            else "local_source_mapping_absent",
        },
    ]


def write_md(rows: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    lines = [
        "# Marine Source Mapping Local Inventory",
        "",
        "Date: 2026-06-30",
        "",
        "This audit checks whether the local marine workspace contains the source",
        "crosswalk needed to use marine readmaps as release-grade GTDB truth.",
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
            "## Inventory",
            "",
            "| Item | Status | Value | Decision |",
            "|---|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(f"| `{row['item']}` | {row['status']} | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "The local marine workspace now contains extracted setup metadata",
            "(`genome_to_id.tsv` and `metadata.tsv`). That improves source-name",
            "auditing but is not, by itself, proof that the truth transfer reaches",
            "the release threshold. Use `audit_marine_setup_metadata_truth_upgrade.py`",
            "for the threshold decision.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    audit = build_audit(rows)
    write_tsv(OUT_TSV, rows, FIELDS)
    write_tsv(OUT_AUDIT, audit, AUDIT_FIELDS)
    write_md(rows, audit)
    for row in audit:
        print(row["metric"], row["value"], row["decision"], sep="\t")


if __name__ == "__main__":
    main()
