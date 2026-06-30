#!/usr/bin/env python3
"""Audit CAMI3 source-readmap truth-transfer scope.

The existing CAMI3 source-readmap scorer reports 79-85% mapped read rows when
all source rows are used as the denominator. This audit separates
GTDB-profiling in-scope source rows from clearly out-of-scope source rows
using cached CAMI taxonomic profiles and the cached per-source mapping table.
It does not rerun raw profiling.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
SOURCE_GENOMES = RESULTS / "cami3_gtdb_source_readmap_source_genomes.tsv"
QUALITY = RESULTS / "cami3_gtdb_source_readmap_quality.tsv"
TAX_PROFILE_ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles")

OUT_SUMMARY = RESULTS / "cami3_source_readmap_scope_summary.tsv"
OUT_CATEGORY = RESULTS / "cami3_source_readmap_unmapped_by_category.tsv"
OUT_TOP = RESULTS / "cami3_source_readmap_top_unmapped_in_scope.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_scope_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_SCOPE_AUDIT.md"

RELEASE_THRESHOLD = 95.0


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


def load_tax_profile(sample: int) -> dict[str, dict[str, str]]:
    path = TAX_PROFILE_ROOT / f"taxonomic_profile_{sample}.txt"
    rows: dict[str, dict[str, str]] = {}
    if not path.exists():
        return rows
    header: list[str] | None = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("@@"):
                header = line[2:].rstrip("\n").split("\t")
                continue
            if line.startswith("@") or not line.strip() or header is None:
                continue
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values + [""] * max(0, len(header) - len(values))))
            genome_id = row.get("_CAMI_genomeID", "")
            if genome_id:
                rows[genome_id] = row
    return rows


def category_for(row: dict[str, str] | None, source_id: str) -> str:
    taxpath = row.get("TAXPATHSN", "") if row else ""
    if "Bacteria" in taxpath or "Archaea" in taxpath:
        return "gtdb_profile_scope"
    if "Homo sapiens" in taxpath or source_id.startswith("hASV"):
        return "human_host"
    if taxpath.startswith("Viruses") or source_id.startswith("vASV"):
        return "viral"
    if "Fungi" in taxpath or source_id.startswith("fASV"):
        return "fungal_eukaryotic"
    if "plasmid" in taxpath or source_id.startswith("pASV") or taxpath.startswith("other entries"):
        return "plasmid_or_other"
    if "Eukaryota" in taxpath:
        return "other_eukaryotic"
    if row is None:
        return "profile_missing"
    return "other_taxonomy"


def species_label(row: dict[str, str] | None) -> str:
    if not row:
        return ""
    parts = [part for part in row.get("TAXPATHSN", "").split("|") if part]
    if not parts:
        return ""
    if row.get("RANK", "") == "strain" and len(parts) >= 2:
        return parts[-2]
    return parts[-1]


def build_rows() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    quality_by_sample = {row["sample"]: row for row in read_tsv(QUALITY)}
    tax_profiles = {str(sample): load_tax_profile(sample) for sample in [0, 1, 2]}
    source_rows = read_tsv(SOURCE_GENOMES)
    category: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    summary: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    top_rows: list[dict[str, object]] = []

    for row in source_rows:
        sample = row["sample"]
        source_id = row["source_genome_id"]
        tax_row = tax_profiles.get(sample, {}).get(source_id)
        cat = category_for(tax_row, source_id)
        read_rows = int(float(row["read_rows"]))
        unmapped = int(float(row["unmapped_read_rows"]))
        mapped = read_rows - unmapped
        category[(sample, cat)]["source_genomes"] += 1
        category[(sample, cat)]["read_rows"] += read_rows
        category[(sample, cat)]["mapped_rows"] += mapped
        category[(sample, cat)]["unmapped_rows"] += unmapped
        summary[sample]["total_read_rows"] += read_rows
        summary[sample]["mapped_rows"] += mapped
        summary[sample]["unmapped_rows"] += unmapped
        if cat == "gtdb_profile_scope":
            summary[sample]["profile_scope_read_rows"] += read_rows
            summary[sample]["profile_scope_mapped_rows"] += mapped
            summary[sample]["profile_scope_unmapped_rows"] += unmapped
            if unmapped > 0:
                top_rows.append(
                    {
                        "sample": sample,
                        "source_genome_id": source_id,
                        "taxid": row.get("taxid", ""),
                        "species_label": species_label(tax_row),
                        "read_rows": read_rows,
                        "mapped_rows": mapped,
                        "unmapped_rows": unmapped,
                        "unmapped_pct_of_profile_scope": 0.0,
                        "primary_mapping_method": row.get("primary_mapping_method", ""),
                    }
                )
        elif unmapped > 0:
            summary[sample]["out_of_scope_unmapped_rows"] += unmapped

    summary_rows: list[dict[str, object]] = []
    for sample in sorted(summary, key=lambda x: int(x)):
        data = summary[sample]
        profile_total = data["profile_scope_read_rows"]
        profile_mapped = data["profile_scope_mapped_rows"]
        profile_unmapped = data["profile_scope_unmapped_rows"]
        profile_pct = 100.0 * profile_mapped / profile_total if profile_total else 0.0
        required = math.ceil(RELEASE_THRESHOLD / 100.0 * profile_total)
        additional_needed = max(0, required - profile_mapped)
        quality = quality_by_sample.get(sample, {})
        summary_rows.append(
            {
                "sample": sample,
                "all_rows_mapped_pct": finite(quality.get("read_rows_mapped_pct")),
                "profile_scope_read_rows": profile_total,
                "profile_scope_mapped_rows": profile_mapped,
                "profile_scope_unmapped_rows": profile_unmapped,
                "profile_scope_mapped_pct": profile_pct,
                "additional_profile_scope_rows_needed_for_95pct": additional_needed,
                "out_of_scope_unmapped_rows": data["out_of_scope_unmapped_rows"],
                "release_threshold_pct": RELEASE_THRESHOLD,
                "profile_scope_release_ready_now": str(profile_pct >= RELEASE_THRESHOLD).lower(),
            }
        )

    top_rows = sorted(
        top_rows,
        key=lambda row: (int(row["sample"]), -int(row["unmapped_rows"])),
    )
    profile_totals = {
        str(row["sample"]): int(row["profile_scope_read_rows"]) for row in summary_rows
    }
    for row in top_rows:
        denom = profile_totals.get(str(row["sample"]), 0)
        row["unmapped_pct_of_profile_scope"] = (
            100.0 * int(row["unmapped_rows"]) / denom if denom else 0.0
        )

    category_rows: list[dict[str, object]] = []
    for (sample, cat), data in sorted(category.items(), key=lambda item: (int(item[0][0]), item[0][1])):
        read_rows = data["read_rows"]
        category_rows.append(
            {
                "sample": sample,
                "category": cat,
                "source_genomes": data["source_genomes"],
                "read_rows": read_rows,
                "mapped_rows": data["mapped_rows"],
                "unmapped_rows": data["unmapped_rows"],
                "mapped_pct": 100.0 * data["mapped_rows"] / read_rows if read_rows else 0.0,
            }
        )

    return summary_rows, category_rows, top_rows


def audit_rows(summary: list[dict[str, object]], top_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    below = [
        row for row in summary if str(row["profile_scope_release_ready_now"]).lower() != "true"
    ]
    needed = sum(int(row["additional_profile_scope_rows_needed_for_95pct"]) for row in summary)
    top_rescue = []
    for sample in sorted({str(row["sample"]) for row in below}, key=int):
        rows = [row for row in top_rows if str(row["sample"]) == sample]
        top_rescue.extend(row["source_genome_id"] for row in rows[:3])
    return [
        {
            "metric": "all_rows_mapping_status",
            "value": ";".join(
                f"sample{row['sample']}={float(row['all_rows_mapped_pct']):.3f}%"
                for row in summary
            ),
            "evidence": "cami3_gtdb_source_readmap_quality.tsv",
            "decision": "all_rows_include_out_of_scope_sources",
        },
        {
            "metric": "profile_scope_mapping_status",
            "value": ";".join(
                f"sample{row['sample']}={float(row['profile_scope_mapped_pct']):.3f}%"
                for row in summary
            ),
            "evidence": "cami3_source_readmap_scope_summary.tsv",
            "decision": "profile_scope_nearly_release_grade",
        },
        {
            "metric": "samples_below_95pct_profile_scope",
            "value": ",".join(str(row["sample"]) for row in below) if below else "none",
            "evidence": "cami3_source_readmap_scope_summary.tsv",
            "decision": "targeted_mapping_needed" if below else "profile_scope_release_ready",
        },
        {
            "metric": "additional_rows_needed_for_95pct",
            "value": needed,
            "evidence": "cami3_source_readmap_scope_summary.tsv",
            "decision": "feasible_if_top_sources_can_be_mapped" if needed else "none_needed",
        },
        {
            "metric": "top_in_scope_sources_to_review",
            "value": ",".join(top_rescue[:12]) if top_rescue else "none",
            "evidence": "cami3_source_readmap_top_unmapped_in_scope.tsv",
            "decision": "manual_or_metadata_mapping_review",
        },
        {
            "metric": "release_upgrade_decision",
            "value": "do_not_promote_yet",
            "evidence": "cami3_source_readmap_scope_audit.tsv",
            "decision": "needs_targeted_in_scope_source_mapping_before_release_grade",
        },
    ]


def write_markdown(
    summary: list[dict[str, object]],
    category_rows: list[dict[str, object]],
    top_rows: list[dict[str, object]],
    audit: list[dict[str, object]],
) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# CAMI3 Source-Readmap Scope Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This cached audit separates CAMI3 source-readmap rows that are in scope",
        "for GTDB species profiling from host, viral, fungal/eukaryotic, plasmid,",
        "and other out-of-scope rows. It does not rerun profiling.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "all_rows_mapping_status",
        "profile_scope_mapping_status",
        "samples_below_95pct_profile_scope",
        "additional_rows_needed_for_95pct",
        "top_in_scope_sources_to_review",
        "release_upgrade_decision",
    ]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Sample | All rows mapped % | Profile-scope mapped % | Rows needed for 95% | Scope ready |",
            "|---:|---:|---:|---:|---|",
        ]
    )
    for row in summary:
        lines.append(
            "| {sample} | {all_pct:.3f} | {scope_pct:.3f} | {needed} | {ready} |".format(
                sample=row["sample"],
                all_pct=float(row["all_rows_mapped_pct"]),
                scope_pct=float(row["profile_scope_mapped_pct"]),
                needed=int(row["additional_profile_scope_rows_needed_for_95pct"]),
                ready=row["profile_scope_release_ready_now"],
            )
        )

    lines.extend(
        [
            "",
            "## Top In-Scope Unmapped Sources",
            "",
            "| Sample | Source | Taxid | Species label | Unmapped rows |",
            "|---:|---|---:|---|---:|",
        ]
    )
    for row in top_rows[:20]:
        lines.append(
            "| {sample} | {source} | {taxid} | {species} | {unmapped} |".format(
                sample=row["sample"],
                source=row["source_genome_id"],
                taxid=row["taxid"],
                species=str(row["species_label"]).replace("|", "/"),
                unmapped=int(row["unmapped_rows"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- CAMI3 should remain diagnostic for now.",
            "- After excluding clear out-of-scope rows, sample 2 crosses 95% mapped",
            "  profile-scope coverage; samples 0 and 1 are close but still below 95%.",
            "- The next CAMI3 release-upgrade task is targeted mapping of the top",
            "  in-scope unmapped sources, not another MinCO call-threshold sweep.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_CATEGORY.relative_to(EXP)}`",
            f"- `{OUT_TOP.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    summary, category_rows, top_rows = build_rows()
    audit = audit_rows(summary, top_rows)
    write_tsv(
        OUT_SUMMARY,
        summary,
        [
            "sample",
            "all_rows_mapped_pct",
            "profile_scope_read_rows",
            "profile_scope_mapped_rows",
            "profile_scope_unmapped_rows",
            "profile_scope_mapped_pct",
            "additional_profile_scope_rows_needed_for_95pct",
            "out_of_scope_unmapped_rows",
            "release_threshold_pct",
            "profile_scope_release_ready_now",
        ],
    )
    write_tsv(
        OUT_CATEGORY,
        category_rows,
        ["sample", "category", "source_genomes", "read_rows", "mapped_rows", "unmapped_rows", "mapped_pct"],
    )
    write_tsv(
        OUT_TOP,
        top_rows,
        [
            "sample",
            "source_genome_id",
            "taxid",
            "species_label",
            "read_rows",
            "mapped_rows",
            "unmapped_rows",
            "unmapped_pct_of_profile_scope",
            "primary_mapping_method",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(summary, category_rows, top_rows, audit)
    print("\n".join(f"{row['sample']}\t{row['profile_scope_mapped_pct']:.3f}\t{row['additional_profile_scope_rows_needed_for_95pct']}" for row in summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
