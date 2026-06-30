#!/usr/bin/env python3
"""Audit whether local marine read maps can upgrade GTDB truth.

CAMI II marine archives contain `reads_mapping.tsv.gz` files with source
genome IDs and source sequence IDs. This audit samples those read maps and
checks whether the local information is sufficient to map unresolved marine
truth rows to GTDB species without external contig-to-assembly/source metadata.
It does not extract full read FASTQs or rerun profilers.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import tarfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

import score_marine_gtdb_taxid_transfer as marine


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

DETAIL = RESULTS / "marine_binomial_transfer_detail.tsv"
OUT_SUMMARY = RESULTS / "marine_source_readmap_feasibility_summary.tsv"
OUT_TOP = RESULTS / "marine_source_readmap_unresolved_top.tsv"
OUT_AUDIT = RESULTS / "marine_source_readmap_feasibility_audit.tsv"
OUT_MD = EXP / "MARINE_SOURCE_READMAP_FEASIBILITY.md"

SAMPLE_ARCHIVES = {
    3: Path("/tmp/cami2_marine_samples3_5_20260625/data/marmgCAMI2_sample_3_reads.tar.gz"),
    4: Path("/tmp/cami2_marine_samples3_5_20260625/data/marmgCAMI2_sample_4_reads.tar.gz"),
    5: Path("/tmp/cami2_marine_samples3_5_20260625/data/marmgCAMI2_sample_5_reads.tar.gz"),
}


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


def source_sequence_id(read_id: str) -> str:
    base = str(read_id or "").rsplit("/", 1)[0]
    if "-" not in base:
        return base
    return base.rsplit("-", 1)[0]


def sequence_id_class(seq_id: str) -> str:
    value = str(seq_id or "")
    if re.match(r"^(NZ_)?[A-Z]{1,4}[0-9]{5,}\\.[0-9]+$", value):
        return "sequence_accession"
    if value.startswith("NODE_"):
        return "assembly_contig_label"
    if value.startswith(("GCA_", "GCF_")):
        return "assembly_accession_like"
    return "other"


def iter_readmap_rows_from_archive(archive: Path):
    with tarfile.open(archive, "r|gz") as tar:
        for member in tar:
            if not member.name.endswith("reads_mapping.tsv.gz"):
                continue
            fileobj = tar.extractfile(member)
            if fileobj is None:
                return
            with gzip.GzipFile(fileobj=fileobj) as handle:
                for raw in handle:
                    line = raw.decode("utf-8", errors="replace").rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) < 4:
                        continue
                    yield {
                        "anonymous_read_id": fields[0],
                        "genome_id": fields[1],
                        "taxid": fields[2],
                        "read_id": fields[3],
                    }
            return


def read_transfer_detail() -> dict[tuple[int, str], dict[str, object]]:
    detail = pd.read_csv(DETAIL, sep="\t")
    out: dict[tuple[int, str], dict[str, object]] = {}
    for row in detail.itertuples(index=False):
        out[(int(row.sample), str(row.ncbi_species_taxid))] = {
            "truth_name": row.truth_name,
            "transfer_method": row.transfer_method,
            "abundance_pct_all": finite(row.abundance_pct_all),
            "taxid_gtdb_species_count": int(row.taxid_gtdb_species_count),
        }
    return out


def gtdb_accession_sets() -> tuple[set[str], set[str]]:
    by_accession, by_core = marine.truth.load_gtdb_metadata(marine.truth.GTDB_METADATA)
    return set(by_accession), set(by_core)


def audit_sample(
    sample: int,
    archive: Path,
    detail_by_key: dict[tuple[int, str], dict[str, object]],
    gtdb_accessions: set[str],
    gtdb_cores: set[str],
    max_rows: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows_read = 0
    genome_ids: set[str] = set()
    taxids: set[str] = set()
    seq_ids: set[str] = set()
    seq_classes: Counter[str] = Counter()
    direct_gtdb_seq_matches = 0
    direct_gtdb_core_matches = 0
    unresolved_rows = 0
    unresolved_genomes: set[str] = set()
    missing_detail_rows = 0
    unresolved_groups: Counter[tuple[str, str, str, str, str]] = Counter()

    for row in iter_readmap_rows_from_archive(archive):
        rows_read += 1
        genome_id = str(row["genome_id"])
        taxid = str(row["taxid"])
        seq_id = source_sequence_id(str(row["read_id"]))
        detail = detail_by_key.get((sample, taxid), {})
        method = str(detail.get("transfer_method", "not_in_gold_detail"))
        truth_name = str(detail.get("truth_name", ""))

        genome_ids.add(genome_id)
        taxids.add(taxid)
        seq_ids.add(seq_id)
        seq_classes[sequence_id_class(seq_id)] += 1
        if seq_id in gtdb_accessions:
            direct_gtdb_seq_matches += 1
        if seq_id.split(".", 1)[0] in gtdb_cores:
            direct_gtdb_core_matches += 1
        if method == "not_in_gold_detail":
            missing_detail_rows += 1
        if method in {"ambiguous_taxid", "unmapped"}:
            unresolved_rows += 1
            unresolved_genomes.add(genome_id)
            unresolved_groups[(genome_id, taxid, truth_name, method, seq_id)] += 1
        if max_rows and rows_read >= max_rows:
            break

    top_rows = [
        {
            "sample": sample,
            "genome_id": genome_id,
            "taxid": taxid,
            "truth_name": truth_name,
            "transfer_method": method,
            "source_sequence_id": seq_id,
            "source_sequence_id_class": sequence_id_class(seq_id),
            "sampled_read_rows": count,
            "taxid_gtdb_species_count": detail_by_key.get((sample, taxid), {}).get(
                "taxid_gtdb_species_count", ""
            ),
            "truth_abundance_pct_all": detail_by_key.get((sample, taxid), {}).get(
                "abundance_pct_all", ""
            ),
        }
        for (genome_id, taxid, truth_name, method, seq_id), count in unresolved_groups.most_common(100)
    ]
    summary = {
        "sample": sample,
        "archive": str(archive),
        "archive_present": str(archive.exists()).lower(),
        "max_rows": max_rows,
        "rows_read": rows_read,
        "distinct_genome_ids": len(genome_ids),
        "distinct_taxids": len(taxids),
        "distinct_source_sequence_ids": len(seq_ids),
        "sequence_accession_rows": seq_classes["sequence_accession"],
        "assembly_contig_label_rows": seq_classes["assembly_contig_label"],
        "assembly_accession_like_rows": seq_classes["assembly_accession_like"],
        "other_sequence_id_rows": seq_classes["other"],
        "direct_gtdb_sequence_accession_row_matches": direct_gtdb_seq_matches,
        "direct_gtdb_core_row_matches": direct_gtdb_core_matches,
        "unresolved_transfer_rows": unresolved_rows,
        "unresolved_distinct_genome_ids": len(unresolved_genomes),
        "not_in_gold_detail_rows": missing_detail_rows,
    }
    return summary, top_rows


def build_audit(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    sampled = sum(int(row["rows_read"]) for row in summary_rows)
    direct = sum(int(row["direct_gtdb_sequence_accession_row_matches"]) for row in summary_rows)
    core = sum(int(row["direct_gtdb_core_row_matches"]) for row in summary_rows)
    unresolved = sum(int(row["unresolved_transfer_rows"]) for row in summary_rows)
    return [
        {
            "metric": "sampled_readmap_rows",
            "value": sampled,
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "source_readmaps_available_in_local_archives",
        },
        {
            "metric": "direct_sequence_to_gtdb_matches",
            "value": f"sequence_accession_rows={direct};core_rows={core}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "no_direct_contig_to_gtdb_assembly_mapping",
        },
        {
            "metric": "unresolved_transfer_rows_sampled",
            "value": unresolved,
            "evidence": str(OUT_TOP.relative_to(EXP)),
            "decision": "readmaps_cover_unresolved_taxids_but_need_source_mapping",
        },
        {
            "metric": "promotion_decision",
            "value": "do_not_promote_marine_source_readmap_from_local_ids",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "contig_or_otu_to_assembly_mapping_missing",
        },
    ]


def write_markdown(summary_rows: list[dict[str, object]], audit_rows: list[dict[str, object]]) -> None:
    lines = [
        "# Marine Source Readmap Feasibility",
        "",
        "Date: 2026-06-29",
        "",
        "This audit samples source read maps embedded in local marine archives.",
        "It checks whether local read IDs can directly resolve unresolved marine",
        "truth rows to GTDB assemblies. It does not extract FASTQ files.",
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
            "## Samples",
            "",
            "| Sample | Rows Read | Genome IDs | Source Seq IDs | Unresolved Rows | Direct GTDB Matches |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary_rows:
        lines.append(
            "| {sample} | {rows_read} | {distinct_genome_ids} | {distinct_source_sequence_ids} | {unresolved_transfer_rows} | {direct_gtdb_sequence_accession_row_matches} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "The read maps are useful, but the local IDs are sequence/contig IDs or",
            "CAMI genome IDs. They do not directly match GTDB assembly accessions in",
            "the local metadata. Marine still needs contig/OTU-to-assembly metadata,",
            "source FASTA provenance, or another clean holdout before release use.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_TOP.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-rows",
        type=int,
        default=2_000_000,
        help="Rows to sample per archive; use 0 for full read maps.",
    )
    args = parser.parse_args()

    detail_by_key = read_transfer_detail()
    gtdb_accessions, gtdb_cores = gtdb_accession_sets()
    summary_rows: list[dict[str, object]] = []
    top_rows: list[dict[str, object]] = []
    for sample, archive in SAMPLE_ARCHIVES.items():
        summary, top = audit_sample(
            sample,
            archive,
            detail_by_key,
            gtdb_accessions,
            gtdb_cores,
            args.max_rows,
        )
        summary_rows.append(summary)
        top_rows.extend(top)
    audit_rows = build_audit(summary_rows)
    write_tsv(
        OUT_SUMMARY,
        summary_rows,
        [
            "sample",
            "archive",
            "archive_present",
            "max_rows",
            "rows_read",
            "distinct_genome_ids",
            "distinct_taxids",
            "distinct_source_sequence_ids",
            "sequence_accession_rows",
            "assembly_contig_label_rows",
            "assembly_accession_like_rows",
            "other_sequence_id_rows",
            "direct_gtdb_sequence_accession_row_matches",
            "direct_gtdb_core_row_matches",
            "unresolved_transfer_rows",
            "unresolved_distinct_genome_ids",
            "not_in_gold_detail_rows",
        ],
    )
    write_tsv(
        OUT_TOP,
        top_rows,
        [
            "sample",
            "genome_id",
            "taxid",
            "truth_name",
            "transfer_method",
            "source_sequence_id",
            "source_sequence_id_class",
            "sampled_read_rows",
            "taxid_gtdb_species_count",
            "truth_abundance_pct_all",
        ],
    )
    write_tsv(OUT_AUDIT, audit_rows, ["metric", "value", "evidence", "decision"])
    write_markdown(summary_rows, audit_rows)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
