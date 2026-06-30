#!/usr/bin/env python3
"""Audit whether strainmadness source coverage can form a GTDB holdout panel."""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
TRUTH_ROOT = Path("/mnt/new3T/minco_cami2_strain_20260621/short_read")
GTDB_META = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]
SAMPLES = [0, 1, 2]


def gtdb_species(taxonomy: str) -> str:
    for part in taxonomy.split(";"):
        if part.startswith("s__") and part != "s__":
            return part
    return ""


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def unique_taxid_maps() -> tuple[dict[str, str], dict[str, str], dict[str, int]]:
    ncbi_taxid: defaultdict[str, set[str]] = defaultdict(set)
    ncbi_species_taxid: defaultdict[str, set[str]] = defaultdict(set)
    for path in GTDB_META:
        with gzip.open(path, "rt", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                species = gtdb_species(row["gtdb_taxonomy"])
                if not species:
                    continue
                for key, mapping in [
                    ("ncbi_taxid", ncbi_taxid),
                    ("ncbi_species_taxid", ncbi_species_taxid),
                ]:
                    value = (row.get(key) or "").strip()
                    if value and value != "none":
                        mapping[value].add(species)
    stats = {
        "unique_ncbi_taxid": sum(1 for value in ncbi_taxid.values() if len(value) == 1),
        "ambiguous_ncbi_taxid": sum(1 for value in ncbi_taxid.values() if len(value) > 1),
        "unique_ncbi_species_taxid": sum(1 for value in ncbi_species_taxid.values() if len(value) == 1),
        "ambiguous_ncbi_species_taxid": sum(1 for value in ncbi_species_taxid.values() if len(value) > 1),
    }
    unique_taxid = {key: next(iter(value)) for key, value in ncbi_taxid.items() if len(value) == 1}
    unique_species_taxid = {
        key: next(iter(value)) for key, value in ncbi_species_taxid.items() if len(value) == 1
    }
    return unique_taxid, unique_species_taxid, stats


def load_metadata() -> dict[str, dict[str, str]]:
    with (TRUTH_ROOT / "metadata.tsv").open(newline="") as handle:
        return {row["genome_ID"]: row for row in csv.DictReader(handle, delimiter="\t")}


def map_taxid(
    taxid: str,
    unique_taxid: dict[str, str],
    unique_species_taxid: dict[str, str],
) -> tuple[str, str]:
    if taxid in unique_species_taxid:
        return "ncbi_species_taxid", unique_species_taxid[taxid]
    if taxid in unique_taxid:
        return "ncbi_taxid", unique_taxid[taxid]
    return "", ""


def source_id(raw_id: str) -> str:
    return raw_id[:-6] if raw_id.endswith(".fasta") else raw_id


def main() -> int:
    unique_taxid, unique_species_taxid, map_stats = unique_taxid_maps()
    metadata = load_metadata()
    summary_rows: list[dict[str, object]] = []
    top_mapped_rows: list[dict[str, object]] = []
    top_unmapped_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        total = 0.0
        positive = 0
        mapped = 0.0
        novelty_mass: Counter[str] = Counter()
        map_method_mass: Counter[str] = Counter()
        top_mapped: list[tuple[float, dict[str, object]]] = []
        top_unmapped: list[tuple[float, dict[str, object]]] = []

        with (TRUTH_ROOT / f"coverage_new{sample}.tsv").open(newline="") as handle:
            for line in handle:
                raw_genome, raw_value = line.rstrip("\n").split("\t")[:2]
                abundance = float(raw_value)
                if abundance <= 0.0:
                    continue
                genome = source_id(raw_genome)
                positive += 1
                total += abundance
                row = metadata.get(genome)
                if not row:
                    novelty_mass["missing_metadata"] += abundance
                    top_unmapped.append(
                        (
                            abundance,
                            {
                                "sample": sample,
                                "genome": genome,
                                "abundance": abundance,
                                "ncbi_id": "",
                                "novelty_category": "missing_metadata",
                            },
                        )
                    )
                    continue
                novelty = row.get("novelty_category", "")
                novelty_mass[novelty] += abundance
                ncbi_id = row.get("NCBI_ID", "")
                map_source, species = map_taxid(ncbi_id, unique_taxid, unique_species_taxid)
                if species:
                    mapped += abundance
                    map_method_mass[map_source] += abundance
                    top_mapped.append(
                        (
                            abundance,
                            {
                                "sample": sample,
                                "genome": genome,
                                "abundance": abundance,
                                "ncbi_id": ncbi_id,
                                "novelty_category": novelty,
                                "map_source": map_source,
                                "gtdb_species": species,
                            },
                        )
                    )
                else:
                    top_unmapped.append(
                        (
                            abundance,
                            {
                                "sample": sample,
                                "genome": genome,
                                "abundance": abundance,
                                "ncbi_id": ncbi_id,
                                "novelty_category": novelty,
                            },
                        )
                    )

        summary_row: dict[str, object] = {
            "sample": sample,
            "positive_sources": positive,
            "total_abundance": total,
            "mapped_abundance": mapped,
            "mapped_abundance_pct": (mapped / total * 100.0) if total else 0.0,
            "release_grade_threshold_pct": 95.0,
            "release_grade_ready": (mapped / total * 100.0) >= 95.0 if total else False,
            **map_stats,
        }
        for key, value in sorted(novelty_mass.items()):
            summary_row[f"{key}_abundance_pct"] = value / total * 100.0 if total else 0.0
        for key, value in sorted(map_method_mass.items()):
            summary_row[f"{key}_abundance_pct"] = value / total * 100.0 if total else 0.0
        summary_rows.append(summary_row)
        top_mapped_rows.extend(
            row for _value, row in sorted(top_mapped, key=lambda item: item[0], reverse=True)[:20]
        )
        top_unmapped_rows.extend(
            row for _value, row in sorted(top_unmapped, key=lambda item: item[0], reverse=True)[:20]
        )

    write_tsv(
        RESULTS / "strain_gtdb_transfer_feasibility_summary.tsv",
        summary_rows,
        [
            "sample",
            "positive_sources",
            "total_abundance",
            "mapped_abundance",
            "mapped_abundance_pct",
            "release_grade_threshold_pct",
            "release_grade_ready",
            "unique_ncbi_taxid",
            "ambiguous_ncbi_taxid",
            "unique_ncbi_species_taxid",
            "ambiguous_ncbi_species_taxid",
            "new_genus_abundance_pct",
            "new_order_abundance_pct",
            "new_strain_abundance_pct",
            "ncbi_species_taxid_abundance_pct",
            "ncbi_taxid_abundance_pct",
        ],
    )
    write_tsv(
        RESULTS / "strain_gtdb_transfer_feasibility_top_mapped.tsv",
        top_mapped_rows,
        [
            "sample",
            "genome",
            "abundance",
            "ncbi_id",
            "novelty_category",
            "map_source",
            "gtdb_species",
        ],
    )
    write_tsv(
        RESULTS / "strain_gtdb_transfer_feasibility_top_unmapped.tsv",
        top_unmapped_rows,
        ["sample", "genome", "abundance", "ncbi_id", "novelty_category"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
