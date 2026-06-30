#!/usr/bin/env python3
"""Audit whether CAMI II plant source abundance can form a GTDB holdout panel."""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
TRUTH_ROOT = Path("/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read")
GTDB_META = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]
SAMPLES = [3, 4, 5]


def gtdb_species(taxonomy: str) -> str:
    for part in taxonomy.split(";"):
        if part.startswith("s__") and part != "s__":
            return part
    return ""


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


def load_camitax() -> dict[str, dict[str, str]]:
    return {row["Genome"]: row for row in read_tsv(TRUTH_ROOT / "camitax.tsv")}


def load_otu_genome() -> dict[str, str]:
    out: dict[str, str] = {}
    with (TRUTH_ROOT / "genome_to_id.tsv").open(newline="") as handle:
        for line in handle:
            otu, source_path = line.rstrip("\n").split("\t")[:2]
            genome = Path(source_path).name
            if genome.endswith(".fasta"):
                genome = genome[:-6]
            out[otu] = genome
    return out


def candidate_taxids(genome: str, camitax: dict[str, dict[str, str]]) -> list[tuple[str, str, str, str]]:
    candidates: list[tuple[str, str, str, str]] = []
    row = camitax.get(genome)
    if row and row.get("taxID"):
        candidates.append(
            (
                "camitax_taxID",
                row["taxID"],
                row.get("taxName", ""),
                row.get("taxLvl", ""),
            )
        )
    match = re.match(r"^(\d+)\.", genome)
    if match:
        candidates.append(("filename_taxID", match.group(1), "", "filename"))
    return candidates


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


def main() -> int:
    unique_taxid, unique_species_taxid, map_stats = unique_taxid_maps()
    camitax = load_camitax()
    otu_genome = load_otu_genome()

    summary_rows: list[dict[str, object]] = []
    top_mapped_rows: list[dict[str, object]] = []
    top_unmapped_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        total = 0.0
        positive = 0
        mapped = 0.0
        source_mass: Counter[str] = Counter()
        map_method_mass: Counter[str] = Counter()
        top_mapped: list[tuple[float, dict[str, object]]] = []
        top_unmapped: list[tuple[float, dict[str, object]]] = []

        with (TRUTH_ROOT / f"abundance{sample}.tsv").open(newline="") as handle:
            for line in handle:
                otu, raw_value = line.rstrip("\n").split("\t")[:2]
                abundance = float(raw_value)
                if abundance <= 0.0:
                    continue
                positive += 1
                total += abundance
                genome = otu_genome.get(otu, "")
                tax_candidates = candidate_taxids(genome, camitax)
                if camitax.get(genome):
                    source_mass[f"camitax_{camitax[genome].get('taxLvl', 'unknown')}"] += abundance
                elif re.match(r"^\d+\.", genome):
                    source_mass["filename_taxid"] += abundance
                else:
                    source_mass["no_taxid"] += abundance

                hit: tuple[str, str, str, str] | None = None
                for source, taxid, _name, _level in tax_candidates:
                    map_source, species = map_taxid(taxid, unique_taxid, unique_species_taxid)
                    if species:
                        hit = (source, map_source, taxid, species)
                        break
                if hit:
                    source, map_source, taxid, species = hit
                    mapped += abundance
                    map_method_mass[f"{source}->{map_source}"] += abundance
                    top_mapped.append(
                        (
                            abundance,
                            {
                                "sample": sample,
                                "otu": otu,
                                "genome": genome,
                                "abundance": abundance,
                                "candidate_source": source,
                                "taxid": taxid,
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
                                "otu": otu,
                                "genome": genome,
                                "abundance": abundance,
                                "candidate_taxids": ";".join(
                                    f"{source}:{taxid}:{level}:{name}"
                                    for source, taxid, name, level in tax_candidates
                                ),
                            },
                        )
                    )

        row: dict[str, object] = {
            "sample": sample,
            "positive_sources": positive,
            "total_abundance": total,
            "mapped_abundance": mapped,
            "mapped_abundance_pct": (mapped / total * 100.0) if total else 0.0,
            "release_grade_threshold_pct": 95.0,
            "release_grade_ready": (mapped / total * 100.0) >= 95.0 if total else False,
            **map_stats,
        }
        for key, value in sorted(source_mass.items()):
            row[f"{key}_abundance_pct"] = value / total * 100.0 if total else 0.0
        for key, value in sorted(map_method_mass.items()):
            row[f"{key}_abundance_pct"] = value / total * 100.0 if total else 0.0
        summary_rows.append(row)
        top_mapped_rows.extend(
            row for _value, row in sorted(top_mapped, key=lambda item: item[0], reverse=True)[:20]
        )
        top_unmapped_rows.extend(
            row for _value, row in sorted(top_unmapped, key=lambda item: item[0], reverse=True)[:20]
        )

    summary_fields = [
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
        "camitax_class_abundance_pct",
        "camitax_family_abundance_pct",
        "camitax_genus_abundance_pct",
        "camitax_species_abundance_pct",
        "filename_taxid_abundance_pct",
        "no_taxid_abundance_pct",
        "camitax_taxID->ncbi_species_taxid_abundance_pct",
        "filename_taxID->ncbi_species_taxid_abundance_pct",
        "filename_taxID->ncbi_taxid_abundance_pct",
    ]
    write_tsv(RESULTS / "plant_gtdb_transfer_feasibility_summary.tsv", summary_rows, summary_fields)
    write_tsv(
        RESULTS / "plant_gtdb_transfer_feasibility_top_mapped.tsv",
        top_mapped_rows,
        [
            "sample",
            "otu",
            "genome",
            "abundance",
            "candidate_source",
            "taxid",
            "map_source",
            "gtdb_species",
        ],
    )
    write_tsv(
        RESULTS / "plant_gtdb_transfer_feasibility_top_unmapped.tsv",
        top_unmapped_rows,
        ["sample", "otu", "genome", "abundance", "candidate_taxids"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
