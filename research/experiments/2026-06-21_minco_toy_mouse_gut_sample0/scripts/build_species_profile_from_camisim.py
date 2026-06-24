#!/usr/bin/env python3
"""Build a species-level CAMI-like profile from CAMISIM distributions."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Sequence


def load_nodes(path: Path) -> tuple[Dict[str, str], Dict[str, str]]:
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 3:
                taxid = parts[0]
                parent[taxid] = parts[1]
                rank[taxid] = parts[2]
    return parent, rank


def load_names(path: Path) -> Dict[str, str]:
    names: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[parts[0]] = parts[1]
    return names


def lineage_to_rank(taxid: str, parent: Dict[str, str], rank: Dict[str, str], target_rank: str) -> str:
    cur = str(taxid)
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        if rank.get(cur) == target_rank:
            return cur
        nxt = parent.get(cur)
        if not nxt or nxt == cur:
            break
        cur = nxt
    return ""


def lineage_path(taxid: str, parent: Dict[str, str], names: Dict[str, str]) -> tuple[str, str]:
    values = []
    cur = str(taxid)
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        values.append(cur)
        nxt = parent.get(cur)
        if not nxt or nxt == cur:
            break
        cur = nxt
    values.reverse()
    if values and values[0] == "1":
        values = values[1:]
    return "|".join(values), "|".join(names.get(t, t) for t in values)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--distribution", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--nodes", type=Path, required=True)
    ap.add_argument("--names", type=Path, required=True)
    ap.add_argument("--sample-id", default="mouse0")
    ap.add_argument("--out-profile", type=Path, required=True)
    ap.add_argument("--out-summary", type=Path, required=True)
    args = ap.parse_args(argv)

    parent, rank = load_nodes(args.nodes)
    names = load_names(args.names)

    genome_to_taxid: Dict[str, str] = {}
    with args.metadata.open() as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx_genome = header.index("genome_ID")
        idx_taxid = header.index("NCBI_ID")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            genome_to_taxid[parts[idx_genome]] = parts[idx_taxid]

    species_abundance: Dict[str, float] = defaultdict(float)
    unmapped = []
    raw_total = 0.0
    positive_genomes = 0
    with args.distribution.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            genome_id, abundance_s = line.rstrip("\n").split("\t")[:2]
            abundance = float(abundance_s)
            if abundance <= 0.0:
                continue
            positive_genomes += 1
            raw_total += abundance
            taxid = genome_to_taxid.get(genome_id, "")
            species = lineage_to_rank(taxid, parent, rank, "species") if taxid else ""
            if species:
                species_abundance[species] += abundance
            else:
                unmapped.append((genome_id, taxid, abundance))

    if raw_total <= 0.0:
        raise SystemExit("distribution has no positive abundance")

    args.out_profile.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    with args.out_profile.open("w") as out:
        out.write(f"@SampleID:{args.sample_id}\n")
        out.write("@Version:0.9.1\n")
        out.write("@Ranks:superkingdom|phylum|class|order|family|genus|species\n\n")
        out.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n")
        for species, abundance in sorted(species_abundance.items(), key=lambda item: (-item[1], item[0])):
            percentage = 100.0 * abundance / raw_total
            taxpath, taxpathsn = lineage_path(species, parent, names)
            out.write(f"{species}\tspecies\t{taxpath}\t{taxpathsn}\t{percentage:.10f}\t\t\n")

    with args.out_summary.open("w") as out:
        out.write("metric\tvalue\n")
        out.write(f"positive_genomes\t{positive_genomes}\n")
        out.write(f"raw_abundance_sum\t{raw_total:.10f}\n")
        out.write(f"species_count\t{len(species_abundance)}\n")
        out.write(f"unmapped_positive_genomes\t{len(unmapped)}\n")
        out.write(f"profile_percentage_sum\t{sum(100.0 * v / raw_total for v in species_abundance.values()):.10f}\n")
        for genome_id, taxid, abundance in unmapped:
            out.write(f"unmapped\t{genome_id}:{taxid}:{abundance}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
