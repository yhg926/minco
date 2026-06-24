#!/usr/bin/env python3
"""Join the L. crispatus readwise trace to CAMISIM read origins."""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_TRACE = Path(
    "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/"
    "l_crispatus_s2000_marker_split_readwise_assignment_trace.tsv"
)
DEFAULT_TAR = Path("/mnt/new3T/minco_cami2_toymouse_20260621/2017.12.29_11.37.26_sample_0_reads.tar")
DEFAULT_TAR_MEMBER = "2017.12.29_11.37.26_sample_0/reads/reads_mapping.tsv.gz"
DEFAULT_SOURCE = Path("/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/positive_source_genomes.tsv")
DEFAULT_NODES = Path("/mnt/new3T/gtdbr220/gtdbr226/taxdump/nodes.dmp")
DEFAULT_NAMES = Path("/mnt/new3T/gtdbr220/gtdbr226/taxdump/names.dmp")
TRUE_LCRI_ACCESSION = "GCF_002218965.1"


def parse_dmp_fields(line: str) -> list[str]:
    return [field.strip() for field in line.rstrip("\n").split("|")[:-1]]


def load_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            fields = parse_dmp_fields(line)
            if len(fields) >= 4 and fields[3] == "scientific name":
                names[fields[0]] = fields[1]
    return names


def load_nodes(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    parent: dict[str, str] = {}
    rank: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            fields = parse_dmp_fields(line)
            if len(fields) >= 3:
                parent[fields[0]] = fields[1]
                rank[fields[0]] = fields[2]
    return parent, rank


def species_taxid(taxid: str, parent: dict[str, str], rank: dict[str, str]) -> str:
    cur = str(taxid)
    seen = set()
    while cur and cur not in seen:
        if rank.get(cur) == "species":
            return cur
        seen.add(cur)
        nxt = parent.get(cur, "")
        if not nxt or nxt == cur:
            break
        cur = nxt
    return str(taxid)


def load_sources(path: Path, names: dict[str, str], parent: dict[str, str], rank: dict[str, str]) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            source_taxid = str(row.get("source_taxid", ""))
            sp_taxid = species_taxid(source_taxid, parent, rank)
            row["source_name"] = names.get(source_taxid, "")
            row["source_species_taxid"] = sp_taxid
            row["source_species_name"] = names.get(sp_taxid, names.get(source_taxid, ""))
            out[str(row["genome_id"])] = row
    return out


def load_trace(path: Path) -> tuple[list[dict[str, str]], set[str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows, {row["read_id"] for row in rows}


def load_mapping_for_reads(tar_path: Path, member_name: str, wanted: set[str]) -> dict[str, dict[str, str]]:
    found: dict[str, dict[str, str]] = {}
    with tarfile.open(tar_path, "r") as tf:
        member = tf.getmember(member_name)
        raw = tf.extractfile(member)
        if raw is None:
            raise RuntimeError(f"cannot extract {member_name} from {tar_path}")
        with gzip.GzipFile(fileobj=raw) as gz:
            text = (line.decode("utf-8").rstrip("\n") for line in gz)
            header = next(text).lstrip("#").split("\t")
            for line in text:
                fields = line.split("\t")
                if len(fields) != len(header):
                    continue
                rec = dict(zip(header, fields))
                rid = rec.get("anonymous_read_id", "")
                if rid in wanted:
                    found[rid] = rec
                    if len(found) == len(wanted):
                        break
    return found


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trace", type=Path, default=DEFAULT_TRACE)
    parser.add_argument("--reads-tar", type=Path, default=DEFAULT_TAR)
    parser.add_argument("--reads-member", default=DEFAULT_TAR_MEMBER)
    parser.add_argument("--source-genomes", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--nodes", type=Path, default=DEFAULT_NODES)
    parser.add_argument("--names", type=Path, default=DEFAULT_NAMES)
    parser.add_argument("--out-prefix", type=Path, default=DEFAULT_TRACE.with_suffix(""))
    args = parser.parse_args()

    names = load_names(args.names)
    parent, rank = load_nodes(args.nodes)
    sources = load_sources(args.source_genomes, names, parent, rank)
    trace_rows, wanted = load_trace(args.trace)
    mapping = load_mapping_for_reads(args.reads_tar, args.reads_member, wanted)

    joined: list[dict[str, str]] = []
    for row in trace_rows:
        rec = dict(row)
        truth = mapping.get(row["read_id"], {})
        genome_id = truth.get("genome_id", "")
        source = sources.get(genome_id, {})
        rec.update(
            {
                "genome_id": genome_id,
                "mapping_tax_id": truth.get("tax_id", ""),
                "orig_read_id": truth.get("read_id", ""),
                "source_taxid": source.get("source_taxid", ""),
                "source_name": source.get("source_name", ""),
                "source_species_taxid": source.get("source_species_taxid", ""),
                "source_species_name": source.get("source_species_name", ""),
                "source_accession": source.get("source_accession", ""),
                "source_basename": source.get("source_basename", ""),
                "is_true_l_crispatus_source": str(source.get("source_accession", "") == TRUE_LCRI_ACCESSION),
            }
        )
        joined.append(rec)

    out_prefix = args.out_prefix
    write_tsv(
        out_prefix.with_name(out_prefix.name + "_with_truth.tsv"),
        joined,
        list(joined[0].keys()) if joined else [],
    )

    total_events = len(joined)
    total_reads = len({row["read_id"] for row in joined})

    by_genome: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in joined:
        by_genome[row["genome_id"]].append(row)
        key = row["source_species_taxid"] or row["source_taxid"] or row["mapping_tax_id"] or "NA"
        by_species[key].append(row)

    genome_rows: list[dict[str, object]] = []
    for genome_id, rows in by_genome.items():
        source_accessions = sorted({r["source_accession"] for r in rows if r["source_accession"]})
        genome_rows.append(
            {
                "genome_id": genome_id,
                "event_n": len(rows),
                "unique_read_n": len({r["read_id"] for r in rows}),
                "event_fraction": len(rows) / total_events if total_events else 0.0,
                "read_fraction": len({r["read_id"] for r in rows}) / total_reads if total_reads else 0.0,
                "mapping_tax_id": rows[0]["mapping_tax_id"],
                "source_taxid": rows[0]["source_taxid"],
                "source_name": rows[0]["source_name"],
                "source_species_taxid": rows[0]["source_species_taxid"],
                "source_species_name": rows[0]["source_species_name"],
                "source_accession": ",".join(source_accessions),
                "source_basename": rows[0]["source_basename"],
                "is_true_l_crispatus_source": rows[0]["is_true_l_crispatus_source"],
            }
        )
    genome_rows.sort(key=lambda r: (-int(r["event_n"]), str(r["genome_id"])))
    write_tsv(out_prefix.with_name(out_prefix.name + "_origin_by_genome.tsv"), genome_rows, list(genome_rows[0].keys()))

    species_rows: list[dict[str, object]] = []
    for species_id, rows in by_species.items():
        source_accessions = sorted({r["source_accession"] for r in rows if r["source_accession"]})
        genomes = sorted({r["genome_id"] for r in rows if r["genome_id"]})
        species_rows.append(
            {
                "source_species_taxid": species_id,
                "source_species_name": rows[0]["source_species_name"] or names.get(species_id, ""),
                "event_n": len(rows),
                "unique_read_n": len({r["read_id"] for r in rows}),
                "event_fraction": len(rows) / total_events if total_events else 0.0,
                "read_fraction": len({r["read_id"] for r in rows}) / total_reads if total_reads else 0.0,
                "genome_n": len(genomes),
                "genome_ids": ",".join(genomes),
                "source_accessions": ",".join(source_accessions),
                "contains_true_l_crispatus_source": str(any(r["source_accession"] == TRUE_LCRI_ACCESSION for r in rows)),
            }
        )
    species_rows.sort(key=lambda r: (-int(r["event_n"]), str(r["source_species_taxid"])))
    write_tsv(out_prefix.with_name(out_prefix.name + "_origin_by_species.tsv"), species_rows, list(species_rows[0].keys()))

    ctx_rows: list[dict[str, object]] = []
    by_ctx: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in joined:
        by_ctx[row["ref_begin"]].append(row)
    for ref_begin, rows in by_ctx.items():
        species_counter = Counter(r["source_species_name"] or r["source_name"] or r["genome_id"] for r in rows)
        genome_counter = Counter(r["genome_id"] for r in rows)
        true_rows = [r for r in rows if r["source_accession"] == TRUE_LCRI_ACCESSION]
        best_rows = [r for r in rows if r["diff"] == r["best_diff"]]
        ctx_rows.append(
            {
                "ref_begin": ref_begin,
                "best_diff": min(int(r["best_diff"]) for r in rows),
                "event_n": len(rows),
                "unique_read_n": len({r["read_id"] for r in rows}),
                "true_l_crispatus_event_n": len(true_rows),
                "true_l_crispatus_read_n": len({r["read_id"] for r in true_rows}),
                "best_event_n": len(best_rows),
                "best_read_n": len({r["read_id"] for r in best_rows}),
                "top_species": species_counter.most_common(1)[0][0],
                "top_species_events": species_counter.most_common(1)[0][1],
                "top_genome": genome_counter.most_common(1)[0][0],
                "top_genome_events": genome_counter.most_common(1)[0][1],
                "all_species": ";".join(f"{k}:{v}" for k, v in species_counter.most_common()),
                "all_genomes": ";".join(f"{k}:{v}" for k, v in genome_counter.most_common()),
            }
        )
    ctx_rows.sort(key=lambda r: (-int(r["event_n"]), str(r["ref_begin"])))
    write_tsv(out_prefix.with_name(out_prefix.name + "_final_context_origin_summary.tsv"), ctx_rows, list(ctx_rows[0].keys()))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
