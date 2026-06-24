#!/usr/bin/env python3
"""Build a CAMISIM sample-aware taxmap from GTDB calls to source NCBI species.

The input skani table compares called GTDB representatives against the positive
source genomes for one CAMISIM sample. For called representatives close enough
to a source genome, this script emits taxmap override rows that map the called
accession to the source genome's NCBI species taxid. The overrides can be
prepended to the normal GTDB taxmap for benchmark scoring.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def load_nodes(path: Path) -> tuple[Dict[str, str], Dict[str, str]]:
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 3:
                parent[parts[0]] = parts[1]
                rank[parts[0]] = parts[2]
    return parent, rank


def load_names(path: Path) -> Dict[str, str]:
    names: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[parts[0]] = parts[1]
    return names


def species_taxid(taxid: str, parent: Dict[str, str], rank: Dict[str, str]) -> str:
    cur = str(taxid)
    seen = set()
    while cur and cur not in seen:
        seen.add(cur)
        if rank.get(cur) == "species":
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
    if values and values[0] == "131567":
        values = values[1:]
    return "|".join(values), "|".join(names.get(t, t) for t in values)


def extract_acc(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-genomes", type=Path, required=True)
    ap.add_argument("--skani", type=Path, required=True)
    ap.add_argument("--base-taxmap", type=Path, required=True)
    ap.add_argument("--nodes", type=Path, required=True)
    ap.add_argument("--names", type=Path, required=True)
    ap.add_argument("--min-ani", type=float, default=95.0)
    ap.add_argument("--min-af", type=float, default=50.0)
    ap.add_argument("--af-mode", choices=["min", "max"], default="min")
    ap.add_argument("--out-overrides", type=Path, required=True)
    ap.add_argument("--out-taxmap", type=Path, required=True)
    ap.add_argument("--out-summary", type=Path, required=True)
    args = ap.parse_args(argv)

    parent, rank = load_nodes(args.nodes)
    names = load_names(args.names)

    source = pd.read_csv(args.source_genomes, sep="\t", keep_default_na=False)
    source["source_species_taxid"] = [
        species_taxid(str(t), parent, rank) for t in source["source_taxid"].astype(str)
    ]
    source["source_taxpath"], source["source_taxpathsn"] = zip(
        *[lineage_path(t, parent, names) if t else ("", "") for t in source["source_species_taxid"]]
    )

    by_acc = source.set_index("source_accession", drop=False).to_dict("index")

    skani = pd.read_csv(args.skani, sep="\t")
    for col in ["ANI", "Align_fraction_ref", "Align_fraction_query"]:
        skani[col] = pd.to_numeric(skani[col], errors="coerce").fillna(0.0)
    skani["called_accession"] = skani["Query_file"].map(extract_acc)
    skani["source_accession"] = skani["Ref_file"].map(extract_acc)
    skani["af_gate"] = (
        skani[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
        if args.af_mode == "min"
        else skani[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
    )
    skani = skani.loc[skani["called_accession"].astype(bool) & skani["source_accession"].astype(bool)].copy()
    best = (
        skani.sort_values(
            ["called_accession", "ANI", "af_gate", "Align_fraction_query", "Align_fraction_ref"],
            ascending=[True, False, False, False, False],
        )
        .groupby("called_accession", as_index=False)
        .head(1)
        .copy()
    )
    best["passes"] = (best["ANI"] >= args.min_ani) & (best["af_gate"] >= args.min_af)

    rows = []
    for _, rec in best.loc[best["passes"]].iterrows():
        src = by_acc.get(str(rec["source_accession"]))
        if not src:
            continue
        taxid = str(src.get("source_species_taxid", ""))
        taxpath = str(src.get("source_taxpath", ""))
        taxpathsn = str(src.get("source_taxpathsn", ""))
        if not taxid or not taxpathsn:
            continue
        rows.append(
            {
                "ref_key": rec["called_accession"],
                "taxid": taxid,
                "rank": "species",
                "taxpath": taxpath,
                "taxpathsn": taxpathsn,
                "_CAMI_genomeID": src.get("genome_id", ""),
                "_CAMI_OTU": "source_aware",
                "source_accession": rec["source_accession"],
                "source_taxid": src.get("source_taxid", ""),
                "ANI_to_source": float(rec["ANI"]),
                "Align_fraction_ref": float(rec["Align_fraction_ref"]),
                "Align_fraction_query": float(rec["Align_fraction_query"]),
                "af_gate": float(rec["af_gate"]),
            }
        )

    overrides = pd.DataFrame(rows).sort_values(["ref_key", "ANI_to_source"], ascending=[True, False])
    args.out_overrides.parent.mkdir(parents=True, exist_ok=True)
    args.out_taxmap.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    overrides.to_csv(args.out_overrides, sep="\t", index=False)

    with args.out_taxmap.open("w") as out:
        for _, rec in overrides.iterrows():
            out.write(
                f"{rec['ref_key']}\t{rec['taxid']}\tspecies\t{rec['taxpath']}\t{rec['taxpathsn']}"
                f"\t{rec['_CAMI_genomeID']}\t{rec['_CAMI_OTU']}\n"
            )
        with args.base_taxmap.open() as fh:
            for line in fh:
                out.write(line)

    with args.out_summary.open("w") as out:
        out.write("metric\tvalue\n")
        out.write(f"selected_called_with_skani_hit\t{best['called_accession'].nunique()}\n")
        out.write(f"override_count\t{len(overrides)}\n")
        out.write(f"min_ani\t{args.min_ani}\n")
        out.write(f"min_af\t{args.min_af}\n")
        out.write(f"af_mode\t{args.af_mode}\n")
        out.write(f"best_hit_ani_ge_95\t{int((best['ANI'] >= 95.0).sum())}\n")
        out.write(f"best_hit_ani_ge_94\t{int((best['ANI'] >= 94.0).sum())}\n")

    print(f"source-aware overrides: {len(overrides)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
