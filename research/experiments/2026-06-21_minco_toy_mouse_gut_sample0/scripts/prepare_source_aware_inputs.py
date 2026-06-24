#!/usr/bin/env python3
"""Prepare source-aware CAMISIM benchmark inputs for one sample."""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path
from typing import Sequence

import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def extract_acc(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--setup-dir", type=Path, required=True)
    ap.add_argument("--unique-minco", type=Path, required=True)
    ap.add_argument("--split-minco", type=Path, required=True)
    ap.add_argument("--sylph-profile", type=Path, required=True)
    ap.add_argument("--gtdb-genome-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    dist = {}
    positive = []
    with (args.setup_dir / "distributions/distribution_0.txt").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            genome_id, abundance_s = line.rstrip("\n").split("\t")[:2]
            abundance = float(abundance_s)
            dist[genome_id] = abundance
            if abundance > 0.0:
                positive.append(genome_id)

    meta = {}
    with (args.setup_dir / "internal/meta_data.tsv").open() as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx_genome = header.index("genome_ID")
        idx_taxid = header.index("NCBI_ID")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            meta[parts[idx_genome]] = parts[idx_taxid]

    loc = {}
    with (args.setup_dir / "internal/genome_locations.tsv").open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                loc[parts[0]] = parts[1]

    source_rows = []
    with (args.outdir / "positive_source_tar_paths.list").open("w") as out:
        for genome_id in positive:
            basename = Path(loc[genome_id]).name
            acc = extract_acc(basename)
            out.write(f"source_genomes/{basename}\n")
            source_rows.append(
                {
                    "genome_id": genome_id,
                    "abundance": dist[genome_id],
                    "source_taxid": meta.get(genome_id, ""),
                    "source_accession": acc,
                    "source_basename": basename,
                    "tar_path": f"source_genomes/{basename}",
                    "local_path": str(args.outdir / "source_positive/source_genomes" / basename),
                }
            )
    pd.DataFrame(source_rows).to_csv(args.outdir / "positive_source_genomes.tsv", sep="\t", index=False)

    selected = set()
    ref_path_by_acc = {}
    for path in [args.unique_minco, args.split_minco]:
        df = pd.read_csv(path, sep="\t")
        for col in ["ANI", "XnY_ctx", "Real_min_align_fraction", "Ref_breadth"]:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        direct = (df["XnY_ctx"] >= 10) & (df["ANI"] >= 0.94) & (df["Real_min_align_fraction"] >= 0.05)
        relaxed = (df["XnY_ctx"] >= 1) & (df["ANI"] >= 0.90) & (df["Ref_breadth"] >= 0.005)
        for ref in df.loc[direct | relaxed, "Ref"].astype(str):
            acc = extract_acc(ref)
            if not acc:
                continue
            selected.add(acc)
            if Path(ref).exists():
                ref_path_by_acc.setdefault(acc, ref)

    sylph = pd.read_csv(args.sylph_profile, sep="\t")
    for ref in sylph["Genome_file"].astype(str):
        acc = extract_acc(ref)
        if acc:
            selected.add(acc)

    for acc in sorted(selected):
        if acc in ref_path_by_acc:
            continue
        patterns = [
            str(args.gtdb_genome_dir / f"{acc}*.fna.gz"),
            str(args.gtdb_genome_dir / f"{acc}*.fa.gz"),
            str(args.gtdb_genome_dir / f"{acc}*.fna"),
            str(args.gtdb_genome_dir / f"{acc}*.fa"),
        ]
        hits = []
        for pattern in patterns:
            hits.extend(glob.glob(pattern))
        if hits:
            ref_path_by_acc[acc] = sorted(hits)[0]

    called_rows = [
        {"called_accession": acc, "called_ref_path": ref_path_by_acc.get(acc, "")}
        for acc in sorted(selected)
    ]
    pd.DataFrame(called_rows).to_csv(args.outdir / "selected_called_refs.tsv", sep="\t", index=False)
    with (args.outdir / "selected_called_ref_paths.list").open("w") as out:
        for row in called_rows:
            if row["called_ref_path"]:
                out.write(row["called_ref_path"] + "\n")

    missing = [row for row in called_rows if not row["called_ref_path"]]
    print(f"positive_sources\t{len(source_rows)}")
    print(f"selected_called_refs\t{len(called_rows)}")
    print(f"missing_called_ref_paths\t{len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
