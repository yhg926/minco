#!/usr/bin/env python3
"""Compute pyani-style ANIm for source-aware called-reference pairs."""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Sequence

import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def extract_acc(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def fasta_length(path: Path) -> int:
    total = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def copy_fasta(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".gz":
        with gzip.open(src, "rt") as inp, dst.open("w") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(src, dst)


def parse_filtered_delta(path: Path) -> tuple[int, int]:
    aln_length = 0
    sim_errors = 0
    with path.open() as fh:
        for raw in fh:
            fields = raw.strip().split()
            if not fields or fields[0] == "NUCMER" or fields[0].startswith(">"):
                continue
            if len(fields) == 7:
                aln_length += abs(int(fields[1]) - int(fields[0]))
                sim_errors += int(fields[4])
    return aln_length, sim_errors


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overrides", type=Path, required=True)
    ap.add_argument("--selected-refs", type=Path, required=True)
    ap.add_argument("--positive-sources", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--taxids", nargs="+", required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    fasta_dir = args.outdir / "fasta"
    delta_dir = args.outdir / "delta"
    fasta_dir.mkdir(exist_ok=True)
    delta_dir.mkdir(exist_ok=True)

    overrides = pd.read_csv(args.overrides, sep="\t", keep_default_na=False)
    selected_refs = pd.read_csv(args.selected_refs, sep="\t", keep_default_na=False)
    sources = pd.read_csv(args.positive_sources, sep="\t", keep_default_na=False)
    summary = pd.read_csv(args.summary, sep="\t", keep_default_na=False)

    ref_path_by_acc: Dict[str, Path] = {
        str(row.called_accession): Path(str(row.called_ref_path))
        for row in selected_refs.itertuples(index=False)
        if str(row.called_ref_path)
    }
    source_path_by_acc: Dict[str, Path] = {
        str(row.source_accession): Path(str(row.local_path))
        for row in sources.itertuples(index=False)
        if str(row.source_accession)
    }
    source_name_by_acc: Dict[str, str] = {
        str(row.source_accession): str(row.source_basename)
        for row in sources.itertuples(index=False)
        if str(row.source_accession)
    }

    taxid_set = {str(t) for t in args.taxids}
    pairs = overrides.loc[overrides["taxid"].astype(str).isin(taxid_set)].copy()
    if pairs.empty:
        raise SystemExit("no matching override pairs")

    summary["taxid"] = summary["taxid"].astype(str)
    summary["sylph_acc"] = summary["sylph_ref"].map(extract_acc)
    summary["minco_best_split_acc"] = summary["best_split_ref"].map(extract_acc)
    sylph_acc_by_taxid = {
        row.taxid: row.sylph_acc for row in summary.itertuples(index=False)
    }
    minco_best_by_taxid = {
        row.taxid: row.minco_best_split_acc for row in summary.itertuples(index=False)
    }

    rows = []
    length_cache: Dict[Path, int] = {}
    for pair in pairs.itertuples(index=False):
        taxid = str(pair.taxid)
        ref_acc = str(pair.ref_key)
        source_acc = str(pair.source_accession)
        source_path = source_path_by_acc.get(source_acc)
        ref_path = ref_path_by_acc.get(ref_acc)
        if source_path is None or ref_path is None:
            continue

        source_fa = fasta_dir / f"src_{source_acc}.fa"
        ref_fa = fasta_dir / f"ref_{ref_acc}.fa"
        if not source_fa.exists():
            copy_fasta(source_path, source_fa)
        if not ref_fa.exists():
            copy_fasta(ref_path, ref_fa)

        prefix = delta_dir / f"src_{source_acc}_vs_ref_{ref_acc}"
        delta = Path(str(prefix) + ".delta")
        filt = Path(str(prefix) + ".filter")
        if not filt.exists():
            subprocess.run(
                ["nucmer", "--mum", "-p", str(prefix), str(source_fa), str(ref_fa)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            with filt.open("w") as out:
                subprocess.run(["delta-filter", "-1", str(delta)], check=True, stdout=out)

        aln_len, sim_errors = parse_filtered_delta(filt)
        if source_fa not in length_cache:
            length_cache[source_fa] = fasta_length(source_fa)
        if ref_fa not in length_cache:
            length_cache[ref_fa] = fasta_length(ref_fa)
        source_len = length_cache[source_fa]
        ref_len = length_cache[ref_fa]
        anim_ani = 1.0 - (float(sim_errors) / float(aln_len)) if aln_len else 0.0
        role = []
        if ref_acc == sylph_acc_by_taxid.get(taxid):
            role.append("sylph_best")
        if ref_acc == minco_best_by_taxid.get(taxid):
            role.append("minco_best_split")
        rows.append(
            {
                "taxid": taxid,
                "species_name": str(pair.taxpathsn).split("|")[-1],
                "source_accession": source_acc,
                "source_basename": source_name_by_acc.get(source_acc, ""),
                "called_ref": ref_acc,
                "role": ",".join(role) if role else "other_sourceaware_ref",
                "ANIm_ANI": anim_ani,
                "ANIm_ANI_percent": anim_ani * 100.0,
                "ANIm_aligned_bp": aln_len,
                "ANIm_similarity_errors": sim_errors,
                "ANIm_source_aligned_fraction": aln_len / source_len if source_len else 0.0,
                "ANIm_ref_aligned_fraction": aln_len / ref_len if ref_len else 0.0,
                "source_length": source_len,
                "ref_length": ref_len,
                "skani_ANI_percent": float(pair.ANI_to_source),
                "skani_refAF_percent": float(pair.Align_fraction_ref),
                "skani_queryAF_percent": float(pair.Align_fraction_query),
            }
        )

    out = pd.DataFrame(rows).sort_values(["taxid", "role", "ANIm_ANI"], ascending=[True, True, False])
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
