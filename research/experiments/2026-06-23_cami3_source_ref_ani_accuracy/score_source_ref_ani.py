#!/usr/bin/env python3
"""Compare read-based ANI against CAMI3 source-genome-to-reference ANI."""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
import subprocess
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
TMP_DIR = Path("/tmp/cami3_source_ref_ani_20260623")
SOURCE_EXTRACT_DIR = TMP_DIR / "source_extract"
SKANI_TMP_DIR = TMP_DIR / "skani"

SOURCE_TAR = Path("/mnt/new3T/minco_cami3_toygut_20260620/source_genomes.tar.gz")
SOURCE_TARGETS = Path("/mnt/new3T/minco_cami3_toygut_20260620/source_targets/source_genomes")

READ_MAPPINGS = {
    0: Path("/mnt/new3T/minco_cami3_toygut_20260620/sample_0_reads_mapping.tsv.gz"),
    1: Path("/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_1_reads/reads_mapping.tsv.gz"),
    2: Path("/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_2_reads/reads_mapping.tsv.gz"),
}

SELECTED_INPUTS = [
    (
        ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue/selected_species.tsv",
        {
            "sylph": "sylph_adjusted",
            "minco_ctxmarker_active_robust_depth_intragenus_rescue": "minco_old_ctxmarker",
        },
    ),
    (
        ROOT
        / "research/experiments/2026-06-23_cami3_toygut_coden15_formula_af/coden15_cami3_selected_species.tsv",
        {"coden15_formula_af_gate_robust_rescue": "minco_coden15_formula_af"},
    ),
]

REF_DIRS = [
    Path("/mnt/new3T/gtdbr220/GTDBr226_genomes"),
    Path("/mnt/new3T/gtdbr220/gtdb232/r232add_genomes"),
]
REF_PATH_LISTS = [
    ROOT
    / "research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse/gtdb232_199924_available_abs_paths.list",
    ROOT
    / "research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse/gtdb232_200709_paths.resolved.list",
]

ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def extract_accession(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def selected_calls() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path, method_map in SELECTED_INPUTS:
        with path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                method = method_map.get(row.get("method", ""))
                if method is None:
                    continue
                sample = int(float(row["sample_id"]))
                accession = extract_accession(row.get("accession", ""))
                taxid = str(row.get("ncbi_species_taxid", "") or "")
                pred_ani = safe_float(row.get("pred_ani"))
                if not accession or not taxid or not math.isfinite(pred_ani):
                    continue
                rows.append(
                    {
                        "sample_id": sample,
                        "method": method,
                        "source_selected_file": str(path),
                        "ncbi_species_taxid": taxid,
                        "ncbi_species": row.get("ncbi_species", ""),
                        "gtdb_species": row.get("gtdb_species", ""),
                        "ref_accession": accession,
                        "read_ani": pred_ani,
                        "support": safe_float(row.get("support")),
                        "pred_abundance_raw": safe_float(row.get("pred_abundance_raw")),
                    }
                )
    out = pd.DataFrame(rows)
    out = out.drop_duplicates(
        ["sample_id", "method", "ncbi_species_taxid", "gtdb_species", "ref_accession"]
    )
    return out.sort_values(["sample_id", "method", "ncbi_species_taxid", "ref_accession"])


def resolve_ref_paths(accessions: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    missing = set(accessions)
    for list_path in REF_PATH_LISTS:
        if not list_path.exists() or not missing:
            continue
        with list_path.open() as fh:
            for line in fh:
                if not missing:
                    break
                path = line.strip()
                acc = extract_accession(path)
                if acc in missing:
                    out[acc] = path
                    missing.discard(acc)
    for ref_dir in REF_DIRS:
        if not ref_dir.exists() or not missing:
            continue
        with os.scandir(ref_dir) as entries:
            for entry in entries:
                if not missing:
                    break
                acc = extract_accession(entry.name)
                if acc in missing and entry.is_file():
                    out[acc] = entry.path
                    missing.discard(acc)
    if missing:
        with (OUT_DIR / "missing_ref_accessions.txt").open("w") as fh:
            for acc in sorted(missing):
                fh.write(f"{acc}\n")
    return out


def contig_from_read_id(read_id: str) -> str:
    read_id = read_id.strip()
    if "-" not in read_id:
        return read_id.split("/", 1)[0]
    return read_id.rsplit("-", 1)[0]


def source_counts_and_contigs(selected: pd.DataFrame):
    wanted_by_sample = {
        sample: set(sub["ncbi_species_taxid"].astype(str))
        for sample, sub in selected.groupby("sample_id")
    }
    counts: dict[tuple[int, str, str], int] = Counter()
    contigs_by_genome: dict[str, set[str]] = defaultdict(set)

    for sample, path in READ_MAPPINGS.items():
        wanted = wanted_by_sample.get(sample, set())
        if not wanted:
            continue
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                _read_name, genome_id, taxid, read_id = fields[:4]
                if taxid not in wanted:
                    continue
                counts[(sample, taxid, genome_id)] += 1
                if len(contigs_by_genome[genome_id]) < 500:
                    contigs_by_genome[genome_id].add(contig_from_read_id(read_id))

    count_rows = [
        {
            "sample_id": sample,
            "ncbi_species_taxid": taxid,
            "source_genome_id": genome_id,
            "read_rows": count,
            "distinct_contigs_sampled": len(contigs_by_genome.get(genome_id, set())),
        }
        for (sample, taxid, genome_id), count in counts.items()
    ]
    count_df = pd.DataFrame(count_rows).sort_values(
        ["sample_id", "ncbi_species_taxid", "read_rows"], ascending=[True, True, False]
    )
    count_df.to_csv(OUT_DIR / "source_genome_read_counts.tsv", sep="\t", index=False)
    return counts, contigs_by_genome


def fasta_headers(path: Path):
    with path.open("rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip().split()[0]


def scan_source_targets(
    unresolved: set[str],
    contig_to_genomes: dict[str, set[str]],
) -> dict[str, str]:
    out: dict[str, str] = {}
    if not SOURCE_TARGETS.exists():
        return out
    for path in SOURCE_TARGETS.glob("*.fna"):
        basename = path.name
        for genome_id in list(unresolved):
            if basename.startswith(f"{genome_id}__") or basename.startswith(f"{genome_id}_"):
                out[genome_id] = str(path)
                unresolved.discard(genome_id)
        if not unresolved:
            break
        for header in fasta_headers(path):
            for genome_id in contig_to_genomes.get(header, set()):
                if genome_id in unresolved:
                    out[genome_id] = str(path)
                    unresolved.discard(genome_id)
            if not unresolved:
                break
    return out


def scan_source_tar_members(
    unresolved: set[str],
    contig_to_genomes: dict[str, set[str]],
) -> dict[str, str]:
    out: dict[str, str] = {}
    if not unresolved:
        return out
    with tarfile.open(SOURCE_TAR, "r|gz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith((".fa", ".fna", ".fasta")):
                continue
            basename = Path(member.name).name
            direct = [
                genome_id
                for genome_id in unresolved
                if basename.startswith(f"{genome_id}__") or basename.startswith(f"{genome_id}_")
            ]
            if direct:
                for genome_id in direct:
                    out[genome_id] = member.name
                    unresolved.discard(genome_id)
                if not unresolved:
                    break
                continue
            fh = tar.extractfile(member)
            if fh is None:
                continue
            matched: set[str] = set()
            for raw in fh:
                if not raw.startswith(b">"):
                    continue
                header = raw[1:].decode("utf-8", "replace").strip().split()[0]
                for genome_id in contig_to_genomes.get(header, set()):
                    if genome_id in unresolved:
                        matched.add(genome_id)
                if matched:
                    for genome_id in matched:
                        out[genome_id] = member.name
                        unresolved.discard(genome_id)
                    break
            if not unresolved:
                break
    return out


def extract_members(member_by_genome: dict[str, str]) -> dict[str, str]:
    SOURCE_EXTRACT_DIR.mkdir(parents=True, exist_ok=True)
    members = sorted({m for m in member_by_genome.values() if not m.startswith("/")})
    if members:
        subprocess.run(
            ["tar", "-xzf", str(SOURCE_TAR), "-C", str(SOURCE_EXTRACT_DIR), *members],
            check=True,
        )
    out: dict[str, str] = {}
    for genome_id, member_or_path in member_by_genome.items():
        if member_or_path.startswith("/"):
            out[genome_id] = member_or_path
        else:
            out[genome_id] = str(SOURCE_EXTRACT_DIR / member_or_path)
    return out


def resolve_source_paths(contigs_by_genome: dict[str, set[str]]) -> dict[str, str]:
    cache = OUT_DIR / "source_genome_paths.tsv"
    if cache.exists():
        rows = pd.read_csv(cache, sep="\t")
        return dict(zip(rows["source_genome_id"].astype(str), rows["source_path"].astype(str)))

    wanted = set(contigs_by_genome)
    contig_to_genomes: dict[str, set[str]] = defaultdict(set)
    for genome_id, contigs in contigs_by_genome.items():
        for contig in contigs:
            contig_to_genomes[contig].add(genome_id)

    unresolved = set(wanted)
    member_or_path: dict[str, str] = {}
    member_or_path.update(scan_source_targets(unresolved, contig_to_genomes))
    member_or_path.update(scan_source_tar_members(unresolved, contig_to_genomes))
    source_paths = extract_members(member_or_path)

    rows = []
    for genome_id in sorted(wanted):
        rows.append(
            {
                "source_genome_id": genome_id,
                "source_path": source_paths.get(genome_id, ""),
                "source_resolved": bool(source_paths.get(genome_id, "")),
                "sampled_contigs": ",".join(sorted(contigs_by_genome[genome_id])[:50]),
                "sampled_contig_count": len(contigs_by_genome[genome_id]),
            }
        )
    pd.DataFrame(rows).to_csv(cache, sep="\t", index=False)
    return source_paths


def run_skani(pairs: set[tuple[str, str]]) -> pd.DataFrame:
    cache = OUT_DIR / "source_ref_skani_pairs.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    SKANI_TMP_DIR.mkdir(parents=True, exist_ok=True)
    by_source: dict[str, set[str]] = defaultdict(set)
    for source, ref in pairs:
        by_source[source].add(ref)

    rows: list[dict[str, object]] = []
    for idx, (source, refs) in enumerate(sorted(by_source.items()), start=1):
        if not Path(source).exists():
            continue
        refs = sorted(ref for ref in refs if Path(ref).exists())
        if not refs:
            continue
        out_path = SKANI_TMP_DIR / f"skani_{idx:05d}.tsv"
        cmd = [
            "skani",
            "dist",
            str(source),
            *refs,
            "--min-af",
            "0",
            "-t",
            "4",
            "-o",
            str(out_path),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if not out_path.exists() or out_path.stat().st_size == 0:
            continue
        part = pd.read_csv(out_path, sep="\t")
        for _, rec in part.iterrows():
            rows.append(
                {
                    "source_path": str(rec["Query_file"]),
                    "ref_path": str(rec["Ref_file"]),
                    "source_ref_ani": safe_float(rec["ANI"]) / 100.0,
                    "skani_ref_af": safe_float(rec["Align_fraction_ref"]) / 100.0,
                    "skani_query_af": safe_float(rec["Align_fraction_query"]) / 100.0,
                    "skani_ref_name": rec.get("Ref_name", ""),
                    "skani_query_name": rec.get("Query_name", ""),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(cache, sep="\t", index=False)
    return out


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (method, source_choice), sub in detail.groupby(["method", "source_choice"]):
        sub = sub.loc[sub["source_ref_ani"].notna() & sub["read_ani"].notna()].copy()
        if sub.empty:
            continue
        err = sub["read_ani"].astype(float) - sub["source_ref_ani"].astype(float)
        rows.append(
            {
                "method": method,
                "source_choice": source_choice,
                "n_calls": len(sub),
                "pearson": sub["read_ani"].corr(sub["source_ref_ani"], method="pearson"),
                "spearman": sub["read_ani"].corr(sub["source_ref_ani"], method="spearman"),
                "mae": err.abs().mean(),
                "mean_error": err.mean(),
                "median_abs_error": err.abs().median(),
                "p95_abs_error": err.abs().quantile(0.95),
                "read_ani_mean": sub["read_ani"].mean(),
                "source_ref_ani_mean": sub["source_ref_ani"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(["source_choice", "mae", "method"])


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    selected = selected_calls()
    selected.to_csv(OUT_DIR / "selected_calls_input.tsv", sep="\t", index=False)

    ref_paths = resolve_ref_paths(set(selected["ref_accession"].astype(str)))
    selected["ref_path"] = selected["ref_accession"].astype(str).map(ref_paths).fillna("")

    counts, contigs_by_genome = source_counts_and_contigs(selected)
    source_paths = resolve_source_paths(contigs_by_genome)

    count_lookup = defaultdict(list)
    for (sample, taxid, genome_id), read_rows in counts.items():
        count_lookup[(sample, taxid)].append((genome_id, read_rows, source_paths.get(genome_id, "")))

    candidate_rows = []
    pairs: set[tuple[str, str]] = set()
    for call_id, call in selected.reset_index(drop=True).iterrows():
        sample = int(call["sample_id"])
        taxid = str(call["ncbi_species_taxid"])
        ref_path = str(call["ref_path"])
        candidates = count_lookup.get((sample, taxid), [])
        for genome_id, read_rows, source_path in candidates:
            row = {
                "call_id": call_id,
                "sample_id": sample,
                "method": call["method"],
                "ncbi_species_taxid": taxid,
                "ncbi_species": call["ncbi_species"],
                "gtdb_species": call["gtdb_species"],
                "ref_accession": call["ref_accession"],
                "ref_path": ref_path,
                "read_ani": call["read_ani"],
                "support": call["support"],
                "source_genome_id": genome_id,
                "source_read_rows": read_rows,
                "source_path": source_path,
                "source_path_resolved": bool(source_path),
                "ref_path_resolved": bool(ref_path),
            }
            candidate_rows.append(row)
            if source_path and ref_path:
                pairs.add((source_path, ref_path))
    candidates = pd.DataFrame(candidate_rows)
    candidates.to_csv(OUT_DIR / "call_source_candidates.tsv", sep="\t", index=False)

    skani = run_skani(pairs)
    skani_lookup = {
        (str(row["source_path"]), str(row["ref_path"])): row for _, row in skani.iterrows()
    }

    annotated = []
    for _, row in candidates.iterrows():
        rec = dict(row)
        hit = skani_lookup.get((str(row["source_path"]), str(row["ref_path"])))
        if hit is None:
            rec.update(
                {
                    "source_ref_ani": math.nan,
                    "skani_ref_af": math.nan,
                    "skani_query_af": math.nan,
                }
            )
        else:
            rec.update(
                {
                    "source_ref_ani": safe_float(hit["source_ref_ani"]),
                    "skani_ref_af": safe_float(hit["skani_ref_af"]),
                    "skani_query_af": safe_float(hit["skani_query_af"]),
                }
            )
        annotated.append(rec)
    candidate_detail = pd.DataFrame(annotated)
    candidate_detail.to_csv(OUT_DIR / "call_source_candidates_with_skani.tsv", sep="\t", index=False)

    chosen_rows = []
    for call_id, sub in candidate_detail.groupby("call_id"):
        sub = sub.loc[sub["source_ref_ani"].notna()].copy()
        if sub.empty:
            continue
        major = sub.sort_values(["source_read_rows", "source_ref_ani"], ascending=False).iloc[0]
        best = sub.sort_values(["source_ref_ani", "source_read_rows"], ascending=False).iloc[0]
        for choice, rec in [("read_major_source", major), ("best_ani_source", best)]:
            out = rec.to_dict()
            out["source_choice"] = choice
            out["ani_error"] = safe_float(out["read_ani"]) - safe_float(out["source_ref_ani"])
            out["ani_abs_error"] = abs(out["ani_error"])
            chosen_rows.append(out)
    detail = pd.DataFrame(chosen_rows)
    detail.to_csv(OUT_DIR / "read_vs_source_ref_ani_detail.tsv", sep="\t", index=False)

    summary = summarize(detail)
    summary.to_csv(OUT_DIR / "summary.tsv", sep="\t", index=False)
    print(summary.to_csv(sep="\t", index=False), end="")

    diagnostics = {
        "selected_calls": len(selected),
        "selected_calls_with_ref_path": int((selected["ref_path"].astype(str) != "").sum()),
        "source_genomes_needed": len(contigs_by_genome),
        "source_genomes_resolved": sum(1 for v in source_paths.values() if v),
        "source_ref_pairs_requested": len(pairs),
        "source_ref_pairs_scored": len(skani),
        "calls_with_any_scored_source": detail["call_id"].nunique() if not detail.empty else 0,
    }
    pd.DataFrame([diagnostics]).to_csv(OUT_DIR / "diagnostics.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
