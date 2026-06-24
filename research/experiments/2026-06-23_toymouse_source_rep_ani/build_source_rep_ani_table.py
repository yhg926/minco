#!/usr/bin/env python3
"""Toy Mouse sample0 source genome vs GTDB representative ANI table."""

from __future__ import annotations

import csv
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
TMP_DIR = Path("/tmp/toymouse_source_rep_anim_20260623")
FASTA_DIR = TMP_DIR / "fasta"
DELTA_DIR = TMP_DIR / "delta"

TRUTH_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
ABUND_EXP = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"
C15_EXP = ROOT / "research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse"
SOURCE_TRUTH = TRUTH_EXP / "mouse0_gtdb_source_species_ground_truth.tsv"
MINCO_PROFILE = Path(
    os.environ.get(
        "MINCO_PROFILE",
        str(OUT_DIR / "toymouse_sample0_ctxmarker_rawani_report.tsv"),
    )
)
MINCO_PROFILE_LABEL = os.environ.get("MINCO_PROFILE_LABEL", "ctx-marker raw-ANI-report")
MINCO_CODEN_LEN = int(os.environ.get("MINCO_CODEN_LEN", "11"))
MINCO_STORAGE = os.environ.get("MINCO_STORAGE", "ctxobj64")
MINCO_REFERENCE_KIND = os.environ.get("MINCO_REFERENCE_KIND", "coden11 s2000 ctx-marker")
OUT_TABLE = os.environ.get("OUT_TABLE", "toymouse_sample0_source_rep_ani.tsv")
OUT_SUMMARY = os.environ.get("OUT_SUMMARY", "summary.tsv")
SYLPH_PROFILE = Path("/mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/profile.tsv")
REF_PATH_LISTS = [
    C15_EXP / "gtdb232_199924_available_abs_paths.list",
    C15_EXP / "gtdb232_200709_paths.resolved.list",
]
GTDB_METADATA = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]

ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
EPSILON = 1e-8


def extract_accession(text: object) -> str:
    match = ACC_RE.search(str(text))
    return match.group(1) if match else ""


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def load_representatives() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for path in GTDB_METADATA:
        with open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                acc = extract_accession(row.get("accession", ""))
                genbank = extract_accession(row.get("ncbi_genbank_assembly_accession", ""))
                rec = {
                    "metadata_accession": row.get("accession", ""),
                    "accession": acc,
                    "genbank_accession": genbank,
                    "gtdb_genome_representative": row.get("gtdb_genome_representative", ""),
                    "gtdb_representative": row.get("gtdb_representative", ""),
                    "gtdb_taxonomy": row.get("gtdb_taxonomy", ""),
                }
                for key in [acc, genbank, row.get("accession", "")]:
                    if key:
                        out.setdefault(str(key), rec)
    return out


def load_ref_paths() -> dict[str, str]:
    out: dict[str, str] = {}
    for path in REF_PATH_LISTS:
        if not path.exists():
            continue
        with path.open() as fh:
            for line in fh:
                fasta = line.strip()
                acc = extract_accession(fasta)
                if acc and acc not in out:
                    out[acc] = fasta
    return out


def copy_fasta(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size > 0:
        return
    if src.suffix == ".gz":
        with gzip.open(src, "rt") as inp, dst.open("w") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(src, dst)


def fasta_length(path: Path) -> int:
    total = 0
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def parse_filtered_delta(path: Path) -> tuple[int, int]:
    aln_len = 0
    sim_errors = 0
    with path.open() as fh:
        for raw in fh:
            fields = raw.strip().split()
            if not fields or fields[0] == "NUCMER" or fields[0].startswith(">"):
                continue
            if len(fields) == 7:
                aln_len += abs(int(fields[1]) - int(fields[0])) + 1
                sim_errors += int(fields[4])
    return aln_len, sim_errors


def anim_pair(source_path: Path, ref_path: Path, source_acc: str, rep_acc: str) -> dict[str, float]:
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    DELTA_DIR.mkdir(parents=True, exist_ok=True)
    src_fa = FASTA_DIR / f"src_{source_acc}.fa"
    ref_fa = FASTA_DIR / f"rep_{rep_acc}.fa"
    copy_fasta(source_path, src_fa)
    copy_fasta(ref_path, ref_fa)

    prefix = DELTA_DIR / f"src_{source_acc}_vs_rep_{rep_acc}"
    delta = Path(str(prefix) + ".delta")
    filt = Path(str(prefix) + ".filter")
    if not filt.exists() or filt.stat().st_size == 0:
        subprocess.run(
            ["nucmer", "--mum", "-p", str(prefix), str(src_fa), str(ref_fa)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        with filt.open("w") as out:
            subprocess.run(["delta-filter", "-1", str(delta)], check=True, stdout=out)

    aln_len, sim_errors = parse_filtered_delta(filt)
    source_len = fasta_length(src_fa)
    ref_len = fasta_length(ref_fa)
    ani = 1.0 - (sim_errors / aln_len) if aln_len else math.nan
    return {
        "ANIm_ANI": ani,
        "ANIm_aligned_bp": aln_len,
        "ANIm_similarity_errors": sim_errors,
        "ANIm_source_aligned_fraction": aln_len / source_len if source_len else math.nan,
        "ANIm_ref_aligned_fraction": aln_len / ref_len if ref_len else math.nan,
        "source_length": source_len,
        "rep_length": ref_len,
    }


def minco_naive_ani(row: pd.Series) -> float:
    xny = safe_float(row.get("XnY_ctx"))
    ndiff = safe_float(row.get("N_diff_obj"))
    nsection = safe_float(row.get("N_diff_obj_section"))
    if not math.isfinite(xny) or xny <= 0.0:
        return math.nan
    ratio = (nsection + EPSILON) / (ndiff + EPSILON)
    dist0 = ndiff / (xny + ndiff) if xny + ndiff > 0.0 else 0.0
    final_dist = 1.0 - math.pow(1.0 - dist0, ratio)
    return 1.0 - final_dist * 0.1544286 if final_dist > 0.0 else 1.0


def load_minco_by_acc() -> dict[str, pd.Series]:
    rows = pd.read_csv(MINCO_PROFILE, sep="\t")
    rows["rep_accession"] = rows["Ref"].map(extract_accession)
    if "Ref_annotation" in rows:
        missing = rows["rep_accession"] == ""
        rows.loc[missing, "rep_accession"] = rows.loc[missing, "Ref_annotation"].map(
            extract_accession
        )
    rows["minco_naive_ANI_calc"] = [minco_naive_ani(row) for _, row in rows.iterrows()]
    rows = rows.sort_values(["rep_accession", "XnY_ctx"], ascending=[True, False])
    return {
        str(row.rep_accession): row
        for row in rows.itertuples(index=False)
        if str(row.rep_accession)
    }


def load_sylph_by_acc() -> dict[str, pd.Series]:
    rows = pd.read_csv(SYLPH_PROFILE, sep="\t")
    rows["rep_accession"] = rows["Genome_file"].map(extract_accession)
    rows = rows.sort_values(["rep_accession", "Eff_cov"], ascending=[True, False])
    return {
        str(row.rep_accession): row
        for row in rows.itertuples(index=False)
        if str(row.rep_accession)
    }


def series_get(row: object, key: str) -> object:
    if row is None:
        return ""
    if hasattr(row, "_asdict"):
        return row._asdict().get(key, "")
    return getattr(row, key, "")


def main() -> int:
    reps = load_representatives()
    ref_paths = load_ref_paths()
    minco = load_minco_by_acc()
    sylph = load_sylph_by_acc()
    sources = pd.read_csv(SOURCE_TRUTH, sep="\t", keep_default_na=False)

    rows = []
    for source in sources.itertuples(index=False):
        source_acc = str(source.source_accession)
        if not source_acc:
            continue
        rec = reps.get(source_acc) or reps.get(str(source.metadata_accession))
        rep_meta = rec.get("gtdb_genome_representative", "") if rec else ""
        rep_acc = extract_accession(rep_meta)
        rep_path = ref_paths.get(rep_acc, "")
        source_path = str(source.local_path)

        out = {
            "genome_id": source.genome_id,
            "source_accession": source_acc,
            "source_basename": source.source_basename,
            "source_abundance": safe_float(source.abundance),
            "source_taxid": source.source_taxid,
            "gtdb_species": source.gtdb_species,
            "gtdb_representative_meta": rep_meta,
            "gtdb_representative_accession": rep_acc,
            "source_path": source_path,
            "representative_path": rep_path,
            "source_is_representative": source_acc == rep_acc,
            "minco_profile_label": MINCO_PROFILE_LABEL,
            "minco_coden_len": MINCO_CODEN_LEN,
            "minco_storage": MINCO_STORAGE,
            "minco_reference_kind": MINCO_REFERENCE_KIND,
        }

        if source_path and rep_path:
            out.update(anim_pair(Path(source_path), Path(rep_path), source_acc, rep_acc))
        else:
            out.update(
                {
                    "ANIm_ANI": math.nan,
                    "ANIm_aligned_bp": math.nan,
                    "ANIm_similarity_errors": math.nan,
                    "ANIm_source_aligned_fraction": math.nan,
                    "ANIm_ref_aligned_fraction": math.nan,
                    "source_length": math.nan,
                    "rep_length": math.nan,
                }
            )

        minco_row = minco.get(rep_acc)
        out.update(
            {
                "minco_ref_present": minco_row is not None,
                "minco_ANI_col": safe_float(series_get(minco_row, "ANI")),
                "minco_naive_ANI_calc": safe_float(series_get(minco_row, "minco_naive_ANI_calc")),
                "minco_aafANI_Ref_zip_aaf_ani": safe_float(series_get(minco_row, "Ref_zip_aaf_ani")),
                "minco_XnY_ctx": safe_float(series_get(minco_row, "XnY_ctx")),
                "minco_Ref_breadth": safe_float(series_get(minco_row, "Ref_breadth")),
            }
        )

        sylph_row = sylph.get(rep_acc)
        out.update(
            {
                "sylph_ref_present": sylph_row is not None,
                "sylph_Adjusted_ANI": safe_float(series_get(sylph_row, "Adjusted_ANI")) / 100.0
                if sylph_row is not None
                else math.nan,
                "sylph_Naive_ANI": safe_float(series_get(sylph_row, "Naive_ANI")) / 100.0
                if sylph_row is not None
                else math.nan,
                "sylph_Eff_cov": safe_float(series_get(sylph_row, "Eff_cov")),
                "sylph_Taxonomic_abundance": safe_float(
                    series_get(sylph_row, "Taxonomic_abundance")
                ),
            }
        )
        rows.append(out)

    out_df = pd.DataFrame(rows).sort_values("source_abundance", ascending=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_DIR / OUT_TABLE, sep="\t", index=False)

    estimator_cols = [
        "minco_ANI_col",
        "minco_naive_ANI_calc",
        "minco_aafANI_Ref_zip_aaf_ani",
        "sylph_Adjusted_ANI",
        "sylph_Naive_ANI",
    ]
    metric_rows = []
    common_mask = out_df["ANIm_ANI"].notna()
    for col in estimator_cols:
        common_mask &= out_df[col].notna()

    for label, base_mask in [
        ("all_available", out_df["ANIm_ANI"].notna()),
        ("common_minco_sylph", common_mask),
    ]:
        for col in estimator_cols:
            sub = out_df.loc[base_mask & out_df[col].notna()].copy()
            err = sub[col] - sub["ANIm_ANI"]
            metric_rows.append(
                {
                    "comparison_set": label,
                    "estimator": col,
                    "n": len(sub),
                    "pearson": sub[col].corr(sub["ANIm_ANI"], method="pearson"),
                    "spearman": sub[col].corr(sub["ANIm_ANI"], method="spearman"),
                    "mae": err.abs().mean(),
                    "mean_error": err.mean(),
                    "median_abs_error": err.abs().median(),
                    "estimator_mean": sub[col].mean(),
                    "ANIm_mean": sub["ANIm_ANI"].mean(),
                }
            )
    metrics = pd.DataFrame(metric_rows).sort_values(["comparison_set", "mae"])
    metrics.to_csv(OUT_DIR / OUT_SUMMARY, sep="\t", index=False)

    print(metrics.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
