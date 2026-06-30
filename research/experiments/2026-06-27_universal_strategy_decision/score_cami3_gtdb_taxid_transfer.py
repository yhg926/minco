#!/usr/bin/env python3
"""Score CAMI3 Toy Human Gut in a systematic GTDB-species transfer namespace.

The transfer is conservative: a CAMI/NCBI species taxid is converted to GTDB
only when GTDB metadata maps that taxid to exactly one GTDB species. Ambiguous
or unmapped truth mass is reported and excluded from the transferred truth set.
This is useful evidence, but not as strong as source-genome accession truth.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))
import build_and_score_gtdb_ground_truth as truth  # noqa: E402


CAMI3_ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620")
RUN = Path("/tmp/cami3_toy_human_gut_20260626/run")

SAMPLE_PATHS = {
    0: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_0.txt",
        "minco": RUN / "sample0_universal_autoexact.tsv",
        "sylph": CAMI3_ROOT / "sylph_sample0/profile.tsv",
    },
    1: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_1.txt",
        "minco": RUN / "sample1_universal_autoexact.tsv",
        "sylph": Path("/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv"),
    },
    2: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_2.txt",
        "minco": RUN / "sample2_universal_autoexact_tail_p025.tsv",
        "sylph": RUN / "sylph_sample2/profile.tsv",
    },
    3: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_3.txt",
        "minco": RUN / "sample3_universal_autoexact_tail_p025.tsv",
        "sylph": RUN / "sylph_sample3/profile.tsv",
    },
    4: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_4.txt",
        "minco": RUN / "sample4_universal_autoexact_guarded.tsv",
        "sylph": RUN / "sylph_sample4/profile.tsv",
    },
    5: {
        "truth": CAMI3_ROOT / "taxonomic_profiles/taxonomic_profile_5.txt",
        "minco": RUN / "sample5_universal_autoexact.tsv",
        "sylph": RUN / "sylph_sample5/profile.tsv",
    },
}


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").replace([math.inf, -math.inf], pd.NA).fillna(0.0)


def read_cami_species_truth(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("@") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            taxid, rank, taxpath, taxpathsn, pct = fields[:5]
            if rank != "species":
                continue
            if not (taxpath.startswith("131567|2") or taxpathsn.startswith("cellular organisms|Bacteria")):
                continue
            try:
                abundance_pct = float(pct)
            except ValueError:
                continue
            if abundance_pct <= 0:
                continue
            rows.append(
                {
                    "ncbi_species_taxid": str(taxid),
                    "truth_name": taxpathsn.split("|")[-1],
                    "abundance_pct_all": abundance_pct,
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["ncbi_species_taxid", "truth_name", "abundance_pct_all"])
    return (
        df.groupby("ncbi_species_taxid", as_index=False)
        .agg(truth_name=("truth_name", "first"), abundance_pct_all=("abundance_pct_all", "sum"))
    )


def build_taxid_transfer(by_accession: dict[str, dict[str, str]]) -> tuple[dict[str, str], dict[str, list[str]]]:
    species_by_taxid: dict[str, set[str]] = defaultdict(set)
    for rec in by_accession.values():
        gtdb_species = rec.get("gtdb_species", "")
        if not gtdb_species:
            continue
        for key in ["ncbi_species_taxid", "ncbi_taxid"]:
            taxid = str(rec.get(key, "")).strip()
            if taxid:
                species_by_taxid[taxid].add(gtdb_species)
    unique = {taxid: next(iter(species)) for taxid, species in species_by_taxid.items() if len(species) == 1}
    ambiguous = {taxid: sorted(species) for taxid, species in species_by_taxid.items() if len(species) > 1}
    return unique, ambiguous


def transfer_truth_to_gtdb(
    sample_id: int,
    truth_rows: pd.DataFrame,
    taxid_to_gtdb: dict[str, str],
    ambiguous_taxids: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    rows: list[dict[str, object]] = []
    ambiguous_mass = 0.0
    unmapped_mass = 0.0
    ambiguous_n = 0
    unmapped_n = 0
    for row in truth_rows.itertuples(index=False):
        taxid = str(row.ncbi_species_taxid)
        abundance = float(row.abundance_pct_all)
        gtdb_species = taxid_to_gtdb.get(taxid, "")
        if gtdb_species:
            rows.append(
                {
                    "sample": sample_id,
                    "ncbi_species_taxid": taxid,
                    "truth_name": row.truth_name,
                    "gtdb_species": gtdb_species,
                    "abundance_pct_all": abundance,
                }
            )
        elif taxid in ambiguous_taxids:
            ambiguous_n += 1
            ambiguous_mass += abundance
        else:
            unmapped_n += 1
            unmapped_mass += abundance
    mapped = pd.DataFrame(rows)
    if not mapped.empty:
        mapped = (
            mapped.groupby("gtdb_species", as_index=False)
            .agg(
                sample=("sample", "first"),
                abundance_pct_all=("abundance_pct_all", "sum"),
                source_taxids=("ncbi_species_taxid", lambda x: ",".join(sorted(set(map(str, x))))),
                truth_names=("truth_name", lambda x: ";".join(sorted(set(map(str, x))))),
            )
        )
        total_mapped = float(mapped["abundance_pct_all"].sum())
        mapped["truth_abundance"] = mapped["abundance_pct_all"] / total_mapped if total_mapped else 0.0
    else:
        mapped = pd.DataFrame(columns=["gtdb_species", "truth_abundance"])
        total_mapped = 0.0
    quality = {
        "sample": sample_id,
        "truth_taxids_total": len(truth_rows),
        "truth_taxids_mapped": int(len(rows)),
        "truth_taxids_ambiguous": ambiguous_n,
        "truth_taxids_unmapped": unmapped_n,
        "truth_mass_mapped_pct_all": total_mapped,
        "truth_mass_ambiguous_pct_all": ambiguous_mass,
        "truth_mass_unmapped_pct_all": unmapped_mass,
    }
    return mapped, quality


def collapse_prediction(rows: Iterable[tuple[str, float, float]]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=["gtdb_species", "abundance_raw", "ani"])
    if df.empty:
        return pd.DataFrame(columns=["gtdb_species", "pred_abundance", "pred_ani"])
    df = df.loc[df["gtdb_species"].astype(bool)].copy()
    if df.empty:
        return pd.DataFrame(columns=["gtdb_species", "pred_abundance", "pred_ani"])
    grouped = (
        df.groupby("gtdb_species", as_index=False)
        .agg(pred_abundance_raw=("abundance_raw", "sum"), pred_ani=("ani", "max"))
    )
    total = float(grouped["pred_abundance_raw"].sum())
    grouped["pred_abundance"] = grouped["pred_abundance_raw"] / total if total else 0.0
    return grouped


def load_minco_predictions(path: Path, taxid_to_gtdb: dict[str, str]) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t")
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[call].copy()
    rows = []
    unmapped = 0
    for row in selected.itertuples(index=False):
        taxid = str(getattr(row, "taxid", ""))
        gtdb_species = taxid_to_gtdb.get(taxid, "")
        if not gtdb_species:
            unmapped += 1
            continue
        abundance = float(getattr(row, "calibrated_abundance", 0.0) or 0.0)
        ani = float(getattr(row, "s_Ref_zip_aaf_ani_max", 0.0) or 0.0)
        rows.append((gtdb_species, abundance, ani))
    return collapse_prediction(rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
    }


def load_sylph_predictions(path: Path, by_accession, by_core) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t")
    rows = []
    unmapped = 0
    for row in raw.itertuples(index=False):
        accession = truth.extract_accession(str(getattr(row, "Genome_file", "")))
        rec, _method, _key = truth.lookup_accession(accession, by_accession, by_core)
        if rec is None or not rec.get("gtdb_species"):
            unmapped += 1
            continue
        abundance = float(getattr(row, "Taxonomic_abundance", 0.0) or 0.0) / 100.0
        ani = float(getattr(row, "Adjusted_ANI", 0.0) or 0.0) / 100.0
        rows.append((rec["gtdb_species"], abundance, ani))
    return collapse_prediction(rows), {
        "pred_rows_called": int(len(raw)),
        "pred_rows_unmapped": unmapped,
    }


def score_prediction(sample_id: int, method: str, pred: pd.DataFrame, truth_df: pd.DataFrame, extra: dict[str, object]) -> dict[str, object]:
    truth_species = set(truth_df["gtdb_species"].astype(str))
    pred_species = set(pred["gtdb_species"].astype(str))
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species
    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0

    truth_abund = dict(zip(truth_df["gtdb_species"].astype(str), numeric(truth_df["truth_abundance"])))
    pred_abund = dict(zip(pred["gtdb_species"].astype(str), numeric(pred["pred_abundance"])))
    union = sorted(truth_species | pred_species)
    y_true = pd.Series([truth_abund.get(species, 0.0) for species in union], dtype=float)
    y_pred = pd.Series([pred_abund.get(species, 0.0) for species in union], dtype=float)
    truth_order = sorted(truth_species)
    y_true_truth = pd.Series([truth_abund.get(species, 0.0) for species in truth_order], dtype=float)
    y_pred_truth = pd.Series([pred_abund.get(species, 0.0) for species in truth_order], dtype=float)

    return {
        "sample": sample_id,
        "method": method,
        "truth_gtdb_species": len(truth_species),
        "pred_gtdb_species": len(pred_species),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "L1_union_pp": float((y_pred - y_true).abs().sum() * 100.0),
        "L1_truth_only_pp": float((y_pred_truth - y_true_truth).abs().sum() * 100.0),
        "Pearson_union": y_pred.corr(y_true, method="pearson") if len(union) > 1 else float("nan"),
        "Pearson_truth_only": y_pred_truth.corr(y_true_truth, method="pearson") if len(truth_order) > 1 else float("nan"),
        **extra,
    }


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, sub in scores.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique()))),
                "mean_F1": sub["F1"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_L1_truth_only_pp": sub["L1_truth_only_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values("mean_F1", ascending=False)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    taxid_to_gtdb, ambiguous_taxids = build_taxid_transfer(by_accession)
    scores = []
    qualities = []
    truth_rows_out = []

    for sample_id, paths in SAMPLE_PATHS.items():
        for label, path in paths.items():
            if not Path(path).exists():
                raise SystemExit(f"missing {label} for sample {sample_id}: {path}")
        cami_truth = read_cami_species_truth(paths["truth"])
        gtdb_truth, quality = transfer_truth_to_gtdb(sample_id, cami_truth, taxid_to_gtdb, ambiguous_taxids)
        qualities.append(quality)
        truth_rows_out.extend(gtdb_truth.to_dict(orient="records"))

        minco_pred, minco_extra = load_minco_predictions(paths["minco"], taxid_to_gtdb)
        sylph_pred, sylph_extra = load_sylph_predictions(paths["sylph"], by_accession, by_core)
        scores.append(score_prediction(sample_id, "minco_universal_autoexact_gtdb_transfer", minco_pred, gtdb_truth, minco_extra))
        scores.append(score_prediction(sample_id, "sylph_gtdb_transfer", sylph_pred, gtdb_truth, sylph_extra))

    score_df = pd.DataFrame(scores)
    quality_df = pd.DataFrame(qualities)
    summary_df = summarize(score_df)
    pd.DataFrame(truth_rows_out).to_csv(RESULTS / "cami3_gtdb_taxid_transfer_truth.tsv", sep="\t", index=False)
    quality_df.to_csv(RESULTS / "cami3_gtdb_taxid_transfer_quality.tsv", sep="\t", index=False)
    score_df.to_csv(RESULTS / "cami3_gtdb_taxid_transfer_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / "cami3_gtdb_taxid_transfer_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
