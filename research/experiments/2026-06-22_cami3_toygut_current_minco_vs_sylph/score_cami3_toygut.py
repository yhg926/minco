#!/usr/bin/env python3
"""Score current MinCO S2000 markerdb against Sylph on CAMI3 ToyGut samples."""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(OLD_EXP))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402


EXP_DIR = Path(__file__).resolve().parent
RUN_DIR = Path("/tmp/cami3_toygut_current_minco_vs_sylph_20260622")
CAMI3_DIR = Path("/mnt/new3T/minco_cami3_toygut_20260620")

EFFECTIVE_CTX_LENGTH = 24
ANI_THRESHOLD = 0.95
ACTIVE_CTX_MIN = 15.0
ACTIVE_RELIABLE_ZTP_AF_FLOOR = 0.40
ACTIVE_MEAN_DEPTH_MIN = 3.0
ACTIVE_VMR_MIN = 50.0
ACTIVE_DELTA_MAX = 0.03

SAMPLES = {
    0: {
        "truth": CAMI3_DIR / "taxonomic_profiles/taxonomic_profile_0.txt",
        "minco": RUN_DIR / "minco_s0_current_product0.tsv",
        "sylph": CAMI3_DIR / "sylph_sample0/profile.tsv",
    },
    1: {
        "truth": CAMI3_DIR / "taxonomic_profiles/taxonomic_profile_1.txt",
        "minco": RUN_DIR / "minco_s1_current_product0.tsv",
        "sylph": RUN_DIR / "sylph_sample1/profile.tsv",
    },
    2: {
        "truth": CAMI3_DIR / "taxonomic_profiles/taxonomic_profile_2.txt",
        "minco": RUN_DIR / "minco_s2_current_product0.tsv",
        "sylph": RUN_DIR / "sylph_sample2/profile.tsv",
    },
}


def numeric(rows: pd.DataFrame, col: str) -> pd.Series:
    if col not in rows.columns:
        return pd.Series([0.0] * len(rows), index=rows.index)
    return pd.to_numeric(rows[col], errors="coerce").fillna(0.0)


def score_sets(predicted: Iterable[str], gold: Iterable[str]):
    pred = {str(x) for x in predicted if str(x)}
    truth_set = {str(x) for x in gold if str(x)}
    tp = pred & truth_set
    fp = pred - truth_set
    fn = truth_set - pred
    precision = len(tp) / len(pred) if pred else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def ztp_lambda_from_pos_mean(mean_pos: float) -> float:
    if not math.isfinite(mean_pos) or mean_pos <= 0.0:
        return 0.0
    if mean_pos <= 1.0 + 1e-10:
        return 1e-10

    lo = 1e-10
    hi = max(2.0, mean_pos * 2.0)

    def cond_mean(lam: float) -> float:
        return lam / (1.0 - math.exp(-lam))

    while cond_mean(hi) < mean_pos and hi < 1e6:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if cond_mean(mid) < mean_pos:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ztp_adjusted_af(breadth: float, mean_pos: float) -> float:
    if not math.isfinite(breadth) or breadth <= 0.0:
        return 0.0
    lam = ztp_lambda_from_pos_mean(mean_pos)
    p_nonzero = 1.0 - math.exp(-lam) if lam > 0.0 else 0.0
    if p_nonzero <= 0.0:
        return 1.0
    return min(1.0, breadth / p_nonzero)


def lookup_rec(accession: str, by_accession, by_core):
    rec, method, key = truth.lookup_accession(accession, by_accession, by_core)
    return rec, method, key


def add_metadata(rows: pd.DataFrame, accession_col: str, by_accession, by_core) -> pd.DataFrame:
    mapped = [lookup_rec(acc, by_accession, by_core) for acc in rows[accession_col].astype(str)]
    rows["metadata_mapping_method"] = [method for _, method, _ in mapped]
    rows["metadata_mapping_key"] = [key for _, _, key in mapped]
    rows["gtdb_species"] = [rec["gtdb_species"] if rec else "" for rec, _, _ in mapped]
    rows["ncbi_species_taxid"] = [
        (rec.get("ncbi_species_taxid") or rec.get("ncbi_taxid") or "") if rec else ""
        for rec, _, _ in mapped
    ]
    rows["ncbi_taxid"] = [rec.get("ncbi_taxid", "") if rec else "" for rec, _, _ in mapped]
    rows["ncbi_species"] = [rec.get("ncbi_species", "") if rec else "" for rec, _, _ in mapped]
    rows["ncbi_organism_name"] = [rec.get("ncbi_organism_name", "") if rec else "" for rec, _, _ in mapped]
    return rows


def add_reliable_ztp_and_active_gate(rows: pd.DataFrame) -> pd.DataFrame:
    rows = rows.copy()
    rows = truth.add_naive_ani(rows)
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_depth_variance",
        "Normalized_abundance_depth",
    ]:
        rows[col] = numeric(rows, col)

    rows["Reliable_ztp_af"] = [
        ztp_adjusted_af(float(b), float(m))
        for b, m in zip(rows["Reliable_Ref_breadth"], rows["Reliable_Ref_hit_mean_depth"])
    ]
    rows["ANI_from_Reliable_ztp_af"] = [
        1.0 + math.log(max(float(af), 1e-300)) / EFFECTIVE_CTX_LENGTH
        for af in rows["Reliable_ztp_af"]
    ]
    rows["ANI_AF_delta"] = rows["ANI_naive_calc"] - rows["ANI_from_Reliable_ztp_af"]
    rows["Reliable_depth_vmr"] = [
        float(var) / float(mean) if float(mean) > 0.0 else 0.0
        for mean, var in zip(rows["Reliable_Ref_hit_mean_depth"], rows["Reliable_Ref_hit_depth_variance"])
    ]
    rows["active_delta_trigger"] = (
        (rows["Reliable_Ref_hit_mean_depth"] > ACTIVE_MEAN_DEPTH_MIN)
        & (rows["Reliable_depth_vmr"] > ACTIVE_VMR_MIN)
    )
    rows["active_delta_pass"] = (
        ~rows["active_delta_trigger"] | (rows["ANI_AF_delta"] < ACTIVE_DELTA_MAX)
    )
    rows["active_gate_pass"] = (
        (rows["XnY_ctx"] >= ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= ACTIVE_RELIABLE_ZTP_AF_FLOOR)
        & rows["active_delta_pass"]
        & rows["ncbi_species_taxid"].astype(bool)
    )
    return rows


def load_minco(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows["accession"] = rows["Ref"].map(truth.extract_accession)
    if "Ref_annotation" in rows.columns:
        missing = rows["accession"] == ""
        rows.loc[missing, "accession"] = rows.loc[missing, "Ref_annotation"].map(truth.extract_accession)
    rows = add_metadata(rows, "accession", by_accession, by_core)
    return add_reliable_ztp_and_active_gate(rows)


def load_sylph(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows["accession"] = rows["Genome_file"].map(truth.extract_accession)
    rows = add_metadata(rows, "accession", by_accession, by_core)
    rows["Taxonomic_abundance"] = numeric(rows, "Taxonomic_abundance")
    rows["Adjusted_ANI"] = numeric(rows, "Adjusted_ANI")
    return rows


def load_truth(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("@") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            taxid, rank, taxpath, taxpathsn, pct = fields[:5]
            if rank != "species":
                continue
            try:
                abundance_pct = float(pct)
            except ValueError:
                continue
            if abundance_pct <= 0.0:
                continue
            if not (taxpath.startswith("131567|2") or taxpathsn.startswith("cellular organisms|Bacteria")):
                continue
            rows.append(
                {
                    "ncbi_species_taxid": taxid,
                    "truth_name": taxpathsn.split("|")[-1],
                    "truth_abundance_pct_all": abundance_pct,
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    grouped = (
        out.groupby("ncbi_species_taxid", as_index=False)
        .agg(
            truth_name=("truth_name", "first"),
            truth_abundance_pct_all=("truth_abundance_pct_all", "sum"),
        )
    )
    total = grouped["truth_abundance_pct_all"].sum()
    grouped["truth_abundance_bacterial_norm"] = (
        grouped["truth_abundance_pct_all"] / total if total else 0.0
    )
    return grouped


def select_predictions(tool: str, rows: pd.DataFrame) -> pd.DataFrame:
    if tool == "minco_current":
        selected = rows.loc[rows["active_gate_pass"]].copy()
        selected["pred_abundance_raw"] = numeric(selected, "Normalized_abundance_depth")
        selected["pred_ani"] = numeric(selected, "ANI_naive_calc")
        selected["support"] = numeric(selected, "XnY_ctx")
    elif tool == "sylph":
        selected = rows.loc[rows["ncbi_species_taxid"].astype(bool)].copy()
        selected["pred_abundance_raw"] = numeric(selected, "Taxonomic_abundance") / 100.0
        selected["pred_ani"] = numeric(selected, "Adjusted_ANI") / 100.0
        selected["support"] = numeric(selected, "Eff_cov")
    else:
        raise ValueError(tool)
    return selected


def collapse_predictions(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame(
            columns=[
                "ncbi_species_taxid",
                "ncbi_species",
                "gtdb_species",
                "accession",
                "pred_abundance_raw",
                "pred_ani",
                "support",
            ]
        )
    selected = selected.copy()
    selected["pred_abundance_raw"] = numeric(selected, "pred_abundance_raw")
    selected["pred_ani"] = numeric(selected, "pred_ani")
    selected["support"] = numeric(selected, "support")
    selected = selected.sort_values(
        ["ncbi_species_taxid", "pred_abundance_raw", "pred_ani", "support"],
        ascending=[True, False, False, False],
    )
    return (
        selected.groupby("ncbi_species_taxid", as_index=False)
        .agg(
            ncbi_species=("ncbi_species", "first"),
            gtdb_species=("gtdb_species", "first"),
            accession=("accession", "first"),
            pred_abundance_raw=("pred_abundance_raw", "max"),
            pred_ani=("pred_ani", "max"),
            support=("support", "max"),
        )
    )


def abundance_metrics(collapsed: pd.DataFrame, truth_rows: pd.DataFrame):
    truth_abund = dict(
        zip(
            truth_rows["ncbi_species_taxid"].astype(str),
            truth_rows["truth_abundance_bacterial_norm"].astype(float),
        )
    )
    pred_raw = dict(
        zip(
            collapsed["ncbi_species_taxid"].astype(str),
            collapsed["pred_abundance_raw"].astype(float),
        )
    )
    out = []
    for renorm in [False, True]:
        pred = pred_raw.copy()
        if renorm:
            total = sum(pred.values())
            if total > 0:
                pred = {k: v / total for k, v in pred.items()}
        y_true = []
        y_pred = []
        for taxid, abundance in truth_abund.items():
            y_true.append(float(abundance))
            y_pred.append(float(pred.get(taxid, 0.0)))
        true_s = pd.Series(y_true, dtype=float)
        pred_s = pd.Series(y_pred, dtype=float)
        out.append(
            {
                "renorm_pred": renorm,
                "truth_taxa": len(y_true),
                "pred_sum_on_truth": sum(y_pred),
                "pred_sum_all": sum(pred.values()),
                "truth_sum": sum(y_true),
                "pearson": pred_s.corr(true_s, method="pearson"),
                "spearman": pred_s.corr(true_s, method="spearman"),
                "mae_pct_points": (pred_s - true_s).abs().mean() * 100.0,
                "l1_pct_points": (pred_s - true_s).abs().sum() * 100.0,
            }
        )
    return out


def detail_rows(sample_id: int, method: str, kind: str, taxids: set[str], collapsed: pd.DataFrame, truth_rows: pd.DataFrame):
    truth_by_taxid = truth_rows.set_index("ncbi_species_taxid").to_dict("index")
    pred_by_taxid = collapsed.set_index("ncbi_species_taxid").to_dict("index")
    rows = []
    for taxid in sorted(taxids):
        truth_rec = truth_by_taxid.get(taxid, {})
        pred_rec = pred_by_taxid.get(taxid, {})
        rows.append(
            {
                "sample_id": sample_id,
                "method": method,
                "kind": kind,
                "ncbi_species_taxid": taxid,
                "truth_name": truth_rec.get("truth_name", ""),
                "truth_abundance_pct_all": truth_rec.get("truth_abundance_pct_all", ""),
                "ncbi_species": pred_rec.get("ncbi_species", ""),
                "gtdb_species": pred_rec.get("gtdb_species", ""),
                "accession": pred_rec.get("accession", ""),
                "pred_abundance_raw": pred_rec.get("pred_abundance_raw", ""),
                "pred_ani": pred_rec.get("pred_ani", ""),
                "support": pred_rec.get("support", ""),
            }
        )
    return rows


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    score_rows = []
    abundance_rows = []
    details = []
    selected_rows = []

    for sample_id, paths in SAMPLES.items():
        truth_rows = load_truth(paths["truth"])
        gold = set(truth_rows["ncbi_species_taxid"].astype(str))

        for method, tool in [("minco_current", "minco"), ("sylph", "sylph")]:
            raw = load_minco(paths["minco"], by_accession, by_core) if tool == "minco" else load_sylph(paths["sylph"], by_accession, by_core)
            selected = select_predictions(method, raw)
            collapsed = collapse_predictions(selected)
            pred = set(collapsed["ncbi_species_taxid"].astype(str))
            tp, fp, fn, precision, recall, f1 = score_sets(pred, gold)
            score_rows.append(
                {
                    "sample_id": sample_id,
                    "method": method,
                    "truth_bacterial_species": len(gold),
                    "truth_bacterial_abundance_pct_all": truth_rows["truth_abundance_pct_all"].sum(),
                    "pred_species": len(pred),
                    "TP": len(tp),
                    "FP": len(fp),
                    "FN": len(fn),
                    "precision": precision,
                    "recall": recall,
                    "F1": f1,
                    "raw_selected_rows": len(selected),
                    "mapped_selected_rows": int(selected["ncbi_species_taxid"].astype(bool).sum()) if "ncbi_species_taxid" in selected else len(selected),
                }
            )
            for row in abundance_metrics(collapsed, truth_rows):
                abundance_rows.append({"sample_id": sample_id, "method": method, **row})
            details.extend(detail_rows(sample_id, method, "FP", fp, collapsed, truth_rows))
            details.extend(detail_rows(sample_id, method, "FN", fn, collapsed, truth_rows))
            selected_rows.append(collapsed.assign(sample_id=sample_id, method=method))

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    details_df = pd.DataFrame(details)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()

    score_path = RUN_DIR / "cami3_toygut_minco_current_vs_sylph_scores.tsv"
    abundance_path = RUN_DIR / "cami3_toygut_minco_current_vs_sylph_abundance.tsv"
    details_path = RUN_DIR / "cami3_toygut_minco_current_vs_sylph_details.tsv"
    selected_path = RUN_DIR / "cami3_toygut_minco_current_vs_sylph_selected_species.tsv"

    score_df.to_csv(score_path, sep="\t", index=False)
    abundance_df.to_csv(abundance_path, sep="\t", index=False)
    details_df.to_csv(details_path, sep="\t", index=False)
    selected_df.to_csv(selected_path, sep="\t", index=False)

    print(score_df.to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {abundance_path}")
    print(f"wrote {details_path}")
    print(f"wrote {selected_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
