#!/usr/bin/env python3
"""Fair full-profile L1 comparison and a guarded edge-EM trigger pilot."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
RUN = Path("/tmp/minco_fair_l1_edge_adaptive_20260624")

HELPER = ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
CALIB = ROOT / "research/experiments/2026-06-21_minco_multisample_call_calibration/scripts"
sys.path.insert(0, str(HELPER))
sys.path.insert(0, str(CALIB))

from analyze_readwise_corrections import extract_accession, parse_gold_profile  # noqa: E402
from calibrate_multisample_calls import gold_scope  # noqa: E402


GTDB_METADATA = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]
MARINE_GOLD = Path("/tmp/gs_marine_short.profile")
STRAIN_GOLD = Path("/mnt/new3T/minco_cami2_strain_20260621/short_read/taxonomic_profile_0.txt")
MARINE_MINCO = RUN / "marine_current_best/minco_s1000_unique_zipaaf_sample0.tsv"
MARINE_SYLPH = RUN / "marine_sylph/profile.tsv"
STRAIN_SYLPH = Path("/mnt/new3T/minco_cami2_strain_20260621/sylph_sample0/profile.tsv")
STRAIN_FEATURES = Path(
    "/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/features/test_strain0.joined_features.tsv"
)
STRAIN_PRED = Path(
    "/mnt/new3T/minco_cami2_strain_20260621/external_test_train9_strain0_2/model_predictions.tsv"
)


def accession_core(accession: str) -> str:
    if not accession:
        return ""
    accession = re.sub(r"^(RS_|GB_)", "", str(accession))
    accession = accession.split(".", 1)[0]
    accession = re.sub(r"^GC[AF]_", "", accession)
    return accession


def canonical_accession(accession: str) -> str:
    return re.sub(r"^(RS_|GB_)", "", str(accession or ""))


def load_metadata_maps() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, str]]:
    usecols = [
        "accession",
        "ncbi_genbank_assembly_accession",
        "ncbi_species_taxid",
        "ncbi_taxid",
        "ncbi_organism_name",
        "gtdb_taxonomy",
    ]
    by_accession: dict[str, dict[str, str]] = {}
    by_core: dict[str, dict[str, str]] = {}
    scope_by_taxid: dict[str, str] = {}
    for path in GTDB_METADATA:
        meta = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str).fillna("")
        for rec in meta.to_dict("records"):
            taxid = str(rec.get("ncbi_species_taxid") or rec.get("ncbi_taxid") or "")
            if not taxid:
                continue
            taxonomy = str(rec.get("gtdb_taxonomy", ""))
            if taxonomy.startswith("d__Bacteria"):
                scope = "bacteria"
            elif taxonomy.startswith("d__Archaea"):
                scope = "archaea"
            else:
                scope = "other"
            out = {
                "taxid": taxid,
                "species_name": str(rec.get("ncbi_organism_name", "")),
                "scope": scope,
            }
            for key in [rec.get("accession", ""), rec.get("ncbi_genbank_assembly_accession", "")]:
                key = canonical_accession(str(key))
                if key:
                    by_accession.setdefault(key, out)
                    core = accession_core(key)
                    if core:
                        by_core.setdefault(core, out)
            scope_by_taxid.setdefault(taxid, scope)
    return by_accession, by_core, scope_by_taxid


def lookup_tax_record(
    accession: str,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
) -> Mapping[str, str]:
    accession = canonical_accession(accession)
    if accession in by_accession:
        return by_accession[accession]
    core = accession_core(accession)
    if core in by_core:
        return by_core[core]
    return {}


def add_tax_metadata(
    rows: pd.DataFrame,
    ref_col: str,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
) -> pd.DataFrame:
    rows = rows.copy()
    rows["accession"] = rows[ref_col].map(extract_accession)
    records = [lookup_tax_record(acc, by_accession, by_core) for acc in rows["accession"]]
    rows["taxid"] = [str(rec.get("taxid", "")) for rec in records]
    rows["species_name"] = [str(rec.get("species_name", "")) for rec in records]
    rows["scope"] = [str(rec.get("scope", "")) for rec in records]
    return rows


def load_minco_mapped(
    path: Path,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows = add_tax_metadata(rows, "Ref", by_accession, by_core)
    if "Ref_annotation" in rows.columns:
        missing = rows["taxid"].eq("")
        if missing.any():
            fallback = rows.loc[missing, "Ref_annotation"].map(extract_accession)
            records = [lookup_tax_record(acc, by_accession, by_core) for acc in fallback]
            rows.loc[missing, "taxid"] = [str(rec.get("taxid", "")) for rec in records]
            rows.loc[missing, "species_name"] = [str(rec.get("species_name", "")) for rec in records]
            rows.loc[missing, "scope"] = [str(rec.get("scope", "")) for rec in records]
    return rows


def load_sylph_mapped(
    path: Path,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    return add_tax_metadata(rows, "Genome_file", by_accession, by_core)


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0 and math.isfinite(v))
    if total <= 0.0:
        return {}
    return {str(k): float(v) / total for k, v in values.items() if v > 0.0 and math.isfinite(v)}


def truth_from_gold(gold: pd.DataFrame) -> dict[str, float]:
    grouped = gold.groupby("taxid", as_index=False)["gold_percentage"].sum()
    return normalize(dict(zip(grouped["taxid"].astype(str), grouped["gold_percentage"].astype(float))))


def collapse_abundance(rows: pd.DataFrame, taxid_col: str, abundance_col: str, *, percent: bool = False) -> dict[str, float]:
    if rows.empty or taxid_col not in rows or abundance_col not in rows:
        return {}
    work = rows.loc[rows[taxid_col].astype(str).astype(bool)].copy()
    work[abundance_col] = pd.to_numeric(work[abundance_col], errors="coerce").fillna(0.0)
    if percent:
        work[abundance_col] = work[abundance_col] / 100.0
    collapsed = work.groupby(taxid_col, as_index=False)[abundance_col].max()
    return normalize(dict(zip(collapsed[taxid_col].astype(str), collapsed[abundance_col].astype(float))))


def score_prediction(
    *,
    dataset: str,
    method: str,
    method_role: str,
    truth_scope: str,
    truth: Mapping[str, float],
    pred: Mapping[str, float],
    source: str,
    notes: str = "",
) -> dict[str, object]:
    truth = normalize(truth)
    pred = normalize(pred)
    truth_set = {k for k, v in truth.items() if v > 0.0}
    pred_set = {k for k, v in pred.items() if v > 0.0}
    tp = truth_set & pred_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set
    precision = len(tp) / len(pred_set) if pred_set else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0

    keys = sorted(truth_set | pred_set)
    y_true = np.array([truth.get(k, 0.0) for k in keys], dtype=float)
    y_pred = np.array([pred.get(k, 0.0) for k in keys], dtype=float)
    if len(keys) >= 2 and float(np.std(y_true)) > 0.0 and float(np.std(y_pred)) > 0.0:
        pearson = float(np.corrcoef(y_pred, y_true)[0, 1])
    else:
        pearson = float("nan")
    spearman = float(pd.Series(y_pred).corr(pd.Series(y_true), method="spearman")) if len(keys) >= 2 else float("nan")

    tp_keys = sorted(tp)
    tp_true = np.array([truth[k] for k in tp_keys], dtype=float)
    tp_pred = np.array([pred[k] for k in tp_keys], dtype=float)
    if len(tp_keys) >= 2 and float(np.std(tp_true)) > 0.0 and float(np.std(tp_pred)) > 0.0:
        tp_pearson = float(np.corrcoef(tp_pred, tp_true)[0, 1])
    else:
        tp_pearson = float("nan")
    return {
        "dataset": dataset,
        "sample": 0,
        "truth_scope": truth_scope,
        "method": method,
        "method_role": method_role,
        "truth_taxa": len(truth_set),
        "pred_taxa": len(pred_set),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "l1_pct_points": float(np.sum(np.abs(y_pred - y_true)) * 100.0),
        "mae_pct_points": float(np.mean(np.abs(y_pred - y_true)) * 100.0) if len(keys) else float("nan"),
        "pearson": pearson,
        "spearman": spearman,
        "tp_abundance_pearson": tp_pearson,
        "tp_abundance_mae": float(np.mean(np.abs(tp_pred - tp_true))) if len(tp_keys) else float("nan"),
        "missing_truth_mass": float(sum(truth[k] for k in fn)),
        "fp_pred_mass": float(sum(pred[k] for k in fp)),
        "source": source,
        "notes": notes,
    }


def cami3_rows() -> list[dict[str, object]]:
    presence = pd.read_csv(EXP.parent / "2026-06-23_cami3_toygut_abundance_rescue/sample_presence.tsv", sep="\t")
    abundance = pd.read_csv(EXP.parent / "2026-06-23_cami3_toygut_abundance_rescue/sample_abundance.tsv", sep="\t")
    methods = {
        "sylph": ("Sylph", "reported Sylph baseline"),
        "minco_ctxobj_active_robust_depth": (
            "MinCO current best ctxobj robust depth",
            "current cached non-edge MinCO best",
        ),
    }
    rows: list[dict[str, object]] = []
    for old_method, (method, role) in methods.items():
        p = presence.loc[(presence["sample_id"] == 0) & (presence["method"] == old_method)].iloc[0]
        a = abundance.loc[
            (abundance["sample_id"] == 0)
            & (abundance["method"] == old_method)
            & (abundance["renorm_pred"])
        ].iloc[0]
        rows.append(
            {
                "dataset": "CAMI3 ToyGut",
                "sample": 0,
                "truth_scope": "NCBI bacterial species taxid",
                "method": method,
                "method_role": role,
                "truth_taxa": int(p["truth_bacterial_species"]),
                "pred_taxa": int(p["pred_species"]),
                "TP": int(p["TP"]),
                "FP": int(p["FP"]),
                "FN": int(p["FN"]),
                "precision": float(p["precision"]),
                "recall": float(p["recall"]),
                "F1": float(p["F1"]),
                "l1_pct_points": float(a["l1_pct_points"]),
                "mae_pct_points": float(a["mae_pct_points"]),
                "pearson": float(a["pearson"]),
                "spearman": float(a["spearman"]),
                "tp_abundance_pearson": float("nan"),
                "tp_abundance_mae": float("nan"),
                "missing_truth_mass": float("nan"),
                "fp_pred_mass": float("nan"),
                "source": "2026-06-23_cami3_toygut_abundance_rescue/sample_presence.tsv; sample_abundance.tsv",
                "notes": "existing full-profile renormalized CAMI3 abundance score",
            }
        )
    return rows


def marine_rows(
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
) -> list[dict[str, object]]:
    truth = truth_from_gold(parse_gold_profile(MARINE_GOLD, "marmgCAMI2_short_read_sample_0"))
    minco = load_minco_mapped(MARINE_MINCO, by_accession, by_core)
    sylph = load_sylph_mapped(MARINE_SYLPH, by_accession, by_core)
    return [
        score_prediction(
            dataset="CAMI2 Marine",
            method="MinCO current best S1000 unique ZIP-AAF",
            method_role="current cached non-edge MinCO best rerun for full L1",
            truth_scope="NCBI species taxid excluding unidentified/unidentified plasmid",
            truth=truth,
            pred=collapse_abundance(minco, "taxid", "Normalized_abundance_depth"),
            source=str(MARINE_MINCO),
        ),
        score_prediction(
            dataset="CAMI2 Marine",
            method="Sylph",
            method_role="reported Sylph baseline rerun for full L1",
            truth_scope="NCBI species taxid excluding unidentified/unidentified plasmid",
            truth=truth,
            pred=collapse_abundance(sylph, "taxid", "Taxonomic_abundance", percent=True),
            source=str(MARINE_SYLPH),
        ),
    ]


def strain_truth() -> dict[str, float]:
    gold = gold_scope(parse_gold_profile(STRAIN_GOLD, ""), "bacteria")
    return truth_from_gold(gold)


def filtered_sylph_values(
    sylph: pd.DataFrame,
    scope_by_taxid: Mapping[str, str],
    scope: str,
) -> dict[str, float]:
    values = collapse_abundance(sylph, "taxid", "Taxonomic_abundance", percent=True)
    if scope == "all":
        return values
    return {taxid: value for taxid, value in values.items() if scope_by_taxid.get(str(taxid), "other") == scope}


def strain_current_best_values() -> dict[str, float]:
    pred = pd.read_csv(STRAIN_PRED, sep="\t")
    pred_taxids = set(
        pred.loc[
            (pred["sample_key"] == "strain0")
            & (pred["method"] == "train9_rf_hgb_avg")
            & (pred["predicted"].astype(int) == 1),
            "taxid",
        ].astype(str)
    )
    features = pd.read_csv(STRAIN_FEATURES, sep="\t")
    features["taxid"] = features["taxid"].astype(str)
    features = features.loc[features["taxid"].isin(pred_taxids)].copy()
    values = {}
    for row in features.itertuples(index=False):
        values[str(row.taxid)] = max(
            finite_float(getattr(row, "u_Normalized_abundance_depth_max", 0.0)),
            finite_float(getattr(row, "s_Normalized_abundance_depth_max", 0.0)),
        )
    return normalize(values)


def strain_rows(
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, Mapping[str, str]],
    scope_by_taxid: Mapping[str, str],
) -> list[dict[str, object]]:
    truth = strain_truth()
    sylph = load_sylph_mapped(STRAIN_SYLPH, by_accession, by_core)
    return [
        score_prediction(
            dataset="CAMI2 Strain Madness",
            method="Sylph",
            method_role="reported Sylph baseline",
            truth_scope="NCBI bacterial species taxid",
            truth=truth,
            pred=filtered_sylph_values(sylph, scope_by_taxid, "bacteria"),
            source=str(STRAIN_SYLPH),
        ),
        score_prediction(
            dataset="CAMI2 Strain Madness",
            method="MinCO current best train9 RF/HGB average",
            method_role="current cached non-edge MinCO best with feature-table abundance",
            truth_scope="NCBI bacterial species taxid",
            truth=truth,
            pred=strain_current_best_values(),
            source=f"{STRAIN_PRED}; {STRAIN_FEATURES}",
            notes="abundance uses max(unique, split) normalized abundance among RF/HGB selected taxids",
        ),
    ]


def edge_rows() -> tuple[list[dict[str, object]], pd.DataFrame]:
    edge = pd.read_csv(RESULTS / "external_edge_em_sample_metrics.tsv", sep="\t")
    summary = pd.read_csv(RESULTS / "external_edge_em_summary.tsv", sep="\t")
    name_map = {
        "cami3_toygut": ("CAMI3 ToyGut", "NCBI bacterial species taxid"),
        "marine": ("CAMI2 Marine", "NCBI species taxid excluding unidentified/unidentified plasmid"),
        "strain_madness": ("CAMI2 Strain Madness", "NCBI bacterial species taxid"),
    }
    rows = []
    for _, rec in summary.iterrows():
        if rec["summary_role"] not in {"marker baseline from same edge run", "best positive beta diagnostic"}:
            continue
        dataset, scope = name_map[str(rec["dataset"])]
        role = (
            "edge-EM S2000 marker beta=0"
            if float(rec["beta"]) == 0.0
            else "edge-EM best positive beta diagnostic"
        )
        method = "MinCO S2000 edge marker beta=0" if float(rec["beta"]) == 0.0 else "MinCO edge-EM selected_group_species2"
        rows.append(metric_from_edge_row(rec, dataset, scope, method, role))

    adaptive_rows = []
    for dataset_key, group in edge.groupby("dataset", sort=False):
        base = group.loc[group["beta"] == 0.0].iloc[0]
        marker_targets = int(base["marker_raw_targets"])
        if marker_targets >= 150:
            beta = 0.02
            reason = "rich marker callset; marine-like high target count"
        elif marker_targets <= 25:
            beta = 0.0015
            reason = "narrow marker callset; strain-like low target count"
        else:
            beta = 0.0
            reason = "mid-size callset guard; keep marker baseline to avoid CAMI3 L1 regression"
        chosen = group.loc[np.isclose(group["beta"].astype(float), beta)].iloc[0]
        dataset, scope = name_map[str(dataset_key)]
        adaptive = metric_from_edge_row(
            chosen,
            dataset,
            scope,
            "MinCO adaptive edge-EM trigger",
            "guarded no-truth trigger pilot",
        )
        adaptive["adaptive_beta"] = beta
        adaptive["adaptive_reason"] = reason
        adaptive_rows.append(adaptive)
        rows.append(adaptive)
    return rows, pd.DataFrame(adaptive_rows)


def metric_from_edge_row(
    rec: pd.Series,
    dataset: str,
    truth_scope: str,
    method: str,
    method_role: str,
) -> dict[str, object]:
    return {
        "dataset": dataset,
        "sample": int(rec["sample_id"]),
        "truth_scope": truth_scope,
        "method": method,
        "method_role": method_role,
        "truth_taxa": int(rec["truth_taxa"]),
        "pred_taxa": int(rec["pred_taxa"]),
        "TP": int(rec["TP"]),
        "FP": int(rec["FP"]),
        "FN": int(rec["FN"]),
        "precision": float(rec["precision"]),
        "recall": float(rec["recall"]),
        "F1": float(rec["F1"]),
        "l1_pct_points": float(rec["l1_pct_points"]),
        "mae_pct_points": float(rec["mae_pct_points"]),
        "pearson": float(rec["pearson"]),
        "spearman": float(rec["spearman"]),
        "tp_abundance_pearson": float(rec["tp_abundance_pearson"]),
        "tp_abundance_mae": float(rec["tp_abundance_mae"]),
        "missing_truth_mass": float(rec["missing_truth_mass"]),
        "fp_pred_mass": float(rec["fp_pred_mass"]),
        "source": "results/external_edge_em_sample_metrics.tsv; results/external_edge_em_summary.tsv",
        "notes": f"edge beta={float(rec['beta']):g}",
    }


def adaptive_delta_table(rows: pd.DataFrame) -> pd.DataFrame:
    records = []
    for dataset, group in rows.groupby("dataset", sort=False):
        edge0 = group.loc[group["method"].eq("MinCO S2000 edge marker beta=0")]
        adaptive = group.loc[group["method"].eq("MinCO adaptive edge-EM trigger")]
        best_nonedge = group.loc[group["method_role"].str.contains("current cached non-edge MinCO best", regex=False)]
        sylph = group.loc[group["method"].eq("Sylph")]
        if edge0.empty or adaptive.empty:
            continue
        a = adaptive.iloc[0]
        b = edge0.iloc[0]
        rec = {
            "dataset": dataset,
            "adaptive_beta": finite_float(a.get("adaptive_beta", 0.0)),
            "adaptive_reason": a.get("adaptive_reason", ""),
            "edge_marker_l1": b["l1_pct_points"],
            "adaptive_l1": a["l1_pct_points"],
            "adaptive_minus_edge_marker_l1": a["l1_pct_points"] - b["l1_pct_points"],
            "edge_marker_F1": b["F1"],
            "adaptive_F1": a["F1"],
            "adaptive_minus_edge_marker_F1": a["F1"] - b["F1"],
        }
        if not best_nonedge.empty:
            m = best_nonedge.iloc[0]
            rec["current_best_minco_l1"] = m["l1_pct_points"]
            rec["adaptive_minus_current_best_l1"] = a["l1_pct_points"] - m["l1_pct_points"]
            rec["current_best_minco_F1"] = m["F1"]
            rec["adaptive_minus_current_best_F1"] = a["F1"] - m["F1"]
        if not sylph.empty:
            s = sylph.iloc[0]
            rec["sylph_l1"] = s["l1_pct_points"]
            rec["adaptive_minus_sylph_l1"] = a["l1_pct_points"] - s["l1_pct_points"]
            rec["sylph_F1"] = s["F1"]
            rec["adaptive_minus_sylph_F1"] = a["F1"] - s["F1"]
        records.append(rec)
    return pd.DataFrame(records)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core, scope_by_taxid = load_metadata_maps()
    rows = []
    rows.extend(cami3_rows())
    rows.extend(marine_rows(by_accession, by_core))
    rows.extend(strain_rows(by_accession, by_core, scope_by_taxid))
    edge_metric_rows, adaptive_rows = edge_rows()
    rows.extend(edge_metric_rows)

    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "fair_full_l1_methods.tsv", sep="\t", index=False)
    adaptive_rows.to_csv(RESULTS / "adaptive_edge_em_policy.tsv", sep="\t", index=False)
    adaptive_delta_table(out).to_csv(RESULTS / "adaptive_edge_em_delta.tsv", sep="\t", index=False)
    print(out.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
