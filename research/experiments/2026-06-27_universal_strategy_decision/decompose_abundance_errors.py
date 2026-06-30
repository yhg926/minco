#!/usr/bin/env python3
"""Decompose MinCO/Sylph abundance L1 errors from cached scored panels.

The release-decision scorers already wrote GTDB-species truth tables. This
diagnostic reuses those truth tables and the existing profile outputs so it can
be rerun quickly. For MinCO calibrated profiles, `species_name` is the GTDB
species label supplied by the species taxmap. For Sylph, reference accessions
are mapped to GTDB species through the same helpers used by the score scripts.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import score_cami3_gtdb_source_readmap as cami3
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp
import score_cami3_gtdb_taxid_transfer as taxid_score


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
TOYMOUSE_SCORER = (
    ROOT
    / "research/experiments/2026-06-26_cami2_toymouse_current_default/"
    "score_sample7_current_default.py"
)

PANELS = {
    "cami2_toy_mouse_gut": {
        "samples": [5, 6, 7],
        "truth": lambda sample: RESULTS / f"mouse{sample}_gtdb_species_profile.tsv",
        "truth_abundance_col": "relative_abundance",
        "minco_method": "minco_current_code_refresh",
        "minco": lambda sample: Path(f"/tmp/minco_current_code_toymouse_refresh_20260627/sample{sample}_current.tsv"),
        "minco_collapse": "max",
        "sylph_method": "sylph_gtdb_profile",
        "sylph": lambda sample: Path(f"/tmp/cami2_toymouse_samples5_7_20260625/run/sylph_sample{sample}/profile.tsv"),
        "sylph_collapse": "max",
        "official_scores": RESULTS / "toymouse_current_refresh_compare.tsv",
        "official_minco_method": "minco_current_code_refresh",
        "official_sylph_method": "sylph_gtdb_profile",
    },
    "hmp_airskin_gtdb_source_abundance": {
        "samples": [0, 1, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 28],
        "truth": lambda _sample: RESULTS / "hmp_current_refresh_r232_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28_truth.tsv",
        "truth_abundance_col": "truth_abundance",
        "minco_method": "minco_current_default_gtdb_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28",
        "minco": lambda sample: Path(
            {
                0: "/tmp/minco_current_code_hmp_airskin0_20260627/minco_sample0_current_default.tsv",
                1: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin1_20260627/minco_sample1_current_default.tsv",
                3: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin3_20260627/minco_sample3_current_default.tsv",
                4: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin4_20260627/minco_sample4_current_default.tsv",
                5: "/tmp/minco_current_code_hmp_airskin5_20260627/minco_sample5_current_default.tsv",
                6: "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_current_exact.tsv",
                7: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin7_20260627/minco_sample7_current_default.tsv",
                9: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin9_20260627/minco_sample9_current_default.tsv",
                10: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin10_20260627/minco_sample10_current_default.tsv",
                11: "/tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_current_exact.tsv",
                13: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin13_20260627/minco_sample13_current_default.tsv",
                14: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin14_20260627/minco_sample14_current_default.tsv",
                15: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin15_20260627/minco_sample15_current_default.tsv",
                16: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin16_20260627/minco_sample16_current_default.tsv",
                17: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin17_20260627/minco_sample17_current_default.tsv",
                18: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin18_20260627/minco_sample18_current_default.tsv",
                19: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin19_20260627/minco_sample19_current_default.tsv",
                20: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin20_20260627/minco_sample20_current_default.tsv",
                21: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin21_20260627/minco_sample21_current_default.tsv",
                22: "/tmp/minco_current_code_hmp_airskin22_20260627/minco_sample22_current_default.tsv",
                23: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin23_20260627/minco_sample23_current_default.tsv",
                24: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin24_20260627/minco_sample24_current_default.tsv",
                25: "/mnt/new3T/minco_release_holdouts_20260628/minco_current_code_hmp_airskin25_20260627/minco_sample25_current_default.tsv",
                28: "/tmp/minco_current_code_hmp_airskin28_20260627/minco_sample28_current_default.tsv",
            }[sample]
        ),
        "minco_collapse": "sum",
        "sylph_method": "sylph_gtdb_source_abundance",
        "sylph": lambda sample: Path(
            {
                0: "/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample0/profile.tsv",
                1: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample1/profile.tsv",
                3: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample3/profile.tsv",
                4: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample4/profile.tsv",
                5: "/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample5/profile.tsv",
                6: "/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample6/profile.tsv",
                7: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample7/profile.tsv",
                9: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample9/profile.tsv",
                10: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample10/profile.tsv",
                11: "/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample11/profile.tsv",
                13: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample13/profile.tsv",
                14: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample14/profile.tsv",
                15: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample15/profile.tsv",
                16: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample16/profile.tsv",
                17: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample17/profile.tsv",
                18: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample18/profile.tsv",
                19: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample19/profile.tsv",
                20: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample20/profile.tsv",
                21: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample21/profile.tsv",
                22: "/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample22/profile.tsv",
                23: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample23/profile.tsv",
                24: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample24/profile.tsv",
                25: "/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232/sylph_sample25/profile.tsv",
                28: "/tmp/cami2_hmp_airskin_20260625/run_sylph_r232/sylph_sample28/profile.chunked.tsv",
            }[sample]
        ),
        "sylph_collapse": "sum",
        "official_scores": RESULTS / "hmp_current_refresh_r232_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28_scores.tsv",
        "official_minco_method": "minco_current_default_gtdb_source_abundance_sample19_0_1_3_4_5_6_7_9_10_11_13_14_15_16_17_18_20_21_22_23_24_25_28",
        "official_sylph_method": "sylph_gtdb_source_abundance",
    },
    "hmp_gastrooral_gtdb_source_abundance": {
        "samples": [0, 6],
        "truth": lambda _sample: RESULTS / "hmp_gastrooral_r232_source_abundance_truth.tsv",
        "truth_abundance_col": "truth_abundance",
        "minco_method": "minco_current_raw_default_gastrooral_source_abundance",
        "minco": lambda sample: Path(
            {
                0: "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample0_current_default.tsv",
                6: "/tmp/minco_current_code_hmp_gastrooral_20260627/minco_sample6_current_default.tsv",
            }[int(sample)]
        ),
        "minco_collapse": "sum",
        "sylph_method": "sylph_gtdb_r232_gastrooral_source_abundance",
        "sylph": lambda sample: hmp_gastro.SAMPLES[int(sample)]["sylph_r232"],
        "sylph_collapse": "sum",
        "official_scores": RESULTS / "hmp_gastrooral_raw_default_r232_source_abundance_scores.tsv",
        "official_minco_method": "minco_current_raw_default_gastrooral_source_abundance",
        "official_sylph_method": "sylph_gtdb_r232_gastrooral_source_abundance",
    },
    "cami3_toy_human_gut_gtdb_source_readmap": {
        "samples": [0, 1, 2],
        "truth": lambda _sample: RESULTS / "cami3_gtdb_source_readmap_truth.tsv",
        "truth_abundance_col": "truth_abundance",
        "minco_method": "minco_universal_autoexact_gtdb_source_readmap",
        "minco": lambda sample: Path(
            {
                0: "/tmp/cami3_toy_human_gut_20260626/run/sample0_universal_autoexact.tsv",
                1: "/tmp/cami3_toy_human_gut_20260626/run/sample1_universal_autoexact.tsv",
                2: "/tmp/cami3_toy_human_gut_20260626/run/sample2_universal_autoexact_tail_p025.tsv",
            }[sample]
        ),
        "minco_collapse": "sum",
        "sylph_method": "sylph_gtdb_source_readmap",
        "sylph": lambda sample: Path(
            {
                0: "/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv",
                1: "/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv",
                2: "/tmp/cami3_toy_human_gut_20260626/run/sylph_sample2/profile.tsv",
            }[sample]
        ),
        "sylph_collapse": "sum",
        "official_scores": RESULTS / "cami3_gtdb_source_readmap_scores.tsv",
        "official_minco_method": "minco_universal_autoexact_gtdb_source_readmap",
        "official_sylph_method": "sylph_gtdb_source_readmap",
    },
}


def load_toy_scorer():
    spec = importlib.util.spec_from_file_location("toy_current_scorer", TOYMOUSE_SCORER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import scorer {TOYMOUSE_SCORER}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(k): max(0.0, finite_float(v)) for k, v in values.items() if str(k)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def collapse_values(df: pd.DataFrame, species_col: str, abundance_col: str, collapse: str) -> dict[str, float]:
    if df.empty:
        return {}
    work = df.loc[df[species_col].fillna("").astype(str).astype(bool)].copy()
    work[abundance_col] = pd.to_numeric(work[abundance_col], errors="coerce").fillna(0.0)
    if collapse == "max":
        grouped = work.groupby(species_col)[abundance_col].max()
    elif collapse == "sum":
        grouped = work.groupby(species_col)[abundance_col].sum()
    else:
        raise ValueError(f"unsupported collapse rule: {collapse}")
    return normalize({str(k): finite_float(v) for k, v in grouped.to_dict().items()})


def load_truth(path: Path, sample: int, abundance_col: str) -> dict[str, float]:
    df = pd.read_csv(path, sep="\t")
    if "sample" in df.columns:
        df = df.loc[pd.to_numeric(df["sample"], errors="coerce").fillna(-1).astype(int).eq(sample)]
    return collapse_values(df, "gtdb_species", abundance_col, "sum")


def load_minco_pred(path: Path, collapse: str) -> dict[str, float]:
    usecols = ["species_name", "calibrated_call", "calibrated_abundance"]
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in usecols, low_memory=False)
    call = df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    df = df.loc[call].copy()
    return collapse_values(df, "species_name", "calibrated_abundance", collapse)


def load_toy_sylph_pred(path: Path, scorer_mod, by_accession, by_core, collapse: str) -> dict[str, float]:
    rows = scorer_mod.score.load_sylph_rows(path, by_accession, by_core)
    rows = rows.loc[rows["active_gate_pass"]].copy()
    return collapse_values(rows, "gtdb_species", "Normalized_abundance_depth", collapse)


def load_sylph_pred(path: Path, by_accession, by_core) -> dict[str, float]:
    pred, _extra = taxid_score.load_sylph_predictions(path, by_accession, by_core)
    if pred.empty:
        return {}
    return normalize(dict(zip(pred["gtdb_species"].astype(str), pred["pred_abundance"])))


def pred_from_collapsed(df: pd.DataFrame) -> dict[str, float]:
    if df.empty:
        return {}
    return normalize(dict(zip(df["gtdb_species"].astype(str), df["pred_abundance"])))


def decompose(
    panel: str,
    sample: str,
    method: str,
    truth: Mapping[str, float],
    pred: Mapping[str, float],
    truth_source: str,
    pred_source: str,
    prediction_rule: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    truth_norm = normalize(truth)
    pred_norm = normalize(pred)
    truth_species = set(truth_norm)
    pred_species = set(pred_norm)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species
    species_rows: list[dict[str, object]] = []
    matched_abs = 0.0
    missing = 0.0
    extra = 0.0
    for species in sorted(truth_species | pred_species):
        t = truth_norm.get(species, 0.0)
        p = pred_norm.get(species, 0.0)
        if species in tp:
            status = "TP"
            matched_abs += abs(p - t)
        elif species in fn:
            status = "FN"
            missing += t
        else:
            status = "FP"
            extra += p
        species_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "method": method,
                "gtdb_species": species,
                "status": status,
                "truth_abundance": t,
                "pred_abundance": p,
                "signed_error_pp": (p - t) * 100.0,
                "abs_error_pp": abs(p - t) * 100.0,
                "truth_source": truth_source,
                "pred_source": pred_source,
                "prediction_rule": prediction_rule,
            }
        )

    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    truth_order = sorted(truth_species)
    union_order = sorted(truth_species | pred_species)
    y_true_truth = pd.Series([truth_norm.get(species, 0.0) for species in truth_order], dtype=float)
    y_pred_truth = pd.Series([pred_norm.get(species, 0.0) for species in truth_order], dtype=float)
    y_true_union = pd.Series([truth_norm.get(species, 0.0) for species in union_order], dtype=float)
    y_pred_union = pd.Series([pred_norm.get(species, 0.0) for species in union_order], dtype=float)
    row = {
        "panel": panel,
        "sample": sample,
        "method": method,
        "truth_species": len(truth_species),
        "pred_species": len(pred_species),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "matched_abs_error_pp": matched_abs * 100.0,
        "missing_truth_mass_pp": missing * 100.0,
        "extra_pred_mass_pp": extra * 100.0,
        "union_L1_pp": (matched_abs + missing + extra) * 100.0,
        "truth_only_L1_pp": (matched_abs + missing) * 100.0,
        "truth_mass_detected_pct": sum(truth_norm.get(species, 0.0) for species in tp) * 100.0,
        "pred_mass_on_truth_pct": sum(pred_norm.get(species, 0.0) for species in tp) * 100.0,
        "pred_mass_extra_pct": extra * 100.0,
        "Pearson_truth_only": y_pred_truth.corr(y_true_truth, method="pearson")
        if len(truth_order) > 1
        else float("nan"),
        "Pearson_union": y_pred_union.corr(y_true_union, method="pearson")
        if len(union_order) > 1
        else float("nan"),
        "truth_source": truth_source,
        "pred_source": pred_source,
        "prediction_rule": prediction_rule,
    }
    return row, species_rows


def official_validation(panel: str, rows: pd.DataFrame) -> pd.DataFrame:
    cfg = PANELS[panel]
    score_path = Path(cfg["official_scores"])
    if not score_path.exists():
        return pd.DataFrame()
    official = pd.read_csv(score_path, sep="\t")
    if "sample" not in official.columns:
        return pd.DataFrame()
    method_map = {
        cfg["minco_method"]: cfg["official_minco_method"],
        cfg["sylph_method"]: cfg["official_sylph_method"],
    }
    out = []
    for method, official_method in method_map.items():
        for sample in cfg["samples"]:
            ours = rows.loc[
                rows["panel"].eq(panel)
                & rows["method"].eq(method)
                & rows["sample"].astype(str).eq(str(sample))
            ]
            theirs = official.loc[
                official["method"].astype(str).eq(str(official_method))
                & official["sample"].astype(str).eq(str(sample))
            ]
            if ours.empty or theirs.empty:
                continue
            ours_row = ours.iloc[0]
            their_row = theirs.iloc[0]
            official_l1_col = "L1_union_pp" if "L1_union_pp" in theirs.columns else "l1_pct_points"
            official_pearson_col = (
                "Pearson_union" if "Pearson_union" in theirs.columns else "pearson"
            )
            out.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "method": method,
                    "official_method": official_method,
                    "TP_delta": int(ours_row["TP"]) - int(their_row["TP"]),
                    "FP_delta": int(ours_row["FP"]) - int(their_row["FP"]),
                    "FN_delta": int(ours_row["FN"]) - int(their_row["FN"]),
                    "F1_delta": finite_float(ours_row["F1"]) - finite_float(their_row["F1"]),
                    "official_L1_col": official_l1_col,
                    "L1_delta_pp": (
                        finite_float(ours_row["union_L1_pp" if official_l1_col == "L1_union_pp" else "truth_only_L1_pp"])
                        - finite_float(their_row[official_l1_col])
                    ),
                    "official_Pearson_col": official_pearson_col,
                    "Pearson_delta": (
                        finite_float(ours_row["Pearson_union" if official_pearson_col == "Pearson_union" else "Pearson_truth_only"])
                        - finite_float(their_row[official_pearson_col])
                    ),
                }
            )
    return pd.DataFrame(out)


def summarize(rows: pd.DataFrame) -> pd.DataFrame:
    out: list[dict[str, object]] = []
    for (panel, method), sub in rows.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        out.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique(), key=str))),
                "mean_F1": sub["F1"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_union_L1_pp": sub["union_L1_pp"].mean(),
                "mean_truth_only_L1_pp": sub["truth_only_L1_pp"].mean(),
                "mean_matched_abs_error_pp": sub["matched_abs_error_pp"].mean(),
                "mean_missing_truth_mass_pp": sub["missing_truth_mass_pp"].mean(),
                "mean_extra_pred_mass_pp": sub["extra_pred_mass_pp"].mean(),
                "mean_truth_mass_detected_pct": sub["truth_mass_detected_pct"].mean(),
                "mean_pred_mass_on_truth_pct": sub["pred_mass_on_truth_pct"].mean(),
                "mean_pred_mass_extra_pct": sub["pred_mass_extra_pct"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
            }
        )
    return pd.DataFrame(out)


def compare_methods(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, cfg in PANELS.items():
        panel_rows = summary.loc[summary["panel"].eq(panel)].set_index("method")
        minco_method = cfg["minco_method"]
        sylph_method = cfg["sylph_method"]
        if minco_method not in panel_rows.index or sylph_method not in panel_rows.index:
            continue
        m = panel_rows.loc[minco_method]
        s = panel_rows.loc[sylph_method]
        component_deltas = {
            "matched_abs_error": m["mean_matched_abs_error_pp"] - s["mean_matched_abs_error_pp"],
            "missing_truth_mass": m["mean_missing_truth_mass_pp"] - s["mean_missing_truth_mass_pp"],
            "extra_pred_mass": m["mean_extra_pred_mass_pp"] - s["mean_extra_pred_mass_pp"],
        }
        rows.append(
            {
                "panel": panel,
                "minco_method": minco_method,
                "sylph_method": sylph_method,
                "minco_minus_sylph_union_L1_pp": m["mean_union_L1_pp"] - s["mean_union_L1_pp"],
                "minco_minus_sylph_truth_only_L1_pp": m["mean_truth_only_L1_pp"] - s["mean_truth_only_L1_pp"],
                "minco_minus_sylph_matched_abs_error_pp": component_deltas["matched_abs_error"],
                "minco_minus_sylph_missing_truth_mass_pp": component_deltas["missing_truth_mass"],
                "minco_minus_sylph_extra_pred_mass_pp": component_deltas["extra_pred_mass"],
                "minco_minus_sylph_pooled_F1": m["pooled_F1"] - s["pooled_F1"],
                "main_gap_component": max(component_deltas, key=lambda key: abs(component_deltas[key])),
            }
        )
    return pd.DataFrame(rows)


def top_species(species_rows: pd.DataFrame, n: int = 25) -> pd.DataFrame:
    parts = []
    for _keys, sub in species_rows.groupby(["panel", "sample", "method"], sort=True):
        top = sub.sort_values("abs_error_pp", ascending=False).head(n).copy()
        top["rank_abs_error"] = range(1, len(top) + 1)
        parts.append(top)
    return pd.concat(parts, ignore_index=True, sort=False) if parts else pd.DataFrame()


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    toy_mod = load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    hmp_taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
    hmp_gastro_taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)
    cami3_wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _map_diag = cami3.build_transfer_maps()
    cami3_taxmap = cami3.parse_species_taxmap(cami3.TAXMAP)
    rows: list[dict[str, object]] = []
    species_rows: list[dict[str, object]] = []

    for panel, cfg in PANELS.items():
        for sample in cfg["samples"]:
            truth = load_truth(Path(cfg["truth"](sample)), int(sample), str(cfg["truth_abundance_col"]))
            minco_path = Path(cfg["minco"](sample))
            sylph_path = Path(cfg["sylph"](sample))
            if panel == "hmp_airskin_gtdb_source_abundance":
                hmp.ensure_sample_record(int(sample))
                hmp_paths = hmp.SAMPLES[int(sample)]
                best_ref_species, _diag = hmp.best_raw_ref_species_by_taxid(
                    {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                    hmp_taxmap,
                    by_accession,
                    by_core,
                )
                minco_collapsed, _extra = hmp.load_minco_predictions(
                    minco_path,
                    best_ref_species,
                    by_accession,
                    by_core,
                )
                minco_pred = pred_from_collapsed(minco_collapsed)
                minco_rule = "official HMP MinCO mapping: best reference accession, fallback raw best ref; summed per GTDB species; renormalized"
            elif panel == "hmp_gastrooral_gtdb_source_abundance":
                hmp_paths = hmp_gastro.SAMPLES[int(sample)]
                best_ref_species, _diag = hmp_gastro.best_raw_ref_species_by_taxid(
                    {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                    hmp_gastro_taxmap,
                    by_accession,
                    by_core,
                )
                minco_collapsed, _extra = hmp_gastro.load_minco_predictions(
                    minco_path,
                    best_ref_species,
                    by_accession,
                    by_core,
                )
                minco_pred = pred_from_collapsed(minco_collapsed)
                minco_rule = (
                    "official HMP gastrooral MinCO mapping: best reference "
                    "accession, fallback raw best ref; summed per GTDB species; "
                    "renormalized"
                )
            elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
                best_ref_species, _diag = cami3.best_raw_ref_species_by_taxid(
                    cami3.RAW_TABLES[int(sample)],
                    cami3_taxmap,
                    by_accession,
                    by_core,
                )
                minco_collapsed, _extra = cami3.load_minco_predictions(
                    minco_path,
                    cami3_taxid_to_species,
                    cami3_name_to_species,
                    best_ref_species,
                    by_accession,
                    by_core,
                )
                minco_pred = pred_from_collapsed(minco_collapsed)
                minco_rule = "official CAMI3 MinCO mapping: best reference accession, fallback raw best ref/taxid/name; summed per GTDB species; renormalized"
            else:
                minco_pred = load_minco_pred(minco_path, str(cfg["minco_collapse"]))
                minco_rule = f"calibrated MinCO called rows; species_name; {cfg['minco_collapse']} per GTDB species; renormalized"
            methods = [
                (
                    cfg["minco_method"],
                    minco_pred,
                    minco_path,
                    minco_rule,
                )
            ]
            if panel == "cami2_toy_mouse_gut":
                sylph_pred = load_toy_sylph_pred(
                    sylph_path,
                    toy_mod,
                    by_accession,
                    by_core,
                    str(cfg["sylph_collapse"]),
                )
                sylph_rule = (
                    f"Sylph rows mapped by accession; {cfg['sylph_collapse']} per "
                    "GTDB species; renormalized to match Toy Mouse scorer"
                )
            else:
                sylph_pred = load_sylph_pred(sylph_path, by_accession, by_core)
                sylph_rule = "Sylph rows mapped by accession; summed per GTDB species; renormalized"
            methods.append((cfg["sylph_method"], sylph_pred, sylph_path, sylph_rule))
            for method, pred, pred_path, rule in methods:
                row, detail = decompose(
                    panel,
                    str(sample),
                    str(method),
                    truth,
                    pred,
                    str(Path(cfg["truth"](sample))),
                    str(pred_path),
                    rule,
                )
                rows.append(row)
                species_rows.extend(detail)

    decomp = pd.DataFrame(rows)
    species = pd.DataFrame(species_rows)
    summary = summarize(decomp)
    delta = compare_methods(summary)
    top = top_species(species)
    validations = pd.concat(
        [official_validation(panel, decomp) for panel in PANELS],
        ignore_index=True,
        sort=False,
    )

    decomp.to_csv(RESULTS / "abundance_error_decomposition.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "abundance_error_decomposition_summary.tsv", sep="\t", index=False)
    delta.to_csv(RESULTS / "abundance_error_decomposition_delta.tsv", sep="\t", index=False)
    top.to_csv(RESULTS / "abundance_error_top_species.tsv", sep="\t", index=False)
    validations.to_csv(RESULTS / "abundance_error_decomposition_validation.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("\nDELTA")
    print(delta.to_string(index=False))
    print("\nVALIDATION")
    print(validations.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
