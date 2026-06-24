#!/usr/bin/env python3
"""Score S2000 ctx+obj markerdb on CAMI II Toy Mouse sample0."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(OLD_EXP))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402


RUN_DIR = Path("/tmp/gtdb232_s2000_ctxobj_marker_20260622")
CTXOBJ_TSV = RUN_DIR / "toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.tsv"
OLD_CTX_TSV = OLD_EXP / "toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv"
PROFILE_TSV = OLD_EXP / "mouse0_gtdb_species_profile.tsv"

EFFECTIVE_CTX_LENGTH = 24
ANI_THRESHOLD = 0.95
ACTIVE_CTX_MIN = 15.0
ACTIVE_RELIABLE_ZTP_AF_FLOOR = 0.40
ACTIVE_MEAN_DEPTH_MIN = 3.0
ACTIVE_VMR_MIN = 50.0
ACTIVE_DELTA_MAX = 0.03


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


def load_rows(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = truth.load_minco_with_gtdb_species(path, by_accession, by_core)
    rows = truth.add_naive_ani(rows)
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_depth_variance",
        "Normalized_abundance_depth",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
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
        (float(var) / float(mean)) if float(mean) > 0.0 else 0.0
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
        & rows["gtdb_species"].astype(bool)
    )
    return rows


def append_details(detail_rows, sample: str, kind: str, species_set: set[str], rows: pd.DataFrame, gold_abundance: dict[str, float]) -> None:
    for species in sorted(species_set):
        rec = {
            "method": sample,
            "kind": kind,
            "gtdb_species": species,
            "gold_abundance": gold_abundance.get(species, ""),
            "accession": "",
            "Ref": "",
            "ANI_naive_calc": "",
            "XnY_ctx": "",
            "Ref_zip_af": "",
            "Reliable_Ref_breadth": "",
            "Reliable_ztp_af": "",
            "Reliable_Ref_hit_mean_depth": "",
            "Reliable_depth_vmr": "",
            "ANI_AF_delta": "",
            "Normalized_abundance_depth": "",
            "active_gate_pass": "",
        }
        examples = rows.loc[rows["gtdb_species"] == species].copy()
        if not examples.empty:
            examples["ANI_naive_calc"] = numeric(examples, "ANI_naive_calc")
            examples["XnY_ctx"] = numeric(examples, "XnY_ctx")
            best = examples.sort_values(["active_gate_pass", "ANI_naive_calc", "XnY_ctx"], ascending=False).iloc[0]
            for field in [
                "accession",
                "Ref",
                "ANI_naive_calc",
                "XnY_ctx",
                "Ref_zip_af",
                "Reliable_Ref_breadth",
                "Reliable_ztp_af",
                "Reliable_Ref_hit_mean_depth",
                "Reliable_depth_vmr",
                "ANI_AF_delta",
                "Normalized_abundance_depth",
                "active_gate_pass",
            ]:
                rec[field] = best.get(field, "")
        detail_rows.append(rec)


def abundance_metrics(selected: pd.DataFrame, profile: pd.DataFrame):
    pred = (
        selected.loc[selected["gtdb_species"].astype(bool)]
        .groupby("gtdb_species")["Normalized_abundance_depth"]
        .max()
        .to_dict()
    )
    gold = dict(zip(profile["gtdb_species"].astype(str), profile["relative_abundance"].astype(float)))
    rows = []
    for renorm in [False, True]:
        values = pred.copy()
        if renorm:
            total = sum(values.values())
            if total > 0.0:
                values = {k: v / total for k, v in values.items()}
        y_true = []
        y_pred = []
        for species, abundance in gold.items():
            y_true.append(float(abundance))
            y_pred.append(float(values.get(species, 0.0)))
        true_s = pd.Series(y_true, dtype=float)
        pred_s = pd.Series(y_pred, dtype=float)
        rows.append(
            {
                "renorm_pred": renorm,
                "truth_species": len(y_true),
                "pred_sum_on_truth": sum(y_pred),
                "pred_sum_all": sum(values.values()),
                "pearson": pred_s.corr(true_s, method="pearson"),
                "spearman": pred_s.corr(true_s, method="spearman"),
                "mae_pct_points": (pred_s - true_s).abs().mean() * 100.0,
                "l1_pct_points": (pred_s - true_s).abs().sum() * 100.0,
            }
        )
    return rows


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    profile = pd.read_csv(PROFILE_TSV, sep="\t")
    gold = set(profile["gtdb_species"].astype(str))
    gold_abundance = dict(zip(profile["gtdb_species"].astype(str), profile["relative_abundance"].astype(float)))

    score_rows = []
    detail_rows = []
    abundance_rows = []
    for method, path in [
        ("ctx_only_current_best_recheck", OLD_CTX_TSV),
        ("ctxobj_markerdb_product0_active", CTXOBJ_TSV),
    ]:
        rows = load_rows(path, by_accession, by_core)
        selected = rows.loc[rows["active_gate_pass"]].copy()
        tp, fp, fn, precision, recall, f1 = score_sets(selected["gtdb_species"], gold)
        score_rows.append(
            {
                "method": method,
                "pred_taxa": len(tp | fp),
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "precision": precision,
                "recall": recall,
                "F1": f1,
                "selected_rows": len(selected),
            }
        )
        append_details(detail_rows, method, "FP", fp, rows, gold_abundance)
        append_details(detail_rows, method, "FN", fn, rows, gold_abundance)
        for rec in abundance_metrics(selected, profile):
            abundance_rows.append({"method": method, **rec})

    score_path = RUN_DIR / "ctxobj_markerdb_gtdb_scores.tsv"
    detail_path = RUN_DIR / "ctxobj_markerdb_gtdb_details.tsv"
    abundance_path = RUN_DIR / "ctxobj_markerdb_gtdb_abundance.tsv"
    pd.DataFrame(score_rows).to_csv(score_path, sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(detail_path, sep="\t", index=False)
    pd.DataFrame(abundance_rows).to_csv(abundance_path, sep="\t", index=False)
    print(pd.DataFrame(score_rows).to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {detail_path}")
    print(f"wrote {abundance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
