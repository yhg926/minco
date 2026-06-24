#!/usr/bin/env python3
"""Compare MinCO abundance estimators against Sylph on GTDB Toy Mouse truth."""

from __future__ import annotations

import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(MORE_EXP))
sys.path.insert(0, str(OLD_EXP))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import score_more_toymouse_gtdb as more  # noqa: E402


SAMPLES = [0, 1, 2]
EPS = 1e-12


def s1000_full_output_path(sample: int) -> Path:
    return (
        EXP_DIR
        / f"toymouse_sample{sample}_s1000_full_split_naive_product_topfrac_median025.tsv"
    )


def s2000_dedup_full_output_path(sample: int) -> Path:
    return (
        EXP_DIR
        / f"toymouse_sample{sample}_s2000_dedup_full_split_naive_product_topfrac_median025.tsv"
    )


def safe_num(value: object) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return 0.0
    return x if math.isfinite(x) else 0.0


def pearson_spearman(y_true: list[float], y_pred: list[float]) -> tuple[float, float]:
    true_s = pd.Series(y_true, dtype=float)
    pred_s = pd.Series(y_pred, dtype=float)
    return (
        safe_num(pred_s.corr(true_s, method="pearson")),
        safe_num(pred_s.corr(true_s, method="spearman")),
    )


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    if total <= 0.0:
        return {k: 0.0 for k in values}
    return {k: (v / total if v > 0.0 else 0.0) for k, v in values.items()}


def estimate_metrics(
    sample: int,
    method: str,
    estimator: str,
    values: dict[str, float],
    gold: dict[str, float],
    renorm: bool,
) -> dict[str, object]:
    pred = normalize(values) if renorm else dict(values)
    species = sorted(gold)
    y_true = [gold[s] for s in species]
    y_pred = [pred.get(s, 0.0) for s in species]
    pearson, spearman = pearson_spearman(y_true, y_pred)
    pred_species = {s for s, v in pred.items() if v > 0.0}
    gold_species = set(gold)
    tp = pred_species & gold_species
    fp = pred_species - gold_species
    fn = gold_species - pred_species
    return {
        "sample": sample,
        "method": method,
        "estimator": estimator,
        "renorm_pred": renorm,
        "truth_species": len(gold_species),
        "pred_species": len(pred_species),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "missing_truth_mass": sum(gold[s] for s in fn),
        "fp_pred_mass": sum(pred.get(s, 0.0) for s in fp),
        "pred_sum_on_truth": sum(y_pred),
        "pred_sum_all": sum(pred.values()),
        "pearson": pearson,
        "spearman": spearman,
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) / len(species) * 100.0
        if species
        else 0.0,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) * 100.0,
    }


def collapse_rows(
    rows: pd.DataFrame,
    value_fn: Callable[[pd.Series], float],
    collapse: str,
) -> dict[str, float]:
    values: dict[str, float] = defaultdict(float)
    best: dict[str, float] = {}
    for _, row in rows.iterrows():
        species = str(row.get("gtdb_species", ""))
        if not species:
            continue
        value = value_fn(row)
        if not math.isfinite(value) or value <= 0.0:
            continue
        if collapse == "sum":
            values[species] += value
        elif collapse == "max":
            best[species] = max(best.get(species, 0.0), value)
        else:
            raise ValueError(collapse)
    return dict(values if collapse == "sum" else best)


def ztp_lambda_from_pos_mean(mean_pos: float) -> float:
    mean_pos = safe_num(mean_pos)
    if mean_pos <= 1.0:
        return max(0.0, mean_pos - 1.0)
    lo = 1e-10
    hi = max(2.0, mean_pos * 2.0)

    def cond_mean(lam: float) -> float:
        return lam / (1.0 - math.exp(-lam))

    while cond_mean(hi) < mean_pos and hi < 1e6:
        hi *= 2.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if cond_mean(mid) < mean_pos:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def load_truth_profile(sample: int) -> dict[str, float]:
    path = MORE_EXP / f"mouse{sample}_gtdb_species_profile.tsv"
    rows = pd.read_csv(path, sep="\t")
    return {
        str(row["gtdb_species"]): float(row["relative_abundance"])
        for _, row in rows.iterrows()
        if str(row.get("gtdb_species", ""))
    }


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def load_genome_size() -> dict[str, float]:
    out: dict[str, float] = {}
    for path in truth.GTDB_METADATA:
        with open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                size = safe_num(row.get("genome_size", ""))
                for acc in [
                    truth.extract_accession(row.get("accession", "")),
                    truth.extract_accession(row.get("ncbi_genbank_assembly_accession", "")),
                ]:
                    if acc and size > 0.0:
                        out.setdefault(acc, size)
    return out


def add_size(rows: pd.DataFrame, genome_size: dict[str, float]) -> pd.DataFrame:
    rows = rows.copy()
    if "accession" not in rows.columns:
        rows["accession"] = rows["Ref"].map(truth.extract_accession)
    rows["genome_size"] = rows["accession"].map(genome_size).fillna(0.0)
    return rows


def minco_estimators(rows: pd.DataFrame) -> list[tuple[str, dict[str, float]]]:
    selected = rows.loc[rows["active_gate_pass"]].copy()
    selected = selected.loc[selected["gtdb_species"].astype(bool)].copy()
    estimators: list[tuple[str, dict[str, float]]] = []

    formulas: list[tuple[str, Callable[[pd.Series], float]]] = [
        ("current_normalized_depth", lambda r: safe_num(r.get("Normalized_abundance_depth"))),
        ("relative_marker_mean_depth", lambda r: safe_num(r.get("Relative_abundance_depth"))),
        ("ref_mean_depth", lambda r: safe_num(r.get("Ref_mean_depth"))),
        ("reliable_mean_depth", lambda r: safe_num(r.get("Reliable_Ref_mean_depth"))),
        ("ref_hit_mean_depth", lambda r: safe_num(r.get("Ref_hit_mean_depth"))),
        ("reliable_hit_mean_depth", lambda r: safe_num(r.get("Reliable_Ref_hit_mean_depth"))),
        (
            "ref_effcov_mean_over_zipaf",
            lambda r: safe_num(r.get("Ref_mean_depth")) / max(safe_num(r.get("Ref_zip_af")), EPS),
        ),
        (
            "reliable_effcov_mean_over_zipaf",
            lambda r: safe_num(r.get("Reliable_Ref_mean_depth"))
            / max(safe_num(r.get("Reliable_Ref_zip_af")), EPS),
        ),
        (
            "reliable_effcov_mean_over_ztpaf",
            lambda r: safe_num(r.get("Reliable_Ref_mean_depth"))
            / max(safe_num(r.get("Reliable_ztp_af")), EPS),
        ),
        (
            "ztp_lambda_from_reliable_hit_mean",
            lambda r: ztp_lambda_from_pos_mean(safe_num(r.get("Reliable_Ref_hit_mean_depth"))),
        ),
        (
            "reliable_effcov_zipaf_size",
            lambda r: safe_num(r.get("Reliable_Ref_mean_depth"))
            / max(safe_num(r.get("Reliable_Ref_zip_af")), EPS)
            * max(safe_num(r.get("genome_size")), 1.0),
        ),
    ]
    for name, fn in formulas:
        for collapse in ["max", "sum"]:
            estimators.append((f"{name}_{collapse}", collapse_rows(selected, fn, collapse)))

    true_species = set(selected.loc[selected["gtdb_species"].astype(bool), "gtdb_species"].astype(str))
    estimators.append(("detected_truth_oracle", {s: 1.0 for s in true_species}))
    return estimators


def sylph_estimators(rows: pd.DataFrame) -> list[tuple[str, dict[str, float]]]:
    selected = rows.loc[rows["gtdb_species"].astype(bool)].copy()
    formulas: list[tuple[str, Callable[[pd.Series], float]]] = [
        ("reported_taxonomic_abundance", lambda r: safe_num(r.get("Taxonomic_abundance")) / 100.0),
        ("effcov", lambda r: safe_num(r.get("Eff_cov"))),
        (
            "effcov_size",
            lambda r: safe_num(r.get("Eff_cov")) * max(safe_num(r.get("genome_size")), 1.0),
        ),
    ]
    out: list[tuple[str, dict[str, float]]] = []
    for name, fn in formulas:
        for collapse in ["max", "sum"]:
            out.append((f"{name}_{collapse}", collapse_rows(selected, fn, collapse)))
    return out


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    genome_size = load_genome_size()
    all_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        gold = load_truth_profile(sample)
        method_loaders = [
            (
                "minco_s1000_full_current_gate",
                s1000_full_output_path(sample),
                "minco",
            ),
            (
                "minco_s2000_dedup_full_current_gate",
                s2000_dedup_full_output_path(sample),
                "minco",
            ),
            (
                "minco_ctx_marker_current_gate",
                more.minco_output_path(sample, "ctx_only_current_best"),
                "minco",
            ),
            (
                "minco_ctxobj_marker_current_gate",
                more.minco_output_path(sample, "ctxobj_markerdb_product0_active"),
                "minco",
            ),
            ("sylph_gtdb_profile", more.sylph_profile_path(sample), "sylph"),
        ]
        for method, path, kind in method_loaders:
            if not path.exists():
                continue
            if kind == "minco":
                rows = more.load_minco_active_rows(path, by_accession, by_core)
                rows = add_size(rows, genome_size)
                estimates = minco_estimators(rows)
            else:
                rows = more.load_sylph_rows(path, by_accession, by_core)
                rows = add_size(rows, genome_size)
                estimates = sylph_estimators(rows)
            for estimator, values in estimates:
                for renorm in [False, True]:
                    all_rows.append(
                        estimate_metrics(sample, method, estimator, values, gold, renorm)
                    )

    metrics = pd.DataFrame(all_rows)
    metrics.to_csv(EXP_DIR / "abundance_estimator_metrics.tsv", sep="\t", index=False)
    mean_cols = [
        "TP",
        "FP",
        "FN",
        "missing_truth_mass",
        "fp_pred_mass",
        "pred_sum_on_truth",
        "pred_sum_all",
        "pearson",
        "spearman",
        "mae_pct_points",
        "l1_pct_points",
    ]
    mean = (
        metrics.groupby(["method", "estimator", "renorm_pred"], as_index=False)[mean_cols]
        .mean()
        .sort_values(["renorm_pred", "l1_pct_points", "mae_pct_points", "method", "estimator"])
    )
    mean.to_csv(EXP_DIR / "abundance_estimator_mean.tsv", sep="\t", index=False)
    best = mean.loc[mean["renorm_pred"]].head(40)
    print(best.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
