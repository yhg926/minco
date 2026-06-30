#!/usr/bin/env python3
"""Audit automatic output-derived switches for call-filter candidates.

The fixed output-only call-filter sweep found filters that improve mean F1 but
hurt some panels. This diagnostic asks whether a simple automatic switch, based
only on MinCO output columns and not dataset labels or truth, can choose when to
apply a stricter call filter.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import score_cami3_gtdb_source_readmap as cami3
import score_cami3_gtdb_taxid_transfer as taxid_score
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp
import sweep_cross_panel_call_filters as call_sweep


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

CURRENT = "current_calls"
MAX_METHODS = 12


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Audit automatic output-derived switches for call-filter candidates.",
    )
    ap.add_argument(
        "--include-hmp-omitted",
        action="store_true",
        help="Add restored HMP omitted samples 2, 8, and 26 to a separate audit output set.",
    )
    ap.add_argument(
        "--output-prefix",
        default="adaptive_call_filter_switch",
        help="Output prefix under results/. Use a non-default prefix for opt-in extensions.",
    )
    return ap.parse_args()


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def entropy(values: np.ndarray) -> float:
    values = np.maximum(values.astype(float), 0.0)
    total = float(values.sum())
    if total <= 0.0:
        return 0.0
    p = values / total
    p = p[p > 0.0]
    return float(-(p * np.log(p)).sum())


def quantile_thresholds(values: pd.Series) -> list[float]:
    clean = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if clean.empty or clean.nunique() <= 1:
        return []
    thresholds = []
    for q in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        value = float(clean.quantile(q))
        if math.isfinite(value):
            thresholds.append(value)
    return sorted(set(thresholds))


def official_l1(row: pd.Series) -> float:
    if str(row["panel"]) == "cami2_toy_mouse_gut":
        return finite(row["L1_truth_only_pp"])
    return finite(row["L1_union_pp"])


def official_pearson(row: pd.Series) -> float:
    if str(row["panel"]) == "cami2_toy_mouse_gut":
        return finite(row["Pearson_truth_only"])
    return finite(row["Pearson_union"])


def loaders() -> dict[str, object]:
    toy_mod = call_sweep.abundance.decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    _cami3_wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _map_diag = (
        cami3.build_transfer_maps()
    )
    return {
        "by_accession": by_accession,
        "by_core": by_core,
        "hmp_taxmap": hmp.parse_species_taxmap(hmp.TAXMAP),
        "hmp_gastro_taxmap": hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP),
        "cami3_taxid_to_species": cami3_taxid_to_species,
        "cami3_name_to_species": cami3_name_to_species,
        "cami3_taxmap": cami3.parse_species_taxmap(cami3.TAXMAP),
    }


def choose_methods(call_filter_prefix: str = "cross_panel_call_filter") -> list[str]:
    overall = pd.read_csv(RESULTS / f"{call_filter_prefix}_overall.tsv", sep="\t")
    candidates = overall.loc[~overall["method"].eq(CURRENT)].copy()
    candidates = candidates.loc[candidates["panel_count"].eq(overall["panel_count"].max())]
    candidates = candidates.sort_values(
        ["mean_delta_F1", "worsened_samples", "worsened_panels", "mean_delta_L1_pp", "method"],
        ascending=[False, True, True, True, True],
        kind="mergesort",
    )
    return list(candidates["method"].head(MAX_METHODS))


def base_features(
    panel: str,
    sample: int,
    profile: Path,
    calls: pd.DataFrame,
    filters: dict[str, np.ndarray],
) -> dict[str, object]:
    features = call_sweep.feature_arrays(calls)
    n = float(len(calls))
    abundance = features["calibrated_abundance"]
    abundance = np.maximum(abundance, 0.0)
    prob = features["prob"]
    max_xny = features["max_xny"]
    max_breadth = features["max_breadth"]
    max_real = features["max_real_min_af"]
    max_zip = features["max_zip_af"]
    max_depth = features["max_mean_depth"]
    max_cv = features["max_depth_cv"]
    total_abundance = float(abundance.sum())

    def frac(mask: np.ndarray) -> float:
        return float(mask.sum() / n) if n else 0.0

    def mass_frac(mask: np.ndarray) -> float:
        if total_abundance <= 0.0:
            return 0.0
        return float(abundance[mask].sum() / total_abundance)

    low_prob_breadth = (prob < 0.7) & (max_breadth < 0.2)
    low_prob_xny = (prob < 0.7) & (max_xny < 50)
    low_real = max_real < 0.05
    low_xny25 = max_xny < 25
    low_breadth05 = max_breadth < 0.05
    top_abundance = np.sort(abundance)[::-1]
    best_filter = filters.get("drop_prob_lt_0.7_and_breadth_lt_0.2")
    best_dropped = ~best_filter if best_filter is not None else low_prob_breadth
    return {
        "panel": panel,
        "sample": str(sample),
        "profile": str(profile),
        "n_calls": n,
        "abundance_entropy": entropy(abundance),
        "top1_abundance": float(top_abundance[0]) if top_abundance.size else 0.0,
        "top3_abundance": float(top_abundance[:3].sum()) if top_abundance.size else 0.0,
        "probability_mean": float(prob.mean()) if n else 0.0,
        "probability_median": float(np.median(prob)) if n else 0.0,
        "probability_min": float(prob.min()) if n else 0.0,
        "max_xny_median": float(np.median(max_xny)) if n else 0.0,
        "max_xny_p25": float(np.quantile(max_xny, 0.25)) if n else 0.0,
        "max_breadth_median": float(np.median(max_breadth)) if n else 0.0,
        "max_breadth_p25": float(np.quantile(max_breadth, 0.25)) if n else 0.0,
        "max_real_min_af_median": float(np.median(max_real)) if n else 0.0,
        "max_zip_af_median": float(np.median(max_zip)) if n else 0.0,
        "max_mean_depth_median": float(np.median(max_depth)) if n else 0.0,
        "max_depth_cv_mean": float(max_cv.mean()) if n else 0.0,
        "low_prob_breadth_frac": frac(low_prob_breadth),
        "low_prob_breadth_mass_frac": mass_frac(low_prob_breadth),
        "low_prob_xny_frac": frac(low_prob_xny),
        "low_prob_xny_mass_frac": mass_frac(low_prob_xny),
        "low_real_frac": frac(low_real),
        "low_xny25_frac": frac(low_xny25),
        "low_breadth05_frac": frac(low_breadth05),
        "best_mean_filter_drop_frac": frac(best_dropped),
        "best_mean_filter_drop_mass_frac": mass_frac(best_dropped),
    }


def load_score_and_feature_tables(
    methods: list[str],
    include_hmp_omitted: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    load = loaders()
    score_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    method_set = {CURRENT, *methods}

    for panel, sample in call_sweep.sample_records(include_hmp_omitted):
        truth, calls, collapse, path = call_sweep.load_panel_inputs(panel, sample, load)
        filters = call_sweep.threshold_filters(calls)
        feature_rows.append(base_features(panel, sample, path, calls, filters))
        base_abundance = call_sweep.arr(calls, "calibrated_abundance", 0.0)
        total_abundance = float(np.maximum(base_abundance, 0.0).sum())
        for method in sorted(method_set):
            keep = np.asarray(filters[method], dtype=bool)
            pred = call_sweep.abundance.collapse_prediction(
                calls.loc[keep].copy(),
                base_abundance[keep],
                collapse,
            )
            row = taxid_score.score_prediction(sample, method, pred, truth, {})
            row["panel"] = panel
            row["sample"] = str(sample)
            row["collapse_rule"] = collapse
            row["kept_rows"] = int(keep.sum())
            row["dropped_rows"] = int(len(keep) - keep.sum())
            row["dropped_fraction"] = float((len(keep) - keep.sum()) / len(keep)) if len(keep) else 0.0
            dropped_abundance = float(np.maximum(base_abundance[~keep], 0.0).sum())
            row["dropped_abundance_fraction"] = (
                dropped_abundance / total_abundance if total_abundance > 0.0 else 0.0
            )
            row["official_L1_pp"] = official_l1(pd.Series(row))
            row["official_Pearson"] = official_pearson(pd.Series(row))
            score_rows.append(row)

    scores = pd.DataFrame(score_rows)
    current = scores.loc[scores["method"].eq(CURRENT)][
        ["panel", "sample", "F1", "official_L1_pp", "official_Pearson"]
    ].rename(
        columns={
            "F1": "current_F1",
            "official_L1_pp": "current_official_L1_pp",
            "official_Pearson": "current_official_Pearson",
        }
    )
    scores = scores.merge(current, on=["panel", "sample"], how="left")
    scores["delta_F1"] = scores["F1"] - scores["current_F1"]
    scores["delta_official_L1_pp"] = scores["official_L1_pp"] - scores["current_official_L1_pp"]
    scores["delta_official_Pearson"] = scores["official_Pearson"] - scores["current_official_Pearson"]
    return scores, pd.DataFrame(feature_rows)


def feature_table_for_method(
    base: pd.DataFrame,
    scores: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    method_features = scores.loc[scores["method"].eq(method)][
        [
            "panel",
            "sample",
            "kept_rows",
            "dropped_rows",
            "dropped_fraction",
            "dropped_abundance_fraction",
        ]
    ].copy()
    method_features["kept_fraction"] = 1.0 - method_features["dropped_fraction"]
    method_features = method_features.rename(
        columns={
            "kept_rows": "method_kept_rows",
            "dropped_rows": "method_dropped_rows",
            "dropped_fraction": "method_dropped_fraction",
            "dropped_abundance_fraction": "method_dropped_abundance_fraction",
            "kept_fraction": "method_kept_fraction",
        }
    )
    return base.merge(method_features, on=["panel", "sample"], how="left")


def sample_outcomes_for_rule(
    scores: pd.DataFrame,
    features: pd.DataFrame,
    method: str,
    feature: str,
    op: str,
    threshold: float,
) -> pd.DataFrame:
    condition = (
        features[feature].astype(float) >= threshold
        if op == ">="
        else features[feature].astype(float) <= threshold
    )
    switched = set(
        zip(
            features.loc[condition, "panel"].astype(str),
            features.loc[condition, "sample"].astype(str),
        )
    )
    rows = []
    sub = scores.loc[scores["method"].isin([CURRENT, method])]
    for row in sub.itertuples(index=False):
        key = (str(getattr(row, "panel")), str(getattr(row, "sample")))
        use_variant = key in switched
        if use_variant != (getattr(row, "method") == method):
            continue
        rows.append(
            {
                "panel": key[0],
                "sample": key[1],
                "method_used": getattr(row, "method"),
                "switched": use_variant,
                "F1": finite(getattr(row, "F1")),
                "current_F1": finite(getattr(row, "current_F1")),
                "delta_F1": finite(getattr(row, "delta_F1")),
                "official_L1_pp": finite(getattr(row, "official_L1_pp")),
                "current_official_L1_pp": finite(getattr(row, "current_official_L1_pp")),
                "delta_official_L1_pp": finite(getattr(row, "delta_official_L1_pp")),
                "official_Pearson": finite(getattr(row, "official_Pearson")),
                "current_official_Pearson": finite(getattr(row, "current_official_Pearson")),
                "delta_official_Pearson": finite(getattr(row, "delta_official_Pearson")),
            }
        )
    return pd.DataFrame(rows)


def panel_stats(rows: pd.DataFrame) -> dict[str, float | int]:
    panel_delta = (
        rows.groupby("panel", as_index=False)
        .agg(
            panel_delta_F1=("delta_F1", "mean"),
            panel_delta_L1_pp=("delta_official_L1_pp", "mean"),
            panel_delta_Pearson=("delta_official_Pearson", "mean"),
        )
    )
    return {
        "sample_count": int(rows.shape[0]),
        "mean_delta_F1": float(rows["delta_F1"].mean()),
        "min_delta_F1": float(rows["delta_F1"].min()),
        "improved_samples": int((rows["delta_F1"] > 1e-9).sum()),
        "worsened_samples": int((rows["delta_F1"] < -1e-9).sum()),
        "mean_delta_L1_pp": float(rows["delta_official_L1_pp"].mean()),
        "max_worse_L1_pp": float(rows["delta_official_L1_pp"].max()),
        "mean_delta_Pearson": float(rows["delta_official_Pearson"].mean()),
        "mean_panel_delta_F1": float(panel_delta["panel_delta_F1"].mean()),
        "min_panel_delta_F1": float(panel_delta["panel_delta_F1"].min()),
        "improved_panels": int((panel_delta["panel_delta_F1"] > 1e-9).sum()),
        "worsened_panels": int((panel_delta["panel_delta_F1"] < -1e-9).sum()),
        "mean_panel_delta_L1_pp": float(panel_delta["panel_delta_L1_pp"].mean()),
        "max_panel_worse_L1_pp": float(panel_delta["panel_delta_L1_pp"].max()),
        "mean_panel_delta_Pearson": float(panel_delta["panel_delta_Pearson"].mean()),
    }


def evaluate_rule(
    scores: pd.DataFrame,
    base_features: pd.DataFrame,
    method: str,
    feature: str,
    op: str,
    threshold: float,
) -> dict[str, object]:
    features = feature_table_for_method(base_features, scores, method)
    outcomes = sample_outcomes_for_rule(scores, features, method, feature, op, threshold)
    stats = panel_stats(outcomes)
    stats.update(
        {
            "variant_method": method,
            "feature": feature,
            "op": op,
            "threshold": threshold,
            "switched_samples": int(outcomes["switched"].sum()),
        }
    )
    return stats


def build_rules(scores: pd.DataFrame, features: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    rows = []
    base_feature_cols = [
        col
        for col in features.columns
        if col not in {"panel", "sample", "profile"}
        and pd.api.types.is_numeric_dtype(features[col])
        and features[col].nunique(dropna=True) > 1
    ]
    method_feature_cols = [
        "method_kept_rows",
        "method_dropped_rows",
        "method_dropped_fraction",
        "method_dropped_abundance_fraction",
        "method_kept_fraction",
    ]
    panel_values = sorted(features["panel"].astype(str).unique())
    for method in methods:
        method_features = feature_table_for_method(features, scores, method)
        method_scores = scores.loc[scores["method"].eq(method)][
            [
                "panel",
                "sample",
                "delta_F1",
                "delta_official_L1_pp",
                "delta_official_Pearson",
            ]
        ].copy()
        data = method_features.merge(method_scores, on=["panel", "sample"], how="left")
        panel_arr = data["panel"].astype(str).to_numpy()
        delta_f1_method = data["delta_F1"].to_numpy(dtype=float)
        delta_l1_method = data["delta_official_L1_pp"].to_numpy(dtype=float)
        delta_pearson_method = data["delta_official_Pearson"].to_numpy(dtype=float)
        feature_cols = base_feature_cols + [
            col for col in method_feature_cols if method_features[col].nunique(dropna=True) > 1
        ]
        for feature in feature_cols:
            values = data[feature].astype(float).to_numpy()
            for threshold in quantile_thresholds(method_features[feature]):
                for op in [">=", "<="]:
                    condition = values >= threshold if op == ">=" else values <= threshold
                    delta_f1 = np.where(condition, delta_f1_method, 0.0)
                    delta_l1 = np.where(condition, delta_l1_method, 0.0)
                    delta_pearson = np.where(condition, delta_pearson_method, 0.0)
                    row: dict[str, object] = {
                        "variant_method": method,
                        "feature": feature,
                        "op": op,
                        "threshold": threshold,
                        "switched_samples": int(condition.sum()),
                        "sample_count": int(len(condition)),
                        "mean_delta_F1": float(delta_f1.mean()),
                        "min_delta_F1": float(delta_f1.min()),
                        "improved_samples": int((delta_f1 > 1e-9).sum()),
                        "worsened_samples": int((delta_f1 < -1e-9).sum()),
                        "mean_delta_L1_pp": float(delta_l1.mean()),
                        "max_worse_L1_pp": float(delta_l1.max()),
                        "mean_delta_Pearson": float(delta_pearson.mean()),
                    }
                    panel_delta_f1 = []
                    panel_delta_l1 = []
                    panel_delta_pearson = []
                    panel_improved = 0
                    panel_worsened = 0
                    for panel in panel_values:
                        mask = panel_arr == panel
                        pf1 = float(delta_f1[mask].mean())
                        pl1 = float(delta_l1[mask].mean())
                        pp = float(delta_pearson[mask].mean())
                        panel_delta_f1.append(pf1)
                        panel_delta_l1.append(pl1)
                        panel_delta_pearson.append(pp)
                        if pf1 > 1e-9:
                            panel_improved += 1
                        elif pf1 < -1e-9:
                            panel_worsened += 1
                        slug = panel.replace("-", "_")
                        row[f"{slug}_delta_F1"] = pf1
                        row[f"{slug}_delta_L1_pp"] = pl1
                        row[f"{slug}_delta_Pearson"] = pp
                        row[f"{slug}_switched_samples"] = int(condition[mask].sum())
                        row[f"{slug}_improved_samples"] = int((delta_f1[mask] > 1e-9).sum())
                        row[f"{slug}_worsened_samples"] = int((delta_f1[mask] < -1e-9).sum())
                    row["mean_panel_delta_F1"] = float(np.mean(panel_delta_f1))
                    row["min_panel_delta_F1"] = float(np.min(panel_delta_f1))
                    row["improved_panels"] = panel_improved
                    row["worsened_panels"] = panel_worsened
                    row["mean_panel_delta_L1_pp"] = float(np.mean(panel_delta_l1))
                    row["max_panel_worse_L1_pp"] = float(np.max(panel_delta_l1))
                    row["mean_panel_delta_Pearson"] = float(np.mean(panel_delta_pearson))
                    rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["mean_delta_F1", "worsened_samples", "mean_delta_L1_pp", "variant_method"],
        ascending=[False, True, True, True],
        kind="mergesort",
    )


def leave_one_panel_out(
    rules: pd.DataFrame,
    scores: pd.DataFrame,
    features: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    panels = sorted(features["panel"].astype(str).unique())
    for holdout in panels:
        candidates = []
        for rule in rules.itertuples(index=False):
            rule_dict = rule._asdict()
            train_panels = [panel for panel in panels if panel != holdout]
            train_f1 = np.array(
                [finite(rule_dict[f"{panel.replace('-', '_')}_delta_F1"]) for panel in train_panels],
                dtype=float,
            )
            train_l1 = np.array(
                [finite(rule_dict[f"{panel.replace('-', '_')}_delta_L1_pp"]) for panel in train_panels],
                dtype=float,
            )
            train_worsened_samples = sum(
                int(rule_dict[f"{panel.replace('-', '_')}_worsened_samples"])
                for panel in train_panels
            )
            train_improved_samples = sum(
                int(rule_dict[f"{panel.replace('-', '_')}_improved_samples"])
                for panel in train_panels
            )
            holdout_slug = holdout.replace("-", "_")
            holdout_f1 = finite(rule_dict[f"{holdout_slug}_delta_F1"])
            holdout_l1 = finite(rule_dict[f"{holdout_slug}_delta_L1_pp"])
            holdout_pearson = finite(rule_dict[f"{holdout_slug}_delta_Pearson"])
            candidates.append(
                {
                    "holdout_panel": holdout,
                    "variant_method": str(getattr(rule, "variant_method")),
                    "feature": str(getattr(rule, "feature")),
                    "op": str(getattr(rule, "op")),
                    "threshold": finite(getattr(rule, "threshold")),
                    "train_sample_safe": (
                        train_worsened_samples == 0
                        and train_improved_samples > 0
                    ),
                    "train_panel_safe": (
                        int((train_f1 < -1e-9).sum()) == 0
                        and int((train_f1 > 1e-9).sum()) > 0
                    ),
                    "train_mean_delta_F1": float(train_f1.mean()),
                    "train_min_delta_F1": float(train_f1.min()),
                    "train_worsened_samples": train_worsened_samples,
                    "train_worsened_panels": int((train_f1 < -1e-9).sum()),
                    "train_mean_delta_L1_pp": float(train_l1.mean()),
                    "holdout_mean_delta_F1": holdout_f1,
                    "holdout_min_delta_F1": holdout_f1,
                    "holdout_improved_samples": int(holdout_f1 > 1e-9),
                    "holdout_worsened_samples": int(holdout_f1 < -1e-9),
                    "holdout_mean_delta_L1_pp": holdout_l1,
                    "holdout_mean_delta_Pearson": holdout_pearson,
                    "holdout_switched_samples": int(
                        getattr(rule, f"{holdout_slug}_switched_samples")
                    ),
                    "holdout_sample_count": int(
                        features.loc[features["panel"].astype(str).eq(holdout)].shape[0]
                    ),
                }
            )
        cand_df = pd.DataFrame(candidates)
        sample_safe = cand_df.loc[cand_df["train_sample_safe"]].copy()
        panel_safe = cand_df.loc[cand_df["train_panel_safe"]].copy()
        if not sample_safe.empty:
            pool = sample_safe
            pool_name = "strict_train_sample_safe"
        elif not panel_safe.empty:
            pool = panel_safe
            pool_name = "strict_train_panel_safe"
        else:
            pool = cand_df
            pool_name = "best_train_mean"
        selected = pool.sort_values(
            ["train_mean_delta_F1", "train_worsened_samples", "train_mean_delta_L1_pp"],
            ascending=[False, True, True],
            kind="mergesort",
        ).iloc[0].to_dict()
        selected["selection_pool"] = pool_name
        rows.append(selected)
    return pd.DataFrame(rows)


def audit_rows(
    rules: pd.DataFrame,
    lopo: pd.DataFrame,
    methods: list[str],
    output_prefix: str = "adaptive_call_filter_switch",
) -> list[dict[str, object]]:
    nontrivial = rules.loc[
        rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])
    ].copy()
    sample_safe = nontrivial.loc[
        nontrivial["worsened_samples"].eq(0) & nontrivial["improved_samples"].gt(0)
    ]
    panel_safe = nontrivial.loc[
        nontrivial["worsened_panels"].eq(0) & nontrivial["improved_panels"].gt(0)
    ]
    best_mean = nontrivial.sort_values(
        ["mean_delta_F1", "worsened_samples", "mean_delta_L1_pp"],
        ascending=[False, True, True],
        kind="mergesort",
    ).iloc[0]
    best_safe = None
    if not sample_safe.empty:
        best_safe = sample_safe.sort_values(
            ["mean_delta_F1", "mean_delta_L1_pp"],
            ascending=[False, True],
            kind="mergesort",
        ).iloc[0]
    lopo_mean = float(lopo["holdout_mean_delta_F1"].mean())
    lopo_l1 = float(lopo["holdout_mean_delta_L1_pp"].mean())
    lopo_worse = int((lopo["holdout_mean_delta_F1"] < -1e-9).sum())
    lopo_improved = int((lopo["holdout_mean_delta_F1"] > 1e-9).sum())
    if best_safe is not None and lopo_worse == 0 and lopo_improved > 0:
        promotion = "candidate_requires_raw_wrapper_validation"
        decision = "do_not_promote_without_independent_holdout"
    elif best_safe is not None:
        promotion = "reject_as_overfit_by_leave_one_panel_out"
        decision = "current_call_gate_remains_default"
    else:
        promotion = "reject_adaptive_output_call_filter_switch"
        decision = "current_call_gate_remains_default"
    rows = [
        {
            "metric": "candidate_methods",
            "value": len(methods),
            "evidence": ",".join(methods),
            "decision": "top_call_filters_by_mean_F1",
        },
        {
            "metric": "tested_switch_rules",
            "value": int(rules.shape[0]),
            "evidence": "output_feature_thresholds",
            "decision": "truth_used_only_for_scoring",
        },
        {
            "metric": "nontrivial_switch_rules",
            "value": int(nontrivial.shape[0]),
            "evidence": "0 < switched_samples < sample_count",
            "decision": "eligible_for_default_audit",
        },
        {
            "metric": "strict_sample_safe_switches",
            "value": int(sample_safe.shape[0]),
            "evidence": f"{output_prefix}_rules.tsv",
            "decision": "must_generalize_before_default",
        },
        {
            "metric": "strict_panel_safe_switches",
            "value": int(panel_safe.shape[0]),
            "evidence": f"{output_prefix}_rules.tsv",
            "decision": "must_generalize_before_default",
        },
        {
            "metric": "best_mean_f1_switch",
            "value": (
                f"{best_mean['variant_method']};{best_mean['feature']}{best_mean['op']}"
                f"{best_mean['threshold']:.6g};mean_delta_F1={best_mean['mean_delta_F1']:.6f};"
                f"worsened_samples={int(best_mean['worsened_samples'])};"
                f"mean_delta_L1_pp={best_mean['mean_delta_L1_pp']:.6f};"
                f"switched={int(best_mean['switched_samples'])}"
            ),
            "evidence": f"{output_prefix}_rules.tsv",
            "decision": "in_sample_tradeoff_check",
        },
        {
            "metric": "lopo_holdout_panels",
            "value": (
                f"mean_delta_F1={lopo_mean:.6f};improved={lopo_improved};"
                f"worsened={lopo_worse};mean_delta_L1_pp={lopo_l1:.6f}"
            ),
            "evidence": f"{output_prefix}_lopo.tsv",
            "decision": "leave_one_panel_out_validation",
        },
        {
            "metric": "promotion_decision",
            "value": promotion,
            "evidence": f"{output_prefix}_rules.tsv;{output_prefix}_lopo.tsv",
            "decision": decision,
        },
    ]
    if best_safe is not None:
        rows.insert(
            6,
            {
                "metric": "best_sample_safe_switch",
                "value": (
                    f"{best_safe['variant_method']};{best_safe['feature']}{best_safe['op']}"
                    f"{best_safe['threshold']:.6g};mean_delta_F1={best_safe['mean_delta_F1']:.6f};"
                    f"mean_delta_L1_pp={best_safe['mean_delta_L1_pp']:.6f};"
                    f"switched={int(best_safe['switched_samples'])}"
                ),
                "evidence": f"{output_prefix}_rules.tsv",
                "decision": "requires_lopo_validation",
            },
        )
    return rows


def compact_rules_for_artifact(rules: pd.DataFrame, lopo: pd.DataFrame) -> pd.DataFrame:
    nontrivial = rules.loc[
        rules["switched_samples"].gt(0) & rules["switched_samples"].lt(rules["sample_count"])
    ].copy()
    sample_safe = nontrivial.loc[
        nontrivial["worsened_samples"].eq(0) & nontrivial["improved_samples"].gt(0)
    ].sort_values(
        ["mean_delta_F1", "mean_delta_L1_pp"],
        ascending=[False, True],
        kind="mergesort",
    )
    panel_safe = nontrivial.loc[
        nontrivial["worsened_panels"].eq(0) & nontrivial["improved_panels"].gt(0)
    ].sort_values(
        ["mean_panel_delta_F1", "mean_delta_L1_pp"],
        ascending=[False, True],
        kind="mergesort",
    )
    selected = []
    key_cols = ["variant_method", "feature", "op", "threshold"]
    for item in lopo.itertuples(index=False):
        mask = (
            rules["variant_method"].eq(str(getattr(item, "variant_method")))
            & rules["feature"].eq(str(getattr(item, "feature")))
            & rules["op"].eq(str(getattr(item, "op")))
            & np.isclose(rules["threshold"].astype(float), finite(getattr(item, "threshold")))
        )
        selected.append(rules.loc[mask])
    frames = [rules.head(200), sample_safe.head(200), panel_safe.head(200), *selected]
    compact = pd.concat(frames, ignore_index=True)
    return compact.drop_duplicates(subset=key_cols).sort_values(
        ["mean_delta_F1", "worsened_samples", "mean_delta_L1_pp", "variant_method"],
        ascending=[False, True, True, True],
        kind="mergesort",
    )


def main() -> int:
    args = parse_args()
    output_prefix = args.output_prefix
    call_filter_prefix = "cross_panel_call_filter"
    if args.include_hmp_omitted:
        if output_prefix == "adaptive_call_filter_switch":
            output_prefix = "adaptive_call_filter_switch_with_hmp_omitted"
        call_filter_prefix = "cross_panel_call_filter_with_hmp_omitted"
    RESULTS.mkdir(parents=True, exist_ok=True)
    methods = choose_methods(call_filter_prefix)
    scores, features = load_score_and_feature_tables(methods, args.include_hmp_omitted)
    rules = build_rules(scores, features, methods)
    lopo = leave_one_panel_out(rules, scores, features)
    audit = audit_rows(rules, lopo, methods, output_prefix)

    scores.to_csv(RESULTS / f"{output_prefix}_scores.tsv", sep="\t", index=False)
    features.sort_values(["panel", "sample"], kind="mergesort").to_csv(
        RESULTS / f"{output_prefix}_features.tsv",
        sep="\t",
        index=False,
    )
    compact_rules_for_artifact(rules, lopo).to_csv(
        RESULTS / f"{output_prefix}_rules.tsv",
        sep="\t",
        index=False,
    )
    lopo.to_csv(RESULTS / f"{output_prefix}_lopo.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / f"{output_prefix}_audit.tsv",
        audit,
        ["metric", "value", "evidence", "decision"],
    )

    print(pd.DataFrame(audit).to_string(index=False))
    print("\nLOPO")
    print(lopo.to_string(index=False))
    print("\nTOP RULES")
    print(rules.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
