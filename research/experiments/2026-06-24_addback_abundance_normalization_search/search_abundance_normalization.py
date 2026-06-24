#!/usr/bin/env python3
"""Search abundance normalization for MinCO marker plus ctx+obj add-back calls."""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
PREV = ROOT / "2026-06-23_threshold_combo_search_all_metrics/search_threshold_combos.py"
OUT = EXP / "results"


def load_prev_module():
    spec = importlib.util.spec_from_file_location("threshold_combo_prev", PREV)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import previous scorer: {PREV}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


prev = load_prev_module()


@dataclass(frozen=True)
class CallSetDef:
    name: str
    marker_combo: object
    ctxobj_combo: object | None
    add_abs_min: float
    add_rel_min: float | None
    include_ctxobj_overlap: bool = False


@dataclass(frozen=True)
class AbundanceRule:
    rule_id: str
    callset: str
    mode: str
    marker_power: float
    add_power: float
    marker_cap_mult: float
    add_cap_mult: float
    add_scale: float
    add_mass_fraction: float
    uniform_alpha: float
    floor_fraction: float


def combo(
    combo_id: str,
    ani_min: float,
    xny_min: float,
    af_mode: str,
    af_value: float,
    delta_mean_min: float,
    delta_vmr_min: float,
    delta_max: float,
    median_cutoff: float,
    abundance_af_exp: float,
    fallback_marker_cutoff: int,
    fallback_xny_min: float,
    fallback_scale: float,
    rescue_mode: str,
):
    return prev.Combo(
        combo_id,
        ani_min,
        xny_min,
        af_mode,
        af_value,
        delta_mean_min,
        delta_vmr_min,
        delta_max,
        median_cutoff,
        abundance_af_exp,
        fallback_marker_cutoff,
        fallback_xny_min,
        fallback_scale,
        rescue_mode,
    )


MARKER_L1 = combo("marker_l1", 0.96, 10.0, "fixed", 0.32, 4.0, 20.0, 0.025, 20.0, 1.0, 300, 700.0, 0.4, "moderate")
MARKER_F1 = combo("marker_f1", 0.95, 10.0, "fixed", 0.30, 3.0, 10.0, 0.030, 10.0, 1.0, 200, 700.0, 0.4, "strict")
CTXOBJ_F1 = combo("ctxobj_f1", 0.96, 8.0, "formula", 27.0, 2.0, 50.0, 0.030, 10.0, 1.05, 100, 700.0, 0.4, "moderate")
CTXOBJ_L1 = combo("ctxobj_l1", 0.955, 8.0, "formula", 24.0, 2.0, 20.0, 0.025, 20.0, 1.0, 500, 1000.0, 0.4, "moderate")
CTXOBJ_CAMI = combo("ctxobj_cami", 0.96, 5.0, "formula", 27.0, 5.0, 50.0, 0.015, 20.0, 1.05, 200, 500.0, 0.4, "moderate")


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0}


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 2:
        return float("nan")
    return float(pd.Series(x, dtype=float).corr(pd.Series(y, dtype=float), method="pearson"))


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) < 2:
        return float("nan")
    return float(pd.Series(x, dtype=float).corr(pd.Series(y, dtype=float), method="spearman"))


def transformed(raw: dict[str, float], power: float, cap_mult: float, floor_fraction: float) -> dict[str, float]:
    vals = {k: float(v) for k, v in raw.items() if float(v) > 0.0 and math.isfinite(float(v))}
    if not vals:
        return {}
    arr = np.array(list(vals.values()), dtype=float)
    if cap_mult > 0.0 and len(arr):
        cap = float(np.median(arr) * cap_mult)
        if cap > 0.0:
            vals = {k: min(v, cap) for k, v in vals.items()}
    if floor_fraction > 0.0 and len(arr):
        floor = float(np.median(arr) * floor_fraction)
        if floor > 0.0:
            vals = {k: v + floor for k, v in vals.items()}
    if power != 1.0:
        vals = {k: v**power for k, v in vals.items()}
    return vals


def build_values(marker_raw: dict[str, float], add_raw: dict[str, float], rule: AbundanceRule) -> dict[str, float]:
    marker = transformed(marker_raw, rule.marker_power, rule.marker_cap_mult, rule.floor_fraction)
    add = transformed(add_raw, rule.add_power, rule.add_cap_mult, rule.floor_fraction)
    add = {k: v * rule.add_scale for k, v in add.items()}

    if rule.mode == "linear":
        values = dict(marker)
        for key, value in add.items():
            if key not in values:
                values[key] = value
        pred = normalize(values)
    elif rule.mode == "linear_sum":
        values = dict(marker)
        for key, value in add.items():
            values[key] = values.get(key, 0.0) + value
        pred = normalize(values)
    elif rule.mode == "group_mass":
        marker_norm = normalize(marker)
        add_norm = normalize(add)
        if marker_norm and add_norm:
            beta = rule.add_mass_fraction
            pred = {k: (1.0 - beta) * v for k, v in marker_norm.items()}
            for key, value in add_norm.items():
                pred[key] = pred.get(key, 0.0) + beta * value
        else:
            pred = normalize({**marker, **add})
    else:
        raise ValueError(rule.mode)

    if rule.uniform_alpha > 0.0 and pred:
        uniform = 1.0 / len(pred)
        pred = {k: (1.0 - rule.uniform_alpha) * v + rule.uniform_alpha * uniform for k, v in pred.items()}
    return pred


def score_prediction(sample, pred: dict[str, float], method: str) -> dict[str, object]:
    truth = sample.truth
    pred_set = set(pred)
    truth_set = set(truth)
    tp = pred_set & truth_set
    fp = pred_set - truth_set
    fn = truth_set - pred_set
    precision = len(tp) / len(pred_set) if pred_set else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    keys = sorted(truth)
    y_true = [truth[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    selected_truth_mass = sum(truth[k] for k in tp)
    # The existing benchmark L1 is calculated over truth taxa only. FP abundance
    # is not a separate L1 row, but it steals normalization mass. If a perfect
    # estimator assigned exact truth mass to selected true species and all
    # remaining mass to non-truth predictions, this is the lower bound.
    oracle_l1 = (1.0 - selected_truth_mass) * 100.0
    return {
        "method": method,
        "dataset": sample.dataset,
        "sample_id": sample.sample_id,
        "truth_taxa": len(truth_set),
        "pred_taxa": len(pred_set),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "selected_truth_mass": selected_truth_mass,
        "oracle_l1_pct_points": oracle_l1,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) * 100.0,
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_pred, y_true)) / max(len(keys), 1) * 100.0,
        "pearson": pearson(y_pred, y_true),
        "spearman": spearman(y_pred, y_true),
    }


def raw_values_for_sample(sample, combo_obj):
    return prev.values_for_sample(sample, combo_obj)[0]


def make_callsets(marker_samples, ctxobj_samples) -> dict[tuple[str, str, int], tuple[dict[str, float], dict[str, float]]]:
    ctx_by_key = {(s.dataset, s.sample_id): s for s in ctxobj_samples}
    callsets = [
        CallSetDef("marker_l1_only", MARKER_L1, None, 0.0, None),
        CallSetDef("marker_f1_only", MARKER_F1, None, 0.0, None),
        CallSetDef("ctxobj_f1_only", CTXOBJ_F1, None, 0.0, None),
        CallSetDef("marker_l1_ctxobj_f1_all", MARKER_L1, CTXOBJ_F1, 0.0, None),
        CallSetDef("marker_l1_ctxobj_l1_all", MARKER_L1, CTXOBJ_L1, 0.0, None),
        CallSetDef("marker_l1_ctxobj_f1_rel0002", MARKER_L1, CTXOBJ_F1, 0.0, 0.002),
        CallSetDef("marker_f1_ctxobj_f1_all", MARKER_F1, CTXOBJ_F1, 0.0, None),
        CallSetDef("marker_l1_ctxobj_cami_all", MARKER_L1, CTXOBJ_CAMI, 0.0, None),
        CallSetDef("marker_l1_ctxobj_f1_blend", MARKER_L1, CTXOBJ_F1, 0.0, None, True),
        CallSetDef("marker_l1_ctxobj_l1_blend", MARKER_L1, CTXOBJ_L1, 0.0, None, True),
        CallSetDef("marker_l1_ctxobj_cami_blend", MARKER_L1, CTXOBJ_CAMI, 0.0, None, True),
        CallSetDef("marker_f1_ctxobj_f1_blend", MARKER_F1, CTXOBJ_F1, 0.0, None, True),
    ]
    raw: dict[tuple[str, str, int], tuple[dict[str, float], dict[str, float]]] = {}
    for cs in callsets:
        for sample in marker_samples:
            marker_raw = raw_values_for_sample(sample, cs.marker_combo)
            add_raw: dict[str, float] = {}
            if cs.ctxobj_combo is not None:
                ctx_sample = ctx_by_key[(sample.dataset, sample.sample_id)]
                ctx_raw = raw_values_for_sample(ctx_sample, cs.ctxobj_combo)
                marker_max = max(marker_raw.values()) if marker_raw else 0.0
                for key, value in ctx_raw.items():
                    if key in marker_raw and not cs.include_ctxobj_overlap:
                        continue
                    if value < cs.add_abs_min:
                        continue
                    if cs.add_rel_min is not None and (marker_max <= 0.0 or value < cs.add_rel_min * marker_max):
                        continue
                    add_raw[key] = value
            elif cs.name.startswith("ctxobj_"):
                ctx_sample = ctx_by_key[(sample.dataset, sample.sample_id)]
                marker_raw = raw_values_for_sample(ctx_sample, CTXOBJ_F1)
            raw[(cs.name, sample.dataset, sample.sample_id)] = (marker_raw, add_raw)
    return raw


def build_rules(callsets: list[str], max_rules: int | None) -> list[AbundanceRule]:
    rules: list[AbundanceRule] = []
    powers = [0.4, 0.5, 0.65, 0.8, 1.0, 1.2]
    caps = [0.0, 5.0, 10.0, 20.0, 50.0]
    add_scales = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    betas = [0.01, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30]
    uniform_alphas = [0.0, 0.01, 0.03, 0.05, 0.10]
    floors = [0.0, 0.01, 0.03, 0.05]

    for callset in callsets:
        for marker_power in powers:
            for add_power in [0.5, 0.8, 1.0, 1.2]:
                for marker_cap in caps:
                    for add_cap in [0.0, 5.0, 10.0, 20.0]:
                        for add_scale in add_scales:
                            for uniform_alpha in uniform_alphas:
                                for floor in floors:
                                    rules.append(
                                        AbundanceRule(
                                            f"{callset}_linear_mp{marker_power}_ap{add_power}_mc{marker_cap}_ac{add_cap}_as{add_scale}_u{uniform_alpha}_f{floor}",
                                            callset,
                                            "linear",
                                            marker_power,
                                            add_power,
                                            marker_cap,
                                            add_cap,
                                            add_scale,
                                            0.0,
                                            uniform_alpha,
                                            floor,
                                        )
                                    )
                                    rules.append(
                                        AbundanceRule(
                                            f"{callset}_sum_mp{marker_power}_ap{add_power}_mc{marker_cap}_ac{add_cap}_as{add_scale}_u{uniform_alpha}_f{floor}",
                                            callset,
                                            "linear_sum",
                                            marker_power,
                                            add_power,
                                            marker_cap,
                                            add_cap,
                                            add_scale,
                                            0.0,
                                            uniform_alpha,
                                            floor,
                                        )
                                    )
                        for beta in betas:
                            for uniform_alpha in uniform_alphas:
                                for floor in floors:
                                    rules.append(
                                        AbundanceRule(
                                            f"{callset}_group_mp{marker_power}_ap{add_power}_mc{marker_cap}_ac{add_cap}_b{beta}_u{uniform_alpha}_f{floor}",
                                            callset,
                                            "group_mass",
                                            marker_power,
                                            add_power,
                                            marker_cap,
                                            add_cap,
                                            1.0,
                                            beta,
                                            uniform_alpha,
                                            floor,
                                        )
                                    )
    if max_rules and len(rules) > max_rules:
        # Deterministic coarse sample across the whole grid, not just the first
        # callset/prefix.
        ordered = sorted(rules, key=lambda r: r.rule_id)
        idx = np.linspace(0, len(ordered) - 1, max_rules, dtype=int)
        rules = [ordered[int(i)] for i in idx]
    return rules


def summarize(sample_rows: pd.DataFrame, rule: AbundanceRule) -> dict[str, object]:
    rec = rule.__dict__.copy()
    for prefix, sub in [
        ("all", sample_rows),
        ("mouse", sample_rows.loc[sample_rows["dataset"] == "mouse_gtdb"]),
        ("cami3", sample_rows.loc[sample_rows["dataset"] == "cami3_ncbi"]),
    ]:
        rec[f"{prefix}_F1"] = float(sub["F1"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_L1"] = float(sub["l1_pct_points"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_oracle_L1"] = float(sub["oracle_l1_pct_points"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_Pearson"] = float(sub["pearson"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_TP"] = float(sub["TP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FP"] = float(sub["FP"].mean()) if len(sub) else float("nan")
        rec[f"{prefix}_FN"] = float(sub["FN"].mean()) if len(sub) else float("nan")
    return rec


def score_sylph_baseline() -> pd.DataFrame:
    path = ROOT / "2026-06-23_threshold_combo_search_all_metrics/results/sylph_baseline.tsv"
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rules", type=int, default=60000)
    args = parser.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    marker_samples = prev.load_mouse_samples() + prev.load_cami3_samples()
    ctxobj_samples = prev.load_mouse_ctxobj_samples() + prev.load_cami3_ctxobj_samples()
    raw = make_callsets(marker_samples, ctxobj_samples)
    callsets = sorted({key[0] for key in raw})
    rules = build_rules(callsets, args.max_rules)

    summary_rows: list[dict[str, object]] = []
    sample_rows_to_save: list[dict[str, object]] = []
    for i, rule in enumerate(rules, start=1):
        sample_rows: list[dict[str, object]] = []
        for sample in marker_samples:
            marker_raw, add_raw = raw[(rule.callset, sample.dataset, sample.sample_id)]
            pred = build_values(marker_raw, add_raw, rule)
            sample_rows.append(score_prediction(sample, pred, rule.rule_id))
        sample_df = pd.DataFrame(sample_rows)
        summary_rows.append(summarize(sample_df, rule))
        if i <= 100:
            sample_rows_to_save.extend(sample_rows)
        if i % 5000 == 0:
            print(f"scored {i}/{len(rules)} abundance rules", file=sys.stderr)

    summary = pd.DataFrame(summary_rows).sort_values(["all_L1", "all_F1"], ascending=[True, False])
    summary.to_csv(OUT / "abundance_rule_summary.tsv", sep="\t", index=False)
    summary.head(200).to_csv(OUT / "abundance_rule_top200_by_l1.tsv", sep="\t", index=False)
    summary.sort_values(["all_F1", "all_L1"], ascending=[False, True]).head(200).to_csv(
        OUT / "abundance_rule_top200_by_f1.tsv", sep="\t", index=False
    )
    pd.DataFrame(sample_rows_to_save).to_csv(OUT / "abundance_rule_sample_metrics_subset.tsv", sep="\t", index=False)

    oracle_rows = []
    for callset in callsets:
        for sample in marker_samples:
            marker_raw, add_raw = raw[(callset, sample.dataset, sample.sample_id)]
            pred = normalize({**marker_raw, **add_raw})
            oracle_rows.append(score_prediction(sample, pred, callset))
    oracle_df = pd.DataFrame(oracle_rows)
    oracle_summary = (
        oracle_df.groupby(["method", "dataset"], as_index=False)[
            ["F1", "l1_pct_points", "oracle_l1_pct_points", "selected_truth_mass", "TP", "FP", "FN"]
        ]
        .mean()
    )
    oracle_summary.to_csv(OUT / "callset_oracle_bounds.tsv", sep="\t", index=False)
    score_sylph_baseline().to_csv(OUT / "sylph_baseline.tsv", sep="\t", index=False)

    print(summary.head(30).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
