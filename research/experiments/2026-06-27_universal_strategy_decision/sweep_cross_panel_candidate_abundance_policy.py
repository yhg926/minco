#!/usr/bin/env python3
"""Stress-test nonzero abundance for the selected cross-panel candidate rule.

The selected candidate rule improves F1 when rescued rows carry zero mass. This
diagnostic keeps that call set fixed and assigns candidate-only abundance from
output/raw-table evidence, then scores Toy Mouse, HMP, and CAMI3 panels with
their existing official L1 convention. It does not change the default wrapper.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_missed_truth_candidates as missed
import decompose_abundance_errors as decomp
import sweep_cross_panel_candidate_rescue as cross
import sweep_hmp_raw_candidate_rescue as hmp_rescue


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
SCORES = RESULTS / "cross_panel_candidate_abundance_policy_scores.tsv"
SUMMARY = RESULTS / "cross_panel_candidate_abundance_policy_summary.tsv"
OVERALL = RESULTS / "cross_panel_candidate_abundance_policy_overall.tsv"
AUDIT = RESULTS / "cross_panel_candidate_abundance_policy_audit.tsv"

PANELS = cross.PANELS
BEST_ANI_MIN = 0.90
BEST_XNY_MIN = 100.0
BEST_BREADTH_MIN = 0.01
BEST_REAL_AF_MIN = 0.70
ALPHAS = [0.05, 0.10, 0.25, 0.50, 1.0, 2.0, 5.0, 10.0]
CAPS = [0.001, 0.003, 0.005, 0.01, 0.02, 0.05]


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(key): max(0.0, finite(value)) for key, value in values.items() if str(key)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([0.0] * len(df), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(0.0)


def selected_rule_mask(candidates: pd.DataFrame) -> pd.Series:
    if candidates.empty:
        return pd.Series([], dtype=bool)
    return (
        (numeric(candidates, "ANI") >= BEST_ANI_MIN)
        & (numeric(candidates, "XnY") >= BEST_XNY_MIN)
        & (numeric(candidates, "breadth") >= BEST_BREADTH_MIN)
        & (numeric(candidates, "real_af") >= BEST_REAL_AF_MIN)
    )


def hmp_candidates(panel: str, sample: int, current_species: set[str]) -> pd.DataFrame:
    cache = hmp_rescue.cache_path(panel, sample)
    if not cache.exists():
        raise SystemExit(f"missing HMP raw candidate cache: {cache}")
    raw = pd.read_csv(cache, sep="\t", low_memory=False)
    out = pd.DataFrame(
        {
            "gtdb_species": raw["gtdb_species"].astype(str),
            "ANI": numeric(raw, "ANI"),
            "XnY": numeric(raw, "XnY_ctx"),
            "breadth": numeric(raw, "Ref_breadth"),
            "real_af": numeric(raw, "Real_min_align_fraction"),
            "candidate_norm_depth": numeric(raw, "Normalized_abundance_depth"),
            "candidate_raw_abundance": numeric(raw, "Relative_abundance_depth"),
            "source": "hmp_raw_best_row_cache",
        }
    )
    out = out.loc[out["gtdb_species"].astype(str).astype(bool)]
    out = out.loc[~out["gtdb_species"].isin(current_species)].copy()
    return out


def emitted_pred_and_candidates(
    panel: str,
    sample: int,
    profile_path: Path,
    collapse: str,
    mapper: missed.CandidateMapper,
) -> tuple[dict[str, float], pd.DataFrame]:
    mapped = mapper.map_profile(panel, sample, profile_path)
    called = mapped.loc[
        mapped["called"] & mapped["mapped_gtdb_species"].astype(str).astype(bool)
    ].copy()
    pred = cross.collapse_abundance(called, "mapped_gtdb_species", "calibrated_abundance", collapse)
    best = missed.candidate_best_rows(mapped)
    if best.empty:
        return pred, pd.DataFrame()
    current_species = set(pred)
    candidates = best.loc[
        best["mapped_gtdb_species"].astype(str).astype(bool)
        & ~best["mapped_gtdb_species"].astype(str).isin(current_species)
    ].copy()
    out = pd.DataFrame(
        {
            "gtdb_species": candidates["mapped_gtdb_species"].astype(str),
            "ANI": numeric(candidates, "best_zip_aaf_ani"),
            "XnY": numeric(candidates, "best_XnY"),
            "breadth": numeric(candidates, "best_breadth"),
            "real_af": numeric(candidates, "best_real_min_af"),
            "candidate_norm_depth": pd.concat(
                [
                    numeric(candidates, "s_Normalized_abundance_depth_max"),
                    numeric(candidates, "u_Normalized_abundance_depth_max"),
                ],
                axis=1,
            ).max(axis=1),
            "candidate_raw_abundance": numeric(candidates, "calibrated_abundance_raw"),
            "source": "emitted_profile",
        }
    )
    return pred, out


def candidate_mass(candidates: pd.DataFrame, method: str, alpha: float, cap: float | None) -> dict[str, float]:
    if candidates.empty:
        return {}
    if method == "normalized_depth":
        base = numeric(candidates, "candidate_norm_depth")
    elif method == "raw_abundance_scaled":
        raw = numeric(candidates, "candidate_raw_abundance")
        denom = max(float(raw.sum()), 1.0)
        base = raw / denom
    elif method == "xny_share":
        xny = numeric(candidates, "XnY")
        denom = max(float(xny.sum()), 1.0)
        base = xny / denom
    else:
        raise ValueError(method)
    values = (base * alpha).clip(lower=0.0)
    if cap is not None:
        values = values.clip(upper=cap)
    work = candidates[["gtdb_species"]].copy()
    work["mass"] = values
    return work.groupby("gtdb_species", sort=False)["mass"].sum().to_dict()


def method_grid() -> list[dict[str, object]]:
    rows = [{"method": "zero_mass", "source": "zero", "alpha": 0.0, "cap": ""}]
    for source in ["normalized_depth", "raw_abundance_scaled", "xny_share"]:
        for alpha in ALPHAS:
            rows.append({"method": f"{source}_alpha{alpha:g}", "source": source, "alpha": alpha, "cap": ""})
        for cap in CAPS:
            rows.append({"method": f"{source}_cap{cap:g}", "source": source, "alpha": 1.0, "cap": cap})
        for alpha in [0.25, 0.50, 1.0, 2.0]:
            for cap in CAPS:
                rows.append(
                    {
                        "method": f"{source}_alpha{alpha:g}_cap{cap:g}",
                        "source": source,
                        "alpha": alpha,
                        "cap": cap,
                    }
                )
    return rows


def score_with_pred(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    pred_abundance: Mapping[str, float],
    call_species: set[str],
    added_mass: float,
) -> dict[str, object]:
    truth_norm = normalize(truth)
    pred_norm = normalize(pred_abundance)
    truth_species = set(truth_norm)
    pred_species = set(call_species)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species
    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    union = sorted(truth_species | pred_species)
    truth_order = sorted(truth_species)
    union_true = [truth_norm.get(species, 0.0) for species in union]
    union_pred = [pred_norm.get(species, 0.0) for species in union]
    truth_true = [truth_norm.get(species, 0.0) for species in truth_order]
    truth_pred = [pred_norm.get(species, 0.0) for species in truth_order]
    l1_union = sum(abs(left - right) for left, right in zip(union_pred, union_true)) * 100.0
    l1_truth = sum(abs(left - right) for left, right in zip(truth_pred, truth_true)) * 100.0
    l1_kind = cross.official_l1_kind(panel)
    return {
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
        "L1_union_pp": l1_union,
        "L1_truth_only_pp": l1_truth,
        "official_L1_kind": l1_kind,
        "official_L1_pp": l1_union if l1_kind == "union" else l1_truth,
        "Pearson_union": cross.pearson(union_pred, union_true),
        "Pearson_truth_only": cross.pearson(truth_pred, truth_true),
        "added_mass_before_norm": added_mass,
    }


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (panel, method), sub in scores.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_official_L1_pp": sub["official_L1_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_added_mass_before_norm": sub["added_mass_before_norm"].mean(),
            }
        )
    return pd.DataFrame(rows)


def compare_to_zero(scores: pd.DataFrame) -> pd.DataFrame:
    zero = {
        (str(row.panel), int(row.sample)): row
        for row in scores.loc[scores["method"].eq("zero_mass")].itertuples(index=False)
    }
    rows = []
    for method, sub in scores.loc[~scores["method"].eq("zero_mass")].groupby("method", sort=True):
        deltas = []
        for row in sub.itertuples(index=False):
            base = zero[(str(row.panel), int(row.sample))]
            deltas.append(
                {
                    "panel": row.panel,
                    "sample": int(row.sample),
                    "delta_L1_pp": finite(row.official_L1_pp) - finite(base.official_L1_pp),
                    "delta_Pearson_union": finite(row.Pearson_union) - finite(base.Pearson_union),
                    "added_mass_before_norm": finite(row.added_mass_before_norm),
                }
            )
        delta_df = pd.DataFrame(deltas)
        panel_delta = delta_df.groupby("panel")["delta_L1_pp"].mean()
        rows.append(
            {
                "method": method,
                "mean_delta_L1_pp": delta_df["delta_L1_pp"].mean(),
                "max_worse_L1_pp": delta_df["delta_L1_pp"].max(),
                "mean_delta_Pearson_union": delta_df["delta_Pearson_union"].mean(),
                "improved_samples_L1": int((delta_df["delta_L1_pp"] < -1e-12).sum()),
                "worsened_samples_L1": int((delta_df["delta_L1_pp"] > 1e-12).sum()),
                "unchanged_samples_L1": int((delta_df["delta_L1_pp"].abs() <= 1e-12).sum()),
                "improved_panels_L1": int((panel_delta < -1e-12).sum()),
                "worsened_panels_L1": int((panel_delta > 1e-12).sum()),
                "mean_added_mass_before_norm": delta_df["added_mass_before_norm"].mean(),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mapper = missed.CandidateMapper()
    payloads = []
    for panel in PANELS:
        cfg = decomp.PANELS[panel]
        for sample_raw in cfg["samples"]:
            sample = int(sample_raw)
            truth = decomp.load_truth(Path(cfg["truth"](sample)), sample, str(cfg["truth_abundance_col"]))
            profile = Path(cfg["minco"](sample))
            if panel in {"hmp_airskin_gtdb_source_abundance", "hmp_gastrooral_gtdb_source_abundance"}:
                current_pred = cross.load_hmp_current_pred(panel, profile, str(cfg["minco_collapse"]), mapper)
                candidates = hmp_candidates(panel, sample, set(current_pred))
            else:
                current_pred, candidates = emitted_pred_and_candidates(
                    panel, sample, profile, str(cfg["minco_collapse"]), mapper
                )
            selected = candidates.loc[selected_rule_mask(candidates)].copy()
            payloads.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "truth": truth,
                    "current_pred": current_pred,
                    "selected": selected,
                }
            )

    score_rows = []
    methods = method_grid()
    for payload in payloads:
        panel = str(payload["panel"])
        sample = int(payload["sample"])
        truth = payload["truth"]
        current_pred = dict(payload["current_pred"])
        selected = payload["selected"]
        rescued = set(selected["gtdb_species"].astype(str)) if not selected.empty else set()
        call_species = set(current_pred) | rescued
        for spec in methods:
            if spec["source"] == "zero":
                pred = current_pred
                added_mass = 0.0
            else:
                added = candidate_mass(
                    selected,
                    str(spec["source"]),
                    float(spec["alpha"]),
                    None if spec["cap"] == "" else float(spec["cap"]),
                )
                pred = dict(current_pred)
                pred.update({species: mass for species, mass in added.items() if species not in pred})
                added_mass = sum(added.values())
            score = score_with_pred(panel, sample, str(spec["method"]), truth, pred, call_species, added_mass)
            score["selected_candidates"] = len(rescued)
            score_rows.append(score)

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    deltas = compare_to_zero(scores)
    overall = summary.groupby("method", as_index=False).agg(
        panels=("panel", lambda values: ",".join(sorted(set(map(str, values))))),
        mean_panel_pooled_F1=("pooled_F1", "mean"),
        mean_panel_official_L1_pp=("mean_official_L1_pp", "mean"),
        mean_panel_Pearson_union=("mean_Pearson_union", "mean"),
        mean_panel_added_mass_before_norm=("mean_added_mass_before_norm", "mean"),
    )
    overall = overall.merge(deltas, on="method", how="left")
    candidates = overall.loc[~overall["method"].eq("zero_mass")].copy()
    best = candidates.sort_values(
        ["mean_delta_L1_pp", "worsened_samples_L1", "worsened_panels_L1"],
        ascending=[True, True, True],
    ).head(1)
    sample_safe = candidates.loc[candidates["worsened_samples_L1"].fillna(0).astype(int).eq(0)].copy()
    panel_safe = candidates.loc[candidates["worsened_panels_L1"].fillna(0).astype(int).eq(0)].copy()
    best_sample_safe = sample_safe.sort_values(["mean_delta_L1_pp", "mean_delta_Pearson_union"], ascending=[True, False]).head(1)
    audit_rows = [
        {
            "metric": "profiles_scored",
            "value": len(payloads),
            "evidence": ",".join(PANELS),
            "decision": "cached_cross_panel_candidate_abundance_sweep",
        },
        {
            "metric": "selected_call_rule",
            "value": "ANI>=0.90;XnY>=100;breadth>=0.01;real_af>=0.70",
            "evidence": "same rule as cross_panel_candidate_rescue best sample-safe zero-mass rule",
            "decision": "fixed_call_set",
        },
        {
            "metric": "tested_nonzero_policies",
            "value": len(methods) - 1,
            "evidence": "normalized-depth, raw-abundance-scaled, and XnY-share candidate-only mass with alpha/cap grid",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "best_mean_L1_policy",
            "value": "" if best.empty else best["method"].iloc[0],
            "evidence": ""
            if best.empty
            else (
                f"mean_delta_L1={best['mean_delta_L1_pp'].iloc[0]};"
                f"worsened_samples={best['worsened_samples_L1'].iloc[0]};"
                f"worsened_panels={best['worsened_panels_L1'].iloc[0]};"
                f"mean_added_mass={best['mean_panel_added_mass_before_norm'].iloc[0]}"
            ),
            "decision": "best_cross_panel_diagnostic_policy",
        },
        {
            "metric": "sample_safe_nonzero_policies",
            "value": len(sample_safe),
            "evidence": "nonzero policies with no per-sample official-L1 increase versus zero-mass selected call set",
            "decision": "safety_screen",
        },
        {
            "metric": "panel_safe_nonzero_policies",
            "value": len(panel_safe),
            "evidence": "nonzero policies with no panel-mean official-L1 increase versus zero-mass selected call set",
            "decision": "safety_screen",
        },
        {
            "metric": "best_sample_safe_policy",
            "value": "" if best_sample_safe.empty else best_sample_safe["method"].iloc[0],
            "evidence": ""
            if best_sample_safe.empty
            else (
                f"mean_delta_L1={best_sample_safe['mean_delta_L1_pp'].iloc[0]};"
                f"mean_added_mass={best_sample_safe['mean_panel_added_mass_before_norm'].iloc[0]}"
            ),
            "decision": "candidate_requires_wrapper_validation" if not best_sample_safe.empty else "no_sample_safe_nonzero_policy",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "cached fixed-call abundance sweep only; no wrapper default change",
            "decision": "do_not_promote",
        },
    ]

    scores.to_csv(SCORES, sep="\t", index=False)
    summary.to_csv(SUMMARY, sep="\t", index=False)
    overall.to_csv(OVERALL, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(AUDIT, sep="\t", index=False)
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
