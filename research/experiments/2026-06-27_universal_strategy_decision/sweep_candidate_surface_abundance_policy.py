#!/usr/bin/env python3
"""Sweep abundance policies for integrated candidate-surface rows.

The integrated candidate surface currently appends zero-mass rows. This script
keeps the call set fixed, assigns mass only to those appended rows, normalizes
the profile, and scores the result on cached HMP source-abundance panels.
It is diagnostic only; truth is used only for scoring.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp
import sweep_hmp_raw_candidate_rescue as raw_rescue


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
APPLY = RESULTS / "candidate_surface_integrated_hmp_apply.tsv"
SCORES = RESULTS / "candidate_surface_abundance_policy_scores.tsv"
SUMMARY = RESULTS / "candidate_surface_abundance_policy_summary.tsv"
OVERALL = RESULTS / "candidate_surface_abundance_policy_overall.tsv"
AUDIT = RESULTS / "candidate_surface_abundance_policy_audit.tsv"

PANELS = [
    "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance",
]
PROFILE_USECOLS = {
    "calibrated_call",
    "calibrated_abundance",
    "candidate_surface_added",
    "candidate_surface_accession",
    "s_best_accession",
    "u_best_accession",
    "s_Normalized_abundance_depth_max",
    "u_Normalized_abundance_depth_max",
    "s_Ref_mean_depth_max",
    "u_Ref_mean_depth_max",
    "s_XnY_ctx_max",
    "u_XnY_ctx_max",
}
ALPHAS = [0.05, 0.10, 0.25, 0.50, 1.0, 2.0, 5.0, 10.0]
CAPS = [0.001, 0.003, 0.005, 0.01, 0.02, 0.05]


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def apply_rows() -> pd.DataFrame:
    if not APPLY.exists():
        raise SystemExit(f"missing integrated apply table: {APPLY}")
    rows = pd.read_csv(APPLY, sep="\t")
    return rows.loc[rows["status"].astype(str).eq("replayed")].copy()


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(key): max(0.0, finite(value)) for key, value in values.items() if str(key)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def row_species(panel: str, row: object, mapper: raw_audit.RawMapper) -> str:
    accessions = [
        getattr(row, "candidate_surface_accession", ""),
        getattr(row, "s_best_accession", ""),
        getattr(row, "u_best_accession", ""),
    ]
    for accession in accessions:
        species = raw_rescue.gtdb_from_accession(panel, accession, mapper)
        if species:
            return species
    return ""


def profile_rows(panel: str, path: Path, mapper: raw_audit.RawMapper) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in PROFILE_USECOLS, low_memory=False)
    for col in PROFILE_USECOLS - {"calibrated_call", "candidate_surface_added", "candidate_surface_accession", "s_best_accession", "u_best_accession"}:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
        else:
            df[col] = 0.0
    if "candidate_surface_added" not in df.columns:
        df["candidate_surface_added"] = False
    called = bool_series(df["calibrated_call"])
    work = df.loc[called].copy()
    work["candidate_surface_added_bool"] = bool_series(work["candidate_surface_added"])
    species = []
    for row in work.itertuples(index=False):
        species.append(row_species(panel, row, mapper))
    work["gtdb_species"] = species
    return work.loc[work["gtdb_species"].astype(str).astype(bool)].copy()


def collapsed_abundance(work: pd.DataFrame, masses: pd.Series) -> dict[str, float]:
    tmp = work[["gtdb_species"]].copy()
    tmp["mass"] = pd.to_numeric(masses, errors="coerce").fillna(0.0).clip(lower=0.0)
    grouped = tmp.groupby("gtdb_species", sort=False)["mass"].sum()
    return normalize(grouped.to_dict())


def candidate_mass(work: pd.DataFrame, mode: str, alpha: float = 1.0, cap: float | None = None) -> pd.Series:
    added = work["candidate_surface_added_bool"]
    out = pd.Series(0.0, index=work.index, dtype=float)
    if not added.any():
        return out
    if mode == "normalized_depth":
        base = pd.concat(
            [
                work["s_Normalized_abundance_depth_max"].astype(float),
                work["u_Normalized_abundance_depth_max"].astype(float),
            ],
            axis=1,
        ).max(axis=1)
    elif mode == "xny_share":
        xny = pd.concat(
            [
                work["s_XnY_ctx_max"].astype(float),
                work["u_XnY_ctx_max"].astype(float),
            ],
            axis=1,
        ).max(axis=1)
        denom = max(float(xny.sum()), 1.0)
        base = xny / denom
    elif mode == "depth_share":
        depth = pd.concat(
            [
                work["s_Ref_mean_depth_max"].astype(float),
                work["u_Ref_mean_depth_max"].astype(float),
            ],
            axis=1,
        ).max(axis=1)
        denom = max(float(depth.sum()), 1.0)
        base = depth / denom
    else:
        raise ValueError(mode)
    values = (base * alpha).clip(lower=0.0)
    if cap is not None:
        values = values.clip(upper=cap)
    out.loc[added] = values.loc[added]
    return out


def method_grid() -> list[dict[str, object]]:
    rows = [{"method": "zero_mass", "mode": "zero", "alpha": 0.0, "cap": ""}]
    for mode in ["normalized_depth", "xny_share", "depth_share"]:
        for alpha in ALPHAS:
            rows.append({"method": f"{mode}_alpha{alpha:g}", "mode": mode, "alpha": alpha, "cap": ""})
        for cap in CAPS:
            rows.append({"method": f"{mode}_cap{cap:g}", "mode": mode, "alpha": 1.0, "cap": cap})
        for alpha in [0.25, 0.50, 1.0, 2.0]:
            for cap in CAPS:
                rows.append(
                    {
                        "method": f"{mode}_alpha{alpha:g}_cap{cap:g}",
                        "mode": mode,
                        "alpha": alpha,
                        "cap": cap,
                    }
                )
    return rows


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
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
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
                    "delta_L1_union_pp": finite(row.L1_union_pp) - finite(base.L1_union_pp),
                    "delta_Pearson_union": finite(row.Pearson_union) - finite(base.Pearson_union),
                    "added_mass_before_norm": finite(row.added_mass_before_norm),
                }
            )
        delta_df = pd.DataFrame(deltas)
        rows.append(
            {
                "method": method,
                "mean_delta_L1_union_pp": delta_df["delta_L1_union_pp"].mean(),
                "max_worse_L1_union_pp": delta_df["delta_L1_union_pp"].max(),
                "mean_delta_Pearson_union": delta_df["delta_Pearson_union"].mean(),
                "improved_samples_L1": int((delta_df["delta_L1_union_pp"] < -1e-12).sum()),
                "worsened_samples_L1": int((delta_df["delta_L1_union_pp"] > 1e-12).sum()),
                "unchanged_samples_L1": int((delta_df["delta_L1_union_pp"].abs() <= 1e-12).sum()),
                "mean_added_mass_before_norm": delta_df["added_mass_before_norm"].mean(),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    mapper = raw_audit.RawMapper()
    profiles: list[dict[str, object]] = []
    for row in apply_rows().itertuples(index=False):
        panel = str(row.panel)
        if panel not in PANELS:
            continue
        sample = int(row.sample)
        path = Path(str(row.out))
        if not path.exists():
            print(f"warning: missing profile {path}", file=sys.stderr)
            continue
        cfg = decomp.PANELS[panel]
        profiles.append(
            {
                "panel": panel,
                "sample": sample,
                "truth": decomp.load_truth(Path(cfg["truth"](sample)), sample, str(cfg["truth_abundance_col"])),
                "work": profile_rows(panel, path, mapper),
            }
        )

    score_rows = []
    methods = method_grid()
    for payload in profiles:
        panel = str(payload["panel"])
        sample = int(payload["sample"])
        truth = payload["truth"]
        work = payload["work"]
        called = set(work["gtdb_species"].astype(str))
        base = work["calibrated_abundance"].astype(float)
        added = work["candidate_surface_added_bool"]
        for spec in methods:
            if spec["mode"] == "zero":
                masses = base.copy()
                added_mass = 0.0
            else:
                add = candidate_mass(
                    work,
                    str(spec["mode"]),
                    float(spec["alpha"]),
                    None if spec["cap"] == "" else float(spec["cap"]),
                )
                masses = base.copy()
                masses.loc[added] = add.loc[added]
                added_mass = float(add.loc[added].sum())
            pred = collapsed_abundance(work, masses)
            score = raw_rescue.score(
                panel,
                sample,
                str(spec["method"]),
                truth,
                called,
                pred,
                "candidate-surface abundance policy; fixed calls",
            )
            score["added_rows"] = int(added.sum())
            score["added_mass_before_norm"] = added_mass
            score_rows.append(score)

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    deltas = compare_to_zero(scores)
    overall = summary.groupby("method", as_index=False).agg(
        panels=("panel", lambda values: ",".join(sorted(set(map(str, values))))),
        mean_panel_pooled_F1=("pooled_F1", "mean"),
        mean_panel_L1_union_pp=("mean_L1_union_pp", "mean"),
        mean_panel_Pearson_union=("mean_Pearson_union", "mean"),
        mean_added_mass_before_norm=("mean_added_mass_before_norm", "mean"),
    )
    overall = overall.merge(deltas, on="method", how="left")
    if "mean_added_mass_before_norm_x" in overall.columns:
        overall = overall.rename(columns={"mean_added_mass_before_norm_x": "mean_panel_added_mass_before_norm"})
    if "mean_added_mass_before_norm_y" in overall.columns:
        overall = overall.rename(columns={"mean_added_mass_before_norm_y": "mean_sample_added_mass_before_norm"})
    nonzero = overall.loc[~overall["method"].eq("zero_mass")].copy()
    best = nonzero.sort_values(
        ["mean_delta_L1_union_pp", "worsened_samples_L1", "mean_delta_Pearson_union"],
        ascending=[True, True, False],
    ).head(1)
    sample_safe = nonzero.loc[nonzero["worsened_samples_L1"].fillna(0).astype(int).eq(0)].copy()
    best_safe = sample_safe.sort_values(
        ["mean_delta_L1_union_pp", "mean_delta_Pearson_union"],
        ascending=[True, False],
    ).head(1)
    audit_rows = [
        {
            "metric": "profiles_scored",
            "value": len(profiles),
            "evidence": str(APPLY),
            "decision": "cached_integrated_surface_abundance_sweep",
        },
        {
            "metric": "tested_nonzero_policies",
            "value": len(methods) - 1,
            "evidence": "normalized-depth, XnY-share, and depth-share candidate-only mass with alpha/cap grid",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "best_mean_L1_policy",
            "value": "" if best.empty else best["method"].iloc[0],
            "evidence": ""
            if best.empty
            else (
                f"mean_delta_L1={best['mean_delta_L1_union_pp'].iloc[0]};"
                f"worsened_samples={best['worsened_samples_L1'].iloc[0]};"
                f"mean_added_mass={best['mean_panel_added_mass_before_norm'].iloc[0]}"
            ),
            "decision": "best_hmp_diagnostic_policy",
        },
        {
            "metric": "sample_safe_nonzero_policies",
            "value": len(sample_safe),
            "evidence": "nonzero policies with no per-sample L1 increase versus zero-mass surface",
            "decision": "safety_screen",
        },
        {
            "metric": "best_sample_safe_policy",
            "value": "" if best_safe.empty else best_safe["method"].iloc[0],
            "evidence": ""
            if best_safe.empty
            else (
                f"mean_delta_L1={best_safe['mean_delta_L1_union_pp'].iloc[0]};"
                f"mean_added_mass={best_safe['mean_panel_added_mass_before_norm'].iloc[0]}"
            ),
            "decision": "candidate_requires_cross_panel_validation" if not best_safe.empty else "no_sample_safe_nonzero_policy",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "HMP-only fixed-call abundance sweep; no default change",
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
