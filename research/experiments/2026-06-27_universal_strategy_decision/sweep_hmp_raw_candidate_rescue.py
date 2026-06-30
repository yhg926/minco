#!/usr/bin/env python3
"""Sweep HMP raw-side-channel rescue rules.

The previous audit showed that high-abundance HMP false negatives missing from
emitted profile rows are present in raw unique/split/exact evidence. This
script tests whether simple raw-table thresholds can add those candidates
without damaging F1. It is diagnostic only: truth is used for scoring and no
wrapper default is changed here.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
TMP_SCORES = Path("/tmp/minco_hmp_raw_candidate_rescue_scores.tsv")
TMP_SUMMARY = Path("/tmp/minco_hmp_raw_candidate_rescue_panel_summary.tsv")
CACHE_DIR = Path("/tmp/minco_hmp_raw_candidate_cache")

PANELS = [
    "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance",
]

ANI_MINS = [0.90, 0.92, 0.93, 0.94, 0.95, 0.97]
XNY_MINS = [100, 300, 500, 700, 900]
BREADTH_MINS = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.50]
REAL_AF_MINS = [0.10, 0.30, 0.50, 0.70, 0.90]

RAW_ABUNDANCE_COLS = [
    "Normalized_effective_abundance_depth",
    "Effective_abundance_depth",
    "Normalized_abundance_depth",
    "Relative_abundance_depth",
    "Ref_mean_depth",
]

RAW_USECOLS = {
    "Ref",
    "Ref_annotation",
    "ANI",
    "Distance",
    "XnY_ctx",
    "Real_min_align_fraction",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_cv",
    "Ref_zero_fraction",
    "Read_match_fraction",
    "Block_match_fraction",
    "Relative_abundance_depth",
    "Normalized_abundance_depth",
    "Effective_abundance_depth",
    "Normalized_effective_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
    "Default_call",
    "Default_call_rule",
}

RAW_NUMERIC_COLS = sorted((set(raw_audit.RAW_NUMERIC_COLS) | set(RAW_ABUNDANCE_COLS)) - {"Default_call"})

PROFILE_USECOLS = {
    "calibrated_call",
    "calibrated_abundance",
    "u_best_accession",
    "s_best_accession",
}


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(key): max(0.0, finite_float(value)) for key, value in values.items() if str(key)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def raw_mass(row: object) -> float:
    for col in RAW_ABUNDANCE_COLS:
        value = finite_float(getattr(row, col, 0.0))
        if value > 0.0:
            return value
    return 0.0


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def score(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    call_species: set[str],
    pred_abundance: Mapping[str, float],
    rule: str,
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
    y_true = pd.Series([truth_norm.get(species, 0.0) for species in union], dtype=float)
    y_pred = pd.Series([pred_norm.get(species, 0.0) for species in union], dtype=float)
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
        "L1_union_pp": float((y_pred - y_true).abs().sum() * 100.0),
        "Pearson_union": y_pred.corr(y_true, method="pearson") if len(union) > 1 else float("nan"),
        "added_species": len(pred_species - set(pred_norm)),
        "rule": rule,
    }


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
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
                "mean_F1": sub["F1"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_pred_species": sub["pred_species"].mean(),
                "rule": sub["rule"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def compare_to_current(scores: pd.DataFrame) -> pd.DataFrame:
    current = scores.loc[scores["method"].eq("current_default")].copy()
    current_key = {
        (row.panel, int(row.sample)): row
        for row in current.itertuples(index=False)
    }
    rows = []
    for method, sub in scores.loc[~scores["method"].eq("current_default")].groupby("method", sort=True):
        deltas = []
        for row in sub.itertuples(index=False):
            base = current_key[(row.panel, int(row.sample))]
            deltas.append(
                {
                    "panel": row.panel,
                    "sample": int(row.sample),
                    "delta_F1": finite_float(row.F1) - finite_float(base.F1),
                    "delta_L1_union_pp": finite_float(row.L1_union_pp) - finite_float(base.L1_union_pp),
                    "delta_Pearson_union": finite_float(row.Pearson_union) - finite_float(base.Pearson_union),
                    "TP_delta": int(row.TP) - int(base.TP),
                    "FP_delta": int(row.FP) - int(base.FP),
                    "FN_delta": int(row.FN) - int(base.FN),
                }
            )
        delta_df = pd.DataFrame(deltas)
        rows.append(
            {
                "method": method,
                "mean_delta_F1": delta_df["delta_F1"].mean(),
                "min_delta_F1": delta_df["delta_F1"].min(),
                "mean_delta_L1_union_pp": delta_df["delta_L1_union_pp"].mean(),
                "mean_delta_Pearson_union": delta_df["delta_Pearson_union"].mean(),
                "total_TP_delta": int(delta_df["TP_delta"].sum()),
                "total_FP_delta": int(delta_df["FP_delta"].sum()),
                "total_FN_delta": int(delta_df["FN_delta"].sum()),
                "improved_samples_F1": int((delta_df["delta_F1"] > 1e-12).sum()),
                "worsened_samples_F1": int((delta_df["delta_F1"] < -1e-12).sum()),
                "unchanged_samples_F1": int((delta_df["delta_F1"].abs() <= 1e-12).sum()),
            }
        )
    return pd.DataFrame(rows)


def extract_accessions(df: pd.DataFrame) -> pd.Series:
    if "Ref" in df.columns:
        accession = df["Ref"].map(raw_audit.hmp.extract_accession).astype(str)
    else:
        accession = pd.Series([""] * len(df), index=df.index, dtype=str)
    if "Ref_annotation" in df.columns:
        missing = accession.eq("")
        if missing.any():
            accession.loc[missing] = df.loc[missing, "Ref_annotation"].map(
                raw_audit.hmp.extract_accession
            )
    return accession


def map_accessions(
    accessions: pd.Series,
    panel: str,
    raw_mapper: raw_audit.RawMapper,
) -> tuple[pd.Series, pd.Series]:
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        mapper = raw_audit.hmp_gastro.gtdb_from_accession
    else:
        mapper = raw_audit.hmp.gtdb_from_accession
    cache = {
        accession: mapper(accession, raw_mapper.by_accession, raw_mapper.by_core)
        for accession in sorted(set(accessions.astype(str)))
        if accession
    }
    species = accessions.map(lambda accession: cache.get(str(accession), ("", ""))[0])
    method = accessions.map(lambda accession: cache.get(str(accession), ("", ""))[1])
    return species.astype(str), method.astype(str)


def gtdb_from_accession(panel: str, accession: object, raw_mapper: raw_audit.RawMapper) -> str:
    accession_s = str(accession or "").strip()
    if not accession_s or accession_s.lower() == "nan":
        return ""
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        species, _method = raw_audit.hmp_gastro.gtdb_from_accession(
            accession_s, raw_mapper.by_accession, raw_mapper.by_core
        )
    else:
        species, _method = raw_audit.hmp.gtdb_from_accession(
            accession_s, raw_mapper.by_accession, raw_mapper.by_core
        )
    return str(species or "")


def load_current_pred_gtdb(
    panel: str,
    profile_path: Path,
    collapse: str,
    raw_mapper: raw_audit.RawMapper,
) -> dict[str, float]:
    df = pd.read_csv(
        profile_path,
        sep="\t",
        usecols=lambda col: col in PROFILE_USECOLS,
        low_memory=False,
    )
    called = df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = df.loc[called].copy()
    rows = []
    unmapped = 0
    for row in selected.itertuples(index=False):
        species = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species = gtdb_from_accession(panel, getattr(row, col), raw_mapper)
                if species:
                    break
        if not species:
            unmapped += 1
            continue
        rows.append(
            {
                "gtdb_species": species,
                "abundance": finite_float(getattr(row, "calibrated_abundance", 0.0)),
            }
        )
    if unmapped:
        log(f"warning: {profile_path} has {unmapped} called rows without GTDB accession mapping")
    if not rows:
        return {}
    work = pd.DataFrame(rows)
    if collapse == "max":
        grouped = work.groupby("gtdb_species")["abundance"].max()
    elif collapse == "sum":
        grouped = work.groupby("gtdb_species")["abundance"].sum()
    else:
        raise ValueError(f"unsupported collapse rule: {collapse}")
    return normalize({str(k): finite_float(v) for k, v in grouped.to_dict().items()})


def slim_best_raw_rows(
    panel: str,
    raw_paths: Mapping[str, Path],
    raw_mapper: raw_audit.RawMapper,
) -> pd.DataFrame:
    frames = []
    for mode, path in raw_paths.items():
        if not path.exists():
            continue
        log(f"loading {panel} {mode}: {path}")
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=lambda col: col in RAW_USECOLS,
            low_memory=False,
        )
        if df.empty:
            continue
        df["accession"] = extract_accessions(df)
        df["gtdb_species"], df["gtdb_mapping_method"] = map_accessions(
            df["accession"], panel, raw_mapper
        )
        df["raw_mode"] = mode
        df["raw_path"] = str(path)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    raw = pd.concat(frames, ignore_index=True, sort=False)
    for col in RAW_NUMERIC_COLS:
        if col in raw.columns:
            raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
        else:
            raw[col] = 0.0
    return raw_audit.best_raw_rows(raw)


def cache_path(panel: str, sample: int) -> Path:
    return CACHE_DIR / f"{panel}.sample{sample}.best_raw_rows.tsv"


def load_raw_best(
    panel: str,
    sample: int,
    profile_path: Path,
    raw_mapper: raw_audit.RawMapper,
) -> pd.DataFrame:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = cache_path(panel, sample)
    if cached.exists():
        log(f"using cache {cached}")
        return pd.read_csv(cached, sep="\t", low_memory=False)
    paths = raw_audit.raw_paths_from_profile(panel, sample, profile_path)
    raw_best = slim_best_raw_rows(panel, paths, raw_mapper)
    raw_best.to_csv(cached, sep="\t", index=False)
    log(f"cached {len(raw_best)} best raw rows at {cached}")
    return raw_best


def raw_abundance_for_calls(raw_best: pd.DataFrame, call_species: set[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    if raw_best.empty:
        return out
    for row in raw_best.itertuples(index=False):
        species = str(getattr(row, "gtdb_species", ""))
        if species in call_species:
            out[species] = raw_mass(row)
    return out


def prefilter_rescue_candidates(raw_best: pd.DataFrame, called: set[str]) -> pd.DataFrame:
    if raw_best.empty:
        return raw_best
    mask = (
        raw_best["gtdb_species"].astype(str).astype(bool)
        & ~raw_best["gtdb_species"].astype(str).isin(called)
        & (pd.to_numeric(raw_best["ANI"], errors="coerce").fillna(0.0) >= min(ANI_MINS))
        & (pd.to_numeric(raw_best["XnY_ctx"], errors="coerce").fillna(0.0) >= min(XNY_MINS))
        & (pd.to_numeric(raw_best["Ref_breadth"], errors="coerce").fillna(0.0) >= min(BREADTH_MINS))
        & (
            pd.to_numeric(raw_best["Real_min_align_fraction"], errors="coerce").fillna(0.0)
            >= min(REAL_AF_MINS)
        )
    )
    return raw_best.loc[mask].copy()


def candidates_for_rule(
    raw_candidates: pd.DataFrame,
    ani_min: float,
    xny_min: float,
    breadth_min: float,
    real_af_min: float,
) -> set[str]:
    if raw_candidates.empty:
        return set()
    mask = (
        (pd.to_numeric(raw_candidates["ANI"], errors="coerce").fillna(0.0) >= ani_min)
        & (pd.to_numeric(raw_candidates["XnY_ctx"], errors="coerce").fillna(0.0) >= xny_min)
        & (pd.to_numeric(raw_candidates["Ref_breadth"], errors="coerce").fillna(0.0) >= breadth_min)
        & (
            pd.to_numeric(raw_candidates["Real_min_align_fraction"], errors="coerce").fillna(0.0)
            >= real_af_min
        )
    )
    return set(raw_candidates.loc[mask, "gtdb_species"].astype(str))


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    log("initializing HMP raw accession mapper")
    raw_mapper = raw_audit.RawMapper()
    log("initialized HMP raw accession mapper")
    sample_payloads = []
    for panel in PANELS:
        cfg = decomp.PANELS[panel]
        for sample in cfg["samples"]:
            sample_id = int(sample)
            log(f"preparing {panel} sample{sample_id}")
            truth = decomp.load_truth(
                Path(cfg["truth"](sample_id)),
                sample_id,
                str(cfg["truth_abundance_col"]),
            )
            profile_path = Path(cfg["minco"](sample_id))
            current_pred = load_current_pred_gtdb(
                panel,
                profile_path,
                str(cfg["minco_collapse"]),
                raw_mapper,
            )
            raw_best = load_raw_best(panel, sample_id, profile_path, raw_mapper)
            current_calls = set(current_pred)
            raw_current_abund = raw_abundance_for_calls(raw_best, current_calls)
            raw_candidates = prefilter_rescue_candidates(raw_best, current_calls)
            log(
                f"candidate rows for {panel} sample{sample_id}: "
                f"{len(raw_candidates)} after loose rescue prefilter"
            )
            sample_payloads.append(
                {
                    "panel": panel,
                    "sample": sample_id,
                    "truth": truth,
                    "current_pred": current_pred,
                    "raw_candidates": raw_candidates,
                    "raw_current_abund": raw_current_abund,
                }
            )

    score_rows: list[dict[str, object]] = []
    for payload in sample_payloads:
        panel = str(payload["panel"])
        sample = int(payload["sample"])
        truth = payload["truth"]
        current_pred = payload["current_pred"]
        current_calls = set(current_pred)
        score_rows.append(
            score(
                panel,
                sample,
                "current_default",
                truth,
                current_calls,
                current_pred,
                "current calibrated calls and abundance",
            )
        )
        raw_current_abund = payload["raw_current_abund"]
        score_rows.append(
            score(
                panel,
                sample,
                "current_calls_raw_depth_abundance",
                truth,
                current_calls,
                raw_current_abund or current_pred,
                "current calls; raw best-row abundance where available",
            )
        )

    for ani_min in ANI_MINS:
        for xny_min in XNY_MINS:
            for breadth_min in BREADTH_MINS:
                for real_af_min in REAL_AF_MINS:
                    method_base = (
                        f"raw_rescue_ani{ani_min:g}_xny{xny_min:g}_"
                        f"br{breadth_min:g}_af{real_af_min:g}"
                    )
                    rule = (
                        f"raw best GTDB species not currently called; ANI>={ani_min:g};"
                        f"XnY>={xny_min:g};Ref_breadth>={breadth_min:g};"
                        f"Real_min_align_fraction>={real_af_min:g}"
                    )
                    for payload in sample_payloads:
                        panel = str(payload["panel"])
                        sample = int(payload["sample"])
                        truth = payload["truth"]
                        current_pred = payload["current_pred"]
                        raw_candidates = payload["raw_candidates"]
                        raw_current_abund = payload["raw_current_abund"]
                        current_calls = set(current_pred)
                        rescued = candidates_for_rule(
                            raw_candidates,
                            ani_min,
                            xny_min,
                            breadth_min,
                            real_af_min,
                        )
                        calls = current_calls | rescued
                        score_rows.append(
                            score(
                                panel,
                                sample,
                                method_base + "_zero_rescue_mass",
                                truth,
                                calls,
                                current_pred,
                                rule + "; rescued species have zero abundance mass",
                            )
                        )
                        raw_abund = dict(raw_current_abund)
                        raw_abund.update(raw_abundance_for_calls(raw_candidates, rescued))
                        if not raw_abund:
                            raw_abund = current_pred
                        score_rows.append(
                            score(
                                panel,
                                sample,
                                method_base + "_raw_depth_abundance",
                                truth,
                                calls,
                                raw_abund,
                                rule + "; all called species abundance from raw best-row depth",
                            )
                        )

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    deltas = compare_to_current(scores)
    overall = summary.groupby("method", as_index=False).agg(
        panels=("panel", lambda values: ",".join(sorted(set(map(str, values))))),
        mean_panel_pooled_F1=("pooled_F1", "mean"),
        mean_panel_L1_union_pp=("mean_L1_union_pp", "mean"),
        mean_panel_Pearson_union=("mean_Pearson_union", "mean"),
        total_TP=("pooled_TP", "sum"),
        total_FP=("pooled_FP", "sum"),
        total_FN=("pooled_FN", "sum"),
    )
    overall = overall.merge(deltas, on="method", how="left")
    noncurrent = overall.loc[~overall["method"].eq("current_default")].copy()
    best_f1 = noncurrent.sort_values(
        ["mean_delta_F1", "mean_delta_L1_union_pp", "total_FP_delta"],
        ascending=[False, True, True],
    ).head(1)
    safe = noncurrent.loc[
        (noncurrent["worsened_samples_F1"].fillna(0).astype(float).eq(0))
        & (noncurrent["mean_delta_F1"].fillna(0.0).astype(float) > 0.0)
    ].copy()
    best_safe = safe.sort_values(
        ["mean_delta_F1", "mean_delta_L1_union_pp", "total_FP_delta"],
        ascending=[False, True, True],
    ).head(1)
    audit_rows = [
        {
            "metric": "tested_rules",
            "value": len(ANI_MINS) * len(XNY_MINS) * len(BREADTH_MINS) * len(REAL_AF_MINS),
            "evidence": "rules are evaluated with zero-rescue-mass and raw-depth-abundance scoring",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "best_mean_F1_rule",
            "value": best_f1["method"].iloc[0] if not best_f1.empty else "",
            "evidence": ""
            if best_f1.empty
            else (
                f"mean_delta_F1={best_f1['mean_delta_F1'].iloc[0]};"
                f"mean_delta_L1={best_f1['mean_delta_L1_union_pp'].iloc[0]};"
                f"TP_delta={best_f1['total_TP_delta'].iloc[0]};"
                f"FP_delta={best_f1['total_FP_delta'].iloc[0]};"
                f"FN_delta={best_f1['total_FN_delta'].iloc[0]}"
            ),
            "decision": "best_diagnostic_rule",
        },
        {
            "metric": "sample_safe_positive_F1_rules",
            "value": int(len(safe)),
            "evidence": "rules with no per-sample F1 decrease and positive mean F1 delta versus current",
            "decision": "safety_screen",
        },
        {
            "metric": "best_sample_safe_rule",
            "value": best_safe["method"].iloc[0] if not best_safe.empty else "",
            "evidence": ""
            if best_safe.empty
            else (
                f"mean_delta_F1={best_safe['mean_delta_F1'].iloc[0]};"
                f"mean_delta_L1={best_safe['mean_delta_L1_union_pp'].iloc[0]};"
                f"TP_delta={best_safe['total_TP_delta'].iloc[0]};"
                f"FP_delta={best_safe['total_FP_delta'].iloc[0]};"
                f"FN_delta={best_safe['total_FN_delta'].iloc[0]}"
            ),
            "decision": "candidate_requires_non_hmp_validation" if not best_safe.empty else "no_safe_rule",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "HMP-only truth-aware raw-side-channel sweep; needs non-HMP validation and wrapper implementation",
            "decision": "do_not_promote_raw_rescue_yet",
        },
    ]
    audit = pd.DataFrame(audit_rows)
    top_overall = overall.sort_values(
        ["mean_delta_F1", "mean_delta_L1_union_pp", "total_FP_delta"],
        ascending=[False, True, True],
    ).head(100)
    scores.to_csv(TMP_SCORES, sep="\t", index=False)
    summary.to_csv(TMP_SUMMARY, sep="\t", index=False)
    overall.to_csv(RESULTS / "hmp_raw_candidate_rescue_overall.tsv", sep="\t", index=False)
    top_overall.to_csv(RESULTS / "hmp_raw_candidate_rescue_top100.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "hmp_raw_candidate_rescue_audit.tsv", sep="\t", index=False)
    print(audit.to_string(index=False))
    if not best_f1.empty:
        print("\nBEST")
        print(best_f1.to_string(index=False))
    if not best_safe.empty:
        print("\nBEST_SAMPLE_SAFE")
        print(best_safe.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
