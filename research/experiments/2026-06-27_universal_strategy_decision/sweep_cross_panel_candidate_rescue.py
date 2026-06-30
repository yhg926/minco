#!/usr/bin/env python3
"""Cross-panel zero-mass candidate rescue diagnostic.

This follows the HMP-only raw-side-channel rescue sweep with a stricter
question: does a simple support-threshold rescue remain useful when HMP raw
candidates are evaluated alongside Toy Mouse and CAMI3 emitted profile
candidates? Rescued species receive zero abundance mass, so this diagnostic is
about call-set/F1 recovery, not abundance replacement.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import pandas as pd

import audit_missed_truth_candidates as missed
import decompose_abundance_errors as decomp
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp
import sweep_hmp_raw_candidate_rescue as hmp_rescue


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
TMP_SCORES = Path("/tmp/minco_cross_panel_candidate_rescue_scores.tsv")
TMP_PANEL_SUMMARY = Path("/tmp/minco_cross_panel_candidate_rescue_panel_summary.tsv")

PANELS = [
    "cami2_toy_mouse_gut",
    "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance",
    "cami3_toy_human_gut_gtdb_source_readmap",
]

ANI_MINS = [0.90, 0.92, 0.93, 0.94, 0.95, 0.97]
XNY_MINS = [10, 25, 50, 100, 300, 500, 700]
BREADTH_MINS = [0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]
REAL_AF_MINS = [0.10, 0.30, 0.50, 0.70, 0.90]

HMP_RAW_USECOLS = {
    "gtdb_species",
    "ANI",
    "XnY_ctx",
    "Ref_breadth",
    "Real_min_align_fraction",
}

PROFILE_PRED_USECOLS = {
    "calibrated_call",
    "calibrated_abundance",
    "s_best_accession",
    "u_best_accession",
}


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {
        str(key): max(0.0, finite_float(value))
        for key, value in values.items()
        if str(key)
    }
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def collapse_abundance(df: pd.DataFrame, species_col: str, abundance_col: str, collapse: str) -> dict[str, float]:
    if df.empty:
        return {}
    work = df.loc[df[species_col].astype(str).astype(bool)].copy()
    work[abundance_col] = pd.to_numeric(work[abundance_col], errors="coerce").fillna(0.0)
    if collapse == "max":
        grouped = work.groupby(species_col)[abundance_col].max()
    elif collapse == "sum":
        grouped = work.groupby(species_col)[abundance_col].sum()
    else:
        raise ValueError(f"unsupported collapse rule: {collapse}")
    return normalize({str(key): finite_float(value) for key, value in grouped.to_dict().items()})


def hmp_accession_to_species(panel: str, accession: object, mapper: missed.CandidateMapper) -> str:
    accession_s = str(accession or "").strip()
    if not accession_s or accession_s.lower() == "nan":
        return ""
    if panel == "hmp_gastrooral_gtdb_source_abundance":
        species, _method = hmp_gastro.gtdb_from_accession(
            accession_s, mapper.by_accession, mapper.by_core
        )
    else:
        species, _method = hmp.gtdb_from_accession(
            accession_s, mapper.by_accession, mapper.by_core
        )
    return str(species or "")


def load_hmp_current_pred(
    panel: str,
    profile_path: Path,
    collapse: str,
    mapper: missed.CandidateMapper,
) -> dict[str, float]:
    df = pd.read_csv(
        profile_path,
        sep="\t",
        usecols=lambda col: col in PROFILE_PRED_USECOLS,
        low_memory=False,
    )
    called = df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    rows = []
    for row in df.loc[called].itertuples(index=False):
        species = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species = hmp_accession_to_species(panel, getattr(row, col), mapper)
                if species:
                    break
        if species:
            rows.append(
                {
                    "gtdb_species": species,
                    "abundance": finite_float(getattr(row, "calibrated_abundance", 0.0)),
                }
            )
    return collapse_abundance(pd.DataFrame(rows), "gtdb_species", "abundance", collapse)


def load_profile_pred_and_candidates(
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
    pred = collapse_abundance(called, "mapped_gtdb_species", "calibrated_abundance", collapse)
    best = missed.candidate_best_rows(mapped)
    if best.empty:
        return pred, pd.DataFrame(columns=["gtdb_species", "ANI", "XnY", "breadth", "real_af", "source"])
    current_species = set(pred)
    candidates = best.loc[
        best["mapped_gtdb_species"].astype(str).astype(bool)
        & ~best["mapped_gtdb_species"].astype(str).isin(current_species)
    ].copy()
    out = pd.DataFrame(
        {
            "gtdb_species": candidates["mapped_gtdb_species"].astype(str),
            "ANI": pd.to_numeric(candidates["best_zip_aaf_ani"], errors="coerce").fillna(0.0),
            "XnY": pd.to_numeric(candidates["best_XnY"], errors="coerce").fillna(0.0),
            "breadth": pd.to_numeric(candidates["best_breadth"], errors="coerce").fillna(0.0),
            "real_af": pd.to_numeric(candidates["best_real_min_af"], errors="coerce").fillna(0.0),
            "source": "emitted_profile",
        }
    )
    return pred, loose_prefilter(out)


def load_hmp_raw_candidates(panel: str, sample: int, current_species: set[str]) -> pd.DataFrame:
    cache = hmp_rescue.cache_path(panel, sample)
    if not cache.exists():
        raise SystemExit(
            f"missing HMP raw candidate cache {cache}; run sweep_hmp_raw_candidate_rescue.py first"
        )
    raw = pd.read_csv(cache, sep="\t", usecols=lambda col: col in HMP_RAW_USECOLS, low_memory=False)
    out = pd.DataFrame(
        {
            "gtdb_species": raw["gtdb_species"].astype(str),
            "ANI": pd.to_numeric(raw["ANI"], errors="coerce").fillna(0.0),
            "XnY": pd.to_numeric(raw["XnY_ctx"], errors="coerce").fillna(0.0),
            "breadth": pd.to_numeric(raw["Ref_breadth"], errors="coerce").fillna(0.0),
            "real_af": pd.to_numeric(raw["Real_min_align_fraction"], errors="coerce").fillna(0.0),
            "source": "hmp_raw_best_row_cache",
        }
    )
    out = out.loc[out["gtdb_species"].astype(str).astype(bool)]
    out = out.loc[~out["gtdb_species"].astype(str).isin(current_species)].copy()
    return loose_prefilter(out)


def loose_prefilter(candidates: pd.DataFrame) -> pd.DataFrame:
    if candidates.empty:
        return candidates
    mask = (
        (pd.to_numeric(candidates["ANI"], errors="coerce").fillna(0.0) >= min(ANI_MINS))
        & (pd.to_numeric(candidates["XnY"], errors="coerce").fillna(0.0) >= min(XNY_MINS))
        & (pd.to_numeric(candidates["breadth"], errors="coerce").fillna(0.0) >= min(BREADTH_MINS))
        & (pd.to_numeric(candidates["real_af"], errors="coerce").fillna(0.0) >= min(REAL_AF_MINS))
    )
    return candidates.loc[mask].copy()


def rescue_species(candidates: pd.DataFrame, ani_min: float, xny_min: float, breadth_min: float, real_af_min: float) -> set[str]:
    if candidates.empty:
        return set()
    mask = (
        (candidates["ANI"] >= ani_min)
        & (candidates["XnY"] >= xny_min)
        & (candidates["breadth"] >= breadth_min)
        & (candidates["real_af"] >= real_af_min)
    )
    return set(candidates.loc[mask, "gtdb_species"].astype(str))


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) < 2:
        return float("nan")
    return float(pd.Series(left, dtype=float).corr(pd.Series(right, dtype=float), method="pearson"))


def official_l1_kind(panel: str) -> str:
    cfg = decomp.PANELS[panel]
    scores = pd.read_csv(Path(cfg["official_scores"]), sep="\t", nrows=1)
    return "union" if "L1_union_pp" in scores.columns else "truth_only"


def official_l1_from_score(row: pd.Series) -> float:
    for col in ["L1_union_pp", "truth_only_L1_pp", "L1_truth_only_pp", "L1"]:
        if col in row.index:
            return finite_float(row[col])
    return float("nan")


def official_baseline(panel: str, sample: int) -> dict[str, float]:
    cfg = decomp.PANELS[panel]
    scores = pd.read_csv(Path(cfg["official_scores"]), sep="\t")
    method = str(cfg["official_minco_method"])
    hit = scores.loc[
        scores["method"].astype(str).eq(method)
        & pd.to_numeric(scores["sample"], errors="coerce").fillna(-1).astype(int).eq(sample)
    ]
    if hit.empty:
        return {"F1": float("nan"), "L1": float("nan")}
    row = hit.iloc[0]
    return {"F1": finite_float(row.get("F1", 0.0)), "L1": official_l1_from_score(row)}


def score_sample(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    current_pred: Mapping[str, float],
    rescued: set[str],
    rule: str,
) -> dict[str, object]:
    truth_norm = normalize(truth)
    current_norm = normalize(current_pred)
    truth_species = set(truth_norm)
    pred_species = set(current_norm) | set(rescued)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species
    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    union = sorted(truth_species | pred_species)
    truth_order = sorted(truth_species)
    union_true = [truth_norm.get(species, 0.0) for species in union]
    union_pred = [current_norm.get(species, 0.0) for species in union]
    truth_true = [truth_norm.get(species, 0.0) for species in truth_order]
    truth_pred = [current_norm.get(species, 0.0) for species in truth_order]
    l1_union = sum(abs(left - right) for left, right in zip(union_pred, union_true)) * 100.0
    l1_truth = sum(abs(left - right) for left, right in zip(truth_pred, truth_true)) * 100.0
    l1_kind = official_l1_kind(panel)
    return {
        "panel": panel,
        "sample": sample,
        "method": method,
        "truth_species": len(truth_species),
        "pred_species": len(pred_species),
        "rescued_species": len(rescued),
        "rescued_TP": len(set(rescued) & truth_species),
        "rescued_FP": len(set(rescued) - truth_species),
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
        "Pearson_union": pearson(union_pred, union_true),
        "Pearson_truth_only": pearson(truth_pred, truth_true),
        "rule": rule,
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
                "mean_F1": sub["F1"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_official_L1_pp": sub["official_L1_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_rescued_species": sub["rescued_species"].mean(),
                "rescued_TP": int(sub["rescued_TP"].sum()),
                "rescued_FP": int(sub["rescued_FP"].sum()),
                "rule": sub["rule"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def compare_to_current(scores: pd.DataFrame) -> pd.DataFrame:
    current = scores.loc[scores["method"].eq("current_default")]
    base = {
        (str(row.panel), int(row.sample)): row
        for row in current.itertuples(index=False)
    }
    rows = []
    for method, sub in scores.loc[~scores["method"].eq("current_default")].groupby("method", sort=True):
        deltas = []
        for row in sub.itertuples(index=False):
            cur = base[(str(row.panel), int(row.sample))]
            deltas.append(
                {
                    "panel": row.panel,
                    "sample": int(row.sample),
                    "delta_F1": finite_float(row.F1) - finite_float(cur.F1),
                    "delta_L1": finite_float(row.official_L1_pp) - finite_float(cur.official_L1_pp),
                    "TP_delta": int(row.TP) - int(cur.TP),
                    "FP_delta": int(row.FP) - int(cur.FP),
                    "FN_delta": int(row.FN) - int(cur.FN),
                }
            )
        delta_df = pd.DataFrame(deltas)
        rows.append(
            {
                "method": method,
                "mean_delta_F1": delta_df["delta_F1"].mean(),
                "min_delta_F1": delta_df["delta_F1"].min(),
                "mean_delta_L1_pp": delta_df["delta_L1"].mean(),
                "total_TP_delta": int(delta_df["TP_delta"].sum()),
                "total_FP_delta": int(delta_df["FP_delta"].sum()),
                "total_FN_delta": int(delta_df["FN_delta"].sum()),
                "improved_samples_F1": int((delta_df["delta_F1"] > 1e-12).sum()),
                "worsened_samples_F1": int((delta_df["delta_F1"] < -1e-12).sum()),
                "unchanged_samples_F1": int((delta_df["delta_F1"].abs() <= 1e-12).sum()),
                "worsened_panels_F1": int((delta_df.groupby("panel")["delta_F1"].mean() < -1e-12).sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    log("initializing cross-panel candidate mapper")
    mapper = missed.CandidateMapper()
    log("initialized cross-panel candidate mapper")
    payloads = []
    validation_rows = []
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
            if panel in {
                "hmp_airskin_gtdb_source_abundance",
                "hmp_gastrooral_gtdb_source_abundance",
            }:
                current_pred = load_hmp_current_pred(
                    panel, profile_path, str(cfg["minco_collapse"]), mapper
                )
                candidates = load_hmp_raw_candidates(panel, sample_id, set(current_pred))
            else:
                current_pred, candidates = load_profile_pred_and_candidates(
                    panel, sample_id, profile_path, str(cfg["minco_collapse"]), mapper
                )
            log(
                f"prepared {panel} sample{sample_id}: "
                f"current_species={len(current_pred)} candidate_rows={len(candidates)}"
            )
            current_score = score_sample(
                panel,
                sample_id,
                "current_default",
                truth,
                current_pred,
                set(),
                "current calibrated calls and abundance",
            )
            official = official_baseline(panel, sample_id)
            validation_rows.append(
                {
                    "panel": panel,
                    "sample": sample_id,
                    "current_F1": current_score["F1"],
                    "official_F1": official["F1"],
                    "delta_F1": finite_float(current_score["F1"]) - official["F1"],
                    "current_official_L1_pp": current_score["official_L1_pp"],
                    "official_L1_pp": official["L1"],
                    "delta_L1_pp": finite_float(current_score["official_L1_pp"]) - official["L1"],
                    "candidate_rows": int(len(candidates)),
                    "candidate_source": ",".join(sorted(set(candidates["source"].astype(str))))
                    if not candidates.empty
                    else "",
                }
            )
            payloads.append(
                {
                    "panel": panel,
                    "sample": sample_id,
                    "truth": truth,
                    "current_pred": current_pred,
                    "candidates": candidates,
                    "current_score": current_score,
                }
            )

    score_rows = [payload["current_score"] for payload in payloads]
    log(
        "scoring threshold grid: "
        f"{len(ANI_MINS) * len(XNY_MINS) * len(BREADTH_MINS) * len(REAL_AF_MINS)} rules "
        f"x {len(payloads)} samples"
    )
    for ani_min in ANI_MINS:
        log(f"scoring ANI_MIN={ani_min:g}")
        for xny_min in XNY_MINS:
            for breadth_min in BREADTH_MINS:
                for real_af_min in REAL_AF_MINS:
                    method = f"cross_rescue_ani{ani_min:g}_xny{xny_min:g}_br{breadth_min:g}_af{real_af_min:g}_zero_mass"
                    rule = (
                        f"uncalled candidate species; ANI>={ani_min:g};XnY>={xny_min:g};"
                        f"breadth>={breadth_min:g};real_af>={real_af_min:g};zero rescued abundance"
                    )
                    for payload in payloads:
                        rescued = rescue_species(
                            payload["candidates"],
                            ani_min,
                            xny_min,
                            breadth_min,
                            real_af_min,
                        )
                        score_rows.append(
                            score_sample(
                                str(payload["panel"]),
                                int(payload["sample"]),
                                method,
                                payload["truth"],
                                payload["current_pred"],
                                rescued,
                                rule,
                            )
                        )

    scores = pd.DataFrame(score_rows)
    summary = summarize(scores)
    deltas = compare_to_current(scores)
    overall = summary.groupby("method", as_index=False).agg(
        panels=("panel", lambda values: ",".join(sorted(set(map(str, values))))),
        mean_panel_pooled_F1=("pooled_F1", "mean"),
        mean_panel_official_L1_pp=("mean_official_L1_pp", "mean"),
        mean_panel_Pearson_union=("mean_Pearson_union", "mean"),
        total_TP=("pooled_TP", "sum"),
        total_FP=("pooled_FP", "sum"),
        total_FN=("pooled_FN", "sum"),
        total_rescued_TP=("rescued_TP", "sum"),
        total_rescued_FP=("rescued_FP", "sum"),
    )
    overall = overall.merge(deltas, on="method", how="left")
    noncurrent = overall.loc[~overall["method"].eq("current_default")].copy()
    top = noncurrent.sort_values(
        ["mean_delta_F1", "worsened_samples_F1", "total_FP_delta"],
        ascending=[False, True, True],
    ).head(100)
    sample_safe = noncurrent.loc[
        (noncurrent["mean_delta_F1"].fillna(0.0) > 0.0)
        & (noncurrent["worsened_samples_F1"].fillna(0).astype(float).eq(0))
    ].copy()
    panel_safe = noncurrent.loc[
        (noncurrent["mean_delta_F1"].fillna(0.0) > 0.0)
        & (noncurrent["worsened_panels_F1"].fillna(0).astype(float).eq(0))
    ].copy()
    best = top.head(1)
    best_sample_safe = sample_safe.sort_values(
        ["mean_delta_F1", "total_FP_delta"], ascending=[False, True]
    ).head(1)
    validation = pd.DataFrame(validation_rows)
    max_abs_validation = max(
        validation["delta_F1"].abs().max(),
        validation["delta_L1_pp"].abs().max(),
    )
    audit = pd.DataFrame(
        [
            {
                "metric": "evaluated_panels",
                "value": len(PANELS),
                "evidence": ",".join(PANELS),
                "decision": "diagnostic_scope",
            },
            {
                "metric": "tested_rules",
                "value": len(ANI_MINS) * len(XNY_MINS) * len(BREADTH_MINS) * len(REAL_AF_MINS),
                "evidence": "zero-mass rescue threshold grid",
                "decision": "diagnostic_scope",
            },
            {
                "metric": "baseline_validation_max_abs_delta",
                "value": max_abs_validation,
                "evidence": "current baseline F1 and official-L1 compared to official score TSVs",
                "decision": "pass" if max_abs_validation <= 1e-9 else "fail",
            },
            {
                "metric": "best_mean_F1_rule",
                "value": best["method"].iloc[0] if not best.empty else "",
                "evidence": ""
                if best.empty
                else (
                    f"mean_delta_F1={best['mean_delta_F1'].iloc[0]};"
                    f"TP_delta={best['total_TP_delta'].iloc[0]};"
                    f"FP_delta={best['total_FP_delta'].iloc[0]};"
                    f"FN_delta={best['total_FN_delta'].iloc[0]};"
                    f"worsened_samples={best['worsened_samples_F1'].iloc[0]};"
                    f"worsened_panels={best['worsened_panels_F1'].iloc[0]}"
                ),
                "decision": "best_diagnostic_rule",
            },
            {
                "metric": "sample_safe_positive_F1_rules",
                "value": int(len(sample_safe)),
                "evidence": "rules with positive mean F1 delta and no sample-level F1 decrease",
                "decision": "safety_screen",
            },
            {
                "metric": "panel_safe_positive_F1_rules",
                "value": int(len(panel_safe)),
                "evidence": "rules with positive mean F1 delta and no panel-level mean F1 decrease",
                "decision": "safety_screen",
            },
            {
                "metric": "best_sample_safe_rule",
                "value": best_sample_safe["method"].iloc[0] if not best_sample_safe.empty else "",
                "evidence": ""
                if best_sample_safe.empty
                else (
                    f"mean_delta_F1={best_sample_safe['mean_delta_F1'].iloc[0]};"
                    f"TP_delta={best_sample_safe['total_TP_delta'].iloc[0]};"
                    f"FP_delta={best_sample_safe['total_FP_delta'].iloc[0]};"
                    f"FN_delta={best_sample_safe['total_FN_delta'].iloc[0]}"
                ),
                "decision": "candidate_requires_wrapper_validation"
                if not best_sample_safe.empty
                else "no_sample_safe_rule",
            },
            {
                "metric": "promotion_decision",
                "value": "diagnostic_only_not_default",
                "evidence": "cross-panel truth-aware zero-mass rescue; no abundance recovery and no wrapper implementation",
                "decision": "do_not_promote_cross_panel_rescue_yet",
            },
        ]
    )

    scores.to_csv(TMP_SCORES, sep="\t", index=False)
    summary.to_csv(TMP_PANEL_SUMMARY, sep="\t", index=False)
    overall.to_csv(RESULTS / "cross_panel_candidate_rescue_overall.tsv", sep="\t", index=False)
    top.to_csv(RESULTS / "cross_panel_candidate_rescue_top100.tsv", sep="\t", index=False)
    validation.to_csv(RESULTS / "cross_panel_candidate_rescue_validation.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "cross_panel_candidate_rescue_audit.tsv", sep="\t", index=False)
    print(audit.to_string(index=False))
    if not best.empty:
        print("\nBEST")
        print(best.to_string(index=False))
    if not best_sample_safe.empty:
        print("\nBEST_SAMPLE_SAFE")
        print(best_sample_safe.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
