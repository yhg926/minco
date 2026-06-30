#!/usr/bin/env python3
"""Cross-validate the CAMI3 source-readmap loose split rescue.

The source-readmap diagnostic found a loose high-depth split rule that improves
CAMI3 samples0-2. This script applies that exact rule to the older 26-sample
readiness panel without using truth to select rows. It is intended to decide
whether the rule is broad enough to promote or should remain CAMI3-specific.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"
PRIOR = REPO / "research/experiments/2026-06-26_cami2_hmp_unseen_transfer"
TMP_RUN = Path("/tmp/cami3_toy_human_gut_20260626/run")

sys.path.insert(0, str(PRIOR))
sys.path.insert(0, str(REPO))

import score_cami3_autoexact as cami3  # noqa: E402
import score_raw_unique_selector as base  # noqa: E402
import search_high_pextra_tail_rescue as tail  # noqa: E402


RULE = {
    "prob_min": 0.02,
    "s_xny_min": 300.0,
    "s_ani_min": 0.93,
    "s_af_min": 0.30,
    "s_breadth_min": 0.20,
    "s_mean_depth_min": 1.0,
    "topn_per_genus": 1,
}


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(np.full(len(df), default), index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(default)


def topn_by_genus(df: pd.DataFrame, mask: np.ndarray, score: np.ndarray, topn: int) -> np.ndarray:
    idx = np.flatnonzero(mask)
    keep = np.zeros(len(df), dtype=bool)
    if len(idx) == 0 or topn <= 0:
        return keep
    if "genus" in df.columns:
        genus = df["genus"].astype(str).to_numpy()
    else:
        names = df.get("species_name", pd.Series([""] * len(df))).astype(str)
        genus = names.map(base.genus_from_name).to_numpy()
    work = pd.DataFrame({"idx": idx, "genus": genus[idx], "score": score[idx]})
    work = work.sort_values(["genus", "score"], ascending=[True, False])
    keep[work.groupby("genus", as_index=False).head(topn)["idx"].to_numpy(dtype=int)] = True
    return keep


def loose_rescue(df: pd.DataFrame, baseline: np.ndarray) -> np.ndarray:
    prob = numeric(df, "calibrated_probability").to_numpy(dtype=float)
    s_x = numeric(df, "s_XnY_ctx_max").to_numpy(dtype=float)
    s_ani = numeric(df, "s_ANI_max").to_numpy(dtype=float)
    s_af = numeric(df, "s_Real_min_align_fraction_max").to_numpy(dtype=float)
    s_b = numeric(df, "s_Ref_breadth_max").to_numpy(dtype=float)
    s_mean = numeric(df, "s_Ref_mean_depth_max").to_numpy(dtype=float)
    raw = (
        (~baseline)
        & (prob >= RULE["prob_min"])
        & (s_x >= RULE["s_xny_min"])
        & (s_ani >= RULE["s_ani_min"])
        & (s_af >= RULE["s_af_min"])
        & (s_b >= RULE["s_breadth_min"])
        & (s_mean >= RULE["s_mean_depth_min"])
    )
    score = s_mean * np.maximum(s_b, 1e-6) * np.maximum(s_ani, 0.0) * np.maximum(s_x, 1.0)
    return topn_by_genus(df, raw, score, int(RULE["topn_per_genus"]))


def abundance_zip_power(df: pd.DataFrame, pred_taxids: set[str], power: float) -> dict[str, float]:
    if not pred_taxids:
        return {}
    work = df.loc[df["taxid"].astype(str).isin(pred_taxids)].copy()
    s_mean = numeric(work, "s_Ref_mean_depth_max").to_numpy(dtype=float)
    s_zip = numeric(work, "s_Ref_zip_af_max").to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        raw = s_mean / np.maximum(s_zip, 1e-6) ** power
    work["__raw"] = np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)
    return dict(work.groupby(work["taxid"].astype(str))["__raw"].max())


def score_one(
    method: str,
    sample: base.Sample,
    gold: pd.DataFrame,
    df: pd.DataFrame,
    mask: np.ndarray,
    abundance_mode: str,
    added: np.ndarray,
) -> dict[str, object]:
    pred = set(df.loc[mask, "taxid"].astype(str))
    if abundance_mode == "tail_safe":
        work = tail.add_abundance_columns(df)
        mean_raw = numeric(work, "s_Ref_mean_depth_max").to_numpy(dtype=float)
        raw = numeric(work, "__abundance_tail_safe").to_numpy(dtype=float)
        raw = np.where(added, mean_raw, raw)
        work["__final_raw"] = raw
        abundance = dict(
            work.loc[work["taxid"].astype(str).isin(pred)].groupby(work["taxid"].astype(str))["__final_raw"].max()
        )
    elif abundance_mode == "zip_p025":
        abundance = abundance_zip_power(df, pred, 0.25)
    else:
        raise ValueError(f"unsupported abundance mode: {abundance_mode}")
    rec = base.score_pred(method, sample, gold, pred, abundance)
    rec["added"] = int(np.sum(added))
    rec["abundance_mode"] = abundance_mode
    return rec


def evaluate_table(
    sample: base.Sample,
    gold: pd.DataFrame,
    df: pd.DataFrame,
    source: str,
) -> list[dict[str, object]]:
    if "species_name" not in df.columns:
        df = df.copy()
        df["species_name"] = ""
    masks, details = tail.adaptive_masks(df)
    baseline = masks["tail_low_uaf_probability"]
    added = loose_rescue(df, baseline)
    final = baseline | added
    rows = [
        score_one("baseline_tail_low_uaf_probability", sample, gold, df, baseline, "tail_safe", np.zeros(len(df), dtype=bool)),
        score_one("loose_cami3_split_rescue_tail_safe", sample, gold, df, final, "tail_safe", added),
        score_one("loose_cami3_split_rescue_zip_p025", sample, gold, df, final, "zip_p025", added),
    ]
    for row in rows:
        row.update({"source": source, **details, **RULE})
    return rows


def load_cami3_profile(sid: int) -> Path:
    guarded = TMP_RUN / f"sample{sid}_universal_autoexact_guarded.tsv"
    if guarded.exists():
        return guarded
    return TMP_RUN / f"sample{sid}_universal_autoexact.tsv"


def summarize(sample_metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method, group in sample_metrics.groupby("method"):
        rec = {
            "method": method,
            "samples": int(group["sample_key"].nunique()),
            "mean_F1": float(group["F1"].mean()),
            "mean_L1": float(group["L1"].mean()),
            "mean_Pearson": float(group["Pearson"].mean(skipna=True)),
            "mean_FP_plus_FN": float(group["FP_plus_FN"].mean()),
            "min_F1": float(group["F1"].min()),
            "max_L1": float(group["L1"].max()),
            "mean_added": float(group["added"].mean()),
        }
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["mean_F1", "mean_L1"], ascending=[False, True])


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    panel, prob, _gold_counts = base.append_airskin_calibrated(pd.read_csv(base.COMBINED, sep="\t"))
    panel["calibrated_probability"] = prob
    if "species_name" not in panel.columns:
        panel["species_name"] = ""
    for sample in base.base_samples():
        sub = panel.loc[panel["sample_key"].astype(str).eq(sample.sample_key)].copy()
        if sub.empty:
            continue
        gold = base.parse_gold_profile(sample.gold_profile, sample.gold_sample_id, sample.scope)
        rows.extend(evaluate_table(sample, gold, sub, "cached20"))

    for sid in range(20):
        path = load_cami3_profile(sid)
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        gold = base.parse_gold_profile(cami3.TRUTH_DIR / f"taxonomic_profile_{sid}.txt", str(sid), "bacteria")
        sample = base.Sample(
            sample_key=f"cami3_toy_human_gut{sid}",
            dataset="cami3_toy_human_gut",
            scope="bacteria",
            gold_profile=Path(""),
            gold_sample_id=str(sid),
            unique_unfiltered=Path(""),
            split_unfiltered=Path(""),
            sylph_profile=None,
        )
        rows.extend(evaluate_table(sample, gold, df, "cami3_current"))

    sample_metrics = pd.DataFrame(rows)
    summary = summarize(sample_metrics)
    by_source = (
        sample_metrics.groupby(["source", "method"], as_index=False)
        .agg(
            samples=("sample_key", "nunique"),
            mean_F1=("F1", "mean"),
            mean_L1=("L1", "mean"),
            mean_Pearson=("Pearson", "mean"),
            mean_FP_plus_FN=("FP_plus_FN", "mean"),
            mean_added=("added", "mean"),
        )
        .sort_values(["source", "mean_F1", "mean_L1"], ascending=[True, False, True])
    )
    sample_metrics.to_csv(RESULTS / "cami3_loose_split_rescue_crossval_sample_metrics.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "cami3_loose_split_rescue_crossval_summary.tsv", sep="\t", index=False)
    by_source.to_csv(RESULTS / "cami3_loose_split_rescue_crossval_by_source.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print()
    print(by_source.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
