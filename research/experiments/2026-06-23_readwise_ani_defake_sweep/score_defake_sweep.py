#!/usr/bin/env python3
"""Score readwise ANI defake/top-fraction sweeps against GTDB source-rep ANIm."""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
TRUTH = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_source_rep_ani.tsv"
ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
RUN_RE = re.compile(
    r"toymouse_sample0_(?P<coden>c(?:11|15))_adjusted_ani_(?:(?P<filter>product1)_)?p(?P<pct>[0-9]{3})\.tsv$"
)


def accession(text: object) -> str:
    match = ACC_RE.search(str(text))
    return match.group(1) if match else ""


def as_float(series: pd.Series, default: float = math.nan) -> pd.Series:
    if series is None:
        return pd.Series(dtype="float64")
    return pd.to_numeric(series, errors="coerce").fillna(default)


def add_metrics(rows: list[dict[str, object]], meta: dict[str, object], detail: pd.DataFrame, estimator: str) -> None:
    sub = detail.dropna(subset=["ANIm_ANI", estimator]).copy()
    if sub.empty:
        rows.append({**meta, "estimator": estimator, "n": 0})
        return
    err = sub[estimator] - sub["ANIm_ANI"]
    rows.append(
        {
            **meta,
            "estimator": estimator,
            "n": len(sub),
            "pearson": sub[estimator].corr(sub["ANIm_ANI"], method="pearson"),
            "spearman": sub[estimator].corr(sub["ANIm_ANI"], method="spearman"),
            "mae": err.abs().mean(),
            "rmse": math.sqrt((err * err).mean()),
            "mean_error": err.mean(),
            "median_abs_error": err.abs().median(),
            "estimator_mean": sub[estimator].mean(),
            "ANIm_mean": sub["ANIm_ANI"].mean(),
            "high_error_gt_0p02": int((err.abs() > 0.02).sum()),
            "high_error_gt_0p03": int((err.abs() > 0.03).sum()),
        }
    )


def load_run(path: Path) -> tuple[dict[str, object], pd.DataFrame]:
    match = RUN_RE.match(path.name)
    if not match:
        raise ValueError(f"unexpected result filename: {path.name}")
    coden = match.group("coden")
    filter_name = "product1-topfrac-median" if match.group("filter") else "product-topfrac-median"
    pct = int(match.group("pct"))
    top_fraction = (100 - pct) / 100.0

    raw = pd.read_csv(path, sep="\t")
    raw["gtdb_representative_accession"] = raw["Ref"].map(accession)
    if "Ref_annotation" in raw:
        missing = raw["gtdb_representative_accession"] == ""
        raw.loc[missing, "gtdb_representative_accession"] = raw.loc[missing, "Ref_annotation"].map(accession)
    raw["XnY_ctx_num"] = pd.to_numeric(raw["XnY_ctx"], errors="coerce")
    raw = raw.sort_values(["gtdb_representative_accession", "XnY_ctx_num"], ascending=[True, False])
    raw = raw.drop_duplicates(["gtdb_representative_accession"])

    keep_cols = [
        "gtdb_representative_accession",
        "ANI",
        "Ref_zip_aaf_ani",
        "XnY_ctx",
        "N_diff_obj",
        "N_diff_obj_section",
        "Ref_breadth",
        "Ref_hit_mean_depth",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_ctx",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_median_depth",
        "Rejected_ctx",
        "Rejected_diff_ctx",
        "Fake_ctx_fraction",
    ]
    for col in keep_cols:
        if col not in raw.columns:
            raw[col] = math.nan
    run = raw[keep_cols].copy()
    numeric_cols = [c for c in keep_cols if c != "gtdb_representative_accession"]
    for col in numeric_cols:
        run[col] = pd.to_numeric(run[col], errors="coerce")

    meta = {
        "run": path.stem,
        "coden": coden,
        "filter": filter_name,
        "top_fraction": top_fraction,
        "percentile_kept": pct / 100.0,
        "path": str(path),
    }
    return meta, run


def main() -> int:
    truth = pd.read_csv(TRUTH, sep="\t")
    truth = truth.loc[truth["ANIm_ANI"].notna()].copy()
    truth["ANIm_ANI"] = pd.to_numeric(truth["ANIm_ANI"], errors="coerce")

    summary_rows: list[dict[str, object]] = []
    joined_frames: list[pd.DataFrame] = []
    crisp_rows: list[pd.DataFrame] = []

    for path in sorted(RESULTS.glob("toymouse_sample0_c*_adjusted_ani_*.tsv")):
        meta, run = load_run(path)
        detail = truth.merge(run, on="gtdb_representative_accession", how="left")
        detail["ANI_adjusted_naive"] = detail["ANI"]
        detail["ANI_zip_aaf"] = detail["Ref_zip_aaf_ani"]
        median_depth = detail["Reliable_Ref_hit_median_depth"]
        detail["ANI_median2_else_zip"] = detail["ANI_adjusted_naive"].where(
            median_depth >= 2.0, detail["ANI_zip_aaf"]
        )
        detail["ANI_median2_else_max"] = detail["ANI_adjusted_naive"].where(
            median_depth >= 2.0,
            detail[["ANI_adjusted_naive", "ANI_zip_aaf"]].max(axis=1),
        )
        detail["abs_err_adjusted_naive"] = (detail["ANI_adjusted_naive"] - detail["ANIm_ANI"]).abs()
        detail["abs_err_median2_else_zip"] = (detail["ANI_median2_else_zip"] - detail["ANIm_ANI"]).abs()
        detail["abs_err_median2_else_max"] = (detail["ANI_median2_else_max"] - detail["ANIm_ANI"]).abs()

        for key, val in meta.items():
            detail[key] = val
        joined_frames.append(detail)

        for estimator in [
            "ANI_adjusted_naive",
            "ANI_zip_aaf",
            "ANI_median2_else_zip",
            "ANI_median2_else_max",
        ]:
            add_metrics(summary_rows, meta, detail, estimator)

        crisp = detail.loc[
            detail["gtdb_species"].astype(str).str.contains("Lactobacillus crispatus", regex=False)
        ].copy()
        if not crisp.empty:
            crisp_rows.append(crisp)

    if not summary_rows:
        raise SystemExit("no sweep result files found")

    summary = pd.DataFrame(summary_rows).sort_values(["mae", "rmse", "coden", "filter", "top_fraction"])
    summary.to_csv(EXP / "summary.tsv", sep="\t", index=False)

    joined = pd.concat(joined_frames, ignore_index=True)
    joined.to_csv(EXP / "source_positive_joined.tsv", sep="\t", index=False)

    worst = joined.sort_values(["run", "abs_err_median2_else_zip"], ascending=[True, False]).groupby("run").head(10)
    worst.to_csv(EXP / "worst10_by_run.tsv", sep="\t", index=False)

    if crisp_rows:
        crisp_out = pd.concat(crisp_rows, ignore_index=True)
        crisp_out.to_csv(EXP / "crispatus_by_run.tsv", sep="\t", index=False)

    print(summary.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
