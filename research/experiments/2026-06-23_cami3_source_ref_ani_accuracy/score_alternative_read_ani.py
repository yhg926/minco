#!/usr/bin/env python3
"""Score alternate read-ANI columns against cached source/ref skani ANI."""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def accession(text: object) -> str:
    match = ACC_RE.search(str(text))
    return match.group(1) if match else ""


def minco_naive_ani(xny_ctx: float, n_diff_obj: float, n_diff_obj_section: float) -> float:
    if xny_ctx <= 0.0:
        return math.nan
    eps = 1e-8
    ratio = (n_diff_obj_section + eps) / (n_diff_obj + eps)
    dist0 = n_diff_obj / (xny_ctx + n_diff_obj) if xny_ctx + n_diff_obj > 0.0 else 0.0
    final_dist = 1.0 - math.pow(1.0 - dist0, ratio)
    return 1.0 - final_dist * 0.1544286 if final_dist > 0.0 else 1.0


def load_minco_raw(paths: dict[int, str]) -> pd.DataFrame:
    frames = []
    for sample, path in paths.items():
        rows = pd.read_csv(path, sep="\t")
        acc = rows["Ref"].map(accession)
        if "Ref_annotation" in rows:
            missing = acc == ""
            acc.loc[missing] = rows.loc[missing, "Ref_annotation"].map(accession)
        rows["ref_accession"] = acc
        rows["ANI_naive_calc_raw"] = [
            minco_naive_ani(float(x), float(d), float(s))
            for x, d, s in zip(
                pd.to_numeric(rows["XnY_ctx"], errors="coerce").fillna(0.0),
                pd.to_numeric(rows["N_diff_obj"], errors="coerce").fillna(0.0),
                pd.to_numeric(rows["N_diff_obj_section"], errors="coerce").fillna(0.0),
            )
        ]
        keep = rows[
            ["ref_accession", "ANI", "Ref_zip_aaf_ani", "ANI_naive_calc_raw", "XnY_ctx"]
        ].copy()
        keep["sample_id"] = sample
        frames.append(keep)
    out = pd.concat(frames, ignore_index=True)
    out = out.sort_values(["sample_id", "ref_accession", "XnY_ctx"], ascending=[True, True, False])
    return out.drop_duplicates(["sample_id", "ref_accession"])


def load_sylph_raw(paths: dict[int, str]) -> pd.DataFrame:
    frames = []
    for sample, path in paths.items():
        rows = pd.read_csv(path, sep="\t")
        rows["ref_accession"] = rows["Genome_file"].map(accession)
        rows["sample_id"] = sample
        frames.append(rows[["sample_id", "ref_accession", "Adjusted_ANI", "Naive_ANI"]])
    return pd.concat(frames, ignore_index=True)


def add_summary(rows: list[dict[str, object]], method: str, detail: pd.DataFrame, metric: str) -> None:
    tmp = detail.copy()
    tmp["read_ani_alt"] = pd.to_numeric(tmp[metric], errors="coerce")
    for choice, group in tmp.groupby("source_choice"):
        group = group.dropna(subset=["read_ani_alt", "source_ref_ani"])
        err = group["read_ani_alt"] - group["source_ref_ani"]
        rows.append(
            {
                "method": method,
                "source_choice": choice,
                "n_calls": len(group),
                "pearson": group["read_ani_alt"].corr(group["source_ref_ani"]),
                "spearman": group["read_ani_alt"].corr(group["source_ref_ani"], method="spearman"),
                "mae": err.abs().mean(),
                "mean_error": err.mean(),
                "median_abs_error": err.abs().median(),
                "read_ani_mean": group["read_ani_alt"].mean(),
                "source_ref_ani_mean": group["source_ref_ani"].mean(),
            }
        )


def main() -> int:
    detail = pd.read_csv(EXP / "read_vs_source_ref_ani_detail.tsv", sep="\t")

    old = load_minco_raw(
        {
            0: "/tmp/cami3_toygut_abundance_rescue_20260623/minco_s0_ctxmarker_current_binary.tsv",
            1: "/tmp/cami3_toygut_abundance_rescue_20260623/minco_s1_ctxmarker_current_binary.tsv",
            2: "/tmp/cami3_toygut_abundance_rescue_20260623/minco_s2_ctxmarker_current_binary.tsv",
        }
    )
    coden15 = load_minco_raw(
        {
            0: str(
                EXP.parent
                / "2026-06-23_cami3_toygut_coden15_formula_af/cami3_sample0_coden15_ctxmarker_split_naive_product.tsv"
            ),
            1: str(
                EXP.parent
                / "2026-06-23_cami3_toygut_coden15_formula_af/cami3_sample1_coden15_ctxmarker_split_naive_product.tsv"
            ),
            2: str(
                EXP.parent
                / "2026-06-23_cami3_toygut_coden15_formula_af/cami3_sample2_coden15_ctxmarker_split_naive_product.tsv"
            ),
        }
    )
    sylph = load_sylph_raw(
        {
            0: "/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv",
            1: "/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv",
            2: "/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample2/profile.tsv",
        }
    )

    rows: list[dict[str, object]] = []
    for method, raw in [
        ("minco_old_ctxmarker", old),
        ("minco_coden15_formula_af", coden15),
    ]:
        sub = detail.loc[detail["method"] == method].merge(
            raw, on=["sample_id", "ref_accession"], how="left"
        )
        for metric in ["Ref_zip_aaf_ani", "ANI_naive_calc_raw", "ANI"]:
            add_summary(rows, f"{method}_{metric}", sub, metric)

    sub = detail.loc[detail["method"] == "sylph_adjusted"].merge(
        sylph, on=["sample_id", "ref_accession"], how="left"
    )
    for metric in ["Adjusted_ANI", "Naive_ANI"]:
        sub_metric = sub.copy()
        sub_metric[metric] = pd.to_numeric(sub_metric[metric], errors="coerce") / 100.0
        add_summary(rows, f"sylph_{metric}", sub_metric, metric)

    out = pd.DataFrame(rows).sort_values(["source_choice", "mae", "method"])
    out.to_csv(EXP / "alternative_read_ani_summary.tsv", sep="\t", index=False)
    print(out.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
