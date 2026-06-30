#!/usr/bin/env python3
"""Summarize current-wrapper HMP gastrooral table-mode scores.

The original HMP pilot note scored an older calibrated table. This diagnostic
keeps the same scorer/truth namespace but points the calibrated-table input at
current universal-auto-exact table-mode outputs, then summarizes the result for
the universal-strategy decision note.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
WORK = Path("/tmp/cami2_hmp_pilot_20260625/run")

SAMPLE_SUMMARIES = {
    "gastrooral0": WORK / "summary_sample0_current_universal_tablemode.tsv",
    "gastrooral6": WORK / "summary_sample6_current_universal_tablemode.tsv",
}

CURRENT_PROFILE_PATHS = {
    "gastrooral0": WORK / "minco_sample0_universal_autoexact_tablemode_check.tsv",
    "gastrooral6": WORK / "minco_sample6_universal_autoexact_tablemode_current.tsv",
}

METHOD_RENAME = {
    "minco_calibrated_train12": "minco_current_universal_tablemode",
}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def pooled_f1(tp: int, fp: int, fn: int) -> float:
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for method, sub in df.groupby("method", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        rows.append(
            {
                "method": method,
                "samples": ",".join(sub["sample"].astype(str).tolist()),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1(tp, fp, fn),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_FP_plus_FN": sub["FP_plus_FN"].mean(),
                "mean_abundance_l1": sub["abundance_l1"].mean(),
                "mean_tp_abundance_pearson": sub["tp_abundance_pearson"].mean(),
                "mean_tp_abundance_mae": sub["tp_abundance_mae"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values(["mean_F1", "mean_abundance_l1"], ascending=[False, True])


def main() -> int:
    missing = [path for path in [*SAMPLE_SUMMARIES.values(), *CURRENT_PROFILE_PATHS.values()] if not path.exists()]
    if missing:
        raise SystemExit("missing HMP current-universal inputs:\n" + "\n".join(map(str, missing)))

    frames = []
    mapping_rows = []
    for sample, path in SAMPLE_SUMMARIES.items():
        df = pd.read_csv(path, sep="\t")
        df["source_summary"] = str(path)
        df["method_original"] = df["method"]
        df["method"] = df["method"].replace(METHOD_RENAME)
        df["current_profile"] = ""
        df.loc[df["method"].eq("minco_current_universal_tablemode"), "current_profile"] = str(
            CURRENT_PROFILE_PATHS[sample]
        )
        frames.append(df)

        profile = pd.read_csv(CURRENT_PROFILE_PATHS[sample], sep="\t", nrows=1)
        mapping_rows.append(
            {
                "sample": sample,
                "current_profile": str(CURRENT_PROFILE_PATHS[sample]),
                "profile_strategy": profile.get("profile_strategy", pd.Series([""])).iloc[0],
                "adaptive_mode": profile.get("adaptive_mode", pd.Series([""])).iloc[0],
                "raw_unique_fallback": profile.get("raw_unique_fallback", pd.Series([""])).iloc[0],
                "raw_unique95_n": profile.get("raw_unique95_n", pd.Series([""])).iloc[0],
                "raw_unique_to_base_ratio": profile.get("raw_unique_to_base_ratio", pd.Series([""])).iloc[0],
                "joined_base_n": profile.get("joined_base_n", pd.Series([""])).iloc[0],
                "joined_base_median_uaf": profile.get("joined_base_median_uaf", pd.Series([""])).iloc[0],
                "probability_extra_mass_ratio": profile.get("probability_extra_mass_ratio", pd.Series([""])).iloc[0],
                "auto_exact_split_requested": profile.get("auto_exact_split_requested", pd.Series([""])).iloc[0],
                "auto_exact_split_used": profile.get("auto_exact_split_used", pd.Series([""])).iloc[0],
                "auto_exact_split_unavailable_reason": profile.get(
                    "auto_exact_split_unavailable_reason",
                    pd.Series([""]),
                ).iloc[0],
            }
        )

    detail = pd.concat(frames, ignore_index=True, sort=False)
    summary = summarize(detail)
    RESULTS.mkdir(parents=True, exist_ok=True)
    detail.to_csv(RESULTS / "hmp_gastrooral_current_universal_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "hmp_gastrooral_current_universal_summary.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "hmp_gastrooral_current_universal_mapping.tsv",
        mapping_rows,
        [
            "sample",
            "current_profile",
            "profile_strategy",
            "adaptive_mode",
            "raw_unique_fallback",
            "raw_unique95_n",
            "raw_unique_to_base_ratio",
            "joined_base_n",
            "joined_base_median_uaf",
            "probability_extra_mass_ratio",
            "auto_exact_split_requested",
            "auto_exact_split_used",
            "auto_exact_split_unavailable_reason",
        ],
    )
    print(summary.to_string(index=False))
    print("\nMAPPING")
    print(pd.DataFrame(mapping_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
