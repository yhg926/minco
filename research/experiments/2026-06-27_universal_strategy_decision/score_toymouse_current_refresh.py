#!/usr/bin/env python3
"""Score refreshed current-code CAMI II Toy Mouse sample5-7 profiles."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
TOYMOUSE_SCORER = (
    ROOT
    / "research/experiments/2026-06-26_cami2_toymouse_current_default/"
    "score_sample7_current_default.py"
)
REFRESH = Path("/tmp/minco_current_code_toymouse_refresh_20260627")
SYLPH_RUN = Path("/tmp/cami2_toymouse_samples5_7_20260625/run")


spec = importlib.util.spec_from_file_location("toy_current_scorer", TOYMOUSE_SCORER)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import scorer {TOYMOUSE_SCORER}")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
score = mod.score


def minco_path(sample: int) -> Path:
    return REFRESH / f"sample{sample}_current.tsv"


def sylph_path(sample: int) -> Path:
    return SYLPH_RUN / f"sylph_sample{sample}/profile.tsv"


def feature_summary(path: Path) -> dict[str, object]:
    raw = pd.read_csv(path, sep="\t")
    calls = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    return {
        "called": int(calls.sum()),
        "rows": int(len(raw)),
        "adaptive_mode": ";".join(map(str, raw.get("adaptive_mode", pd.Series()).dropna().astype(str).unique())),
        "low_extra_split_rescue_added_n": raw.get("low_extra_split_rescue_added_n", pd.Series([0])).iloc[0]
        if len(raw)
        else 0,
        "reported_ani_present": "reported_ani" in raw.columns,
        "abundance_rule": raw.get("tail_rescue_abundance_rule", pd.Series([""])).iloc[0]
        if len(raw)
        else "",
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core = score.truth.load_gtdb_metadata(score.truth.GTDB_METADATA)
    genome_taxid, genome_location = mod.load_meta()
    score.EXP_DIR = RESULTS
    score.BASE = mod.DATA_BASE
    score.distribution_path = mod.distribution_path

    score_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []

    for sample in [5, 6, 7]:
        source_rows = score.build_gtdb_truth(
            sample, by_accession, by_core, genome_taxid, genome_location
        )
        profile = score.write_truth_files(sample, source_rows)
        gold = {str(row["gtdb_species"]) for row in profile}
        gold_abundance = {
            str(row["gtdb_species"]): float(row["relative_abundance"]) for row in profile
        }

        methods = [
            ("minco_current_code_refresh", minco_path(sample), "minco"),
            ("sylph_gtdb_profile", sylph_path(sample), "sylph"),
        ]
        for method, path, kind in methods:
            if not path.exists():
                raise FileNotFoundError(path)
            if kind == "sylph":
                rows = score.load_sylph_rows(path, by_accession, by_core)
            else:
                rows = mod.score_calibrated_minco(path)
                feature_rows.append({"sample": sample, "method": method, "path": str(path), **feature_summary(path)})
            selected = rows.loc[rows["active_gate_pass"]].copy()
            tp, fp, fn, precision, recall, f1 = score.score_sets(selected["gtdb_species"], gold)
            score_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "path": str(path),
                    "gold_taxa": len(gold),
                    "pred_taxa": len(tp | fp),
                    "TP": len(tp),
                    "FP": len(fp),
                    "FN": len(fn),
                    "precision": precision,
                    "recall": recall,
                    "F1": f1,
                }
            )
            score.append_details(detail_rows, sample, method, "FP", fp, rows, gold_abundance)
            score.append_details(detail_rows, sample, method, "FN", fn, rows, gold_abundance)
            for rec in score.abundance_metrics(selected, profile):
                abundance_rows.append({"sample": sample, "method": method, **rec})

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    compare = score_df.merge(
        abundance_df.loc[abundance_df["renorm_pred"].astype(str).str.lower().eq("true")],
        on=["sample", "method"],
        how="left",
    )
    compare.to_csv(RESULTS / "toymouse_current_refresh_compare.tsv", sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(RESULTS / "toymouse_current_refresh_details.tsv", sep="\t", index=False)
    pd.DataFrame(feature_rows).to_csv(RESULTS / "toymouse_current_refresh_features.tsv", sep="\t", index=False)

    summary_rows: list[dict[str, object]] = []
    for method, sub in compare.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        summary_rows.append(
            {
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique()))),
                "mean_F1": sub["F1"].mean(),
                "mean_L1_pp": sub["l1_pct_points"].mean(),
                "mean_Pearson": sub["pearson"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": pooled_f1,
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("mean_F1", ascending=False)
    summary.to_csv(RESULTS / "toymouse_current_refresh_summary.tsv", sep="\t", index=False)
    print(compare[["sample", "method", "TP", "FP", "FN", "precision", "recall", "F1", "l1_pct_points", "pearson"]].to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
