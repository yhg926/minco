#!/usr/bin/env python3
"""Score coden15 S2000 markerdb on Toy Mouse samples 0-2."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
sys.path.insert(0, str(MORE_EXP))

import score_more_toymouse_gtdb as more  # noqa: E402


EXP_DIR = ROOT / "research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse"
WORK = Path("/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623")


def coden15_output_path(sample: int) -> Path:
    return WORK / f"toymouse_sample{sample}_T15_s2000_ctxmarker_split_naive_product.tsv"


def profile_rows(sample: int) -> list[dict[str, object]]:
    path = MORE_EXP / f"mouse{sample}_gtdb_species_profile.tsv"
    return pd.read_csv(path, sep="\t").to_dict("records")


def coden15_relaxed_mask(rows: pd.DataFrame) -> pd.Series:
    return (
        (rows["XnY_ctx"] >= 20.0)
        & (rows["ANI_naive_calc"] > 0.95)
        & (rows["Reliable_ztp_af"] >= 0.25)
        & rows["active_delta_pass"]
        & rows["gtdb_species"].astype(bool)
    )


def score_method(sample: int, method: str, rows: pd.DataFrame, gold: set[str],
                 gold_abundance: dict[str, float], detail_rows: list[dict[str, object]],
                 abundance_rows: list[dict[str, object]]) -> dict[str, object]:
    selected = rows.loc[rows["active_gate_pass"]].copy()
    tp, fp, fn, precision, recall, f1 = more.score_sets(selected["gtdb_species"], gold)
    more.append_details(detail_rows, sample, method, "FP", fp, rows, gold_abundance)
    more.append_details(detail_rows, sample, method, "FN", fn, rows, gold_abundance)
    for rec in more.abundance_metrics(selected, profile_rows(sample)):
        abundance_rows.append({"sample": sample, "method": method, **rec})
    return {
        "sample": sample,
        "method": method,
        "gold_taxa": len(gold),
        "pred_taxa": len(tp | fp),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "selected_rows": len(selected),
    }


def main() -> int:
    by_accession, by_core = more.truth.load_gtdb_metadata(more.truth.GTDB_METADATA)
    score_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []

    for sample in more.SAMPLES:
        profile = profile_rows(sample)
        gold = {str(row["gtdb_species"]) for row in profile}
        gold_abundance = {
            str(row["gtdb_species"]): float(row["relative_abundance"])
            for row in profile
        }

        old_rows = more.load_minco_active_rows(
            more.minco_output_path(sample, "ctx_only_current_best"),
            by_accession,
            by_core,
        )
        score_rows.append(score_method(
            sample, "ctx_only_current_best", old_rows, gold,
            gold_abundance, detail_rows, abundance_rows,
        ))

        sylph_rows = more.load_sylph_rows(
            more.sylph_profile_path(sample), by_accession, by_core
        )
        score_rows.append(score_method(
            sample, "sylph_gtdb_profile", sylph_rows, gold,
            gold_abundance, detail_rows, abundance_rows,
        ))

        c15_rows = more.load_minco_active_rows(coden15_output_path(sample), by_accession, by_core)
        score_rows.append(score_method(
            sample, "coden15_ctxmarker_active_gate", c15_rows, gold,
            gold_abundance, detail_rows, abundance_rows,
        ))

        relaxed_rows = c15_rows.copy()
        relaxed_rows["active_gate_pass"] = coden15_relaxed_mask(relaxed_rows)
        score_rows.append(score_method(
            sample, "coden15_ctxmarker_relaxed_xny20_ztp025",
            relaxed_rows, gold, gold_abundance, detail_rows, abundance_rows,
        ))

    score_df = pd.DataFrame(score_rows)
    mean_df = (
        score_df.groupby("method", as_index=False)[["precision", "recall", "F1"]]
        .mean()
        .rename(
            columns={
                "precision": "mean_precision",
                "recall": "mean_recall",
                "F1": "mean_F1",
            }
        )
    )
    total_rows = []
    for method, sub in score_df.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        total_rows.append(
            {
                "method": method,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": precision,
                "recall": recall,
                "F1": f1,
            }
        )

    score_path = EXP_DIR / "coden15_three_sample_scores.tsv"
    mean_path = EXP_DIR / "coden15_three_sample_mean_scores.tsv"
    total_path = EXP_DIR / "coden15_three_sample_total_scores.tsv"
    detail_path = EXP_DIR / "coden15_three_sample_details.tsv"
    abundance_path = EXP_DIR / "coden15_three_sample_abundance.tsv"
    score_df.to_csv(score_path, sep="\t", index=False)
    mean_df.to_csv(mean_path, sep="\t", index=False)
    pd.DataFrame(total_rows).to_csv(total_path, sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(detail_path, sep="\t", index=False)
    pd.DataFrame(abundance_rows).to_csv(abundance_path, sep="\t", index=False)
    print(score_df.to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {mean_path}")
    print(f"wrote {total_path}")
    print(f"wrote {abundance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
