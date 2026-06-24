#!/usr/bin/env python3
"""Diagnose simple call filters on the CAMI strain external holdout."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
HELPER_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(HELPER_DIR))

from analyze_readwise_corrections import parse_species_taxmap  # noqa: E402


def invalid_species_name(name: str) -> bool:
    text = str(name).lower().strip()
    generic_markers = [
        "uncultured",
        " bacterium",
        "bacterium ",
        " metagenome",
        "environmental sample",
    ]
    return text == "bacterium" or any(marker in text for marker in generic_markers)


def score_mask(
    df: pd.DataFrame,
    mask: pd.Series,
    gold_counts: Mapping[str, int],
    *,
    split: str,
    method: str,
) -> Dict[str, object]:
    values: List[Dict[str, float]] = []
    for sample_key, sample_df in df.groupby("sample_key"):
        sample_mask = mask.loc[sample_df.index]
        tp = int((sample_mask & sample_df["label"].eq(1)).sum())
        fp = int((sample_mask & sample_df["label"].eq(0)).sum())
        fn = int(gold_counts[str(sample_key)] - tp)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        values.append(
            {
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": precision,
                "recall": recall,
                "F1": f1,
                "FP_plus_FN": fp + fn,
            }
        )
    return {
        "split": split,
        "method": method,
        "samples": len(values),
        "mean_precision": float(np.mean([v["precision"] for v in values])),
        "mean_recall": float(np.mean([v["recall"] for v in values])),
        "mean_F1": float(np.mean([v["F1"] for v in values])),
        "mean_FP": float(np.mean([v["FP"] for v in values])),
        "mean_FN": float(np.mean([v["FN"] for v in values])),
        "mean_FP_plus_FN": float(np.mean([v["FP_plus_FN"] for v in values])),
    }


def add_taxon_quality(df: pd.DataFrame, taxmap: Mapping[str, Mapping[str, str]]) -> pd.DataFrame:
    name_by_taxid = {
        str(rec.get("taxid", "")): str(rec.get("species", "") or rec.get("taxpathsn", "")).split("|")[-1]
        for rec in taxmap.values()
        if str(rec.get("taxid", ""))
    }
    out = df.copy()
    out["taxon_name"] = out["taxid"].astype(str).map(name_by_taxid).fillna("")
    out["valid_taxon_name"] = ~out["taxon_name"].map(invalid_species_name)
    return out


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-dir", type=Path, required=True)
    ap.add_argument("--taxmap", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    taxmap = parse_species_taxmap(args.taxmap)
    summary = pd.read_csv(args.score_dir / "summary.tsv", sep="\t")
    gold_counts = {
        split: summary.loc[summary["split"].eq(split)]
        .drop_duplicates("sample_key")
        .set_index("sample_key")["gold_taxa"]
        .astype(int)
        .to_dict()
        for split in ["train", "test"]
    }

    records: List[Dict[str, object]] = []
    for split in ["train", "test"]:
        features = add_taxon_quality(pd.read_csv(args.score_dir / f"{split}.joined_features.tsv", sep="\t").fillna(0.0), taxmap)
        masks = {
            "u_direct": features["u_direct_call"] >= 0.5,
            "s_direct": features["s_direct_call"] >= 0.5,
            "u_or_s_direct": (features["u_direct_call"] >= 0.5) | (features["s_direct_call"] >= 0.5),
            "u_relaxed": features["u_relaxed_call"] >= 0.5,
            "s_relaxed": features["s_relaxed_call"] >= 0.5,
        }
        for method, mask in masks.items():
            records.append(score_mask(features, mask, gold_counts[split], split=split, method=method))
            records.append(score_mask(features, mask & features["valid_taxon_name"], gold_counts[split], split=split, method=f"{method}_validtax"))

        if split == "test":
            best_simple = masks["u_or_s_direct"] & features["valid_taxon_name"]
            diagnostic_cols = [
                "sample_key",
                "taxid",
                "taxon_name",
                "u_direct_call",
                "s_direct_call",
                "u_relaxed_call",
                "s_relaxed_call",
                "u_XnY_ctx_max",
                "s_XnY_ctx_max",
                "u_ANI_max",
                "s_ANI_max",
                "u_Real_min_align_fraction_max",
                "s_Real_min_align_fraction_max",
                "u_Ref_breadth_max",
                "s_Ref_breadth_max",
                "u_Ref_mean_depth_max",
                "s_Ref_mean_depth_max",
                "u_Ref_depth_cv_min",
                "s_Ref_depth_cv_min",
                "u_Normalized_abundance_depth_max",
                "s_Normalized_abundance_depth_max",
            ]
            features.loc[features["label"].eq(1) & ~best_simple, diagnostic_cols].sort_values(
                ["sample_key", "s_XnY_ctx_max"], ascending=[True, False]
            ).to_csv(args.outdir / "u_or_s_direct_validtax_missed_gold.tsv", sep="\t", index=False)
            features.loc[features["label"].eq(0) & best_simple, diagnostic_cols].sort_values(
                ["sample_key", "s_XnY_ctx_max"], ascending=[True, False]
            ).to_csv(args.outdir / "u_or_s_direct_validtax_false_positives.tsv", sep="\t", index=False)

            fp = features.loc[
                ((features["u_direct_call"] >= 0.5) | (features["s_direct_call"] >= 0.5))
                & features["label"].eq(0)
                & ~features["valid_taxon_name"],
                [
                    "sample_key",
                    "taxid",
                    "taxon_name",
                    "u_direct_call",
                    "s_direct_call",
                    "u_XnY_ctx_max",
                    "s_XnY_ctx_max",
                    "u_ANI_max",
                    "s_ANI_max",
                    "u_Ref_breadth_max",
                    "s_Ref_breadth_max",
                ],
            ].sort_values(["sample_key", "s_XnY_ctx_max"], ascending=[True, False])
            fp.to_csv(args.outdir / "invalid_taxon_direct_false_positives.tsv", sep="\t", index=False)

    predictions = pd.read_csv(args.score_dir / "model_predictions.tsv", sep="\t")
    test_features = add_taxon_quality(
        pd.read_csv(args.score_dir / "test.joined_features.tsv", sep="\t", usecols=["sample_key", "taxid", "label"]).fillna(0.0),
        taxmap,
    )
    predictions = predictions.merge(test_features, on=["sample_key", "taxid", "label"], how="left")
    for method, group in predictions.groupby("method"):
        raw_mask = group["predicted"].astype(bool)
        records.append(score_mask(group, raw_mask, gold_counts["test"], split="test", method=method))
        records.append(
            score_mask(
                group,
                raw_mask & group["valid_taxon_name"].astype(bool),
                gold_counts["test"],
                split="test",
                method=f"{method}_validtax",
            )
        )

    out = pd.DataFrame(records).sort_values(["split", "mean_F1"], ascending=[True, False])
    out.to_csv(args.outdir / "simple_filter_diagnostics.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
