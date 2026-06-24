#!/usr/bin/env python3
"""Score coden15 S2000 markerdb readwise output on Toy Mouse sample0."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SCORE_EXP = ROOT / "research/experiments/2026-06-22_s2000_ctxobj_markerdb_toymouse"
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(SCORE_EXP))

import score_ctxobj_markerdb as active  # noqa: E402


EXP_DIR = ROOT / "research/experiments/2026-06-23_gtdb_s2000_coden15_toymouse"
WORK = Path("/mnt/new3T/gtdbr220/eval_runs/minco_coden15_s2000_toymouse_20260623")

OLD_CTX_TSV = OLD_EXP / "toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv"
CODEN15_TSV = WORK / "toymouse_sample0_T15_s2000_ctxmarker_split_naive_product.tsv"
PROFILE_TSV = OLD_EXP / "mouse0_gtdb_species_profile.tsv"


def main() -> int:
    by_accession, by_core = active.truth.load_gtdb_metadata(active.truth.GTDB_METADATA)
    profile = pd.read_csv(PROFILE_TSV, sep="\t")
    gold = set(profile["gtdb_species"].astype(str))
    gold_abundance = dict(
        zip(profile["gtdb_species"].astype(str), profile["relative_abundance"].astype(float))
    )

    score_rows = []
    detail_rows = []
    abundance_rows = []
    row_paths = [
        ("ctx_only_current_best_recheck", OLD_CTX_TSV),
        ("coden15_ctxmarker_active_gate", CODEN15_TSV),
    ]
    for method, path in row_paths:
        rows = active.load_rows(path, by_accession, by_core)
        selected = rows.loc[rows["active_gate_pass"]].copy()
        tp, fp, fn, precision, recall, f1 = active.score_sets(
            selected["gtdb_species"], gold
        )
        score_rows.append(
            {
                "method": method,
                "input_tsv": str(path),
                "pred_taxa": len(tp | fp),
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "precision": precision,
                "recall": recall,
                "F1": f1,
                "selected_rows": len(selected),
            }
        )
        active.append_details(detail_rows, method, "FP", fp, rows, gold_abundance)
        active.append_details(detail_rows, method, "FN", fn, rows, gold_abundance)
        for rec in active.abundance_metrics(selected, profile):
            abundance_rows.append({"method": method, **rec})

    score_path = EXP_DIR / "coden15_sample0_gtdb_scores.tsv"
    detail_path = EXP_DIR / "coden15_sample0_gtdb_details.tsv"
    abundance_path = EXP_DIR / "coden15_sample0_gtdb_abundance.tsv"
    pd.DataFrame(score_rows).to_csv(score_path, sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(detail_path, sep="\t", index=False)
    pd.DataFrame(abundance_rows).to_csv(abundance_path, sep="\t", index=False)
    print(pd.DataFrame(score_rows).to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {detail_path}")
    print(f"wrote {abundance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
