#!/usr/bin/env python3
"""Grid-search simple species reporting thresholds for minco CAMI tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
PREV_SCRIPT_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(PREV_SCRIPT_DIR))

from analyze_readwise_corrections import (  # noqa: E402
    load_minco,
    load_sylph,
    parse_gold_profile,
    parse_species_taxmap,
    score_taxids,
)


def score_mask(rows: pd.DataFrame, mask: pd.Series, gold_taxids: Set[str]) -> Dict[str, object]:
    selected = rows.loc[mask & rows["taxid"].astype(bool)]
    stats = score_taxids(selected["taxid"].astype(str), gold_taxids)
    stats["selected_rows"] = int(len(selected))
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--minco", type=Path, required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--gold", type=Path, default=Path("/tmp/gs_marine_short.profile"))
    ap.add_argument("--taxmap", type=Path, default=Path("/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv"))
    ap.add_argument("--sylph", type=Path, default=Path("/tmp/sylph_marine_sample0/profile.tsv"))
    ap.add_argument("--sample-id", default="marmgCAMI2_short_read_sample_0")
    ap.add_argument("--ctx-k", type=float, default=11.0)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    gold = parse_gold_profile(args.gold, args.sample_id)
    gold_taxids = set(gold["taxid"].astype(str))
    taxmap = parse_species_taxmap(args.taxmap)
    rows = load_minco(args.minco, taxmap, args.ctx_k)
    sylph = load_sylph(args.sylph, taxmap)
    sylph_stats = score_mask(sylph.assign(XnY_ctx=1, Ref_breadth=1.0), sylph["taxid"].astype(bool), gold_taxids)
    sylph_tp = set(sylph.loc[sylph["taxid"].astype(str).isin(gold_taxids), "taxid"].astype(str))

    for col in ["XnY_ctx", "Ref_breadth", "ANI", "Relative_abundance_depth", "Ref_depth_cv", "Ref_mean_depth"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce").fillna(0.0)

    support_arr = rows["XnY_ctx"].to_numpy(dtype=float)
    breadth_arr = rows["Ref_breadth"].to_numpy(dtype=float)
    ani_arr = rows["ANI"].to_numpy(dtype=float)
    rel_depth_arr = rows["Relative_abundance_depth"].to_numpy(dtype=float)
    depth_cv_arr = rows["Ref_depth_cv"].to_numpy(dtype=float)
    mapped = rows["taxid"].astype(bool).to_numpy()
    taxid_values = rows["taxid"].astype(str).to_numpy()
    unique_taxids = sorted({t for t in taxid_values[mapped] if t})
    taxid_to_code = {t: i for i, t in enumerate(unique_taxids)}
    codes = np.full(len(rows), -1, dtype=np.int32)
    for i, taxid in enumerate(taxid_values):
        if taxid:
            codes[i] = taxid_to_code.get(taxid, -1)
    gold_codes = {taxid_to_code[t] for t in gold_taxids if t in taxid_to_code}
    gold_code_arr = np.array(sorted(gold_codes), dtype=np.int32)
    sylph_codes = {taxid_to_code[t] for t in sylph_tp if t in taxid_to_code}

    def fast_score(mask: np.ndarray) -> Dict[str, object]:
        selected_codes = np.unique(codes[mask & (codes >= 0)])
        if selected_codes.size == 0:
            pred_codes: Set[int] = set()
        else:
            pred_codes = set(int(x) for x in selected_codes)
        tp = len(pred_codes & gold_codes)
        fp = len(pred_codes - gold_codes)
        fn = len(gold_taxids) - tp
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        return {
            "pred_taxa": len(pred_codes),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "F1": f1,
            "selected_rows": int(mask.sum()),
            "recovered_sylph_tp": len(pred_codes & sylph_codes),
        }

    support_cuts = [10, 25, 50, 100, 200, 500, 1000, 2000]
    breadth_cuts = [0.0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3]
    ani_cuts = [0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995]
    rel_depth_cuts = [0.0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
    depth_cv_maxes = [float("inf"), 50.0, 25.0, 10.0, 5.0]

    records: List[Dict[str, object]] = []
    for support in support_cuts:
        support_mask = rows["XnY_ctx"] >= support
        for breadth in breadth_cuts:
            breadth_mask = rows["Ref_breadth"] >= breadth
            for ani in ani_cuts:
                ani_mask = rows["ANI"] >= ani
                for rel_depth in rel_depth_cuts:
                    base = (
                        (support_arr >= support)
                        & (breadth_arr >= breadth)
                        & (ani_arr >= ani)
                        & (rel_depth_arr >= rel_depth)
                        & mapped
                    )
                    for depth_cv_max in depth_cv_maxes:
                        mask = base if np.isinf(depth_cv_max) else base & (depth_cv_arr <= depth_cv_max)
                        stats = fast_score(mask)
                        stats.update(
                            {
                                "label": args.label,
                                "support_cut": support,
                                "breadth_cut": breadth,
                                "ani_cut": ani,
                                "rel_depth_cut": rel_depth,
                                "depth_cv_max": "inf" if np.isinf(depth_cv_max) else depth_cv_max,
                                "FP_plus_FN": stats["FP"] + stats["FN"],
                            }
                        )
                        records.append(stats)

    grid = pd.DataFrame(records)
    grid = grid.sort_values(["F1", "FP_plus_FN", "TP"], ascending=[False, True, False])
    grid.to_csv(args.outdir / f"{args.label}.threshold_grid.tsv", sep="\t", index=False)

    best_f1 = grid.head(20).copy()
    best_fpfh = grid.sort_values(["FP_plus_FN", "F1", "TP"], ascending=[True, False, False]).head(20).copy()
    best = pd.concat(
        [
            pd.DataFrame(
                [
                    {
                        "label": "sylph_default",
                        "pred_taxa": sylph_stats["pred_taxa"],
                        "TP": sylph_stats["TP"],
                        "FP": sylph_stats["FP"],
                        "FN": sylph_stats["FN"],
                        "precision": sylph_stats["precision"],
                        "recall": sylph_stats["recall"],
                        "F1": sylph_stats["F1"],
                        "selected_rows": sylph_stats["selected_rows"],
                    }
                ]
            ),
            best_f1.assign(selection="best_F1"),
            best_fpfh.assign(selection="best_FP_plus_FN"),
        ],
        ignore_index=True,
        sort=False,
    )
    best.to_csv(args.outdir / f"{args.label}.threshold_best.tsv", sep="\t", index=False)
    print(best.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
