#!/usr/bin/env python3
"""Score Toy Mouse sample0 direct readwise calls with the source-aware taxmap."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(
    0,
    str(REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"),
)
sys.path.insert(
    0,
    str(REPO_ROOT / "research/experiments/2026-06-21_minco_multisample_call_calibration/scripts"),
)

from analyze_readwise_corrections import extract_accession, load_minco, parse_gold_profile, parse_species_taxmap
from calibrate_multisample_calls import direct_mask, score_prediction


GOLD_PROFILE = Path("/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/taxonomic_profile_0.txt")
SOURCE_AWARE_TAXMAP = Path(
    "/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/"
    "mouse0_sourceaware_taxmap_ani95_minaf50.tsv"
)

METHODS = [
    (
        "s1000_unique_direct",
        Path(
            "/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/"
            "toymouse_sample0_s1000_gtdb_unique_zip_unfiltered.tsv"
        ),
    ),
    (
        "s1000_split_direct",
        Path(
            "/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/"
            "toymouse_sample0_s1000_gtdb_split_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_unique_direct",
        Path(
            "/tmp/gtdb232_s2000_dedup_marker.qKJofv/"
            "toymouse_sample0_s2000_dedup_ctxmarker_unique_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_split_direct",
        Path(
            "/tmp/gtdb232_s2000_dedup_marker.qKJofv/"
            "toymouse_sample0_s2000_dedup_ctxmarker_split_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_split_naive_unfiltered",
        REPO_ROOT
        / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/"
        "toymouse_sample0_s2000_marker_split_naive_unfiltered.tsv",
    ),
    (
        "s2000_marker_split_zip_poisson_depth_p005",
        REPO_ROOT
        / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/"
        "toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.tsv",
    ),
    (
        "s2000_marker_split_naive_poisson_depth_p005",
        REPO_ROOT
        / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/"
        "toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.tsv",
    ),
    (
        "s2000_marker_split_naive_poisson_product_p005",
        REPO_ROOT
        / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise/"
        "toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.tsv",
    ),
]

TARGET_ACCESSION = "GCF_018987235.1"
EPSILON = 1e-8


def minco_naive_ani(xny_ctx: float, n_diff_obj: float, n_diff_obj_section: float) -> float:
    if xny_ctx == 0:
        return 0.0
    ratio = (n_diff_obj_section + EPSILON) / (n_diff_obj + EPSILON)
    dist0 = n_diff_obj / (xny_ctx + n_diff_obj)
    final_dist = 1.0 - pow(1.0 - dist0, ratio)
    if final_dist <= 0.0:
        naive_dist = 0.0
    else:
        naive_dist = final_dist * 0.1544286
    return max(0.0, min(1.0, 1.0 - naive_dist))


def main() -> int:
    out_dir = Path(__file__).resolve().parent
    taxmap = parse_species_taxmap(SOURCE_AWARE_TAXMAP)
    gold = parse_gold_profile(GOLD_PROFILE, "0")
    gold_taxids = set(gold["taxid"].astype(str))

    score_rows = []
    target_rows = []
    for method, path in METHODS:
        rows = load_minco(path, taxmap, 11.0)
        rows = rows.loc[rows["taxid"].astype(bool)].copy()
        selected = rows.loc[direct_mask(rows), "taxid"].astype(str)
        stats = score_prediction(selected, gold_taxids)
        score_rows.append(
            {
                "method": method,
                "gold_taxa": len(gold_taxids),
                "pred_taxa": stats["pred_taxa"],
                "TP": stats["TP"],
                "FP": stats["FP"],
                "FN": stats["FN"],
                "precision": f"{stats['precision']:.12g}",
                "recall": f"{stats['recall']:.12g}",
                "F1": f"{stats['F1']:.12g}",
            }
        )

        work = rows.copy()
        work["accession"] = work["Ref"].map(extract_accession)
        target = work.loc[work["accession"] == TARGET_ACCESSION]
        if target.empty:
            continue
        rec = target.iloc[0]
        naive_ani = minco_naive_ani(
            float(rec["XnY_ctx"]),
            float(rec["N_diff_obj"]),
            float(rec["N_diff_obj_section"]),
        )
        target_rows.append(
            {
                "method": method,
                "accession": TARGET_ACCESSION,
                "selected_zip_aaf_ANI": f"{float(rec['ANI']):.6f}",
                "readwise_naive_ANI_calc": f"{naive_ani:.12g}",
                "XnY_ctx": int(rec["XnY_ctx"]),
                "N_diff_obj": int(rec["N_diff_obj"]),
                "N_diff_obj_section": int(rec["N_diff_obj_section"]),
                "N_mut2_ctx": int(rec["N_mut2_ctx"]),
                "Ref_zip_af": f"{float(rec.get('Ref_zip_af', 0.0)):.6f}",
                "Ref_breadth": f"{float(rec.get('Ref_breadth', 0.0)):.6f}",
                "Ref_mean_depth": f"{float(rec.get('Ref_mean_depth', 0.0)):.6f}",
                "Normalized_abundance_depth": f"{float(rec.get('Normalized_abundance_depth', 0.0)):.6f}",
                "Ref_zip_aaf_ani": f"{float(rec.get('Ref_zip_aaf_ani', 0.0)):.6f}",
                "Default_call": rec.get("Default_call", ""),
            }
        )

    with (out_dir / "sourceaware_direct_scores.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(score_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(score_rows)

    with (out_dir / "l_crispatus_markerdb_comparison.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(target_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(target_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
