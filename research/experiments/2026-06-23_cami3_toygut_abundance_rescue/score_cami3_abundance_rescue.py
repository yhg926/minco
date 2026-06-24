#!/usr/bin/env python3
"""Score Toy Mouse abundance-rescue rule on CAMI3 ToyGut samples."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
sys.path.insert(0, str(BASE_EXP))

import score_cami3_toygut as base  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
RUN_DIR = Path("/tmp/cami3_toygut_abundance_rescue_20260623")

MEDIAN_DEPTH_CUTOFF = 20.0
AF_EXPONENT = 1.05
RESCUE_ANI_MIN = 0.999
RESCUE_XNY_MIN = 100.0
RESCUE_EFFECTIVE_MIN = 5.0
RESCUE_ZIP_AF_MAX = 0.25
RESCUE_ACTIVE_GENUS_RATIO = 3.0

MINCO_OUTPUTS = {
    "ctxmarker": {
        0: RUN_DIR / "minco_s0_ctxmarker_current_binary.tsv",
        1: RUN_DIR / "minco_s1_ctxmarker_current_binary.tsv",
        2: RUN_DIR / "minco_s2_ctxmarker_current_binary.tsv",
    },
    "ctxobj": {
        0: RUN_DIR / "minco_s0_ctxobj_current_binary.tsv",
        1: RUN_DIR / "minco_s1_ctxobj_current_binary.tsv",
        2: RUN_DIR / "minco_s2_ctxobj_current_binary.tsv",
    },
}

RESCUE_COLUMNS = [
    "sample_id",
    "method",
    "ref_kind",
    "gtdb_species",
    "ncbi_species_taxid",
    "ncbi_species",
    "accession",
    "genus",
    "effective_depth",
    "active_genus_max_effective_depth",
    "XnY_ctx",
    "ANI_naive_calc",
    "Reliable_Ref_zip_af",
    "Reliable_Ref_breadth",
    "Reliable_Ref_hit_mean_depth",
    "Reliable_Ref_hit_median_depth",
]


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def safe_num(value) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return 0.0
    return out if math.isfinite(out) else 0.0


def robust_effective_depth(row: pd.Series) -> float:
    emitted = safe_num(row.get("Effective_abundance_depth"))
    if emitted > 0.0:
        return emitted
    median = safe_num(row.get("Reliable_Ref_hit_median_depth"))
    if median >= MEDIAN_DEPTH_CUTOFF:
        return median
    mean = safe_num(row.get("Reliable_Ref_mean_depth"))
    af = max(safe_num(row.get("Reliable_Ref_zip_af")), 1e-12)
    return mean / (af**AF_EXPONENT) if mean > 0.0 else 0.0


def add_effective_depth(rows: pd.DataFrame) -> pd.DataFrame:
    rows = rows.copy()
    rows["effective_depth"] = rows.apply(robust_effective_depth, axis=1)
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(species_genus)
    return rows


def select_minco_old(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    selected = rows.loc[rows["active_gate_pass"]].copy()
    selected["pred_abundance_raw"] = base.numeric(selected, "Normalized_abundance_depth")
    selected["pred_ani"] = base.numeric(selected, "ANI_naive_calc")
    selected["support"] = base.numeric(selected, "XnY_ctx")
    selected["selection_status"] = "active"
    return selected, []


def select_minco_robust(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    selected = rows.loc[rows["active_gate_pass"]].copy()
    selected["pred_abundance_raw"] = base.numeric(selected, "effective_depth")
    selected["pred_ani"] = base.numeric(selected, "ANI_naive_calc")
    selected["support"] = base.numeric(selected, "XnY_ctx")
    selected["selection_status"] = "active"
    return selected, []


def select_minco_robust_rescue(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    active = rows.loc[rows["active_gate_pass"]].copy()
    active["pred_abundance_raw"] = base.numeric(active, "effective_depth")
    active["pred_ani"] = base.numeric(active, "ANI_naive_calc")
    active["support"] = base.numeric(active, "XnY_ctx")
    active["selection_status"] = "active"

    active_genus_max: dict[str, float] = {}
    for _, row in active.iterrows():
        genus = str(row.get("species_genus", ""))
        value = safe_num(row.get("effective_depth"))
        if genus and value > 0.0:
            active_genus_max[genus] = max(active_genus_max.get(genus, 0.0), value)

    rescue_best: dict[str, pd.Series] = {}
    for _, row in rows.loc[~rows["active_gate_pass"]].iterrows():
        genus = str(row.get("species_genus", ""))
        active_value = active_genus_max.get(genus, 0.0)
        value = safe_num(row.get("effective_depth"))
        if active_value <= 0.0 or value <= 0.0:
            continue
        if safe_num(row.get("ANI_naive_calc")) < RESCUE_ANI_MIN:
            continue
        if safe_num(row.get("XnY_ctx")) < RESCUE_XNY_MIN:
            continue
        if value < RESCUE_EFFECTIVE_MIN:
            continue
        if safe_num(row.get("Reliable_Ref_zip_af")) > RESCUE_ZIP_AF_MAX:
            continue
        if value < RESCUE_ACTIVE_GENUS_RATIO * active_value:
            continue
        old = rescue_best.get(genus)
        if old is None or (
            value,
            safe_num(row.get("XnY_ctx")),
        ) > (
            safe_num(old.get("effective_depth")),
            safe_num(old.get("XnY_ctx")),
        ):
            rescue_best[genus] = row

    rescued_rows = []
    rescue_details = []
    for row in rescue_best.values():
        rec = row.copy()
        rec["pred_abundance_raw"] = safe_num(row.get("effective_depth"))
        rec["pred_ani"] = safe_num(row.get("ANI_naive_calc"))
        rec["support"] = safe_num(row.get("XnY_ctx"))
        rec["selection_status"] = "rescued"
        rescued_rows.append(rec)
        rescue_details.append(
            {
                "gtdb_species": row.get("gtdb_species", ""),
                "ncbi_species_taxid": row.get("ncbi_species_taxid", ""),
                "ncbi_species": row.get("ncbi_species", ""),
                "accession": row.get("accession", ""),
                "genus": row.get("species_genus", ""),
                "effective_depth": safe_num(row.get("effective_depth")),
                "active_genus_max_effective_depth": active_genus_max.get(
                    str(row.get("species_genus", "")), 0.0
                ),
                "XnY_ctx": safe_num(row.get("XnY_ctx")),
                "ANI_naive_calc": safe_num(row.get("ANI_naive_calc")),
                "Reliable_Ref_zip_af": safe_num(row.get("Reliable_Ref_zip_af")),
                "Reliable_Ref_breadth": safe_num(row.get("Reliable_Ref_breadth")),
                "Reliable_Ref_hit_mean_depth": safe_num(row.get("Reliable_Ref_hit_mean_depth")),
                "Reliable_Ref_hit_median_depth": safe_num(
                    row.get("Reliable_Ref_hit_median_depth")
                ),
            }
        )

    if rescued_rows:
        selected = pd.concat([active, pd.DataFrame(rescued_rows)], ignore_index=True)
    else:
        selected = active
    return selected, rescue_details


def select_sylph(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    return base.select_predictions("sylph", rows), []


def score_method(sample_id: int, method: str, selected: pd.DataFrame, truth_rows: pd.DataFrame):
    collapsed = base.collapse_predictions(selected)
    gold = set(truth_rows["ncbi_species_taxid"].astype(str))
    pred = set(collapsed["ncbi_species_taxid"].astype(str))
    tp, fp, fn, precision, recall, f1 = base.score_sets(pred, gold)
    score = {
        "sample_id": sample_id,
        "method": method,
        "truth_bacterial_species": len(gold),
        "pred_species": len(pred),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "raw_selected_rows": len(selected),
        "mapped_selected_rows": int(selected["ncbi_species_taxid"].astype(bool).sum())
        if "ncbi_species_taxid" in selected
        else len(selected),
    }
    abundance = [
        {"sample_id": sample_id, "method": method, **row}
        for row in base.abundance_metrics(collapsed, truth_rows)
    ]
    details = []
    details.extend(base.detail_rows(sample_id, method, "FP", fp, collapsed, truth_rows))
    details.extend(base.detail_rows(sample_id, method, "FN", fn, collapsed, truth_rows))
    selected_out = collapsed.assign(sample_id=sample_id, method=method)
    return score, abundance, details, selected_out


def main() -> int:
    by_accession, by_core = base.truth.load_gtdb_metadata(base.truth.GTDB_METADATA)
    score_rows = []
    abundance_rows = []
    detail_rows = []
    selected_rows = []
    rescue_rows = []

    for sample_id, paths in base.SAMPLES.items():
        truth_rows = base.load_truth(paths["truth"])
        sylph = base.load_sylph(paths["sylph"], by_accession, by_core)

        for ref_kind, sample_paths in MINCO_OUTPUTS.items():
            rows = base.load_minco(sample_paths[sample_id], by_accession, by_core)
            rows = add_effective_depth(rows)
            method_selectors = [
                (f"minco_{ref_kind}_active_old_norm", select_minco_old),
                (f"minco_{ref_kind}_active_robust_depth", select_minco_robust),
                (
                    f"minco_{ref_kind}_active_robust_depth_intragenus_rescue",
                    select_minco_robust_rescue,
                ),
            ]
            for method, selector in method_selectors:
                selected, rescued = selector(rows)
                for rec in rescued:
                    rec["sample_id"] = sample_id
                    rec["method"] = method
                    rec["ref_kind"] = ref_kind
                    rescue_rows.append(rec)
                score, abundance, details, selected_out = score_method(
                    sample_id, method, selected, truth_rows
                )
                score_rows.append(score)
                abundance_rows.extend(abundance)
                detail_rows.extend(details)
                selected_rows.append(selected_out)

        selected, _ = select_sylph(sylph)
        score, abundance, details, selected_out = score_method(
            sample_id, "sylph", selected, truth_rows
        )
        score_rows.append(score)
        abundance_rows.extend(abundance)
        detail_rows.extend(details)
        selected_rows.append(selected_out)

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    detail_df = pd.DataFrame(detail_rows)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()
    rescue_df = pd.DataFrame(rescue_rows, columns=RESCUE_COLUMNS)

    mean_presence = (
        score_df.groupby("method", as_index=False)[
            ["pred_species", "TP", "FP", "FN", "precision", "recall", "F1"]
        ].mean()
    )
    mean_abundance = (
        abundance_df.loc[abundance_df["renorm_pred"]]
        .groupby("method", as_index=False)[
            [
                "pred_sum_on_truth",
                "pred_sum_all",
                "pearson",
                "spearman",
                "mae_pct_points",
                "l1_pct_points",
            ]
        ]
        .mean()
    )

    score_df.to_csv(OUT_DIR / "sample_presence.tsv", sep="\t", index=False)
    abundance_df.to_csv(OUT_DIR / "sample_abundance.tsv", sep="\t", index=False)
    detail_df.to_csv(OUT_DIR / "fp_fn_details.tsv", sep="\t", index=False)
    selected_df.to_csv(OUT_DIR / "selected_species.tsv", sep="\t", index=False)
    rescue_df.to_csv(OUT_DIR / "rescued_species.tsv", sep="\t", index=False)
    mean_presence.to_csv(OUT_DIR / "mean_presence.tsv", sep="\t", index=False)
    mean_abundance.to_csv(OUT_DIR / "mean_abundance_renorm.tsv", sep="\t", index=False)

    merged = mean_presence.merge(mean_abundance, on="method", how="outer")
    merged = merged.sort_values(["l1_pct_points", "F1"], ascending=[True, False])
    merged.to_csv(OUT_DIR / "summary.tsv", sep="\t", index=False)
    print(merged.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
