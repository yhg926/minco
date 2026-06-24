#!/usr/bin/env python3
"""Compare CAMI3 ToyGut coden15 gates against old MinCO and Sylph."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
OLD_RESCUE_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue"
sys.path.insert(0, str(BASE_EXP))

import score_cami3_toygut as base  # noqa: E402


ANI_THRESHOLD = 0.95
FORMULA_EFFECTIVE_CTX_LEN = 24
FORMULA_AF_FLOOR = ANI_THRESHOLD ** FORMULA_EFFECTIVE_CTX_LEN


def output_path(sample: int) -> Path:
    return OUT_DIR / f"cami3_sample{sample}_coden15_ctxmarker_split_naive_product.tsv"


def safe_num(value) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return 0.0
    return out if math.isfinite(out) else 0.0


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def robust_depth(row: pd.Series, exponent: float, median_cutoff: float) -> float:
    median = safe_num(row.get("Reliable_Ref_hit_median_depth"))
    if median >= median_cutoff:
        return median
    mean = safe_num(row.get("Reliable_Ref_mean_depth"))
    af = max(safe_num(row.get("Reliable_Ref_zip_af")), 1e-12)
    return mean / (af**exponent) if mean > 0.0 else 0.0


def add_depth_and_genus(rows: pd.DataFrame, exponent: float, median_cutoff: float) -> pd.DataFrame:
    rows = rows.copy()
    rows["effective_depth_custom"] = [
        robust_depth(row, exponent, median_cutoff) for _, row in rows.iterrows()
    ]
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(species_genus)
    return rows


def coden15_formula_mask(rows: pd.DataFrame) -> pd.Series:
    delta_trigger = (
        (rows["Reliable_Ref_hit_mean_depth"] > 3.0)
        & (rows["Reliable_depth_vmr"] > 50.0)
    )
    delta_pass = ~delta_trigger | (rows["ANI_AF_delta"] < 0.02)
    return (
        (rows["XnY_ctx"] >= 10.0)
        & (rows["ANI_naive_calc"] > ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= FORMULA_AF_FLOOR)
        & delta_pass
        & rows["ncbi_species_taxid"].astype(bool)
    )


def select_robust_rescue(
    rows: pd.DataFrame,
    active_mask: pd.Series,
    value_col: str = "effective_depth_custom",
) -> tuple[pd.DataFrame, list[str]]:
    active = rows.loc[active_mask].copy()
    active["pred_abundance_raw"] = base.numeric(active, value_col)
    active["pred_ani"] = base.numeric(active, "ANI_naive_calc")
    active["support"] = base.numeric(active, "XnY_ctx")
    active["selection_status"] = "active"

    active_genus_max = {}
    for _, row in active.iterrows():
        genus = str(row.get("species_genus", ""))
        value = safe_num(row.get(value_col))
        if genus and value > 0.0:
            active_genus_max[genus] = max(active_genus_max.get(genus, 0.0), value)

    rescue_best: dict[str, pd.Series] = {}
    for _, row in rows.loc[~active_mask].iterrows():
        genus = str(row.get("species_genus", ""))
        active_value = active_genus_max.get(genus, 0.0)
        value = safe_num(row.get(value_col))
        if active_value <= 0.0 or value <= 0.0:
            continue
        if safe_num(row.get("ANI_naive_calc")) < 0.999:
            continue
        if safe_num(row.get("XnY_ctx")) < 100.0:
            continue
        if value < 5.0:
            continue
        if safe_num(row.get("Reliable_Ref_zip_af")) > 0.25:
            continue
        if value < 3.0 * active_value:
            continue
        old = rescue_best.get(genus)
        if old is None or (
            value,
            safe_num(row.get("XnY_ctx")),
        ) > (
            safe_num(old.get(value_col)),
            safe_num(old.get("XnY_ctx")),
        ):
            rescue_best[genus] = row

    rescued_rows = []
    rescued_species = []
    for row in rescue_best.values():
        rec = row.copy()
        rec["pred_abundance_raw"] = safe_num(row.get(value_col))
        rec["pred_ani"] = safe_num(row.get("ANI_naive_calc"))
        rec["support"] = safe_num(row.get("XnY_ctx"))
        rec["selection_status"] = "rescued"
        rescued_rows.append(rec)
        rescued_species.append(str(row.get("gtdb_species", "")))

    if rescued_rows:
        selected = pd.concat([active, pd.DataFrame(rescued_rows)], ignore_index=True)
    else:
        selected = active
    return selected, rescued_species


def score_selected(sample: int, method: str, selected: pd.DataFrame, truth_rows: pd.DataFrame):
    collapsed = base.collapse_predictions(selected)
    gold = set(truth_rows["ncbi_species_taxid"].astype(str))
    pred = set(collapsed["ncbi_species_taxid"].astype(str))
    tp, fp, fn, precision, recall, f1 = base.score_sets(pred, gold)
    score = {
        "sample_id": sample,
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
        {"sample_id": sample, "method": method, **row}
        for row in base.abundance_metrics(collapsed, truth_rows)
    ]
    return score, abundance, collapsed.assign(sample_id=sample, method=method)


def main() -> int:
    by_accession, by_core = base.truth.load_gtdb_metadata(base.truth.GTDB_METADATA)
    score_rows = []
    abundance_rows = []
    selected_rows = []
    rescue_rows = []

    for sample, paths in base.SAMPLES.items():
        truth_rows = base.load_truth(paths["truth"])
        rows = base.load_minco(output_path(sample), by_accession, by_core)

        original_rows = add_depth_and_genus(rows, exponent=1.05, median_cutoff=20.0)
        formula_rows = add_depth_and_genus(rows, exponent=1.0, median_cutoff=10.0)

        configs = [
            (
                "coden15_original_gate_robust_rescue",
                original_rows,
                original_rows["active_gate_pass"],
            ),
            (
                "coden15_formula_af_gate_robust_rescue",
                formula_rows,
                coden15_formula_mask(formula_rows),
            ),
        ]
        for method, use_rows, mask in configs:
            selected, rescued = select_robust_rescue(use_rows, mask)
            score, abundance, collapsed = score_selected(sample, method, selected, truth_rows)
            score_rows.append(score)
            abundance_rows.extend(abundance)
            selected_rows.append(collapsed)
            rescue_rows.append(
                {
                    "sample_id": sample,
                    "method": method,
                    "rescued_species": ",".join(sorted(set(rescued))),
                }
            )

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()
    rescue_df = pd.DataFrame(rescue_rows)

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

    old_summary = pd.read_csv(OLD_RESCUE_EXP / "summary.tsv", sep="\t")
    old_keep = old_summary.loc[
        old_summary["method"].isin(
            [
                "minco_ctxmarker_active_robust_depth_intragenus_rescue",
                "sylph",
            ]
        )
    ].copy()
    rename = {
        "minco_ctxmarker_active_robust_depth_intragenus_rescue": "old_ctxmarker_robust_rescue",
        "sylph": "sylph",
    }
    old_keep["method"] = old_keep["method"].map(rename)

    merged = mean_presence.merge(mean_abundance, on="method", how="outer")
    combined = pd.concat([old_keep, merged], ignore_index=True, sort=False)
    combined = combined.sort_values(["l1_pct_points", "F1"], ascending=[True, False])

    score_df.to_csv(OUT_DIR / "coden15_cami3_sample_presence.tsv", sep="\t", index=False)
    abundance_df.to_csv(OUT_DIR / "coden15_cami3_sample_abundance.tsv", sep="\t", index=False)
    selected_df.to_csv(OUT_DIR / "coden15_cami3_selected_species.tsv", sep="\t", index=False)
    rescue_df.to_csv(OUT_DIR / "coden15_cami3_rescued_species.tsv", sep="\t", index=False)
    mean_presence.to_csv(OUT_DIR / "coden15_cami3_mean_presence.tsv", sep="\t", index=False)
    mean_abundance.to_csv(OUT_DIR / "coden15_cami3_mean_abundance_renorm.tsv", sep="\t", index=False)
    combined.to_csv(OUT_DIR / "summary.tsv", sep="\t", index=False)
    print(combined.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
