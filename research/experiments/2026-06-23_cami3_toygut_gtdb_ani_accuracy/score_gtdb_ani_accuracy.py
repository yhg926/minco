#!/usr/bin/env python3
"""ANI-threshold accuracy on CAMI3 ToyGut using GTDB-mapped truth."""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = Path(__file__).resolve().parent
BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
RESCUE_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue"
C15_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_coden15_formula_af"

sys.path.insert(0, str(BASE_EXP))
sys.path.insert(0, str(RESCUE_EXP))
sys.path.insert(0, str(C15_EXP))

import score_cami3_toygut as base  # noqa: E402
import score_cami3_abundance_rescue as old_rescue  # noqa: E402
import score_cami3_coden15_formula_af as c15_score  # noqa: E402


ANI_THRESHOLD = 0.95


def safe_num(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return 0.0
    return out if math.isfinite(out) else 0.0


def metadata_taxid_map(by_accession: dict[str, dict[str, str]]) -> dict[str, set[str]]:
    taxid_to_species: dict[str, set[str]] = defaultdict(set)
    seen_records: set[tuple[str, str]] = set()
    for rec in by_accession.values():
        key = (str(rec.get("metadata_accession", "")), str(rec.get("gtdb_species", "")))
        if key in seen_records:
            continue
        seen_records.add(key)
        species = str(rec.get("gtdb_species", ""))
        if not species:
            continue
        for field in ("ncbi_species_taxid", "ncbi_taxid"):
            taxid = str(rec.get(field, "") or "")
            if taxid:
                taxid_to_species[taxid].add(species)
    return taxid_to_species


def truth_gtdb_for_sample(
    sample: int,
    taxid_to_species: dict[str, set[str]],
) -> tuple[set[str], set[str], list[dict[str, object]]]:
    truth_rows = base.load_truth(base.SAMPLES[sample]["truth"])
    unique_species: set[str] = set()
    possible_species: set[str] = set()
    audit_rows: list[dict[str, object]] = []
    for _, row in truth_rows.iterrows():
        taxid = str(row["ncbi_species_taxid"])
        species_set = sorted(taxid_to_species.get(taxid, set()))
        possible_species.update(species_set)
        status = (
            "unique"
            if len(species_set) == 1
            else "ambiguous"
            if len(species_set) > 1
            else "unmapped"
        )
        if status == "unique":
            unique_species.add(species_set[0])
        audit_rows.append(
            {
                "sample_id": sample,
                "ncbi_species_taxid": taxid,
                "truth_name": row.get("truth_name", ""),
                "truth_abundance_bacterial_norm": safe_num(
                    row.get("truth_abundance_bacterial_norm")
                ),
                "gtdb_mapping_status": status,
                "gtdb_species_count": len(species_set),
                "gtdb_species": ",".join(species_set),
            }
        )
    return unique_species, possible_species, audit_rows


def collapse_by_gtdb(selected: pd.DataFrame) -> pd.DataFrame:
    if selected.empty:
        return pd.DataFrame(
            columns=[
                "gtdb_species",
                "ncbi_species_taxid",
                "ncbi_species",
                "accession",
                "pred_ani",
                "support",
                "pred_abundance_raw",
                "selection_status",
            ]
        )
    rows = selected.loc[selected["gtdb_species"].astype(bool)].copy()
    if rows.empty:
        return pd.DataFrame(columns=["gtdb_species", "pred_ani", "support"])
    for col in ("pred_ani", "support", "pred_abundance_raw"):
        rows[col] = base.numeric(rows, col)
    sort_cols = ["gtdb_species", "pred_ani", "support", "pred_abundance_raw"]
    rows = rows.sort_values(sort_cols, ascending=[True, False, False, False])
    agg = {
        "ncbi_species_taxid": ("ncbi_species_taxid", "first")
        if "ncbi_species_taxid" in rows
        else ("gtdb_species", "first"),
        "ncbi_species": ("ncbi_species", "first")
        if "ncbi_species" in rows
        else ("gtdb_species", "first"),
        "accession": ("accession", "first") if "accession" in rows else ("gtdb_species", "first"),
        "pred_ani": ("pred_ani", "max"),
        "support": ("support", "max"),
        "pred_abundance_raw": ("pred_abundance_raw", "max"),
    }
    if "selection_status" in rows:
        agg["selection_status"] = ("selection_status", "first")
    return rows.groupby("gtdb_species", as_index=False).agg(**agg)


def score_sets(pred: set[str], truth: set[str]) -> tuple[int, int, int, float, float, float]:
    tp = len(pred & truth)
    fp = len(pred - truth)
    fn = len(truth - pred)
    precision = tp / len(pred) if pred else 0.0
    recall = tp / len(truth) if truth else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def selected_ani_metrics(
    sample: int,
    method: str,
    selected: pd.DataFrame,
    truth_species: set[str],
    possible_truth_species: set[str],
) -> tuple[dict[str, object], pd.DataFrame]:
    collapsed = collapse_by_gtdb(selected)
    pred_species = set(collapsed["gtdb_species"].astype(str)) if not collapsed.empty else set()
    tp, fp, fn, precision, recall, f1 = score_sets(pred_species, truth_species)

    if collapsed.empty:
        metric = {
            "sample_id": sample,
            "method": method,
            "truth_gtdb_species_unique": len(truth_species),
            "truth_gtdb_species_possible": len(possible_truth_species),
            "pred_gtdb_species": 0,
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "precision": precision,
            "recall": recall,
            "F1": f1,
            "pred_in_possible_truth": 0,
            "pred_outside_possible_truth": 0,
            "possible_precision": 0.0,
            "tp_ani_mean": 0.0,
            "tp_ani_median": 0.0,
            "tp_ani_min": 0.0,
            "tp_ani_below_095": 0,
            "fp_ani_mean": 0.0,
            "fp_ani_median": 0.0,
            "fp_ani_min": 0.0,
            "fp_ani_ge_095": 0,
            "outside_possible_ani_mean": 0.0,
            "outside_possible_ani_ge_095": 0,
        }
        return metric, collapsed

    collapsed["is_truth_unique"] = collapsed["gtdb_species"].astype(str).isin(truth_species)
    collapsed["is_truth_possible"] = collapsed["gtdb_species"].astype(str).isin(
        possible_truth_species
    )
    collapsed["ani_ge_095"] = base.numeric(collapsed, "pred_ani") >= ANI_THRESHOLD
    tp_rows = collapsed.loc[collapsed["is_truth_unique"]].copy()
    fp_rows = collapsed.loc[~collapsed["is_truth_unique"]].copy()
    outside_possible_rows = collapsed.loc[~collapsed["is_truth_possible"]].copy()
    pred_in_possible = int(collapsed["is_truth_possible"].sum())
    pred_outside_possible = len(collapsed) - pred_in_possible

    def stat(rows: pd.DataFrame, fn: str) -> float:
        if rows.empty:
            return 0.0
        return safe_num(getattr(base.numeric(rows, "pred_ani"), fn)())

    metric = {
        "sample_id": sample,
        "method": method,
        "truth_gtdb_species_unique": len(truth_species),
        "truth_gtdb_species_possible": len(possible_truth_species),
        "pred_gtdb_species": len(pred_species),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "pred_in_possible_truth": pred_in_possible,
        "pred_outside_possible_truth": pred_outside_possible,
        "possible_precision": pred_in_possible / len(collapsed) if len(collapsed) else 0.0,
        "tp_ani_mean": stat(tp_rows, "mean"),
        "tp_ani_median": stat(tp_rows, "median"),
        "tp_ani_min": stat(tp_rows, "min"),
        "tp_ani_below_095": int((base.numeric(tp_rows, "pred_ani") < ANI_THRESHOLD).sum()),
        "fp_ani_mean": stat(fp_rows, "mean"),
        "fp_ani_median": stat(fp_rows, "median"),
        "fp_ani_min": stat(fp_rows, "min"),
        "fp_ani_ge_095": int((base.numeric(fp_rows, "pred_ani") >= ANI_THRESHOLD).sum()),
        "outside_possible_ani_mean": stat(outside_possible_rows, "mean"),
        "outside_possible_ani_ge_095": int(
            (base.numeric(outside_possible_rows, "pred_ani") >= ANI_THRESHOLD).sum()
        ),
    }
    collapsed.insert(0, "method", method)
    collapsed.insert(0, "sample_id", sample)
    return metric, collapsed


def candidate_species_table(raw: pd.DataFrame, ani_col: str, support_col: str) -> pd.DataFrame:
    rows = raw.loc[raw["gtdb_species"].astype(bool)].copy()
    if rows.empty:
        return pd.DataFrame(columns=["gtdb_species", "pred_ani", "support"])
    rows["pred_ani"] = base.numeric(rows, ani_col)
    rows["support"] = base.numeric(rows, support_col)
    rows = rows.sort_values(
        ["gtdb_species", "pred_ani", "support"], ascending=[True, False, False]
    )
    return rows.groupby("gtdb_species", as_index=False).agg(
        pred_ani=("pred_ani", "max"),
        support=("support", "max"),
        accession=("accession", "first") if "accession" in rows else ("gtdb_species", "first"),
        ncbi_species_taxid=("ncbi_species_taxid", "first")
        if "ncbi_species_taxid" in rows
        else ("gtdb_species", "first"),
        ncbi_species=("ncbi_species", "first")
        if "ncbi_species" in rows
        else ("gtdb_species", "first"),
    )


def candidate_ani95_metrics(
    sample: int,
    method: str,
    raw: pd.DataFrame,
    ani_col: str,
    support_col: str,
    truth_species: set[str],
    possible_truth_species: set[str],
) -> tuple[dict[str, object], pd.DataFrame]:
    table = candidate_species_table(raw, ani_col, support_col)
    ani_pos = set(
        table.loc[base.numeric(table, "pred_ani") >= ANI_THRESHOLD, "gtdb_species"].astype(str)
    )
    tp, fp, fn, precision, recall, f1 = score_sets(ani_pos, truth_species)
    table["is_truth_unique"] = table["gtdb_species"].astype(str).isin(truth_species)
    table["is_truth_possible"] = table["gtdb_species"].astype(str).isin(possible_truth_species)
    table["ani_ge_095"] = base.numeric(table, "pred_ani") >= ANI_THRESHOLD
    ani_pos_possible = set(table.loc[table["ani_ge_095"] & table["is_truth_possible"], "gtdb_species"].astype(str))
    ani_pos_outside_possible = set(table.loc[table["ani_ge_095"] & ~table["is_truth_possible"], "gtdb_species"].astype(str))
    table.insert(0, "method", method)
    table.insert(0, "sample_id", sample)
    return (
        {
            "sample_id": sample,
            "method": method,
            "truth_gtdb_species_unique": len(truth_species),
            "truth_gtdb_species_possible": len(possible_truth_species),
            "candidate_gtdb_species": len(table),
            "candidate_truth_species_seen": int(table["is_truth_unique"].sum()),
            "candidate_possible_species_seen": int(table["is_truth_possible"].sum()),
            "ani_ge_095_species": len(ani_pos),
            "ANI95_TP": tp,
            "ANI95_FP": fp,
            "ANI95_FN": fn,
            "ANI95_precision": precision,
            "ANI95_recall": recall,
            "ANI95_F1": f1,
            "ANI95_possible_hits": len(ani_pos_possible),
            "ANI95_outside_possible": len(ani_pos_outside_possible),
            "ANI95_possible_precision": len(ani_pos_possible) / len(ani_pos)
            if ani_pos
            else 0.0,
        },
        table,
    )


def load_old_selected_and_raw(sample: int, by_accession, by_core) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = base.load_minco(old_rescue.MINCO_OUTPUTS["ctxmarker"][sample], by_accession, by_core)
    raw = old_rescue.add_effective_depth(raw)
    selected, _ = old_rescue.select_minco_robust_rescue(raw)
    return selected, raw


def load_coden15_selected_and_raw(
    sample: int,
    method: str,
    by_accession,
    by_core,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = base.load_minco(c15_score.output_path(sample), by_accession, by_core)
    if method == "coden15_original_gate_robust_rescue":
        rows = c15_score.add_depth_and_genus(raw, exponent=1.05, median_cutoff=20.0)
        mask = rows["active_gate_pass"]
    elif method == "coden15_formula_af_gate_robust_rescue":
        rows = c15_score.add_depth_and_genus(raw, exponent=1.0, median_cutoff=10.0)
        mask = c15_score.coden15_formula_mask(rows)
    else:
        raise ValueError(method)
    selected, _ = c15_score.select_robust_rescue(rows, mask)
    return selected, rows


def load_sylph_selected_and_raw(sample: int, by_accession, by_core) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = base.load_sylph(base.SAMPLES[sample]["sylph"], by_accession, by_core)
    selected = base.select_predictions("sylph", raw)
    return selected, raw


def main() -> int:
    by_accession, by_core = base.truth.load_gtdb_metadata(base.truth.GTDB_METADATA)
    taxid_to_species = metadata_taxid_map(by_accession)

    truth_audit_rows = []
    selected_metric_rows = []
    selected_detail_frames = []
    candidate_metric_rows = []
    candidate_detail_frames = []

    for sample in sorted(base.SAMPLES):
        truth_species, possible_truth_species, audit = truth_gtdb_for_sample(
            sample, taxid_to_species
        )
        truth_audit_rows.extend(audit)

        loaders = [
            (
                "sylph",
                lambda s=sample: load_sylph_selected_and_raw(s, by_accession, by_core),
                "Adjusted_ANI",
                "Eff_cov",
            ),
            (
                "old_ctxmarker_robust_rescue",
                lambda s=sample: load_old_selected_and_raw(s, by_accession, by_core),
                "ANI_naive_calc",
                "XnY_ctx",
            ),
            (
                "coden15_formula_af_gate_robust_rescue",
                lambda s=sample: load_coden15_selected_and_raw(
                    s, "coden15_formula_af_gate_robust_rescue", by_accession, by_core
                ),
                "ANI_naive_calc",
                "XnY_ctx",
            ),
            (
                "coden15_original_gate_robust_rescue",
                lambda s=sample: load_coden15_selected_and_raw(
                    s, "coden15_original_gate_robust_rescue", by_accession, by_core
                ),
                "ANI_naive_calc",
                "XnY_ctx",
            ),
        ]
        for method, loader, ani_col, support_col in loaders:
            selected, raw = loader()
            selected_metric, selected_detail = selected_ani_metrics(
                sample, method, selected, truth_species, possible_truth_species
            )
            selected_metric_rows.append(selected_metric)
            selected_detail_frames.append(selected_detail)

            candidate_metric, candidate_detail = candidate_ani95_metrics(
                sample, method, raw, ani_col, support_col, truth_species, possible_truth_species
            )
            candidate_metric_rows.append(candidate_metric)
            candidate_detail_frames.append(candidate_detail)

    truth_audit = pd.DataFrame(truth_audit_rows)
    selected_metrics = pd.DataFrame(selected_metric_rows)
    candidate_metrics = pd.DataFrame(candidate_metric_rows)
    selected_details = pd.concat(selected_detail_frames, ignore_index=True)
    candidate_details = pd.concat(candidate_detail_frames, ignore_index=True)

    selected_mean = selected_metrics.groupby("method", as_index=False).mean(numeric_only=True)
    candidate_mean = candidate_metrics.groupby("method", as_index=False).mean(numeric_only=True)

    truth_summary = (
        truth_audit.groupby(["sample_id", "gtdb_mapping_status"], as_index=False)
        .agg(
            ncbi_species=("ncbi_species_taxid", "count"),
            abundance=("truth_abundance_bacterial_norm", "sum"),
        )
        .sort_values(["sample_id", "gtdb_mapping_status"])
    )

    truth_audit.to_csv(OUT_DIR / "gtdb_truth_transfer_audit.tsv", sep="\t", index=False)
    truth_summary.to_csv(OUT_DIR / "gtdb_truth_transfer_summary.tsv", sep="\t", index=False)
    selected_metrics.to_csv(OUT_DIR / "selected_call_gtdb_ani_metrics.tsv", sep="\t", index=False)
    selected_mean.to_csv(OUT_DIR / "selected_call_gtdb_ani_mean_metrics.tsv", sep="\t", index=False)
    selected_details.to_csv(OUT_DIR / "selected_call_gtdb_ani_details.tsv", sep="\t", index=False)
    candidate_metrics.to_csv(OUT_DIR / "candidate_gtdb_ani95_metrics.tsv", sep="\t", index=False)
    candidate_mean.to_csv(OUT_DIR / "candidate_gtdb_ani95_mean_metrics.tsv", sep="\t", index=False)
    candidate_details.to_csv(OUT_DIR / "candidate_gtdb_ani95_details.tsv", sep="\t", index=False)

    print("GTDB truth transfer summary")
    print(truth_summary.to_csv(sep="\t", index=False), end="")
    print("\nSelected-call GTDB ANI mean metrics")
    print(selected_mean.to_csv(sep="\t", index=False), end="")
    print("\nCandidate ANI>=0.95 GTDB mean metrics")
    print(candidate_mean.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
