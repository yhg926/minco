#!/usr/bin/env python3
"""Score Toy Mouse sample0 on the S2000 pairwise context markerdb."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import pandas as pd

import build_and_score_gtdb_ground_truth as truth
import score_reliable_truncated_depth_gate as reliable


EXP_DIR = Path(__file__).resolve().parent
RUN_DIR = Path("/tmp/gtdb232_s2000_pairctx_af005.YHvD9c")
ZIP_TSV = RUN_DIR / "toymouse_sample0_pairctx_af005_split_zip_unfiltered.tsv"
PRODUCT_NB_TSV = RUN_DIR / "toymouse_sample0_pairctx_af005_split_naive_product_nb_p005.tsv"
PRODUCT_TOPFRAC_MEDIAN_TSV = RUN_DIR / "toymouse_sample0_pairctx_af005_split_naive_product_topfrac_median025.tsv"
OLD_PRODUCT_TOPFRAC_MEDIAN_TSV = EXP_DIR / "toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv"
SCORE_OUT = RUN_DIR / "pairctx_af005_gtdb_scores.tsv"
DETAIL_OUT = RUN_DIR / "pairctx_af005_gtdb_score_details.tsv"
ABUND_OUT = RUN_DIR / "pairctx_af005_gtdb_abundance.tsv"

ANI_THRESHOLD = 0.95
EFFECTIVE_CTX_LENGTH = 22 + 2
AF_FLOOR = math.exp((ANI_THRESHOLD - 1.0) * EFFECTIVE_CTX_LENGTH)
ACTIVE_CTX_MIN = 15.0
ACTIVE_RELIABLE_ZTP_AF_FLOOR = 0.40
ACTIVE_MEAN_DEPTH_MIN = 3.0
ACTIVE_VMR_MIN = 50.0
ACTIVE_DELTA_MAX = 0.03


def score_sets(predicted: Iterable[str], gold: set[str]) -> tuple[set[str], set[str], set[str], float, float, float]:
    pred = {str(x) for x in predicted if str(x)}
    tp = pred & gold
    fp = pred - gold
    fn = gold - pred
    precision = len(tp) / len(pred) if pred else 0.0
    recall = len(tp) / len(gold) if gold else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def numeric(rows: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(rows[col], errors="coerce").fillna(0.0)


def load_rows(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = truth.load_minco_with_gtdb_species(path, by_accession, by_core)
    rows = truth.add_naive_ani(rows)
    return rows


def add_reliable_adjusted_af(rows: pd.DataFrame, fit_ztnb: bool = False) -> pd.DataFrame:
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_ctx",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_depth_variance",
        "Normalized_abundance_depth",
    ]:
        rows[col] = numeric(rows, col)

    ztp_af = []
    ztnb_af = []
    ztnb_model = []
    for breadth, mean_pos, var_pos, hit_ctx in zip(
        rows["Reliable_Ref_breadth"],
        rows["Reliable_Ref_hit_mean_depth"],
        rows["Reliable_Ref_hit_depth_variance"],
        rows["Reliable_Ref_hit_ctx"],
    ):
        ztp_af.append(reliable.ztp_adjusted_af(float(breadth), float(mean_pos)))
        if fit_ztnb:
            p_nonzero, model = reliable.ztnb_p_nonzero(
                float(mean_pos), float(var_pos), int(hit_ctx)
            )
            ztnb_model.append(model)
            ztnb_af.append(min(1.0, float(breadth) / p_nonzero) if p_nonzero > 0.0 else 0.0)
        else:
            ztnb_model.append("not_fit")
            ztnb_af.append(float("nan"))
    rows["Reliable_ztp_af"] = ztp_af
    rows["Reliable_ztnb_af"] = ztnb_af
    rows["Reliable_depth_model"] = ztnb_model
    return rows


def add_active_gate_columns(rows: pd.DataFrame) -> pd.DataFrame:
    rows = rows.copy()
    rows["ANI_from_Reliable_ztp_af"] = [
        1.0 + math.log(max(float(af), 1e-300)) / EFFECTIVE_CTX_LENGTH
        for af in rows["Reliable_ztp_af"]
    ]
    rows["ANI_AF_delta"] = rows["ANI_naive_calc"] - rows["ANI_from_Reliable_ztp_af"]
    rows["Reliable_depth_vmr"] = [
        (float(var) / float(mean)) if float(mean) > 0.0 else 0.0
        for mean, var in zip(rows["Reliable_Ref_hit_mean_depth"], rows["Reliable_Ref_hit_depth_variance"])
    ]
    rows["active_base_pass"] = (
        (rows["XnY_ctx"] >= ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= ACTIVE_RELIABLE_ZTP_AF_FLOOR)
        & rows["gtdb_species"].astype(bool)
    )
    rows["active_delta_trigger"] = (
        (rows["Reliable_Ref_hit_mean_depth"] > ACTIVE_MEAN_DEPTH_MIN)
        & (rows["Reliable_depth_vmr"] > ACTIVE_VMR_MIN)
    )
    rows["active_delta_pass"] = (
        (~rows["active_delta_trigger"]) | (rows["ANI_AF_delta"] < ACTIVE_DELTA_MAX)
    )
    rows["active_gate_pass"] = rows["active_base_pass"] & rows["active_delta_pass"]
    return rows


def append_details(
    detail_rows,
    method: str,
    kind: str,
    species_set: set[str],
    selected: pd.DataFrame,
    abundance: dict[str, float],
    lookup_rows: pd.DataFrame | None = None,
) -> None:
    for species in sorted(species_set):
        rec = {
            "method": method,
            "kind": kind,
            "gtdb_species": species,
            "gold_abundance": abundance.get(species, ""),
            "accession": "",
            "Ref": "",
            "ANI": "",
            "ANI_naive_calc": "",
            "XnY_ctx": "",
            "Real_min_align_fraction": "",
            "Ref_zip_af": "",
            "Reliable_Ref_zip_af": "",
            "Reliable_ztp_af": "",
            "Reliable_ztnb_af": "",
            "Reliable_Ref_breadth": "",
            "Reliable_Ref_hit_ctx": "",
            "Reliable_Ref_hit_mean_depth": "",
            "Reliable_Ref_hit_depth_variance": "",
            "Reliable_depth_vmr": "",
            "ANI_from_Reliable_ztp_af": "",
            "ANI_AF_delta": "",
            "active_base_pass": "",
            "active_delta_trigger": "",
            "active_delta_pass": "",
            "active_gate_pass": "",
            "Normalized_abundance_depth": "",
        }
        source_rows = selected if lookup_rows is None else lookup_rows
        examples = source_rows.loc[source_rows["gtdb_species"] == species].copy()
        if not examples.empty:
            examples["ANI_naive_calc"] = numeric(examples, "ANI_naive_calc")
            examples["XnY_ctx"] = numeric(examples, "XnY_ctx")
            best = examples.sort_values(["ANI_naive_calc", "XnY_ctx"], ascending=False).iloc[0]
            for field in [
                "accession",
                "Ref",
                "ANI",
                "ANI_naive_calc",
                "XnY_ctx",
                "Real_min_align_fraction",
                "Ref_zip_af",
                "Reliable_Ref_zip_af",
                "Reliable_ztp_af",
                "Reliable_ztnb_af",
                "Reliable_Ref_breadth",
                "Reliable_Ref_hit_ctx",
                "Reliable_Ref_hit_mean_depth",
                "Reliable_Ref_hit_depth_variance",
                "Reliable_depth_vmr",
                "ANI_from_Reliable_ztp_af",
                "ANI_AF_delta",
                "active_base_pass",
                "active_delta_trigger",
                "active_delta_pass",
                "active_gate_pass",
                "Normalized_abundance_depth",
            ]:
                rec[field] = best.get(field, "")
        detail_rows.append(rec)


def abundance_metrics(selected: pd.DataFrame, gold_profile: pd.DataFrame) -> list[dict[str, object]]:
    selected = selected.copy()
    selected["Normalized_abundance_depth"] = numeric(selected, "Normalized_abundance_depth")
    pred = (
        selected.loc[selected["gtdb_species"].astype(bool)]
        .groupby("gtdb_species")["Normalized_abundance_depth"]
        .max()
        .to_dict()
    )
    gold = dict(
        zip(gold_profile["gtdb_species"].astype(str), gold_profile["relative_abundance"].astype(float))
    )
    rows = []
    for renorm in [False, True]:
        values = pred.copy()
        if renorm:
            total = sum(values.values())
            if total > 0.0:
                values = {k: v / total for k, v in values.items()}
        y_true = []
        y_pred = []
        for species, truth_abund in gold.items():
            y_true.append(float(truth_abund))
            y_pred.append(float(values.get(species, 0.0)))
        true_s = pd.Series(y_true)
        pred_s = pd.Series(y_pred)
        mae = (pred_s - true_s).abs().mean()
        l1 = (pred_s - true_s).abs().sum()
        rows.append(
            {
                "view": "gold_zero_missing",
                "renorm": renorm,
                "n": len(y_true),
                "pred_sum": sum(y_pred),
                "truth_sum": sum(y_true),
                "pearson": pred_s.corr(true_s, method="pearson"),
                "spearman": pred_s.corr(true_s, method="spearman"),
                "mae_pct": mae * 100.0,
                "l1_pct": l1 * 100.0,
            }
        )
    return rows


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    profile = pd.read_csv(EXP_DIR / "mouse0_gtdb_species_profile.tsv", sep="\t")
    gold = set(profile["gtdb_species"].astype(str))
    abundance = dict(zip(profile["gtdb_species"].astype(str), profile["relative_abundance"].astype(float)))

    score_rows = []
    detail_rows = []

    zip_rows = load_rows(ZIP_TSV, by_accession, by_core)
    zip_mask = truth.direct_mask(zip_rows, "ANI", 0.94, False)
    zip_selected = zip_rows.loc[zip_mask & zip_rows["gtdb_species"].astype(bool)].copy()
    tp, fp, fn, precision, recall, f1 = score_sets(zip_selected["gtdb_species"], gold)
    score_rows.append(
        {
            "method": "pairctx_af005_split_zip_direct",
            "gate": "XnY>=10; ANI>=0.94; Real_min_af>=0.05",
            "af_floor": "",
            "pred_taxa": len(tp | fp),
            "TP": len(tp),
            "FP": len(fp),
            "FN": len(fn),
            "precision": precision,
            "recall": recall,
            "F1": f1,
            "selected_rows": len(zip_selected),
        }
    )
    append_details(detail_rows, "pairctx_af005_split_zip_direct", "FP", fp, zip_selected, abundance)
    append_details(detail_rows, "pairctx_af005_split_zip_direct", "FN", fn, zip_selected, abundance)

    nb_rows = add_reliable_adjusted_af(load_rows(PRODUCT_NB_TSV, by_accession, by_core))
    base = (
        (nb_rows["XnY_ctx"] >= 10.0)
        & (nb_rows["ANI_naive_calc"] > ANI_THRESHOLD)
        & nb_rows["gtdb_species"].astype(bool)
    )
    variants = {
        "pairctx_af005_product_nb_raw_ref_zip": "Ref_zip_af",
        "pairctx_af005_product_nb_reliable_zip": "Reliable_Ref_zip_af",
        "pairctx_af005_product_nb_reliable_ztp": "Reliable_ztp_af",
    }
    best_selected = None
    for method, col in variants.items():
        selected = nb_rows.loc[base & (nb_rows[col] >= AF_FLOOR)].copy()
        tp, fp, fn, precision, recall, f1 = score_sets(selected["gtdb_species"], gold)
        score_rows.append(
            {
                "method": method,
                "gate": f"XnY>=10; ANI_naive>{ANI_THRESHOLD}; {col}>={AF_FLOOR:.12g}",
                "af_floor": AF_FLOOR,
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
        append_details(detail_rows, method, "FP", fp, selected, abundance)
        append_details(detail_rows, method, "FN", fn, selected, abundance)
        if method == "pairctx_af005_product_nb_reliable_ztp":
            best_selected = selected

    active_inputs = {
        "old_s2000_globalctx_product0_active_recheck": OLD_PRODUCT_TOPFRAC_MEDIAN_TSV,
        "pairctx_af005_product0_active": PRODUCT_TOPFRAC_MEDIAN_TSV,
    }
    active_gate = (
        f"XnY>={ACTIVE_CTX_MIN:g}; ANI_naive>{ANI_THRESHOLD}; "
        f"Reliable_ztp_af>={ACTIVE_RELIABLE_ZTP_AF_FLOOR:g}; "
        f"if hit_mean>{ACTIVE_MEAN_DEPTH_MIN:g} and var/mean>{ACTIVE_VMR_MIN:g} "
        f"then ANI_delta<{ACTIVE_DELTA_MAX:g}; k_eff={EFFECTIVE_CTX_LENGTH:g}"
    )
    for method, path in active_inputs.items():
        active_rows = add_active_gate_columns(
            add_reliable_adjusted_af(load_rows(path, by_accession, by_core))
        )
        selected = active_rows.loc[active_rows["active_gate_pass"]].copy()
        tp, fp, fn, precision, recall, f1 = score_sets(selected["gtdb_species"], gold)
        score_rows.append(
            {
                "method": method,
                "gate": active_gate,
                "af_floor": ACTIVE_RELIABLE_ZTP_AF_FLOOR,
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
        append_details(detail_rows, method, "FP", fp, selected, abundance, active_rows)
        append_details(detail_rows, method, "FN", fn, selected, abundance, active_rows)
        if method == "pairctx_af005_product0_active":
            best_selected = selected

    pd.DataFrame(score_rows).to_csv(SCORE_OUT, sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(DETAIL_OUT, sep="\t", index=False)
    if best_selected is not None:
        pd.DataFrame(abundance_metrics(best_selected, profile)).to_csv(
            ABUND_OUT, sep="\t", index=False
        )

    print(pd.DataFrame(score_rows).to_csv(sep="\t", index=False), end="")
    print(f"wrote {SCORE_OUT}")
    print(f"wrote {DETAIL_OUT}")
    print(f"wrote {ABUND_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
