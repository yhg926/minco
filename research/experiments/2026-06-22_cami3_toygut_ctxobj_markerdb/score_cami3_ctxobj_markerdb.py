#!/usr/bin/env python3
"""Compare ctx-only and ctx+obj MinCO markerdbs against Sylph on CAMI3 ToyGut."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
sys.path.insert(0, str(BASE_EXP))

import score_cami3_toygut as base  # noqa: E402


RUN_DIR = Path("/tmp/cami3_toygut_ctxobj_markerdb_20260622")
OLD_RUN_DIR = Path("/tmp/cami3_toygut_current_minco_vs_sylph_20260622")
OUT_DIR = Path(__file__).resolve().parent

CTXOBJ_SAMPLES = {
    0: RUN_DIR / "minco_s0_ctxobj_product0.tsv",
    1: RUN_DIR / "minco_s1_ctxobj_product0.tsv",
    2: RUN_DIR / "minco_s2_ctxobj_product0.tsv",
}

OLD_CTX_SAMPLES = {
    0: OLD_RUN_DIR / "minco_s0_current_product0.tsv",
    1: OLD_RUN_DIR / "minco_s1_current_product0.tsv",
    2: OLD_RUN_DIR / "minco_s2_current_product0.tsv",
}


def minco_gate(rows: pd.DataFrame, ztp_floor: float, breadth_floor: float | None = None) -> pd.Series:
    gate = (
        (rows["XnY_ctx"] >= base.ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > base.ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= ztp_floor)
        & rows["active_delta_pass"]
        & rows["ncbi_species_taxid"].astype(bool)
    )
    if breadth_floor is not None:
        gate &= rows["Reliable_Ref_breadth"] >= breadth_floor
    return gate


def select_minco(rows: pd.DataFrame, ztp_floor: float, breadth_floor: float | None = None) -> pd.DataFrame:
    selected = rows.loc[minco_gate(rows, ztp_floor, breadth_floor)].copy()
    selected["pred_abundance_raw"] = base.numeric(selected, "Normalized_abundance_depth")
    selected["pred_ani"] = base.numeric(selected, "ANI_naive_calc")
    selected["support"] = base.numeric(selected, "XnY_ctx")
    return selected


def score_method(sample_id: int, method: str, selected: pd.DataFrame, truth_rows: pd.DataFrame):
    collapsed = base.collapse_predictions(selected)
    gold = set(truth_rows["ncbi_species_taxid"].astype(str))
    pred = set(collapsed["ncbi_species_taxid"].astype(str))
    tp, fp, fn, precision, recall, f1 = base.score_sets(pred, gold)
    score = {
        "sample_id": sample_id,
        "method": method,
        "truth_bacterial_species": len(gold),
        "truth_bacterial_abundance_pct_all": truth_rows["truth_abundance_pct_all"].sum(),
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

    for sample_id, paths in base.SAMPLES.items():
        truth_rows = base.load_truth(paths["truth"])

        ctx_only = base.load_minco(OLD_CTX_SAMPLES[sample_id], by_accession, by_core)
        ctxobj = base.load_minco(CTXOBJ_SAMPLES[sample_id], by_accession, by_core)
        sylph = base.load_sylph(paths["sylph"], by_accession, by_core)

        methods = [
            ("minco_ctx_only_current", select_minco(ctx_only, 0.40)),
            ("minco_ctxobj_active", select_minco(ctxobj, 0.40)),
            ("minco_ctxobj_relaxed", select_minco(ctxobj, 0.35, 0.03)),
            ("sylph", base.select_predictions("sylph", sylph)),
        ]
        for method, selected in methods:
            score, abundance, details, selected_out = score_method(
                sample_id, method, selected, truth_rows
            )
            score_rows.append(score)
            abundance_rows.extend(abundance)
            detail_rows.extend(details)
            selected_rows.append(selected_out)

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    details_df = pd.DataFrame(detail_rows)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()

    score_path = RUN_DIR / "cami3_toygut_ctxobj_markerdb_scores.tsv"
    abundance_path = RUN_DIR / "cami3_toygut_ctxobj_markerdb_abundance.tsv"
    details_path = RUN_DIR / "cami3_toygut_ctxobj_markerdb_details.tsv"
    selected_path = RUN_DIR / "cami3_toygut_ctxobj_markerdb_selected_species.tsv"
    summary_path = OUT_DIR / "summary.tsv"
    abundance_summary_path = OUT_DIR / "abundance.tsv"

    score_df.to_csv(score_path, sep="\t", index=False)
    abundance_df.to_csv(abundance_path, sep="\t", index=False)
    details_df.to_csv(details_path, sep="\t", index=False)
    selected_df.to_csv(selected_path, sep="\t", index=False)
    score_df.to_csv(summary_path, sep="\t", index=False)
    abundance_df.to_csv(abundance_summary_path, sep="\t", index=False)

    print(score_df.to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {abundance_path}")
    print(f"wrote {details_path}")
    print(f"wrote {selected_path}")
    print(f"wrote {summary_path}")
    print(f"wrote {abundance_summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
