#!/usr/bin/env python3
"""Score hybrid full-sketch fallback on CAMI3 ToyGut samples."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
CAMI3_BASE_EXP = ROOT / "research/experiments/2026-06-22_cami3_toygut_current_minco_vs_sylph"
CAMI3_RESCUE_EXP = ROOT / "research/experiments/2026-06-23_cami3_toygut_abundance_rescue"

sys.path.insert(0, str(CAMI3_BASE_EXP))
sys.path.insert(0, str(CAMI3_RESCUE_EXP))

import score_cami3_abundance_rescue as cami3_rescue  # noqa: E402
import score_cami3_toygut as base  # noqa: E402


RUN_DIR = Path("/tmp/cami3_toygut_hybrid_full_20260623")
CTX_MARKER_RUN_DIR = Path("/tmp/cami3_toygut_abundance_rescue_20260623")
CTX_MARKER_PSM_PATH = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv")

MARKER_SIZE_CUTOFF = 200
FALLBACK_ABUNDANCE_SCALE = 0.4
HYBRID_STRATEGIES = [
    ("hybrid_baseline_plus_full_xny500_scale04", 500),
    ("hybrid_baseline_plus_full_xny700_scale04", 700),
    ("hybrid_baseline_plus_full_xny1000_scale04", 1000),
]


def marker_path(sample: int) -> Path:
    return CTX_MARKER_RUN_DIR / f"minco_s{sample}_ctxmarker_current_binary.tsv"


def full_path(sample: int) -> Path:
    return RUN_DIR / f"minco_s{sample}_full_product0.tsv"


def load_psmp_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open() as fh:
        for raw in fh:
            if not raw.strip():
                continue
            size_s, ref = raw.rstrip("\n").split("\t")[:2]
            acc = base.truth.extract_accession(ref)
            if not acc:
                continue
            try:
                sizes[acc] = int(float(size_s))
            except ValueError:
                continue
    return sizes


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def add_hybrid_columns(
    rows: pd.DataFrame,
    source_mode: str,
    ctx_marker_sizes: dict[str, int],
) -> pd.DataFrame:
    rows = cami3_rescue.add_effective_depth(rows)
    rows = rows.copy()
    rows["source_mode"] = source_mode
    rows["ctx_marker_size"] = rows["accession"].map(ctx_marker_sizes).fillna(0).astype(int)
    rows["species_genus"] = rows["gtdb_species"].astype(str).map(species_genus)
    rows["pred_depth"] = base.numeric(rows, "effective_depth")
    rows.loc[rows["source_mode"] == "full_fallback", "pred_depth"] *= FALLBACK_ABUNDANCE_SCALE
    return rows


def select_robust_rescue_fast(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    rows = rows.loc[rows["ncbi_species_taxid"].astype(bool)].copy()
    active = rows.loc[rows["active_gate_pass"]].copy()
    marker_taxids = set(
        active.loc[active["source_mode"] == "ctx_marker", "ncbi_species_taxid"].astype(str)
    )
    if marker_taxids:
        active = active.loc[
            ~(
                (active["source_mode"] == "full_fallback")
                & active["ncbi_species_taxid"].astype(str).isin(marker_taxids)
            )
        ].copy()
    active["pred_abundance_raw"] = base.numeric(active, "pred_depth")
    active["pred_ani"] = base.numeric(active, "ANI_naive_calc")
    active["support"] = base.numeric(active, "XnY_ctx")
    active["selection_status"] = "active"
    active_genus_max = active.loc[active["pred_depth"] > 0.0].groupby("species_genus")["pred_depth"].max()

    candidates = rows.loc[
        ~rows["active_gate_pass"]
        & (base.numeric(rows, "ANI_naive_calc") >= cami3_rescue.RESCUE_ANI_MIN)
        & (base.numeric(rows, "XnY_ctx") >= cami3_rescue.RESCUE_XNY_MIN)
        & (base.numeric(rows, "pred_depth") >= cami3_rescue.RESCUE_EFFECTIVE_MIN)
        & (base.numeric(rows, "Reliable_Ref_zip_af") <= cami3_rescue.RESCUE_ZIP_AF_MAX)
    ].copy()

    rescued_rows = []
    rescue_details: list[dict[str, object]] = []
    if not candidates.empty and not active_genus_max.empty:
        candidates["active_genus_value"] = candidates["species_genus"].map(active_genus_max).fillna(0.0)
        candidates = candidates.loc[
            (candidates["active_genus_value"] > 0.0)
            & (candidates["pred_depth"] >= cami3_rescue.RESCUE_ACTIVE_GENUS_RATIO * candidates["active_genus_value"])
        ].copy()
        if not candidates.empty:
            winners = (
                candidates.sort_values(["species_genus", "pred_depth", "XnY_ctx"])
                .groupby("species_genus", as_index=False)
                .tail(1)
            )
            for _, row in winners.iterrows():
                taxid = str(row.get("ncbi_species_taxid", ""))
                if row.get("source_mode", "") == "full_fallback" and taxid in marker_taxids:
                    continue
                rec = row.copy()
                rec["pred_abundance_raw"] = float(row["pred_depth"])
                rec["pred_ani"] = float(row["ANI_naive_calc"])
                rec["support"] = float(row["XnY_ctx"])
                rec["selection_status"] = "rescued"
                rescued_rows.append(rec)
                rescue_details.append(
                    {
                        "gtdb_species": row.get("gtdb_species", ""),
                        "ncbi_species_taxid": row.get("ncbi_species_taxid", ""),
                        "ncbi_species": row.get("ncbi_species", ""),
                        "accession": row.get("accession", ""),
                        "source_mode": row.get("source_mode", ""),
                        "ctx_marker_size": row.get("ctx_marker_size", ""),
                        "effective_depth": row.get("effective_depth", ""),
                        "pred_depth": row.get("pred_depth", ""),
                        "active_genus_max_pred_depth": row.get("active_genus_value", ""),
                        "XnY_ctx": row.get("XnY_ctx", ""),
                        "ANI_naive_calc": row.get("ANI_naive_calc", ""),
                        "Reliable_Ref_zip_af": row.get("Reliable_Ref_zip_af", ""),
                    }
                )

    selected = pd.concat([active, pd.DataFrame(rescued_rows)], ignore_index=True) if rescued_rows else active
    return selected, rescue_details


def score_selected(sample_id: int, method: str, selected: pd.DataFrame, truth_rows: pd.DataFrame):
    return cami3_rescue.score_method(sample_id, method, selected, truth_rows)


def main() -> int:
    ctx_marker_sizes = load_psmp_sizes(CTX_MARKER_PSM_PATH)
    by_accession, by_core = base.truth.load_gtdb_metadata(base.truth.GTDB_METADATA)
    score_rows = []
    abundance_rows = []
    detail_rows = []
    selected_rows = []
    rescue_rows = []

    for sample_id, paths in base.SAMPLES.items():
        truth_rows = base.load_truth(paths["truth"])

        marker = add_hybrid_columns(
            base.load_minco(marker_path(sample_id), by_accession, by_core),
            "ctx_marker",
            ctx_marker_sizes,
        )
        full = add_hybrid_columns(
            base.load_minco(full_path(sample_id), by_accession, by_core),
            "full_fallback",
            ctx_marker_sizes,
        )

        methods = {"old_ctxmarker_robust_rescue_recalc": marker}
        for method_name, fallback_xny_min in HYBRID_STRATEGIES:
            methods[method_name] = pd.concat(
                [
                    marker.copy(),
                    full.loc[
                        (full["ctx_marker_size"] < MARKER_SIZE_CUTOFF)
                        & (full["XnY_ctx"] >= fallback_xny_min)
                    ].copy(),
                ],
                ignore_index=True,
                sort=False,
            ).fillna(0)

        for method, rows in methods.items():
            selected, rescued = select_robust_rescue_fast(rows)
            for rec in rescued:
                rec["sample_id"] = sample_id
                rec["method"] = method
                rescue_rows.append(rec)
            score, abundance, details, selected_out = score_selected(sample_id, method, selected, truth_rows)
            score_rows.append(score)
            abundance_rows.extend(abundance)
            detail_rows.extend(details)
            selected_rows.append(selected_out)

        sylph = base.load_sylph(paths["sylph"], by_accession, by_core)
        selected = base.select_predictions("sylph", sylph)
        score, abundance, details, selected_out = score_selected(sample_id, "sylph", selected, truth_rows)
        score_rows.append(score)
        abundance_rows.extend(abundance)
        detail_rows.extend(details)
        selected_rows.append(selected_out)

    score_df = pd.DataFrame(score_rows)
    abundance_df = pd.DataFrame(abundance_rows)
    detail_df = pd.DataFrame(detail_rows)
    selected_df = pd.concat(selected_rows, ignore_index=True) if selected_rows else pd.DataFrame()
    rescue_df = pd.DataFrame(rescue_rows)

    score_df.to_csv(EXP_DIR / "hybrid_cami3_sample_presence.tsv", sep="\t", index=False)
    abundance_df.to_csv(EXP_DIR / "hybrid_cami3_sample_abundance.tsv", sep="\t", index=False)
    detail_df.to_csv(EXP_DIR / "hybrid_cami3_fp_fn_details.tsv", sep="\t", index=False)
    selected_df.to_csv(EXP_DIR / "hybrid_cami3_selected_species.tsv", sep="\t", index=False)
    rescue_df.to_csv(EXP_DIR / "hybrid_cami3_rescued_species.tsv", sep="\t", index=False)

    mean_presence = (
        score_df.groupby("method", as_index=False)[
            ["pred_species", "TP", "FP", "FN", "precision", "recall", "F1"]
        ].mean()
    )
    mean_abundance = (
        abundance_df.loc[abundance_df["renorm_pred"]]
        .groupby("method", as_index=False)[
            ["pred_sum_on_truth", "pred_sum_all", "pearson", "spearman", "mae_pct_points", "l1_pct_points"]
        ]
        .mean()
    )
    summary = mean_presence.merge(mean_abundance, on="method", how="outer")
    summary = summary.sort_values(["F1", "l1_pct_points"], ascending=[False, True])
    summary.to_csv(EXP_DIR / "hybrid_cami3_summary.tsv", sep="\t", index=False)
    print(summary.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
