#!/usr/bin/env python3
"""Score wrapper-emitted genus-XnY abundance-blend profiles.

This validates the implementation in `scripts/minco_profile_calibrated.py`
against the offline fixed-call sweep. The profiles are regenerated from cached
unique/split tables by `run_abundance_blend_wrapper_validation.sh`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd

import decompose_abundance_errors as decomp
import score_cami3_gtdb_source_readmap as cami3
import score_cami3_gtdb_taxid_transfer as taxid_score
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
VALIDATION_DIR = Path("/tmp/minco_abundance_blend_validation_20260627")
METHOD = "wrapper_genus_xny_blend_a0.25"
OFFLINE_METHOD = "blend_current_genus_realloc_xny_a0.25"


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(k): max(0.0, float(v)) for k, v in values.items() if str(k)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def toy_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / f"mouse{sample}_gtdb_species_profile.tsv", sep="\t")
    return truth.rename(columns={"relative_abundance": "truth_abundance"})[
        ["gtdb_species", "truth_abundance"]
    ]


def toy_prediction(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(
        path,
        sep="\t",
        usecols=lambda col: col
        in {
            "species_name",
            "calibrated_call",
            "calibrated_abundance",
            "abundance_rule",
            "abundance_genus_xny_blend_alpha",
        },
        low_memory=False,
    )
    called = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[called].copy()
    selected["calibrated_abundance"] = pd.to_numeric(
        selected["calibrated_abundance"], errors="coerce"
    ).fillna(0.0)
    grouped = selected.groupby("species_name", as_index=False)["calibrated_abundance"].max()
    pred = normalize(dict(zip(grouped["species_name"].astype(str), grouped["calibrated_abundance"])))
    pred_df = pd.DataFrame(
        {
            "gtdb_species": list(pred),
            "pred_abundance": list(pred.values()),
            "pred_ani": 0.0,
        }
    )
    return pred_df, {
        "pred_rows_called": int(len(selected)),
        "abundance_rule": str(raw["abundance_rule"].dropna().astype(str).iloc[0])
        if "abundance_rule" in raw and raw["abundance_rule"].notna().any()
        else "",
        "abundance_genus_xny_blend_alpha": float(
            pd.to_numeric(raw.get("abundance_genus_xny_blend_alpha", 0.0), errors="coerce")
            .dropna()
            .iloc[0]
        )
        if "abundance_genus_xny_blend_alpha" in raw
        else 0.0,
    }


def hmp_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "hmp_current_refresh_r232_source_abundance_truth.tsv", sep="\t")
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def hmp_gastro_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "hmp_gastrooral_r232_source_abundance_truth.tsv", sep="\t")
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def cami3_truth(sample: int) -> pd.DataFrame:
    truth = pd.read_csv(RESULTS / "cami3_gtdb_source_readmap_truth.tsv", sep="\t")
    return truth.loc[truth["sample"].astype(int).eq(sample), ["gtdb_species", "truth_abundance"]].copy()


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, sub in scores.groupby("panel", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        official_l1_col = "L1_truth_only_pp" if panel == "cami2_toy_mouse_gut" else "L1_union_pp"
        official_pearson_col = "Pearson_truth_only" if panel == "cami2_toy_mouse_gut" else "Pearson_union"
        rows.append(
            {
                "panel": panel,
                "method": METHOD,
                "samples": ",".join(map(str, sorted(sub["sample"].unique(), key=str))),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_L1_truth_only_pp": sub["L1_truth_only_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "mean_Pearson_truth_only": sub["Pearson_truth_only"].mean(),
                "official_L1_pp": sub[official_l1_col].mean(),
                "official_Pearson": sub[official_pearson_col].mean(),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    missing = []
    paths = {
        ("cami2_toy_mouse_gut", 5): VALIDATION_DIR / "toymouse_sample5.tsv",
        ("cami2_toy_mouse_gut", 6): VALIDATION_DIR / "toymouse_sample6.tsv",
        ("cami2_toy_mouse_gut", 7): VALIDATION_DIR / "toymouse_sample7.tsv",
        ("hmp_airskin_gtdb_source_abundance", 6): VALIDATION_DIR / "hmp_airskin_sample6.tsv",
        ("hmp_airskin_gtdb_source_abundance", 11): VALIDATION_DIR / "hmp_airskin_sample11.tsv",
        ("hmp_gastrooral_gtdb_source_abundance", 0): VALIDATION_DIR / "hmp_gastrooral_sample0.tsv",
        ("hmp_gastrooral_gtdb_source_abundance", 6): VALIDATION_DIR / "hmp_gastrooral_sample6.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 0): VALIDATION_DIR / "cami3_sample0.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 1): VALIDATION_DIR / "cami3_sample1.tsv",
        ("cami3_toy_human_gut_gtdb_source_readmap", 2): VALIDATION_DIR / "cami3_sample2.tsv",
    }
    for path in paths.values():
        if not path.exists():
            missing.append(path)
    if missing:
        raise SystemExit("missing wrapper validation profiles:\n" + "\n".join(map(str, missing)))

    toy_mod = decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    hmp_taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
    hmp_gastro_taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)
    _wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _diag = cami3.build_transfer_maps()
    cami3_taxmap = cami3.parse_species_taxmap(cami3.TAXMAP)

    score_rows: list[dict[str, object]] = []
    map_rows: list[dict[str, object]] = []
    for (panel, sample), path in paths.items():
        if panel == "cami2_toy_mouse_gut":
            truth_df = toy_truth(sample)
            pred, extra = toy_prediction(path)
        elif panel == "hmp_airskin_gtdb_source_abundance":
            truth_df = hmp_truth(sample)
            hmp_paths = hmp.SAMPLES[int(sample)]
            best_ref_species, best_diag = hmp.best_raw_ref_species_by_taxid(
                {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                hmp_taxmap,
                by_accession,
                by_core,
            )
            pred, extra = hmp.load_minco_predictions(path, best_ref_species, by_accession, by_core)
            extra.update(best_diag)
        elif panel == "hmp_gastrooral_gtdb_source_abundance":
            truth_df = hmp_gastro_truth(sample)
            hmp_paths = hmp_gastro.SAMPLES[int(sample)]
            best_ref_species, best_diag = hmp_gastro.best_raw_ref_species_by_taxid(
                {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                hmp_gastro_taxmap,
                by_accession,
                by_core,
            )
            pred, extra = hmp_gastro.load_minco_predictions(path, best_ref_species, by_accession, by_core)
            extra.update(best_diag)
        elif panel == "cami3_toy_human_gut_gtdb_source_readmap":
            truth_df = cami3_truth(sample)
            best_ref_species, best_diag = cami3.best_raw_ref_species_by_taxid(
                cami3.RAW_TABLES[int(sample)],
                cami3_taxmap,
                by_accession,
                by_core,
            )
            pred, extra = cami3.load_minco_predictions(
                path,
                cami3_taxid_to_species,
                cami3_name_to_species,
                best_ref_species,
                by_accession,
                by_core,
            )
            extra.update(best_diag)
        else:
            raise AssertionError(panel)
        row = taxid_score.score_prediction(sample, METHOD, pred, truth_df, extra)
        row["panel"] = panel
        row["profile"] = str(path)
        score_rows.append(row)
        map_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "profile": str(path),
                "pred_species": int(len(pred)),
                "pred_abundance_sum": float(pd.to_numeric(pred["pred_abundance"], errors="coerce").sum())
                if not pred.empty
                else 0.0,
            }
        )

    scores = pd.DataFrame(score_rows)
    panel_summary = summarize(scores)
    offline_scores = pd.read_csv(RESULTS / "cross_panel_abundance_variant_scores.tsv", sep="\t")
    validation_rows = []
    for row in scores.itertuples(index=False):
        panel = str(getattr(row, "panel"))
        sample = str(getattr(row, "sample"))
        offline = offline_scores.loc[
            offline_scores["panel"].astype(str).eq(panel)
            & offline_scores["sample"].astype(str).eq(sample)
            & offline_scores["method"].astype(str).eq(OFFLINE_METHOD)
        ]
        if offline.empty:
            continue
        offline_row = offline.iloc[0]
        l1_col = "L1_truth_only_pp" if panel == "cami2_toy_mouse_gut" else "L1_union_pp"
        pearson_col = "Pearson_truth_only" if panel == "cami2_toy_mouse_gut" else "Pearson_union"
        validation_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "offline_method": OFFLINE_METHOD,
                "wrapper_method": METHOD,
                "TP_delta": int(getattr(row, "TP")) - int(offline_row["TP"]),
                "FP_delta": int(getattr(row, "FP")) - int(offline_row["FP"]),
                "FN_delta": int(getattr(row, "FN")) - int(offline_row["FN"]),
                "F1_delta": float(getattr(row, "F1")) - float(offline_row["F1"]),
                "official_L1_delta_pp": float(getattr(row, l1_col)) - float(offline_row[l1_col]),
                "official_Pearson_delta": float(getattr(row, pearson_col)) - float(offline_row[pearson_col]),
            }
        )

    validation = pd.DataFrame(validation_rows)
    scores.to_csv(RESULTS / "abundance_blend_wrapper_validation_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(
        RESULTS / "abundance_blend_wrapper_validation_panel_summary.tsv",
        sep="\t",
        index=False,
    )
    validation.to_csv(
        RESULTS / "abundance_blend_wrapper_validation_vs_offline.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(map_rows).to_csv(
        RESULTS / "abundance_blend_wrapper_validation_mapping.tsv",
        sep="\t",
        index=False,
    )

    print(panel_summary.to_string(index=False))
    print("\nVALIDATION")
    print(validation.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
