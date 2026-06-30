#!/usr/bin/env python3
"""Score wrapper-emitted adaptive call-filter profiles.

This validates the off-by-default `--adaptive-call-filter-switch
lopo-min-xny25` implementation in `scripts/minco_profile_calibrated.py`.
Profiles are regenerated from cached unique/split tables by
`run_adaptive_call_filter_wrapper_validation.sh`, then compared to the offline
adaptive switch audit. The expected method is `min_xny_ge_25` only when the
wrapper's sample-level switch fires; otherwise it is `current_calls`.
"""

from __future__ import annotations

import math
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
VALIDATION_DIR = Path("/tmp/minco_adaptive_call_filter_validation_20260629")
METHOD = "wrapper_adaptive_call_filter_lopo_min_xny25"
CURRENT_METHOD = "current_calls"
CANDIDATE_METHOD = "min_xny_ge_25"


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


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


def profile_metadata(path: Path) -> dict[str, object]:
    cols = {
        "calibrated_call",
        "adaptive_call_filter_switch",
        "adaptive_call_filter_applied",
        "adaptive_call_filter_removed_n",
        "adaptive_call_filter_input_call_n",
        "adaptive_call_filter_output_call_n",
        "adaptive_call_filter_rule",
        "adaptive_call_filter_max_xny_median",
        "adaptive_call_filter_max_xny_median_threshold",
        "adaptive_call_filter_min_xny_threshold",
    }
    raw = pd.read_csv(path, sep="\t", usecols=lambda col: col in cols, low_memory=False)
    called = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    first = raw.iloc[0] if not raw.empty else pd.Series(dtype=object)
    applied = bool_value(first.get("adaptive_call_filter_applied", False))
    expected_method = CANDIDATE_METHOD if applied else CURRENT_METHOD
    return {
        "profile": str(path),
        "adaptive_call_filter_switch": str(first.get("adaptive_call_filter_switch", "")),
        "adaptive_call_filter_applied": applied,
        "adaptive_call_filter_removed_n": int(finite(first.get("adaptive_call_filter_removed_n", 0))),
        "adaptive_call_filter_input_call_n": int(
            finite(first.get("adaptive_call_filter_input_call_n", 0))
        ),
        "adaptive_call_filter_output_call_n": int(
            finite(first.get("adaptive_call_filter_output_call_n", called.sum()))
        ),
        "observed_called_n": int(called.sum()),
        "adaptive_call_filter_max_xny_median": finite(
            first.get("adaptive_call_filter_max_xny_median", 0.0)
        ),
        "adaptive_call_filter_max_xny_median_threshold": finite(
            first.get("adaptive_call_filter_max_xny_median_threshold", 0.0)
        ),
        "adaptive_call_filter_min_xny_threshold": finite(
            first.get("adaptive_call_filter_min_xny_threshold", 0.0)
        ),
        "expected_offline_method": expected_method,
    }


def toy_prediction(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(
        path,
        sep="\t",
        usecols=lambda col: col
        in {
            "species_name",
            "calibrated_call",
            "calibrated_abundance",
            "adaptive_call_filter_switch",
            "adaptive_call_filter_applied",
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
        "adaptive_call_filter_switch": str(
            raw["adaptive_call_filter_switch"].dropna().astype(str).iloc[0]
        )
        if "adaptive_call_filter_switch" in raw and raw["adaptive_call_filter_switch"].notna().any()
        else "",
        "adaptive_call_filter_applied": bool_value(
            raw["adaptive_call_filter_applied"].dropna().astype(str).iloc[0]
        )
        if "adaptive_call_filter_applied" in raw and raw["adaptive_call_filter_applied"].notna().any()
        else False,
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


def official_l1_col(panel: str) -> str:
    return "L1_truth_only_pp" if panel == "cami2_toy_mouse_gut" else "L1_union_pp"


def official_pearson_col(panel: str) -> str:
    return "Pearson_truth_only" if panel == "cami2_toy_mouse_gut" else "Pearson_union"


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, sub in scores.groupby("panel", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        l1_col = official_l1_col(panel)
        pearson_col = official_pearson_col(panel)
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
                "official_L1_pp": sub[l1_col].mean(),
                "official_Pearson": sub[pearson_col].mean(),
                "switched_samples": int(sub["adaptive_call_filter_applied"].astype(bool).sum()),
                "removed_rows": int(sub["adaptive_call_filter_removed_n"].sum()),
            }
        )
    return pd.DataFrame(rows)


def path_map() -> dict[tuple[str, int], Path]:
    return {
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


def score_profiles(paths: dict[tuple[str, int], Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    toy_mod = decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    hmp_taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
    hmp_gastro_taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)
    _wgs_to_species, cami3_taxid_to_species, cami3_name_to_species, _diag = cami3.build_transfer_maps()
    cami3_taxmap = cami3.parse_species_taxmap(cami3.TAXMAP)

    score_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    for (panel, sample), path in paths.items():
        metadata = profile_metadata(path)
        metadata.update({"panel": panel, "sample": str(sample)})
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
        extra.update(metadata)
        row = taxid_score.score_prediction(sample, METHOD, pred, truth_df, extra)
        row["panel"] = panel
        row.update(metadata)
        score_rows.append(row)
        metadata["pred_species"] = int(len(pred))
        metadata["pred_abundance_sum"] = (
            float(pd.to_numeric(pred["pred_abundance"], errors="coerce").sum())
            if not pred.empty
            else 0.0
        )
        metadata_rows.append(metadata)
    return pd.DataFrame(score_rows), pd.DataFrame(metadata_rows)


def compare_to_offline(scores: pd.DataFrame) -> pd.DataFrame:
    offline = pd.read_csv(RESULTS / "adaptive_call_filter_switch_scores.tsv", sep="\t")
    rows: list[dict[str, object]] = []
    for row in scores.itertuples(index=False):
        panel = str(getattr(row, "panel"))
        sample = str(getattr(row, "sample"))
        expected_method = str(getattr(row, "expected_offline_method"))
        match = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(str).eq(sample)
            & offline["method"].astype(str).eq(expected_method)
        ]
        current = offline.loc[
            offline["panel"].astype(str).eq(panel)
            & offline["sample"].astype(str).eq(sample)
            & offline["method"].astype(str).eq(CURRENT_METHOD)
        ]
        if match.empty or current.empty:
            continue
        expected = match.iloc[0]
        current_row = current.iloc[0]
        l1_col = official_l1_col(panel)
        pearson_col = official_pearson_col(panel)
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "expected_offline_method": expected_method,
                "adaptive_call_filter_applied": bool(getattr(row, "adaptive_call_filter_applied")),
                "removed_rows": int(getattr(row, "adaptive_call_filter_removed_n")),
                "wrapper_TP": int(getattr(row, "TP")),
                "wrapper_FP": int(getattr(row, "FP")),
                "wrapper_FN": int(getattr(row, "FN")),
                "TP_delta_vs_expected": int(getattr(row, "TP")) - int(expected["TP"]),
                "FP_delta_vs_expected": int(getattr(row, "FP")) - int(expected["FP"]),
                "FN_delta_vs_expected": int(getattr(row, "FN")) - int(expected["FN"]),
                "F1_delta_vs_expected": finite(getattr(row, "F1")) - finite(expected["F1"]),
                "official_L1_delta_vs_expected_pp": finite(getattr(row, l1_col)) - finite(expected[l1_col]),
                "official_Pearson_delta_vs_expected": finite(getattr(row, pearson_col))
                - finite(expected[pearson_col]),
                "F1_delta_vs_current": finite(getattr(row, "F1")) - finite(current_row["F1"]),
                "official_L1_delta_vs_current_pp": finite(getattr(row, l1_col)) - finite(current_row[l1_col]),
                "official_Pearson_delta_vs_current": finite(getattr(row, pearson_col))
                - finite(current_row[pearson_col]),
                "expected_F1_delta_vs_current": finite(expected["F1"]) - finite(current_row["F1"]),
                "expected_official_L1_delta_vs_current_pp": finite(expected[l1_col])
                - finite(current_row[l1_col]),
                "expected_official_Pearson_delta_vs_current": finite(expected[pearson_col])
                - finite(current_row[pearson_col]),
            }
        )
    return pd.DataFrame(rows)


def audit(validation: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    if validation.empty:
        return pd.DataFrame(
            [
                {
                    "metric": "validation_rows",
                    "value": 0,
                    "evidence": "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
                    "decision": "failed_missing_expected_rows",
                }
            ]
        )
    max_count_delta = int(
        validation[
            ["TP_delta_vs_expected", "FP_delta_vs_expected", "FN_delta_vs_expected"]
        ]
        .abs()
        .max()
        .max()
    )
    max_f1_delta = float(validation["F1_delta_vs_expected"].abs().max())
    max_l1_delta = float(validation["official_L1_delta_vs_expected_pp"].abs().max())
    max_pearson_delta = float(validation["official_Pearson_delta_vs_expected"].abs().max())
    call_matches = (
        max_count_delta == 0
        and max_f1_delta < 1e-9
    )
    rows = [
        {
            "metric": "wrapper_profiles",
            "value": int(scores.shape[0]),
            "evidence": str(VALIDATION_DIR),
            "decision": "cached_table_wrapper_replay",
        },
        {
            "metric": "switched_profiles",
            "value": int(scores["adaptive_call_filter_applied"].astype(bool).sum()),
            "evidence": "adaptive_call_filter_switch=lopo-min-xny25",
            "decision": "candidate_switch_behavior",
        },
        {
            "metric": "removed_rows",
            "value": int(scores["adaptive_call_filter_removed_n"].sum()),
            "evidence": "wrapper profile metadata",
            "decision": "candidate_switch_behavior",
        },
        {
            "metric": "max_abs_call_delta_vs_offline_expected",
            "value": f"counts={max_count_delta};F1={max_f1_delta:.3g}",
            "evidence": "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
            "decision": "wrapper_call_mask_matches_expected" if call_matches else "wrapper_call_mismatch_review",
        },
        {
            "metric": "max_abs_abundance_delta_vs_offline_expected",
            "value": f"L1_pp={max_l1_delta:.3g};Pearson={max_pearson_delta:.3g}",
            "evidence": "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
            "decision": "informational_offline_abundance_not_raw_wrapper_contract",
        },
        {
            "metric": "mean_expected_delta_vs_current",
            "value": (
                f"F1={validation['expected_F1_delta_vs_current'].mean():.6f};"
                f"L1_pp={validation['expected_official_L1_delta_vs_current_pp'].mean():.6f};"
                f"Pearson={validation['expected_official_Pearson_delta_vs_current'].mean():.6f}"
            ),
            "evidence": "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
            "decision": "offline_candidate_effect_on_validation_subset",
        },
        {
            "metric": "mean_wrapper_score_delta_vs_current",
            "value": (
                f"F1={validation['F1_delta_vs_current'].mean():.6f};"
                f"L1_pp={validation['official_L1_delta_vs_current_pp'].mean():.6f};"
                f"Pearson={validation['official_Pearson_delta_vs_current'].mean():.6f}"
            ),
            "evidence": "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
            "decision": "includes_cached_table_abundance_replay_differences",
        },
        {
            "metric": "promotion_decision",
            "value": "candidate_requires_independent_holdout_not_default",
            "evidence": "wrapper replay validates implementation but not generalization",
            "decision": "current_call_gate_remains_default",
        },
    ]
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    paths = path_map()
    missing = [path for path in paths.values() if not path.exists()]
    if missing:
        raise SystemExit("missing wrapper validation profiles:\n" + "\n".join(map(str, missing)))

    scores, metadata = score_profiles(paths)
    panel_summary = summarize(scores)
    validation = compare_to_offline(scores)
    audit_df = audit(validation, scores)

    scores.to_csv(RESULTS / "adaptive_call_filter_wrapper_validation_scores.tsv", sep="\t", index=False)
    panel_summary.to_csv(
        RESULTS / "adaptive_call_filter_wrapper_validation_panel_summary.tsv",
        sep="\t",
        index=False,
    )
    metadata.to_csv(
        RESULTS / "adaptive_call_filter_wrapper_validation_metadata.tsv",
        sep="\t",
        index=False,
    )
    validation.to_csv(
        RESULTS / "adaptive_call_filter_wrapper_validation_vs_offline.tsv",
        sep="\t",
        index=False,
    )
    audit_df.to_csv(
        RESULTS / "adaptive_call_filter_wrapper_validation_audit.tsv",
        sep="\t",
        index=False,
    )

    print(panel_summary.to_string(index=False))
    print("\nVALIDATION")
    print(validation.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
