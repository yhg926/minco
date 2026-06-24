#!/usr/bin/env python3
"""Two-stage diagnostic optimization for coden15 Toy Mouse abundance."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
OLD_ABUND_EXP = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"
MORE_EXP = ROOT / "research/experiments/2026-06-22_cami2_toymouse_more_gtdb"
OLD_TRUTH_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"

sys.path.insert(0, str(OLD_ABUND_EXP))
sys.path.insert(0, str(MORE_EXP))
sys.path.insert(0, str(OLD_TRUTH_EXP))
sys.path.insert(0, str(EXP_DIR))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import evaluate_abundance_estimators as ev  # noqa: E402
import score_coden15_toymouse_three_samples as c15  # noqa: E402
import score_intragenus_winner_rescue as rescue  # noqa: E402
import score_more_toymouse_gtdb as more  # noqa: E402


SAMPLES = [0, 1, 2]
EPS = 1e-12


def safe_num(value: object) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return 0.0
    return x if math.isfinite(x) else 0.0


def species_genus(species: str) -> str:
    name = species[3:] if species.startswith("s__") else species
    return name.split()[0] if name else ""


def load_cached_rows(sample: int, by_accession, by_core) -> pd.DataFrame:
    cache = EXP_DIR / f"coden15_preprocessed_sample{sample}.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    rows = more.load_minco_active_rows(c15.coden15_output_path(sample), by_accession, by_core)
    keep = [
        "gtdb_species",
        "Ref",
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_ztp_af",
        "ANI_AF_delta",
        "Reliable_depth_vmr",
        "Reliable_Ref_breadth",
        "Reliable_Ref_mean_depth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_median_depth",
        "Reliable_Ref_hit_depth_variance",
        "Reliable_Ref_zip_af",
    ]
    out = rows.loc[rows["gtdb_species"].astype(bool), keep].copy()
    out["species_genus"] = out["gtdb_species"].astype(str).map(species_genus)
    for col in keep:
        if col not in {"gtdb_species", "Ref"}:
            out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    out.to_csv(cache, sep="\t", index=False)
    return out


def load_truth(sample: int) -> dict[str, float]:
    return ev.load_truth_profile(sample)


def robust_value(rows: pd.DataFrame, exponent: float, median_cutoff: float) -> pd.Series:
    median = pd.to_numeric(rows["Reliable_Ref_hit_median_depth"], errors="coerce").fillna(0.0)
    mean = pd.to_numeric(rows["Reliable_Ref_mean_depth"], errors="coerce").fillna(0.0)
    af = pd.to_numeric(rows["Reliable_Ref_zip_af"], errors="coerce").fillna(0.0).clip(lower=EPS)
    corrected = mean / (af**exponent)
    return corrected.where(median < median_cutoff, median).where(mean > 0.0, 0.0)


def gate_mask(
    rows: pd.DataFrame,
    xny_min: float,
    af_floor: float,
    vmr_min: float,
    delta_max: float,
    hit_mean_min: float,
    ani_min: float,
) -> pd.Series:
    delta_trigger = (
        (rows["Reliable_Ref_hit_mean_depth"] > hit_mean_min)
        & (rows["Reliable_depth_vmr"] > vmr_min)
    )
    delta_pass = ~delta_trigger | (rows["ANI_AF_delta"] < delta_max)
    return (
        (rows["XnY_ctx"] >= xny_min)
        & (rows["ANI_naive_calc"] > ani_min)
        & (rows["Reliable_ztp_af"] >= af_floor)
        & delta_pass
    )


def values_for_config(
    rows: pd.DataFrame,
    active_mask: pd.Series,
    values: pd.Series,
    rescue_enabled: bool,
    rescue_ratio: float,
    rescue_zip_af_max: float,
    rescue_xny_min: float,
    rescue_effective_min: float,
) -> tuple[dict[str, float], list[str]]:
    work = rows.copy()
    work["value"] = values

    active = work.loc[active_mask & (work["value"] > 0.0)]
    species_values = active.groupby("gtdb_species")["value"].max().to_dict()
    active_genus_max = active.groupby("species_genus")["value"].max().to_dict()
    if not rescue_enabled or not active_genus_max:
        return species_values, []

    nonactive = work.loc[~active_mask & (work["value"] > 0.0)].copy()
    if nonactive.empty:
        return species_values, []

    nonactive["active_genus_value"] = nonactive["species_genus"].map(active_genus_max).fillna(0.0)
    rescue_mask = (
        (nonactive["active_genus_value"] > 0.0)
        & (nonactive["ANI_naive_calc"] >= rescue.RESCUE_ANI_MIN)
        & (nonactive["XnY_ctx"] >= rescue_xny_min)
        & (nonactive["value"] >= rescue_effective_min)
        & (nonactive["Reliable_Ref_zip_af"] <= rescue_zip_af_max)
        & (nonactive["value"] >= rescue_ratio * nonactive["active_genus_value"])
    )
    candidates = nonactive.loc[rescue_mask].copy()
    if candidates.empty:
        return species_values, []
    candidates = candidates.sort_values(["species_genus", "value", "XnY_ctx"], ascending=[True, False, False])
    best = candidates.drop_duplicates("species_genus", keep="first")
    rescued = []
    for _, row in best.iterrows():
        species = str(row["gtdb_species"])
        species_values[species] = max(species_values.get(species, 0.0), safe_num(row["value"]))
        rescued.append(species)
    return species_values, rescued


def score_values(values: dict[str, float], gold: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0)
    pred = {k: v / total for k, v in values.items() if v > 0.0} if total > 0.0 else {}
    keys = sorted(gold)
    y_true = [gold[k] for k in keys]
    y_pred = [pred.get(k, 0.0) for k in keys]
    pearson, spearman = ev.pearson_spearman(y_true, y_pred)
    pred_species = set(pred)
    gold_species = set(gold)
    tp = pred_species & gold_species
    fp = pred_species - gold_species
    fn = gold_species - pred_species
    return {
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "missing_truth_mass": sum(gold[s] for s in fn),
        "fp_pred_mass": sum(pred.get(s, 0.0) for s in fp),
        "pred_sum_on_truth": sum(y_pred),
        "pred_sum_all": sum(pred.values()),
        "pearson": pearson,
        "spearman": spearman,
        "mae_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) / len(keys) * 100.0,
        "l1_pct_points": sum(abs(a - b) for a, b in zip(y_true, y_pred)) * 100.0,
    }


def mean_metrics(sample_metrics: list[dict[str, float]]) -> dict[str, float]:
    df = pd.DataFrame(sample_metrics)
    return df.mean(numeric_only=True).to_dict()


def score_config(
    rows_by_sample: dict[int, pd.DataFrame],
    truth_by_sample: dict[int, dict[str, float]],
    params: dict[str, float | bool],
) -> tuple[dict[str, float], list[dict[str, object]]]:
    sample_rows = []
    for sample, rows in rows_by_sample.items():
        values = robust_value(
            rows,
            exponent=float(params["abundance_exponent"]),
            median_cutoff=float(params["median_cutoff"]),
        )
        active = gate_mask(
            rows,
            xny_min=float(params["xny_min"]),
            af_floor=float(params["af_floor"]),
            vmr_min=float(params["vmr_min"]),
            delta_max=float(params["delta_max"]),
            hit_mean_min=float(params["hit_mean_min"]),
            ani_min=float(params["ani_min"]),
        )
        pred_values, rescued = values_for_config(
            rows,
            active,
            values,
            rescue_enabled=bool(params["rescue_enabled"]),
            rescue_ratio=float(params["rescue_ratio"]),
            rescue_zip_af_max=float(params["rescue_zip_af_max"]),
            rescue_xny_min=float(params["rescue_xny_min"]),
            rescue_effective_min=float(params["rescue_effective_min"]),
        )
        rec = score_values(pred_values, truth_by_sample[sample])
        rec["sample"] = sample
        rec["rescued_species"] = ",".join(sorted(set(rescued)))
        sample_rows.append(rec)
    return mean_metrics(sample_rows), sample_rows


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    rows_by_sample = {
        sample: load_cached_rows(sample, by_accession, by_core)
        for sample in SAMPLES
    }
    truth_by_sample = {sample: load_truth(sample) for sample in SAMPLES}

    default_params = {
        "xny_min": 15.0,
        "af_floor": 0.40,
        "vmr_min": 50.0,
        "delta_max": 0.03,
        "hit_mean_min": 3.0,
        "ani_min": 0.95,
        "abundance_exponent": 1.05,
        "median_cutoff": 20.0,
        "rescue_enabled": True,
        "rescue_ratio": 3.0,
        "rescue_zip_af_max": 0.25,
        "rescue_xny_min": 100.0,
        "rescue_effective_min": 5.0,
    }

    stage1_rows = []
    for xny_min in [10.0, 15.0, 20.0, 30.0, 50.0]:
        for af_floor in [0.30, 0.32, 0.35, 0.37, 0.40]:
            for vmr_min in [5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 50.0]:
                for delta_max in [0.02, 0.03, 0.04]:
                    params = dict(default_params)
                    params.update(
                        {
                            "xny_min": xny_min,
                            "af_floor": af_floor,
                            "vmr_min": vmr_min,
                            "delta_max": delta_max,
                            "abundance_exponent": 1.0,
                        }
                    )
                    mean, _ = score_config(rows_by_sample, truth_by_sample, params)
                    stage1_rows.append({**params, **mean})
    stage1_df = pd.DataFrame(stage1_rows).sort_values("l1_pct_points")
    stage1_df.to_csv(EXP_DIR / "coden15_gate_stage1_grid.tsv", sep="\t", index=False)

    stage2_rows = []
    top_gates = stage1_df.head(20)
    for _, gate in top_gates.iterrows():
        for exponent in [0.75, 0.9, 1.0, 1.05, 1.15]:
            for median_cutoff in [10.0, 20.0, 50.0]:
                for rescue_enabled in [False, True]:
                    for rescue_ratio in [2.0, 3.0, 4.0]:
                        if not rescue_enabled and rescue_ratio != 2.0:
                            continue
                        for rescue_zip_af_max in [0.20, 0.25, 0.30, 0.40]:
                            if not rescue_enabled and rescue_zip_af_max != 0.20:
                                continue
                            for rescue_xny_min in [50.0, 100.0]:
                                if not rescue_enabled and rescue_xny_min != 50.0:
                                    continue
                                params = dict(default_params)
                                params.update(
                                    {
                                        "xny_min": safe_num(gate["xny_min"]),
                                        "af_floor": safe_num(gate["af_floor"]),
                                        "vmr_min": safe_num(gate["vmr_min"]),
                                        "delta_max": safe_num(gate["delta_max"]),
                                        "abundance_exponent": exponent,
                                        "median_cutoff": median_cutoff,
                                        "rescue_enabled": rescue_enabled,
                                        "rescue_ratio": rescue_ratio,
                                        "rescue_zip_af_max": rescue_zip_af_max,
                                        "rescue_xny_min": rescue_xny_min,
                                    }
                                )
                                mean, _ = score_config(rows_by_sample, truth_by_sample, params)
                                stage2_rows.append({**params, **mean})
    stage2_df = pd.DataFrame(stage2_rows).sort_values("l1_pct_points")
    stage2_df.to_csv(EXP_DIR / "coden15_param_stage2_grid.tsv", sep="\t", index=False)

    best = stage2_df.iloc[0].to_dict()
    mean, sample_rows = score_config(rows_by_sample, truth_by_sample, best)
    sample_df = pd.DataFrame(sample_rows)
    for key, value in best.items():
        if key not in sample_df.columns:
            sample_df[key] = value
    sample_df.to_csv(EXP_DIR / "coden15_param_best_sample_metrics.tsv", sep="\t", index=False)

    print("stage1 best")
    print(stage1_df.head(10).to_csv(sep="\t", index=False), end="")
    print("\nstage2 best")
    print(stage2_df.head(20).to_csv(sep="\t", index=False), end="")
    print("\nbest sample metrics")
    print(sample_df.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
