#!/usr/bin/env python3
"""Replay raw-candidate retention rules across cached panels.

This validates the raw-side candidate-retention direction beyond the two newest
negative routes. Rescued species receive zero abundance mass in this diagnostic,
so abundance L1/Pearson stay tied to the current default while F1 changes.

This script is intentionally offline and truth-aware for scoring only. It does
not change MinCO defaults.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Callable, Iterable, Mapping

import pandas as pd

import audit_negative_route_failure_modes as failure
import decompose_abundance_errors as decomp
import score_cami3_gtdb_source_readmap as cami3
import score_hmp_gastrooral_gtdb_source_abundance as hmp_gastro
import score_hmp_gtdb_source_abundance as hmp
import sweep_hmp_raw_candidate_rescue as hmp_raw
import sweep_raw_candidate_retention_negative_routes as neg_sweep


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
CACHE_DIR = Path("/tmp/minco_cross_panel_raw_candidate_retention_cache")

HMP_OMITTED_TRUTH = (
    RESULTS
    / "hmp_current_refresh_r232_source_abundance_sample12_0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_17_18_19_20_21_22_23_24_25_26_27_28_truth.tsv"
)


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def normalize(values: Mapping[str, float]) -> dict[str, float]:
    clean = {str(key): max(0.0, finite(value)) for key, value in values.items() if str(key)}
    total = sum(clean.values())
    if total <= 0.0:
        return clean
    return {key: value / total for key, value in clean.items()}


def f1(tp: int, fp: int, fn: int) -> float:
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) < 2:
        return float("nan")
    return float(pd.Series(left, dtype=float).corr(pd.Series(right, dtype=float), method="pearson"))


def score_sample(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    current_pred: Mapping[str, float],
    call_species: set[str],
    rescued: set[str],
    rule: str,
) -> dict[str, object]:
    truth_norm = normalize(truth)
    pred_norm = normalize(current_pred)
    truth_species = set(truth_norm)
    pred_species = set(call_species) | set(rescued)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species
    union = sorted(truth_species | pred_species)
    y_true = [truth_norm.get(species, 0.0) for species in union]
    y_pred = [pred_norm.get(species, 0.0) for species in union]
    return {
        "panel": panel,
        "sample": sample,
        "method": method,
        "truth_species": len(truth_species),
        "pred_species": len(pred_species),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": len(tp) / len(pred_species) if pred_species else 0.0,
        "recall": len(tp) / len(truth_species) if truth_species else 0.0,
        "F1": f1(len(tp), len(fp), len(fn)),
        "L1_union_pp": sum(abs(a - b) for a, b in zip(y_pred, y_true)) * 100.0,
        "Pearson_union": pearson(y_pred, y_true),
        "rescued_species": len(rescued),
        "rescued_TP": len(set(rescued) & truth_species),
        "rescued_FP": len(set(rescued) - truth_species),
        "rule": rule,
    }


def collapse_profile_by_accession(path: Path, mapper: failure.GtdbMapper, collapse: str = "sum") -> dict[str, float]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    if "calibrated_call" not in raw.columns:
        return {}
    selected = raw.loc[raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})]
    rows: list[dict[str, object]] = []
    for row in selected.itertuples(index=False):
        species = ""
        for col in ["candidate_surface_accession", "s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species, _method = mapper.species_from_accession(getattr(row, col))
                if species:
                    break
        if species:
            rows.append(
                {
                    "gtdb_species": species,
                    "abundance": finite(getattr(row, "calibrated_abundance", 0.0)),
                }
            )
    if not rows:
        return {}
    df = pd.DataFrame(rows)
    if collapse == "max":
        grouped = df.groupby("gtdb_species")["abundance"].max()
    else:
        grouped = df.groupby("gtdb_species")["abundance"].sum()
    return normalize({str(key): finite(value) for key, value in grouped.to_dict().items()})


def current_pred_for_config(config: Mapping[str, object], mapper: failure.GtdbMapper) -> dict[str, float]:
    kind = str(config["kind"])
    panel = str(config["panel"])
    sample = int(config["sample"])
    profile = Path(str(config["profile"]))
    if kind == "toy":
        return decomp.load_minco_pred(profile, str(config.get("collapse", "max")))
    if kind == "hmp_airskin":
        hmp.ensure_sample_record(sample)
        hmp_paths = hmp.SAMPLES[sample]
        taxmap = hmp.parse_species_taxmap(hmp.TAXMAP)
        best_ref_species, _diag = hmp.best_raw_ref_species_by_taxid(
            {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
            taxmap,
            mapper.by_accession,
            mapper.by_core,
        )
        pred, _extra = hmp.load_minco_predictions(
            profile,
            best_ref_species,
            mapper.by_accession,
            mapper.by_core,
        )
        return decomp.pred_from_collapsed(pred)
    if kind == "hmp_gastro":
        hmp_paths = hmp_gastro.SAMPLES[sample]
        taxmap = hmp_gastro.parse_species_taxmap(hmp_gastro.TAXMAP)
        best_ref_species, _diag = hmp_gastro.best_raw_ref_species_by_taxid(
            {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
            taxmap,
            mapper.by_accession,
            mapper.by_core,
        )
        pred, _extra = hmp_gastro.load_minco_predictions(
            profile,
            best_ref_species,
            mapper.by_accession,
            mapper.by_core,
        )
        return decomp.pred_from_collapsed(pred)
    if kind == "cami3":
        _wgs, taxid_to_species, name_to_species, _diag = cami3.build_transfer_maps()
        taxmap = cami3.parse_species_taxmap(cami3.TAXMAP)
        best_ref_species, _diag2 = cami3.best_raw_ref_species_by_taxid(
            cami3.RAW_TABLES[sample],
            taxmap,
            mapper.by_accession,
            mapper.by_core,
        )
        pred, _extra = cami3.load_minco_predictions(
            profile,
            taxid_to_species,
            name_to_species,
            best_ref_species,
            mapper.by_accession,
            mapper.by_core,
        )
        return decomp.pred_from_collapsed(pred)
    if kind in {"marine", "cami3_extension"}:
        return collapse_profile_by_accession(profile, mapper, str(config.get("collapse", "sum")))
    raise AssertionError(f"unsupported kind: {kind} for {panel}")


def truth_for_config(config: Mapping[str, object]) -> dict[str, float]:
    return decomp.load_truth(
        Path(str(config["truth"])),
        int(config["sample"]),
        str(config["truth_abundance_col"]),
    )


def raw_cache_path(panel: str, sample: int) -> Path:
    clean = panel.replace("/", "_")
    return CACHE_DIR / f"{clean}.sample{sample}.best_raw.tsv"


def raw_best_for_config(
    config: Mapping[str, object],
    mapper: failure.GtdbMapper,
    hmp_raw_mapper: hmp_raw.raw_audit.RawMapper,
) -> pd.DataFrame:
    kind = str(config["kind"])
    panel = str(config["panel"])
    sample = int(config["sample"])
    profile = Path(str(config["profile"]))
    if kind == "hmp_airskin":
        return hmp_raw.load_raw_best(
            "hmp_airskin_gtdb_source_abundance",
            sample,
            profile,
            hmp_raw_mapper,
        )
    if kind == "hmp_gastro":
        return hmp_raw.load_raw_best(
            "hmp_gastrooral_gtdb_source_abundance",
            sample,
            profile,
            hmp_raw_mapper,
        )

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = raw_cache_path(panel, sample)
    if cached.is_file():
        return pd.read_csv(cached, sep="\t", low_memory=False)
    raw = failure.best_raw_rows(config["raw_paths"], mapper)  # type: ignore[arg-type]
    raw.to_csv(cached, sep="\t", index=False)
    return raw


def configs() -> list[dict[str, object]]:
    out: list[dict[str, object]] = []

    toy_raw_root = Path("/tmp/cami2_toymouse_current_default_20260626/run")
    toy_cfg = decomp.PANELS["cami2_toy_mouse_gut"]
    for sample in toy_cfg["samples"]:
        sample = int(sample)
        out.append(
            {
                "panel": "cami2_toy_mouse_gut",
                "kind": "toy",
                "sample": sample,
                "truth": toy_cfg["truth"](sample),
                "truth_abundance_col": toy_cfg["truth_abundance_col"],
                "profile": toy_cfg["minco"](sample),
                "collapse": toy_cfg["minco_collapse"],
                "raw_paths": {
                    "unique": toy_raw_root / f"sample{sample}_default/minco.best_diff_unique.unfiltered.tsv",
                    "split": toy_raw_root / f"sample{sample}_default/minco.best_diff_split.unfiltered.tsv",
                    "split_exact": toy_raw_root
                    / f"sample{sample}_default/minco.best_diff_split.exact.unfiltered.tsv",
                },
            }
        )

    hmp_cfg = decomp.PANELS["hmp_airskin_gtdb_source_abundance"]
    for sample in hmp_cfg["samples"]:
        sample = int(sample)
        out.append(
            {
                "panel": "hmp_airskin_gtdb_source_abundance",
                "kind": "hmp_airskin",
                "sample": sample,
                "truth": hmp_cfg["truth"](sample),
                "truth_abundance_col": hmp_cfg["truth_abundance_col"],
                "profile": hmp_cfg["minco"](sample),
                "collapse": hmp_cfg["minco_collapse"],
            }
        )

    for sample in [2, 8, 12, 26, 27]:
        hmp.ensure_sample_record(sample)
        out.append(
            {
                "panel": "hmp_airskin_omitted_gtdb_source_abundance",
                "kind": "hmp_airskin",
                "sample": sample,
                "truth": HMP_OMITTED_TRUTH,
                "truth_abundance_col": "truth_abundance",
                "profile": hmp.SAMPLES[sample]["minco"],
                "collapse": "sum",
            }
        )

    gastro_cfg = decomp.PANELS["hmp_gastrooral_gtdb_source_abundance"]
    for sample in gastro_cfg["samples"]:
        sample = int(sample)
        out.append(
            {
                "panel": "hmp_gastrooral_gtdb_source_abundance",
                "kind": "hmp_gastro",
                "sample": sample,
                "truth": gastro_cfg["truth"](sample),
                "truth_abundance_col": gastro_cfg["truth_abundance_col"],
                "profile": gastro_cfg["minco"](sample),
                "collapse": gastro_cfg["minco_collapse"],
            }
        )

    cami3_cfg = decomp.PANELS["cami3_toy_human_gut_gtdb_source_readmap"]
    for sample in cami3_cfg["samples"]:
        sample = int(sample)
        out.append(
            {
                "panel": "cami3_toy_human_gut_gtdb_source_readmap",
                "kind": "cami3",
                "sample": sample,
                "truth": cami3_cfg["truth"](sample),
                "truth_abundance_col": cami3_cfg["truth_abundance_col"],
                "profile": cami3_cfg["minco"](sample),
                "collapse": cami3_cfg["minco_collapse"],
                "raw_paths": cami3.RAW_TABLES[sample],
            }
        )

    # Completed negative extension routes.
    for route in failure.route_configs():
        panel = str(route["route"])
        for sample in route["samples"]:  # type: ignore[index]
            sample = int(sample)
            out.append(
                {
                    "panel": panel,
                    "kind": "marine" if "marine" in panel else "cami3_extension",
                    "sample": sample,
                    "truth": route["truth"],
                    "truth_abundance_col": "truth_abundance",
                    "profile": route["profile"](sample),  # type: ignore[operator]
                    "collapse": "sum",
                    "raw_paths": route["raw_paths"](sample),  # type: ignore[operator]
                }
            )
    return out


def rule_species(raw: pd.DataFrame, current_calls: set[str], rule: Callable[[pd.DataFrame], pd.Series]) -> set[str]:
    if raw.empty:
        return set()
    work = raw.loc[
        raw["gtdb_species"].fillna("").astype(str).ne("")
        & ~raw["gtdb_species"].astype(str).isin(current_calls)
    ].copy()
    if work.empty:
        return set()
    mask = rule(work).fillna(False)
    return set(work.loc[mask, "gtdb_species"].astype(str))


def summarize(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_base = {
        (str(row.panel), int(row.sample)): row
        for row in scores.loc[scores["method"].eq("current_default")].itertuples(index=False)
    }
    delta_rows = []
    for row in scores.loc[~scores["method"].eq("current_default")].itertuples(index=False):
        base = sample_base[(str(row.panel), int(row.sample))]
        delta_rows.append(
            {
                "panel": str(row.panel),
                "sample": int(row.sample),
                "method": str(row.method),
                "delta_F1": finite(row.F1) - finite(base.F1),
                "delta_L1_union_pp": finite(row.L1_union_pp) - finite(base.L1_union_pp),
                "delta_Pearson_union": finite(row.Pearson_union) - finite(base.Pearson_union),
                "TP_delta": int(row.TP) - int(base.TP),
                "FP_delta": int(row.FP) - int(base.FP),
                "FN_delta": int(row.FN) - int(base.FN),
                "rescued_species": int(row.rescued_species),
                "rescued_TP": int(row.rescued_TP),
                "rescued_FP": int(row.rescued_FP),
            }
        )
    deltas = pd.DataFrame(delta_rows)

    panel_rows = []
    for (panel, method), sub in scores.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        panel_rows.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "pooled_F1": f1(tp, fp, fn),
                "mean_F1": sub["F1"].mean(),
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "rescued_TP": int(sub["rescued_TP"].sum()),
                "rescued_FP": int(sub["rescued_FP"].sum()),
            }
        )
    panels = pd.DataFrame(panel_rows)
    return deltas, panels


def overall(deltas: pd.DataFrame, panels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, sub in deltas.groupby("method", sort=True):
        panel_sub = panels.loc[panels["method"].eq(method)]
        base_panel_sub = panels.loc[panels["method"].eq("current_default")]
        panel_delta = panel_sub.set_index("panel")["pooled_F1"] - base_panel_sub.set_index("panel")[
            "pooled_F1"
        ]
        rows.append(
            {
                "method": method,
                "mean_delta_F1": sub["delta_F1"].mean(),
                "min_delta_F1": sub["delta_F1"].min(),
                "sample_worsen_n": int((sub["delta_F1"] < -1e-12).sum()),
                "sample_improve_n": int((sub["delta_F1"] > 1e-12).sum()),
                "mean_delta_L1_union_pp": sub["delta_L1_union_pp"].mean(),
                "mean_delta_Pearson_union": sub["delta_Pearson_union"].mean(),
                "total_TP_delta": int(sub["TP_delta"].sum()),
                "total_FP_delta": int(sub["FP_delta"].sum()),
                "total_FN_delta": int(sub["FN_delta"].sum()),
                "rescued_TP": int(sub["rescued_TP"].sum()),
                "rescued_FP": int(sub["rescued_FP"].sum()),
                "panel_worsen_n": int((panel_delta < -1e-12).sum()),
                "panel_improve_n": int((panel_delta > 1e-12).sum()),
                "min_panel_delta_pooled_F1": float(panel_delta.min()) if len(panel_delta) else 0.0,
                "mean_panel_delta_pooled_F1": float(panel_delta.mean()) if len(panel_delta) else 0.0,
                "strict_pass": bool(
                    (sub["delta_F1"] >= -1e-12).all()
                    and (panel_delta >= -1e-12).all()
                    and (sub["delta_F1"] > 1e-12).any()
                ),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        [
            "strict_pass",
            "min_panel_delta_pooled_F1",
            "mean_panel_delta_pooled_F1",
            "sample_worsen_n",
            "total_FP_delta",
        ],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    )


def build_audit(overall_df: pd.DataFrame, sample_n: int, panel_n: int) -> list[dict[str, object]]:
    strict = overall_df.loc[overall_df["strict_pass"].astype(bool)] if not overall_df.empty else pd.DataFrame()
    best = overall_df.iloc[0].to_dict() if not overall_df.empty else {}
    return [
        {
            "metric": "scope",
            "value": f"samples={sample_n};panels={panel_n}",
            "evidence": "current cached panels plus completed negative route extensions",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "strict_pass_rules",
            "value": int(len(strict)),
            "evidence": "requires no sample or panel F1 regression",
            "decision": "candidate_found" if not strict.empty else "no_cross_panel_safe_rule",
        },
        {
            "metric": "best_rule",
            "value": best.get("method", ""),
            "evidence": (
                f"mean_delta_F1={best.get('mean_delta_F1', '')};"
                f"min_delta_F1={best.get('min_delta_F1', '')};"
                f"sample_worsen_n={best.get('sample_worsen_n', '')};"
                f"panel_worsen_n={best.get('panel_worsen_n', '')};"
                f"TP_delta={best.get('total_TP_delta', '')};"
                f"FP_delta={best.get('total_FP_delta', '')}"
            ),
            "decision": "review_candidate" if bool(best.get("strict_pass", False)) else "diagnostic_only",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "offline zero-mass raw-candidate retention replay",
            "decision": "implement_opt_in_wrapper_then_validate_runtime_and_independent_holdout"
            if not strict.empty
            else "do_not_implement_default_from_current_rules",
        },
    ]


def main() -> int:
    mapper = failure.GtdbMapper()
    hmp_raw_mapper = hmp_raw.raw_audit.RawMapper()
    rules = neg_sweep.build_rules()
    configs_list = configs()
    score_rows: list[dict[str, object]] = []
    panel_set = {str(config["panel"]) for config in configs_list}

    for index, config in enumerate(configs_list, start=1):
        panel = str(config["panel"])
        sample = int(config["sample"])
        log(f"[{index}/{len(configs_list)}] preparing {panel} sample{sample}")
        truth = truth_for_config(config)
        current_pred = current_pred_for_config(config, mapper)
        current_calls = set(current_pred)
        raw = raw_best_for_config(config, mapper, hmp_raw_mapper)
        score_rows.append(
            score_sample(
                panel,
                sample,
                "current_default",
                truth,
                current_pred,
                current_calls,
                set(),
                "current selected default calls",
            )
        )
        for name, rule in rules:
            if name == "baseline_no_raw_add":
                continue
            rescued = rule_species(raw, current_calls, rule)
            score_rows.append(
                score_sample(
                    panel,
                    sample,
                    name,
                    truth,
                    current_pred,
                    current_calls,
                    rescued,
                    f"zero-mass raw candidate retention rule {name}",
                )
            )

    scores = pd.DataFrame(score_rows)
    deltas, panels = summarize(scores)
    overall_df = overall(deltas, panels)
    audit_rows = build_audit(overall_df, len(configs_list), len(panel_set))

    scores.to_csv(RESULTS / "raw_candidate_retention_cross_panel_sample_scores.tsv", sep="\t", index=False)
    deltas.to_csv(RESULTS / "raw_candidate_retention_cross_panel_sample_deltas.tsv", sep="\t", index=False)
    panels.to_csv(RESULTS / "raw_candidate_retention_cross_panel_panel_summary.tsv", sep="\t", index=False)
    overall_df.to_csv(RESULTS / "raw_candidate_retention_cross_panel_overall.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "raw_candidate_retention_cross_panel_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )

    print(overall_df.head(20).to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
