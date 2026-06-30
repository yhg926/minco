#!/usr/bin/env python3
"""Estimate abundance-allocation headroom under the current MinCO call sets.

This is a diagnostic, not a default strategy. It asks whether a perfect
truth-aware abundance allocator, constrained to the species currently called by
MinCO, could close the abundance gap to Sylph. The oracle rows keep MinCO's call
set for F1 accounting and only change the abundance vector.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping

import pandas as pd

import decompose_abundance_errors as decomp


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def clean_abundance(values: Mapping[str, float]) -> dict[str, float]:
    return {
        str(key): max(0.0, finite_float(value))
        for key, value in values.items()
        if str(key)
    }


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) < 2:
        return float("nan")
    return float(pd.Series(left, dtype=float).corr(pd.Series(right, dtype=float), method="pearson"))


def official_l1_kind(panel: str) -> str:
    cfg = decomp.PANELS[panel]
    scores = pd.read_csv(Path(cfg["official_scores"]), sep="\t", nrows=1)
    return "union" if "L1_union_pp" in scores.columns else "truth_only"


def official_pearson_kind(panel: str) -> str:
    cfg = decomp.PANELS[panel]
    scores = pd.read_csv(Path(cfg["official_scores"]), sep="\t", nrows=1)
    return "union" if "Pearson_union" in scores.columns else "truth_only"


def metric_row(
    panel: str,
    sample: int,
    method: str,
    truth: Mapping[str, float],
    call_species: set[str],
    pred_abundance: Mapping[str, float],
    prediction_rule: str,
    truth_source: str,
    pred_source: str,
) -> dict[str, object]:
    truth_norm = decomp.normalize(truth)
    abundance = clean_abundance(pred_abundance)
    truth_species = set(truth_norm)
    pred_species = set(call_species)
    tp = truth_species & pred_species
    fp = pred_species - truth_species
    fn = truth_species - pred_species

    matched_abs = sum(abs(abundance.get(species, 0.0) - truth_norm.get(species, 0.0)) for species in tp)
    missing = sum(truth_norm.get(species, 0.0) for species in fn)
    extra = sum(abundance.get(species, 0.0) for species in fp)
    union_species = sorted(truth_species | pred_species)
    truth_order = sorted(truth_species)
    union_true = [truth_norm.get(species, 0.0) for species in union_species]
    union_pred = [abundance.get(species, 0.0) for species in union_species]
    truth_true = [truth_norm.get(species, 0.0) for species in truth_order]
    truth_pred = [abundance.get(species, 0.0) for species in truth_order]
    precision = len(tp) / len(pred_species) if pred_species else 0.0
    recall = len(tp) / len(truth_species) if truth_species else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    union_l1 = sum(abs(p - t) for p, t in zip(union_pred, union_true)) * 100.0
    truth_l1 = sum(abs(p - t) for p, t in zip(truth_pred, truth_true)) * 100.0
    l1_kind = official_l1_kind(panel)
    pearson_kind = official_pearson_kind(panel)
    return {
        "panel": panel,
        "sample": sample,
        "method": method,
        "truth_species": len(truth_species),
        "pred_species": len(pred_species),
        "TP": len(tp),
        "FP": len(fp),
        "FN": len(fn),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "union_L1_pp": union_l1,
        "truth_only_L1_pp": truth_l1,
        "official_L1_kind": l1_kind,
        "official_L1_pp": union_l1 if l1_kind == "union" else truth_l1,
        "Pearson_union": pearson(union_pred, union_true),
        "Pearson_truth_only": pearson(truth_pred, truth_true),
        "official_Pearson_kind": pearson_kind,
        "official_Pearson": pearson(union_pred, union_true)
        if pearson_kind == "union"
        else pearson(truth_pred, truth_true),
        "matched_abs_error_pp": matched_abs * 100.0,
        "missing_truth_mass_pp": missing * 100.0,
        "extra_pred_mass_pp": extra * 100.0,
        "truth_mass_detected_pct": sum(truth_norm.get(species, 0.0) for species in tp) * 100.0,
        "pred_mass_total_pct": sum(abundance.values()) * 100.0,
        "pred_mass_on_truth_pct": sum(abundance.get(species, 0.0) for species in tp) * 100.0,
        "pred_mass_extra_pct": extra * 100.0,
        "truth_source": truth_source,
        "pred_source": pred_source,
        "prediction_rule": prediction_rule,
    }


class PanelLoader:
    def __init__(self) -> None:
        self.toy_mod = decomp.load_toy_scorer()
        self.by_accession, self.by_core = self.toy_mod.score.truth.load_gtdb_metadata(
            self.toy_mod.score.truth.GTDB_METADATA
        )
        self.hmp_taxmap = decomp.hmp.parse_species_taxmap(decomp.hmp.TAXMAP)
        self.hmp_gastro_taxmap = decomp.hmp_gastro.parse_species_taxmap(decomp.hmp_gastro.TAXMAP)
        (
            _cami3_wgs_to_species,
            self.cami3_taxid_to_species,
            self.cami3_name_to_species,
            _map_diag,
        ) = decomp.cami3.build_transfer_maps()
        self.cami3_taxmap = decomp.cami3.parse_species_taxmap(decomp.cami3.TAXMAP)

    def minco_pred(self, panel: str, sample: int, path: Path, collapse: str) -> tuple[dict[str, float], str]:
        if panel == "hmp_airskin_gtdb_source_abundance":
            decomp.hmp.ensure_sample_record(sample)
            hmp_paths = decomp.hmp.SAMPLES[sample]
            best_ref_species, _diag = decomp.hmp.best_raw_ref_species_by_taxid(
                {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                self.hmp_taxmap,
                self.by_accession,
                self.by_core,
            )
            collapsed, _extra = decomp.hmp.load_minco_predictions(
                path,
                best_ref_species,
                self.by_accession,
                self.by_core,
            )
            return decomp.pred_from_collapsed(collapsed), "official HMP MinCO GTDB mapping; summed and normalized"
        if panel == "hmp_gastrooral_gtdb_source_abundance":
            hmp_paths = decomp.hmp_gastro.SAMPLES[sample]
            best_ref_species, _diag = decomp.hmp_gastro.best_raw_ref_species_by_taxid(
                {"unique": hmp_paths["unique"], "split": hmp_paths["split"]},
                self.hmp_gastro_taxmap,
                self.by_accession,
                self.by_core,
            )
            collapsed, _extra = decomp.hmp_gastro.load_minco_predictions(
                path,
                best_ref_species,
                self.by_accession,
                self.by_core,
            )
            return decomp.pred_from_collapsed(collapsed), "official HMP gastrooral MinCO GTDB mapping; summed and normalized"
        if panel == "cami3_toy_human_gut_gtdb_source_readmap":
            best_ref_species, _diag = decomp.cami3.best_raw_ref_species_by_taxid(
                decomp.cami3.RAW_TABLES[sample],
                self.cami3_taxmap,
                self.by_accession,
                self.by_core,
            )
            collapsed, _extra = decomp.cami3.load_minco_predictions(
                path,
                self.cami3_taxid_to_species,
                self.cami3_name_to_species,
                best_ref_species,
                self.by_accession,
                self.by_core,
            )
            return decomp.pred_from_collapsed(collapsed), "official CAMI3 MinCO GTDB source-readmap mapping; summed and normalized"
        return (
            decomp.load_minco_pred(path, collapse),
            f"calibrated MinCO called rows; {collapse} per GTDB species; normalized",
        )

    def sylph_pred(self, panel: str, sample: int, path: Path, collapse: str) -> tuple[dict[str, float], str]:
        if panel == "cami2_toy_mouse_gut":
            return (
                decomp.load_toy_sylph_pred(
                    path,
                    self.toy_mod,
                    self.by_accession,
                    self.by_core,
                    collapse,
                ),
                f"Sylph Toy Mouse GTDB mapping; {collapse} per GTDB species; normalized",
            )
        return (
            decomp.load_sylph_pred(path, self.by_accession, self.by_core),
            "Sylph GTDB accession mapping; summed per GTDB species; normalized",
        )


def oracle_vectors(
    truth: Mapping[str, float],
    minco_pred: Mapping[str, float],
) -> dict[str, tuple[set[str], dict[str, float], str]]:
    truth_norm = decomp.normalize(truth)
    calls = set(minco_pred)
    detected = calls & set(truth_norm)
    detected_truth_mass = sum(truth_norm.get(species, 0.0) for species in detected)
    minco_tp_mass = sum(max(0.0, finite_float(minco_pred.get(species, 0.0))) for species in detected)
    raw_truth = {species: truth_norm.get(species, 0.0) for species in detected}
    renorm_truth = (
        {species: truth_norm.get(species, 0.0) / detected_truth_mass for species in detected}
        if detected_truth_mass > 0.0
        else {}
    )
    no_fp_current = (
        {species: max(0.0, finite_float(minco_pred.get(species, 0.0))) / minco_tp_mass for species in detected}
        if minco_tp_mass > 0.0
        else {}
    )
    return {
        "oracle_fixed_calls_truth_mass": (
            calls,
            raw_truth,
            "same MinCO call set; true abundance assigned to true-positive calls; false-positive calls get zero; no renormalization",
        ),
        "oracle_fixed_calls_truth_renorm": (
            calls,
            renorm_truth,
            "same MinCO call set; true-positive truth abundance renormalized over detected truth mass; false-positive calls get zero",
        ),
        "oracle_current_tp_abundance_renorm": (
            calls,
            no_fp_current,
            "same MinCO call set; current MinCO abundance kept on true-positive calls then renormalized; false-positive calls get zero",
        ),
    }


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (panel, method), sub in scores.groupby(["panel", "method"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "panel": panel,
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].unique()))),
                "mean_F1": sub["F1"].mean(),
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "pooled_F1": f1,
                "mean_official_L1_pp": sub["official_L1_pp"].mean(),
                "mean_union_L1_pp": sub["union_L1_pp"].mean(),
                "mean_truth_only_L1_pp": sub["truth_only_L1_pp"].mean(),
                "mean_official_Pearson": sub["official_Pearson"].mean(),
                "mean_truth_mass_detected_pct": sub["truth_mass_detected_pct"].mean(),
                "mean_pred_mass_total_pct": sub["pred_mass_total_pct"].mean(),
                "mean_matched_abs_error_pp": sub["matched_abs_error_pp"].mean(),
                "mean_missing_truth_mass_pp": sub["missing_truth_mass_pp"].mean(),
                "mean_extra_pred_mass_pp": sub["extra_pred_mass_pp"].mean(),
            }
        )
    return pd.DataFrame(rows)


def compare_panel_bounds(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel, cfg in decomp.PANELS.items():
        sub = summary.loc[summary["panel"].eq(panel)].set_index("method")
        minco_method = cfg["minco_method"]
        sylph_method = cfg["sylph_method"]
        required = [
            minco_method,
            sylph_method,
            "oracle_fixed_calls_truth_mass",
            "oracle_fixed_calls_truth_renorm",
            "oracle_current_tp_abundance_renorm",
        ]
        if any(method not in sub.index for method in required):
            continue
        current = sub.loc[minco_method]
        sylph = sub.loc[sylph_method]
        oracle_raw = sub.loc["oracle_fixed_calls_truth_mass"]
        oracle_renorm = sub.loc["oracle_fixed_calls_truth_renorm"]
        oracle_current_tp = sub.loc["oracle_current_tp_abundance_renorm"]
        rows.append(
            {
                "panel": panel,
                "samples": current["samples"],
                "current_minco_L1_pp": current["mean_official_L1_pp"],
                "sylph_L1_pp": sylph["mean_official_L1_pp"],
                "oracle_truth_mass_L1_pp": oracle_raw["mean_official_L1_pp"],
                "oracle_truth_renorm_L1_pp": oracle_renorm["mean_official_L1_pp"],
                "oracle_current_tp_renorm_L1_pp": oracle_current_tp["mean_official_L1_pp"],
                "current_minus_sylph_L1_pp": current["mean_official_L1_pp"] - sylph["mean_official_L1_pp"],
                "current_minus_oracle_truth_renorm_L1_pp": current["mean_official_L1_pp"] - oracle_renorm["mean_official_L1_pp"],
                "oracle_truth_renorm_minus_sylph_L1_pp": oracle_renorm["mean_official_L1_pp"] - sylph["mean_official_L1_pp"],
                "current_pooled_F1": current["pooled_F1"],
                "sylph_pooled_F1": sylph["pooled_F1"],
                "mean_truth_mass_detected_pct": current["mean_truth_mass_detected_pct"],
                "allocation_oracle_beats_sylph": oracle_renorm["mean_official_L1_pp"] < sylph["mean_official_L1_pp"],
            }
        )
    return pd.DataFrame(rows)


def validate_against_decomposition(scores: pd.DataFrame) -> pd.DataFrame:
    path = RESULTS / "abundance_error_decomposition.tsv"
    if not path.exists():
        return pd.DataFrame()
    expected = pd.read_csv(path, sep="\t")
    methods = set()
    for cfg in decomp.PANELS.values():
        methods.add(cfg["minco_method"])
        methods.add(cfg["sylph_method"])
    rows: list[dict[str, object]] = []
    for _, row in scores.loc[scores["method"].isin(methods)].iterrows():
        match = expected.loc[
            expected["panel"].astype(str).eq(str(row["panel"]))
            & expected["sample"].astype(str).eq(str(row["sample"]))
            & expected["method"].astype(str).eq(str(row["method"]))
        ]
        if match.empty:
            continue
        exp = match.iloc[0]
        rows.append(
            {
                "panel": row["panel"],
                "sample": row["sample"],
                "method": row["method"],
                "TP_delta": int(row["TP"]) - int(exp["TP"]),
                "FP_delta": int(row["FP"]) - int(exp["FP"]),
                "FN_delta": int(row["FN"]) - int(exp["FN"]),
                "F1_delta": finite_float(row["F1"]) - finite_float(exp["F1"]),
                "union_L1_delta_pp": finite_float(row["union_L1_pp"]) - finite_float(exp["union_L1_pp"]),
                "truth_only_L1_delta_pp": finite_float(row["truth_only_L1_pp"]) - finite_float(exp["truth_only_L1_pp"]),
            }
        )
    return pd.DataFrame(rows)


def build_audit(scores: pd.DataFrame, bounds: pd.DataFrame, validation: pd.DataFrame) -> pd.DataFrame:
    validation_max = 0.0
    if not validation.empty:
        numeric_cols = [
            "TP_delta",
            "FP_delta",
            "FN_delta",
            "F1_delta",
            "union_L1_delta_pp",
            "truth_only_L1_delta_pp",
        ]
        validation_max = float(validation[numeric_cols].abs().max().max())
    recoverable = bounds["current_minus_oracle_truth_renorm_L1_pp"].mean() if not bounds.empty else float("nan")
    remaining = bounds["oracle_truth_renorm_L1_pp"].mean() if not bounds.empty else float("nan")
    sylph = bounds["sylph_L1_pp"].mean() if not bounds.empty else float("nan")
    beats = int(bounds["allocation_oracle_beats_sylph"].sum()) if not bounds.empty else 0
    panels = int(len(bounds))
    decision = (
        "allocation_model_can_close_some_not_all_sylph_gap"
        if beats and beats < panels
        else "allocation_model_insufficient_vs_sylph"
        if beats == 0
        else "allocation_model_sufficient_under_fixed_calls"
    )
    rows = [
        {
            "metric": "evaluated_panels",
            "value": panels,
            "evidence": ",".join(bounds["panel"].astype(str)) if not bounds.empty else "",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "evaluated_samples",
            "value": int(scores[["panel", "sample"]].drop_duplicates().shape[0]),
            "evidence": "current MinCO/Sylph GTDB abundance panels reused from abundance decomposition",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "baseline_validation_max_abs_delta",
            "value": validation_max,
            "evidence": "current MinCO and Sylph rows compared with abundance_error_decomposition.tsv",
            "decision": "pass" if validation_max <= 1e-9 else "check_before_interpretation",
        },
        {
            "metric": "mean_current_minus_oracle_truth_renorm_L1_pp",
            "value": recoverable,
            "evidence": "mean recoverable official L1 if truth-aware allocation is constrained to current MinCO calls",
            "decision": "allocation_headroom",
        },
        {
            "metric": "mean_oracle_truth_renorm_L1_pp",
            "value": remaining,
            "evidence": "remaining official L1 after perfect truth-aware renormalized allocation over detected true species",
            "decision": "lower_bound_under_current_calls",
        },
        {
            "metric": "mean_sylph_L1_pp",
            "value": sylph,
            "evidence": "same panels and official L1 conventions",
            "decision": "reference",
        },
        {
            "metric": "oracle_truth_renorm_beats_sylph_panels",
            "value": f"{beats}/{panels}",
            "evidence": ";".join(
                f"{row.panel}:{row.allocation_oracle_beats_sylph}"
                for row in bounds.itertuples(index=False)
            )
            if not bounds.empty
            else "",
            "decision": decision,
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": "truth-aware upper bound only; not an implementable default strategy",
            "decision": "keep_current_default_and_prioritize_abundance_allocation_research",
        },
    ]
    return pd.DataFrame(rows)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    loader = PanelLoader()
    rows: list[dict[str, object]] = []
    for panel, cfg in decomp.PANELS.items():
        for sample in cfg["samples"]:
            sample_id = int(sample)
            truth_path = Path(cfg["truth"](sample_id))
            truth = decomp.load_truth(truth_path, sample_id, str(cfg["truth_abundance_col"]))
            minco_path = Path(cfg["minco"](sample_id))
            sylph_path = Path(cfg["sylph"](sample_id))
            minco_pred, minco_rule = loader.minco_pred(
                panel,
                sample_id,
                minco_path,
                str(cfg["minco_collapse"]),
            )
            sylph_pred, sylph_rule = loader.sylph_pred(
                panel,
                sample_id,
                sylph_path,
                str(cfg["sylph_collapse"]),
            )
            rows.append(
                metric_row(
                    panel,
                    sample_id,
                    str(cfg["minco_method"]),
                    truth,
                    set(minco_pred),
                    minco_pred,
                    minco_rule,
                    str(truth_path),
                    str(minco_path),
                )
            )
            rows.append(
                metric_row(
                    panel,
                    sample_id,
                    str(cfg["sylph_method"]),
                    truth,
                    set(sylph_pred),
                    sylph_pred,
                    sylph_rule,
                    str(truth_path),
                    str(sylph_path),
                )
            )
            for method, (call_species, abundance, rule) in oracle_vectors(truth, minco_pred).items():
                rows.append(
                    metric_row(
                        panel,
                        sample_id,
                        method,
                        truth,
                        call_species,
                        abundance,
                        rule,
                        str(truth_path),
                        str(minco_path),
                    )
                )
    scores = pd.DataFrame(rows)
    summary = summarize(scores)
    bounds = compare_panel_bounds(summary)
    validation = validate_against_decomposition(scores)
    audit = build_audit(scores, bounds, validation)
    scores.to_csv(RESULTS / "abundance_oracle_bounds_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "abundance_oracle_bounds_summary.tsv", sep="\t", index=False)
    bounds.to_csv(RESULTS / "abundance_oracle_bounds_panel_delta.tsv", sep="\t", index=False)
    validation.to_csv(RESULTS / "abundance_oracle_bounds_validation.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "abundance_oracle_bounds_audit.tsv", sep="\t", index=False)
    print(bounds.to_string(index=False))
    print("\nAUDIT")
    print(audit.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
