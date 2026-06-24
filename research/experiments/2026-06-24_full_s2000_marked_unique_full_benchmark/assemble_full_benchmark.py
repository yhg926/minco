#!/usr/bin/env python3
"""Assemble the current MinCO/Sylph benchmark comparison from cached results."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "research/experiments/2026-06-24_full_s2000_marked_unique_full_benchmark"
RESULTS = EXP / "results"

HYBRID = ROOT / "research/experiments/2026-06-23_hybrid_full_sketch_unique_marker_fallback"
ABUND = ROOT / "research/experiments/2026-06-24_abundance_feature_calibration"
ANI_SRC = ROOT / "research/experiments/2026-06-23_cami3_source_ref_ani_accuracy"
ANI_SWEEP = ROOT / "research/experiments/2026-06-23_readwise_ani_defake_sweep"
VIRTUAL = ROOT / "research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark"
TOYMOUSE = ROOT / "research/experiments/2026-06-22_minco_sylph_abundance_model"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_row(rows: list[dict[str, str]], column: str, value: str) -> dict[str, str]:
    for row in rows:
        if row.get(column) == value:
            return row
    raise KeyError(f"{value!r} not found in {column!r}")


def f(row: dict[str, str], key: str, default: float | None = None) -> float:
    val = row.get(key, "")
    if val == "":
        if default is None:
            raise KeyError(key)
        return default
    return float(val)


def fmt(value: float | str | int | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return f"{value:.12g}"


def mean2(a: float, b: float) -> float:
    return (a + b) / 2.0


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fields})


def build_f1_abundance() -> list[dict[str, object]]:
    final_rows = read_tsv(ABUND / "results/final_comparison.tsv")
    top_l1 = read_tsv(ABUND / "results/calibration_top200_by_l1.tsv")
    top_f1 = read_tsv(ABUND / "results/calibration_top200_by_f1.tsv")
    mouse_rows = read_tsv(HYBRID / "hybrid_mouse_baseline_plus_mean_metrics.tsv")
    cami_rows = read_tsv(HYBRID / "hybrid_cami3_summary.tsv")
    toy_rows = read_tsv(TOYMOUSE / "intragenus_winner_rescue_mean_metrics.tsv")

    sylph_final = find_row(final_rows, "method", "Sylph")
    sylph_toy = find_row(toy_rows, "method", "sylph_gtdb_profile_reported_taxonomic_abundance")
    sylph_cami = find_row(cami_rows, "method", "sylph")

    current_l1 = find_row(top_l1, "method", "raw_value_sum")
    current_f1 = find_row(top_f1, "method", "raw_value_sum")
    old_mouse = find_row(mouse_rows, "method", "old_ctxmarker_robust_rescue_recalc")
    old_cami = find_row(cami_rows, "method", "old_ctxmarker_robust_rescue_recalc")
    virtual_mouse = find_row(mouse_rows, "method", "hybrid_baseline_plus_full_xny700_scale04")
    virtual_cami = find_row(cami_rows, "method", "hybrid_baseline_plus_full_xny700_scale04")

    rows: list[dict[str, object]] = []

    rows.append(
        {
            "method": "Sylph baseline",
            "strategy_role": "external baseline",
            "refdb_or_signal": "Sylph reported taxonomic profile",
            "benchmark_scope": "mouse_gtdb0-2+cami3_ncbi0-2",
            "sample_count": 6,
            "all_F1": f(sylph_final, "F1"),
            "all_L1_pp": f(sylph_final, "L1"),
            "all_Pearson": mean2(f(sylph_toy, "pearson"), f(sylph_cami, "pearson")),
            "all_TP_mean": mean2(f(sylph_toy, "TP"), f(sylph_cami, "TP")),
            "all_FP_mean": mean2(f(sylph_toy, "FP"), f(sylph_cami, "FP")),
            "all_FN_mean": mean2(f(sylph_toy, "FN"), f(sylph_cami, "FN")),
            "mouse_F1": f(sylph_final, "mouse_F1"),
            "mouse_L1_pp": f(sylph_final, "mouse_L1"),
            "mouse_Pearson": f(sylph_toy, "pearson"),
            "mouse_TP_mean": f(sylph_toy, "TP"),
            "mouse_FP_mean": f(sylph_toy, "FP"),
            "mouse_FN_mean": f(sylph_toy, "FN"),
            "cami3_F1": f(sylph_final, "cami3_F1"),
            "cami3_L1_pp": f(sylph_final, "cami3_L1"),
            "cami3_Pearson": f(sylph_cami, "pearson"),
            "cami3_TP_mean": f(sylph_cami, "TP"),
            "cami3_FP_mean": f(sylph_cami, "FP"),
            "cami3_FN_mean": f(sylph_cami, "FN"),
            "source_file": str(ABUND / "results/final_comparison.tsv"),
            "caveat": "same six samples; Sylph abundance table comes from cached benchmark scorer",
        }
    )

    for row, source_path, method, role, caveat in [
        (
            current_f1,
            ABUND / "results/calibration_top200_by_f1.tsv",
            "MinCO current mixed F1-priority",
            "current MinCO best F1 on six-sample mixed panel",
            "postprocessed raw_value_sum marker_l1_ctxobj_f1_blend; not the virtual S2000 implementation",
        ),
        (
            current_l1,
            ABUND / "results/calibration_top200_by_l1.tsv",
            "MinCO current mixed L1-priority",
            "current MinCO best L1 on six-sample mixed panel",
            "postprocessed raw_value_sum marker_l1_ctxobj_cami_blend; not the virtual S2000 implementation",
        ),
    ]:
        rows.append(
            {
                "method": method,
                "strategy_role": role,
                "refdb_or_signal": row["callset"],
                "benchmark_scope": "mouse_gtdb0-2+cami3_ncbi0-2",
                "sample_count": 6,
                "all_F1": f(row, "all_F1"),
                "all_L1_pp": f(row, "all_L1"),
                "all_Pearson": f(row, "all_Pearson"),
                "all_TP_mean": f(row, "all_TP"),
                "all_FP_mean": f(row, "all_FP"),
                "all_FN_mean": f(row, "all_FN"),
                "mouse_F1": f(row, "mouse_F1"),
                "mouse_L1_pp": f(row, "mouse_L1"),
                "mouse_Pearson": f(row, "mouse_Pearson"),
                "mouse_TP_mean": f(row, "mouse_TP"),
                "mouse_FP_mean": f(row, "mouse_FP"),
                "mouse_FN_mean": f(row, "mouse_FN"),
                "cami3_F1": f(row, "cami3_F1"),
                "cami3_L1_pp": f(row, "cami3_L1"),
                "cami3_Pearson": f(row, "cami3_Pearson"),
                "cami3_TP_mean": f(row, "cami3_TP"),
                "cami3_FP_mean": f(row, "cami3_FP"),
                "cami3_FN_mean": f(row, "cami3_FN"),
                "source_file": str(source_path),
                "caveat": caveat,
            }
        )

    for mouse, cami, method, role, ref_signal, caveat in [
        (
            old_mouse,
            old_cami,
            "MinCO old ctx-marker robust rescue",
            "previous markerdb baseline",
            "coden11 S2000 physical ctx-markerdb",
            "computed from physical ctx-markerdb rows; current robust-depth/intra-genus rescue baseline",
        ),
        (
            virtual_mouse,
            virtual_cami,
            "MinCO full S2000 marked-unique proxy",
            "tested proxy for requested full S2000 marked-unique strategy",
            "full S2000 plus virtual-marker-equivalent marker branch and marker-poor full fallback",
            "proxy row: integrated dual-counter ani is not implemented; marker branch is exact by virtual-index equivalence",
        ),
    ]:
        rows.append(
            {
                "method": method,
                "strategy_role": role,
                "refdb_or_signal": ref_signal,
                "benchmark_scope": "mouse_gtdb0-2+cami3_ncbi0-2",
                "sample_count": 6,
                "all_F1": mean2(f(mouse, "F1"), f(cami, "F1")),
                "all_L1_pp": mean2(f(mouse, "l1_pct_points"), f(cami, "l1_pct_points")),
                "all_Pearson": mean2(f(mouse, "pearson"), f(cami, "pearson")),
                "all_TP_mean": mean2(f(mouse, "TP"), f(cami, "TP")),
                "all_FP_mean": mean2(f(mouse, "FP"), f(cami, "FP")),
                "all_FN_mean": mean2(f(mouse, "FN"), f(cami, "FN")),
                "mouse_F1": f(mouse, "F1"),
                "mouse_L1_pp": f(mouse, "l1_pct_points"),
                "mouse_Pearson": f(mouse, "pearson"),
                "mouse_TP_mean": f(mouse, "TP"),
                "mouse_FP_mean": f(mouse, "FP"),
                "mouse_FN_mean": f(mouse, "FN"),
                "cami3_F1": f(cami, "F1"),
                "cami3_L1_pp": f(cami, "l1_pct_points"),
                "cami3_Pearson": f(cami, "pearson"),
                "cami3_TP_mean": f(cami, "TP"),
                "cami3_FP_mean": f(cami, "FP"),
                "cami3_FN_mean": f(cami, "FN"),
                "source_file": str(HYBRID / "hybrid_mouse_baseline_plus_mean_metrics.tsv"),
                "caveat": caveat,
            }
        )

    return rows


def build_ani() -> list[dict[str, object]]:
    alt = read_tsv(ANI_SRC / "alternative_read_ani_summary.tsv")
    sweep = read_tsv(ANI_SWEEP / "summary.tsv")
    rows: list[dict[str, object]] = []

    methods = [
        "minco_coden15_formula_af_Ref_zip_aaf_ani",
        "sylph_Adjusted_ANI",
        "minco_old_ctxmarker_Ref_zip_aaf_ani",
        "sylph_Naive_ANI",
    ]
    for method in methods:
        row = next(r for r in alt if r["method"] == method and r["source_choice"] == "read_major_source")
        rows.append(
            {
                "dataset": "CAMI3 ToyGut selected calls",
                "truth": "skani source-genome-to-called-reference ANI",
                "method": method,
                "source_choice": row["source_choice"],
                "n": f(row, "n_calls"),
                "pearson": f(row, "pearson"),
                "spearman": f(row, "spearman"),
                "mae": f(row, "mae"),
                "rmse": "",
                "mean_error": f(row, "mean_error"),
                "high_error_gt_0p02": "",
                "high_error_gt_0p03": "",
                "source_file": str(ANI_SRC / "alternative_read_ani_summary.tsv"),
                "caveat": "selected calls with same-species source pairing; read_major_source is the primary non-optimistic pairing",
            }
        )

    best = find_row(sweep, "run", "toymouse_sample0_c15_adjusted_ani_product1_p090")
    rows.append(
        {
            "dataset": "Toy Mouse sample0 source-positive GTDB representatives",
            "truth": "source-to-GTDB-representative ANIm",
            "method": "minco_coden15_product1_p90_ANI_median2_else_max",
            "source_choice": "source-positive",
            "n": f(best, "n"),
            "pearson": f(best, "pearson"),
            "spearman": f(best, "spearman"),
            "mae": f(best, "mae"),
            "rmse": f(best, "rmse"),
            "mean_error": f(best, "mean_error"),
            "high_error_gt_0p02": f(best, "high_error_gt_0p02"),
            "high_error_gt_0p03": f(best, "high_error_gt_0p03"),
            "source_file": str(ANI_SWEEP / "summary.tsv"),
            "caveat": "best MinCO ANI reporter; no Sylph row was generated in this toy-mouse ANIm sweep",
        }
    )

    zip_row = next(
        r
        for r in sweep
        if r["run"] == "toymouse_sample0_c15_adjusted_ani_p075" and r["estimator"] == "ANI_zip_aaf"
    )
    rows.append(
        {
            "dataset": "Toy Mouse sample0 source-positive GTDB representatives",
            "truth": "source-to-GTDB-representative ANIm",
            "method": "minco_coden15_ZIP_AAF_baseline",
            "source_choice": "source-positive",
            "n": f(zip_row, "n"),
            "pearson": f(zip_row, "pearson"),
            "spearman": f(zip_row, "spearman"),
            "mae": f(zip_row, "mae"),
            "rmse": f(zip_row, "rmse"),
            "mean_error": f(zip_row, "mean_error"),
            "high_error_gt_0p02": f(zip_row, "high_error_gt_0p02"),
            "high_error_gt_0p03": f(zip_row, "high_error_gt_0p03"),
            "source_file": str(ANI_SWEEP / "summary.tsv"),
            "caveat": "baseline for same toy-mouse ANI sweep",
        }
    )

    return rows


def build_decision(f1_rows: list[dict[str, object]], ani_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_method = {str(row["method"]): row for row in f1_rows}
    ani_by_method = {str(row["method"]): row for row in ani_rows}
    virtual = by_method["MinCO full S2000 marked-unique proxy"]
    old = by_method["MinCO old ctx-marker robust rescue"]
    sylph = by_method["Sylph baseline"]
    minco_f1 = by_method["MinCO current mixed F1-priority"]
    minco_l1 = by_method["MinCO current mixed L1-priority"]
    minco_ani = ani_by_method["minco_coden15_formula_af_Ref_zip_aaf_ani"]
    sylph_ani = ani_by_method["sylph_Adjusted_ANI"]

    return [
        {
            "priority_metric": "six-sample F1",
            "winner": "Sylph baseline",
            "winner_value": sylph["all_F1"],
            "runner_up": "MinCO current mixed F1-priority",
            "runner_up_value": minco_f1["all_F1"],
            "delta_winner_minus_runner": float(sylph["all_F1"]) - float(minco_f1["all_F1"]),
            "interpretation": "Sylph is still marginally ahead on mixed F1; margin is about 0.0011.",
        },
        {
            "priority_metric": "six-sample abundance L1",
            "winner": "Sylph baseline",
            "winner_value": sylph["all_L1_pp"],
            "runner_up": "MinCO current mixed L1-priority",
            "runner_up_value": minco_l1["all_L1_pp"],
            "delta_winner_minus_runner": float(sylph["all_L1_pp"]) - float(minco_l1["all_L1_pp"]),
            "interpretation": "Lower is better; Sylph remains clearly better by 2.785 pp.",
        },
        {
            "priority_metric": "MinCO-only S2000 marked-unique effect",
            "winner": "MinCO full S2000 marked-unique proxy",
            "winner_value": f"F1={fmt(virtual['all_F1'])}; L1={fmt(virtual['all_L1_pp'])}",
            "runner_up": "MinCO old ctx-marker robust rescue",
            "runner_up_value": f"F1={fmt(old['all_F1'])}; L1={fmt(old['all_L1_pp'])}",
            "delta_winner_minus_runner": f"F1_delta={fmt(float(virtual['all_F1']) - float(old['all_F1']))}; L1_delta={fmt(float(virtual['all_L1_pp']) - float(old['all_L1_pp']))}",
            "interpretation": "Virtual-marked full S2000 proxy is a small MinCO-only gain, driven by Toy Mouse; CAMI3 is neutral.",
        },
        {
            "priority_metric": "CAMI3 source/ref ANI MAE",
            "winner": "MinCO coden15 Ref_zip_aaf_ani",
            "winner_value": minco_ani["mae"],
            "runner_up": "Sylph Adjusted_ANI",
            "runner_up_value": sylph_ani["mae"],
            "delta_winner_minus_runner": float(minco_ani["mae"]) - float(sylph_ani["mae"]),
            "interpretation": "Lower is better; MinCO coden15 ZIP-AAF is slightly better on this source/ref ANI subset.",
        },
    ]


def build_summary(f1_rows: list[dict[str, object]], ani_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_method = {str(row["method"]): row for row in f1_rows}
    ani_by_method = {str(row["method"]): row for row in ani_rows}
    virtual_index = {row["metric"]: row for row in read_tsv(VIRTUAL / "summary.tsv")}

    return [
        {
            "metric": "main_answer",
            "value": "not yet a universal win",
            "unit": "",
            "notes": "Full S2000 marked-unique proxy improves the old MinCO ctx-marker row slightly, but Sylph remains best on six-sample F1 and abundance L1.",
        },
        {
            "metric": "sylph_all_F1",
            "value": by_method["Sylph baseline"]["all_F1"],
            "unit": "F1",
            "notes": "mouse_gtdb0-2+cami3_ncbi0-2",
        },
        {
            "metric": "best_minco_all_F1",
            "value": by_method["MinCO current mixed F1-priority"]["all_F1"],
            "unit": "F1",
            "notes": "current mixed-panel F1-priority row",
        },
        {
            "metric": "sylph_all_L1",
            "value": by_method["Sylph baseline"]["all_L1_pp"],
            "unit": "percentage_points",
            "notes": "lower is better",
        },
        {
            "metric": "best_minco_all_L1",
            "value": by_method["MinCO current mixed L1-priority"]["all_L1_pp"],
            "unit": "percentage_points",
            "notes": "lower is better",
        },
        {
            "metric": "virtual_s2000_proxy_all_F1",
            "value": by_method["MinCO full S2000 marked-unique proxy"]["all_F1"],
            "unit": "F1",
            "notes": "proxy for integrated full S2000 marked-unique strategy",
        },
        {
            "metric": "virtual_s2000_proxy_all_L1",
            "value": by_method["MinCO full S2000 marked-unique proxy"]["all_L1_pp"],
            "unit": "percentage_points",
            "notes": "proxy for integrated full S2000 marked-unique strategy",
        },
        {
            "metric": "virtual_vs_physical_marker_mismatched_refs",
            "value": virtual_index["virtual_vs_physical_mismatched_refs"]["value"],
            "unit": "refs",
            "notes": "exact-marker equivalence prerequisite for proxy interpretation",
        },
        {
            "metric": "refs_marker_lt_500",
            "value": virtual_index["refs_marker_lt_500"]["value"],
            "unit": "refs",
            "notes": "default low-marker warning count in S2000 refdb",
        },
        {
            "metric": "best_minco_cami3_source_ref_ani_mae",
            "value": ani_by_method["minco_coden15_formula_af_Ref_zip_aaf_ani"]["mae"],
            "unit": "ANI",
            "notes": "read_major_source CAMI3 selected-call source/ref comparison",
        },
        {
            "metric": "sylph_adjusted_cami3_source_ref_ani_mae",
            "value": ani_by_method["sylph_Adjusted_ANI"]["mae"],
            "unit": "ANI",
            "notes": "read_major_source CAMI3 selected-call source/ref comparison",
        },
        {
            "metric": "best_minco_toymouse_animm_ani_mae",
            "value": ani_by_method["minco_coden15_product1_p90_ANI_median2_else_max"]["mae"],
            "unit": "ANI",
            "notes": "Toy Mouse sample0 source-positive GTDB representative ANIm sweep",
        },
    ]


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    f1_rows = build_f1_abundance()
    ani_rows = build_ani()
    decision_rows = build_decision(f1_rows, ani_rows)
    summary_rows = build_summary(f1_rows, ani_rows)

    f1_fields = [
        "method",
        "strategy_role",
        "refdb_or_signal",
        "benchmark_scope",
        "sample_count",
        "all_F1",
        "all_L1_pp",
        "all_Pearson",
        "all_TP_mean",
        "all_FP_mean",
        "all_FN_mean",
        "mouse_F1",
        "mouse_L1_pp",
        "mouse_Pearson",
        "mouse_TP_mean",
        "mouse_FP_mean",
        "mouse_FN_mean",
        "cami3_F1",
        "cami3_L1_pp",
        "cami3_Pearson",
        "cami3_TP_mean",
        "cami3_FP_mean",
        "cami3_FN_mean",
        "source_file",
        "caveat",
    ]
    ani_fields = [
        "dataset",
        "truth",
        "method",
        "source_choice",
        "n",
        "pearson",
        "spearman",
        "mae",
        "rmse",
        "mean_error",
        "high_error_gt_0p02",
        "high_error_gt_0p03",
        "source_file",
        "caveat",
    ]
    decision_fields = [
        "priority_metric",
        "winner",
        "winner_value",
        "runner_up",
        "runner_up_value",
        "delta_winner_minus_runner",
        "interpretation",
    ]
    summary_fields = ["metric", "value", "unit", "notes"]

    write_tsv(RESULTS / "f1_abundance_comparison.tsv", f1_rows, f1_fields)
    write_tsv(RESULTS / "ani_accuracy_comparison.tsv", ani_rows, ani_fields)
    write_tsv(RESULTS / "benchmark_decision.tsv", decision_rows, decision_fields)
    write_tsv(EXP / "summary.tsv", summary_rows, summary_fields)

    print(f"wrote {RESULTS / 'f1_abundance_comparison.tsv'}")
    print(f"wrote {RESULTS / 'ani_accuracy_comparison.tsv'}")
    print(f"wrote {RESULTS / 'benchmark_decision.tsv'}")
    print(f"wrote {EXP / 'summary.tsv'}")


if __name__ == "__main__":
    main()
