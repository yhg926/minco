#!/usr/bin/env python3
"""Score exact-hit abundance trigger reruns for formerly missing samples."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Mapping

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
DECISION_EXP = EXP.parents[0] / "2026-06-27_universal_strategy_decision"
sys.path.insert(0, str(DECISION_EXP))

import sweep_cross_panel_abundance_variants as sweep  # noqa: E402


RERUNS = {
    ("cami3_toy_human_gut_gtdb_source_readmap", 0): {
        "baseline": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample0_baseline_table.tsv"),
        "trigger010": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample0_trigger010.tsv"),
    },
    ("cami3_toy_human_gut_gtdb_source_readmap", 1): {
        "baseline": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample1_baseline_table.tsv"),
        "trigger010": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample1_trigger010.tsv"),
    },
    ("cami3_toy_human_gut_gtdb_source_readmap", 2): {
        "baseline": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample2_baseline_table.tsv"),
        "trigger010": Path("/tmp/minco_exact_trigger_validation_20260630/run/cami3_sample2_trigger010.tsv"),
    },
    ("hmp_airskin_gtdb_source_abundance", 13): {
        "baseline": Path("/tmp/minco_exact_trigger_validation_20260630/run/hmp_airskin13_baseline_table.tsv"),
        "trigger010": Path("/tmp/minco_exact_trigger_validation_20260630/run/hmp_airskin13_trigger010.tsv"),
    },
}


def call_mask(df: pd.DataFrame) -> pd.Series:
    return df["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})


def map_accession_to_gtdb(
    accession: object,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    lookup_accession,
) -> str:
    rec, _method, _key = lookup_accession(
        str(accession or ""),
        by_accession,
        by_core,
    )
    return str(rec.get("gtdb_species", "")) if rec else ""


def load_profile_calls(
    path: Path,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    lookup_accession,
) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    calls = raw.loc[call_mask(raw)].copy()
    species = []
    for row in calls.itertuples(index=False):
        mapped = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                mapped = map_accession_to_gtdb(
                    getattr(row, col),
                    by_accession,
                    by_core,
                    lookup_accession,
                )
                if mapped:
                    break
        species.append(mapped)
    calls["gtdb_species"] = species
    return calls.loc[calls["gtdb_species"].astype(str).astype(bool)].copy()


def load_truth(panel: str, sample: int) -> pd.DataFrame:
    if panel == "cami3_toy_human_gut_gtdb_source_readmap":
        return sweep.cami3_truth(sample)
    if panel == "hmp_airskin_gtdb_source_abundance":
        return sweep.hmp_truth(sample)
    raise ValueError(panel)


def summarize_against_current(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel_summary, _overall = sweep.summarize(scores)
    baseline_method = "validated_current_no_exact"
    metric_cols = [
        "current_official_L1_pp",
        "current_official_Pearson",
        "delta_current_L1_pp",
        "delta_current_Pearson",
    ]
    panel_summary = panel_summary.drop(columns=[c for c in metric_cols if c in panel_summary.columns])
    current = panel_summary.loc[panel_summary["method"].eq(baseline_method)][
        ["panel", "official_L1_pp", "official_Pearson"]
    ].rename(
        columns={
            "official_L1_pp": "current_official_L1_pp",
            "official_Pearson": "current_official_Pearson",
        }
    )
    merged = panel_summary.merge(current, on="panel", how="left")
    merged["delta_current_L1_pp"] = merged["official_L1_pp"] - merged["current_official_L1_pp"]
    merged["delta_current_Pearson"] = merged["official_Pearson"] - merged["current_official_Pearson"]

    expected_panel_count = int(merged["panel"].nunique())
    overall_rows: list[dict[str, object]] = []
    for method, sub in merged.groupby("method", sort=True):
        panel_count = int(sub["panel"].nunique())
        improved = int((sub["delta_current_L1_pp"] < -1e-9).sum())
        worsened = int((sub["delta_current_L1_pp"] > 1e-9).sum())
        max_worse = float(sub["delta_current_L1_pp"].max())
        mean_delta = float(sub["delta_current_L1_pp"].mean())
        if method == baseline_method:
            decision = "validated_current"
        elif method == "validated_exact_hit_ge_0.10" and worsened > 0:
            decision = "reject_raw_regression"
        elif panel_count == expected_panel_count and max_worse <= 1e-9 and improved > 0:
            decision = "promote_default_candidate"
        elif mean_delta < 0.0 and worsened > 0:
            decision = "local_tradeoff_not_default"
        else:
            decision = "worse_or_equal"
        overall_rows.append(
            {
                "method": method,
                "panels": ",".join(sorted(sub["panel"].astype(str).unique())),
                "panel_count": panel_count,
                "mean_official_L1_pp": float(sub["official_L1_pp"].mean()),
                "mean_delta_current_L1_pp": mean_delta,
                "max_worse_current_L1_pp": max_worse,
                "improved_panel_count": improved,
                "worsened_panel_count": worsened,
                "mean_official_Pearson": float(sub["official_Pearson"].mean()),
                "mean_delta_current_Pearson": float(sub["delta_current_Pearson"].mean()),
                "decision": decision,
            }
        )
    overall = pd.DataFrame(overall_rows).sort_values(
        ["decision", "mean_delta_current_L1_pp", "max_worse_current_L1_pp", "method"]
    )
    return merged.sort_values(["panel", "official_L1_pp", "method"]), overall


def score_one(
    panel: str,
    sample: int,
    method: str,
    path: Path,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    lookup_accession,
) -> dict[str, object]:
    calls = load_profile_calls(path, by_accession, by_core, lookup_accession)
    collapse = "sum"
    raw_values = sweep.numeric(calls, "calibrated_abundance").to_numpy(dtype=float)
    pred = sweep.collapse_prediction(calls, raw_values, collapse)
    row = sweep.taxid_score.score_prediction(sample, method, pred, load_truth(panel, sample), {})
    row["panel"] = panel
    row["collapse_rule"] = collapse
    row["profile"] = str(path)
    row["called_rows_mapped"] = len(calls)
    row["called_species_mapped"] = calls["gtdb_species"].nunique()
    return row


def official_l1(row: pd.Series) -> float:
    return float(row["L1_truth_only_pp"] if row["panel"] == "cami2_toy_mouse_gut" else row["L1_union_pp"])


def official_pearson(row: pd.Series) -> float:
    return float(row["Pearson_truth_only"] if row["panel"] == "cami2_toy_mouse_gut" else row["Pearson_union"])


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    toy_mod = sweep.decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    lookup_accession = toy_mod.score.truth.lookup_accession

    rerun_rows: list[dict[str, object]] = []
    for (panel, sample), paths in RERUNS.items():
        for label, path in paths.items():
            if not path.exists():
                raise SystemExit(f"missing rerun profile: {path}")
            method = "validated_current_no_exact" if label == "baseline" else "validated_exact_hit_ge_0.10"
            rerun_rows.append(
                score_one(panel, sample, method, path, by_accession, by_core, lookup_accession)
            )
    rerun_scores = pd.DataFrame(rerun_rows)
    rerun_scores.to_csv(RESULTS / "exact_hit_abundance_trigger_rerun_sample_scores.tsv", sep="\t", index=False)

    previous = pd.read_csv(RESULTS / "exact_hit_abundance_trigger_scores.tsv", sep="\t")
    keys = {(panel, sample) for panel, sample in RERUNS}
    keep = ~previous.apply(lambda row: (row["panel"], int(row["sample"])) in keys, axis=1)
    current_prev = previous.loc[keep & previous["method"].eq("current_calibrated_abundance")].copy()
    trigger_prev = previous.loc[keep & previous["method"].eq("exact_hit_abund_ge_0.10")].copy()
    current_prev["method"] = "validated_current_no_exact"
    trigger_prev["method"] = "validated_exact_hit_ge_0.10"
    full_scores = pd.concat([current_prev, trigger_prev, rerun_scores], ignore_index=True, sort=False)
    full_scores.to_csv(RESULTS / "exact_hit_abundance_trigger_validated_scores.tsv", sep="\t", index=False)

    trigger030_prev = previous.loc[keep & previous["method"].eq("exact_hit_abund_ge_0.30")].copy()
    trigger030_prev["method"] = "validated_exact_hit_ge_0.30"
    trigger030_override = rerun_scores.loc[
        rerun_scores["method"].eq("validated_current_no_exact")
    ].copy()
    trigger030_override["method"] = "validated_exact_hit_ge_0.30"
    full_scores = pd.concat(
        [full_scores, trigger030_prev, trigger030_override],
        ignore_index=True,
        sort=False,
    )
    full_scores.to_csv(RESULTS / "exact_hit_abundance_trigger_validated_scores.tsv", sep="\t", index=False)

    panel_summary, overall = summarize_against_current(full_scores)
    panel_with_sylph = sweep.add_external_sylph(panel_summary)
    panel_summary.to_csv(RESULTS / "exact_hit_abundance_trigger_validated_panel_summary.tsv", sep="\t", index=False)
    panel_with_sylph.to_csv(
        RESULTS / "exact_hit_abundance_trigger_validated_panel_summary_with_sylph.tsv",
        sep="\t",
        index=False,
    )
    overall.to_csv(RESULTS / "exact_hit_abundance_trigger_validated_overall.tsv", sep="\t", index=False)

    sample_delta_rows = []
    for (panel, sample), sub in rerun_scores.groupby(["panel", "sample"]):
        base = sub.loc[sub["method"].eq("validated_current_no_exact")].iloc[0]
        trig = sub.loc[sub["method"].eq("validated_exact_hit_ge_0.10")].iloc[0]
        sample_delta_rows.append(
            {
                "panel": panel,
                "sample": sample,
                "baseline_F1": float(base["F1"]),
                "trigger_F1": float(trig["F1"]),
                "delta_F1": float(trig["F1"]) - float(base["F1"]),
                "baseline_official_L1_pp": official_l1(base),
                "trigger_official_L1_pp": official_l1(trig),
                "delta_official_L1_pp": official_l1(trig) - official_l1(base),
                "baseline_official_Pearson": official_pearson(base),
                "trigger_official_Pearson": official_pearson(trig),
                "delta_official_Pearson": official_pearson(trig) - official_pearson(base),
            }
        )
    sample_delta = pd.DataFrame(sample_delta_rows)
    sample_delta.to_csv(
        RESULTS / "exact_hit_abundance_trigger_rerun_sample_deltas.tsv",
        sep="\t",
        index=False,
    )

    print("RERUN SAMPLE DELTAS")
    print(sample_delta.to_string(index=False))
    print("\nVALIDATED OVERALL")
    print(overall.to_string(index=False))
    print("\nVALIDATED PANEL SUMMARY WITH SYLPH")
    show = panel_with_sylph.loc[
        panel_with_sylph["method"].isin(
            {
                "validated_current_no_exact",
                "validated_exact_hit_ge_0.10",
                "validated_exact_hit_ge_0.30",
                "sylph_external_baseline",
            }
        )
    ].copy()
    print(show.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
