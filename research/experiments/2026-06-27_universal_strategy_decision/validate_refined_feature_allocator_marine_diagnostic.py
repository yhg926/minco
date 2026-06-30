#!/usr/bin/env python3
"""Independent marine diagnostic replay for the refined feature allocator.

This is a cached-profile replay on CAMI II marine exact-split profiles that
were not part of the 38-profile refined-guard selection cache. The scorer uses
the existing conservative GTDB taxid-transfer marine truth, which is useful as
diagnostic stress evidence but is not release-grade GTDB truth.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
sys.path.insert(0, str(ROOT))

from scripts import minco_profile_calibrated as wrapper  # noqa: E402
import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402
import score_marine_gtdb_taxid_transfer as marine  # noqa: E402

TRUTH_HELPER = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))
import build_and_score_gtdb_ground_truth as truth  # noqa: E402


EXACT_DIR = Path("/tmp/minco_exactsplit_universal_20260626")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
SAMPLES = {
    0: EXACT_DIR / "marine0_universal_exactsplit_strategy.tsv",
    2: EXACT_DIR / "marine2_universal_exactsplit_strategy.tsv",
}
SWITCH = wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230

OUT_SCORES = RESULTS / "feature_allocator_refined_marine_diagnostic_scores.tsv"
OUT_SUMMARY = RESULTS / "feature_allocator_refined_marine_diagnostic_summary.tsv"
OUT_VALIDATION = RESULTS / "feature_allocator_refined_marine_diagnostic_validation.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_refined_marine_diagnostic_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_REFINED_MARINE_DIAGNOSTIC.md"

METHOD_BASE = "marine_exactsplit_current_gtdb_transfer"
METHOD_REFINED = "marine_exactsplit_refined_allocator_xny230"


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def finite_max_abs(values: Iterable[object]) -> float:
    clean = []
    for value in values:
        out = finite_float(value, float("nan"))
        if math.isfinite(out):
            clean.append(abs(out))
    return max(clean) if clean else 0.0


def raw_values(rows: pd.DataFrame, col: str) -> np.ndarray:
    raw = pd.to_numeric(rows[col], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return np.where(raw.to_numpy(dtype=float) > 0.0, raw.to_numpy(dtype=float), 0.0)


def current_mask(rows: pd.DataFrame) -> np.ndarray:
    return rows["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()


def truth_for_sample(
    sample: int,
    taxid_to_gtdb: Mapping[str, str],
    ambiguous_taxids: Mapping[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    ba_truth, scope = marine.read_marine_ba_truth(sample)
    gtdb_truth, quality = taxid_score.transfer_truth_to_gtdb(
        sample,
        ba_truth,
        dict(taxid_to_gtdb),
        dict(ambiguous_taxids),
    )
    quality.update(scope)
    ba_mass = float(scope["gold_bacteria_archaea_species_mass_pct_all"])
    quality["truth_mass_mapped_pct_bacteria_archaea"] = (
        float(quality["truth_mass_mapped_pct_all"]) / ba_mass * 100.0 if ba_mass else 0.0
    )
    return gtdb_truth, quality


def collapse_profile_prediction(
    rows: pd.DataFrame,
    call: np.ndarray,
    abundance_raw: np.ndarray,
    taxid_to_gtdb: Mapping[str, str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    pred_rows: list[tuple[str, float, float]] = []
    unmapped = 0
    selected = rows.loc[call].copy()
    selected_raw = np.asarray(abundance_raw, dtype=float)[np.flatnonzero(call)]
    for idx, row in enumerate(selected.itertuples(index=False)):
        species = str(taxid_to_gtdb.get(str(getattr(row, "taxid", "")), ""))
        if not species:
            unmapped += 1
            continue
        abundance = finite_float(selected_raw[idx])
        ani = max(
            finite_float(getattr(row, "s_Ref_zip_aaf_ani_max", 0.0)),
            finite_float(getattr(row, "u_Ref_zip_aaf_ani_max", 0.0)),
        )
        pred_rows.append((species, abundance, ani))
    return taxid_score.collapse_prediction(pred_rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_mapped": int(len(selected) - unmapped),
    }


def score_profile(
    sample: int,
    method: str,
    profile: Path,
    rows: pd.DataFrame,
    call: np.ndarray,
    abundance_raw: np.ndarray,
    truth_df: pd.DataFrame,
    taxid_to_gtdb: Mapping[str, str],
    extra: Mapping[str, object] | None = None,
) -> dict[str, object]:
    pred, pred_extra = collapse_profile_prediction(rows, call, abundance_raw, taxid_to_gtdb)
    details = dict(pred_extra)
    details.update(extra or {})
    details.update(
        {
            "sample_key": f"marine{sample}",
            "dataset": "marine",
            "profile": str(profile),
            "truth_namespace": "conservative_GTDB_taxid_transfer",
            "audit_grade": "diagnostic_nonrelease",
        }
    )
    score = taxid_score.score_prediction(sample, method, pred, truth_df, details)
    score["L1"] = finite_float(score.get("L1_union_pp")) / 100.0
    score["Pearson"] = finite_float(score.get("Pearson_union"), float("nan"))
    return score


def evaluate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = [marine.GOLD, TAXMAP, *SAMPLES.values()]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required marine diagnostic inputs:\n" + "\n".join(missing))

    by_accession, _by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    taxid_to_gtdb, ambiguous_taxids = taxid_score.build_taxid_transfer(by_accession)
    allocator_taxmap = wrapper.parse_species_taxmap(TAXMAP)

    score_rows: list[dict[str, object]] = []
    validation_rows: list[dict[str, object]] = []
    quality_rows: list[dict[str, object]] = []
    for sample, profile in SAMPLES.items():
        truth_df, quality = truth_for_sample(sample, taxid_to_gtdb, ambiguous_taxids)
        quality["sample"] = sample
        quality_rows.append(quality)
        rows = pd.read_csv(profile, sep="\t", low_memory=False)
        call = current_mask(rows)
        raw = raw_values(rows, "calibrated_abundance_raw")
        calibrated = raw_values(rows, "calibrated_abundance")
        base_raw_score = score_profile(
            sample,
            METHOD_BASE,
            profile,
            rows,
            call,
            raw,
            truth_df,
            taxid_to_gtdb,
        )
        base_calibrated_score = score_profile(
            sample,
            "marine_exactsplit_current_calibrated_column_gtdb_transfer",
            profile,
            rows,
            call,
            calibrated,
            truth_df,
            taxid_to_gtdb,
        )
        validation_rows.append(
            {
                "sample": sample,
                "F1_delta_raw_vs_calibrated": finite_float(base_raw_score["F1"])
                - finite_float(base_calibrated_score["F1"]),
                "L1_delta_raw_vs_calibrated": finite_float(base_raw_score["L1"])
                - finite_float(base_calibrated_score["L1"]),
                "Pearson_delta_raw_vs_calibrated": finite_float(
                    base_raw_score["Pearson"],
                    float("nan"),
                )
                - finite_float(base_calibrated_score["Pearson"], float("nan")),
                "pred_rows_mapped": base_raw_score.get("pred_rows_mapped", 0),
                "pred_rows_unmapped": base_raw_score.get("pred_rows_unmapped", 0),
                "truth_mass_mapped_pct_bacteria_archaea": quality.get(
                    "truth_mass_mapped_pct_bacteria_archaea",
                    0.0,
                ),
            }
        )
        adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
            rows,
            call,
            raw,
            SWITCH,
            allocator_taxmap,
        )
        refined_score = score_profile(
            sample,
            METHOD_REFINED,
            profile,
            rows,
            call,
            adjusted,
            truth_df,
            taxid_to_gtdb,
            details,
        )
        score_rows.extend([base_raw_score, refined_score])

    scores = pd.DataFrame(score_rows)
    base = scores.loc[scores["method"].eq(METHOD_BASE)].set_index("sample")
    rows_with_delta = []
    for row in scores.to_dict("records"):
        if row["method"] == METHOD_BASE:
            row["delta_F1"] = 0.0
            row["delta_L1"] = 0.0
            row["delta_Pearson"] = 0.0
        else:
            base_row = base.loc[int(row["sample"])]
            row["delta_F1"] = finite_float(row["F1"]) - finite_float(base_row["F1"])
            row["delta_L1"] = finite_float(row["L1"]) - finite_float(base_row["L1"])
            row["delta_Pearson"] = finite_float(row["Pearson"], float("nan")) - finite_float(
                base_row["Pearson"],
                float("nan"),
            )
        rows_with_delta.append(row)
    scores = pd.DataFrame(rows_with_delta)
    validation = pd.DataFrame(validation_rows)
    quality = pd.DataFrame(quality_rows)
    return scores, validation, quality


def summarize(scores: pd.DataFrame, validation: pd.DataFrame, quality: pd.DataFrame) -> pd.DataFrame:
    refined = scores.loc[scores["method"].eq(METHOD_REFINED)].copy()
    base = scores.loc[scores["method"].eq(METHOD_BASE)].copy()
    rows = []
    for method, sub in scores.groupby("method", sort=True):
        rows.append(
            {
                "method": method,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(sub["sample"].nunique()),
                "mean_F1": float(sub["F1"].mean()),
                "mean_L1": float(sub["L1"].mean()),
                "mean_Pearson": float(sub["Pearson"].mean()),
                "mean_delta_F1": float(sub["delta_F1"].mean()),
                "mean_delta_L1": float(sub["delta_L1"].mean()),
                "max_worse_L1": float(sub["delta_L1"].max()),
                "improved_L1_samples": int((sub["delta_L1"] < -1e-12).sum()),
                "worsened_L1_samples": int((sub["delta_L1"] > 1e-12).sum()),
                "switched_samples": int(
                    sub.get("abundance_feature_allocator_applied", pd.Series([], dtype=object))
                    .astype(str)
                    .str.lower()
                    .isin({"true", "1"})
                    .sum()
                ),
            }
        )
    overall = {
        "method": "refined_vs_current",
        "samples": ",".join(map(str, sorted(refined["sample"].astype(int).unique()))),
        "sample_count": int(refined["sample"].nunique()),
        "mean_F1": float(refined["F1"].mean()),
        "mean_L1": float(refined["L1"].mean()),
        "mean_Pearson": float(refined["Pearson"].mean()),
        "mean_delta_F1": float(refined["delta_F1"].mean()),
        "mean_delta_L1": float(refined["delta_L1"].mean()),
        "max_worse_L1": float(refined["delta_L1"].max()),
        "improved_L1_samples": int((refined["delta_L1"] < -1e-12).sum()),
        "worsened_L1_samples": int((refined["delta_L1"] > 1e-12).sum()),
        "switched_samples": int(
            refined["abundance_feature_allocator_applied"].astype(str).str.lower().isin({"true", "1"}).sum()
        ),
        "mean_truth_mass_mapped_pct_bacteria_archaea": float(
            quality["truth_mass_mapped_pct_bacteria_archaea"].mean()
        ),
        "min_truth_mass_mapped_pct_bacteria_archaea": float(
            quality["truth_mass_mapped_pct_bacteria_archaea"].min()
        ),
        "max_validation_L1_delta_raw_vs_calibrated": finite_max_abs(
            validation["L1_delta_raw_vs_calibrated"]
        ),
    }
    rows.append(overall)
    return pd.DataFrame(rows)


def audit_rows(summary: pd.DataFrame, validation: pd.DataFrame) -> list[dict[str, object]]:
    overall = summary.loc[summary["method"].eq("refined_vs_current")].iloc[0]
    max_validation = finite_max_abs(validation["L1_delta_raw_vs_calibrated"])
    worsened = int(overall["worsened_L1_samples"])
    improved = int(overall["improved_L1_samples"])
    min_truth = float(overall["min_truth_mass_mapped_pct_bacteria_archaea"])
    decision = (
        "diagnostic_supports_refined_allocator_holdout_candidate"
        if worsened == 0 and improved > 0
        else "diagnostic_warns_refined_allocator_not_sample_safe"
    )
    return [
        {
            "metric": "diagnostic_scope",
            "value": (
                "samples=marine0,marine2;"
                f"min_truth_mass_mapped_pct_bacteria_archaea={min_truth:.6f}"
            ),
            "evidence": f"{OUT_SCORES.relative_to(EXP)};{OUT_VALIDATION.relative_to(EXP)}",
            "decision": "diagnostic_nonrelease_truth_transfer",
        },
        {
            "metric": "baseline_internal_validation",
            "value": f"max_L1_delta_raw_vs_calibrated={max_validation:.12g}",
            "evidence": str(OUT_VALIDATION.relative_to(EXP)),
            "decision": "pass" if max_validation < 1e-9 else "review",
        },
        {
            "metric": "refined_allocator_marine_effect",
            "value": (
                f"samples={int(overall['sample_count'])};"
                f"switched={int(overall['switched_samples'])};"
                f"improved={improved};worsened={worsened};"
                f"mean_L1_delta={float(overall['mean_delta_L1']):.9f};"
                f"max_worse_L1={float(overall['max_worse_L1']):.9f};"
                f"mean_F1_delta={float(overall['mean_delta_F1']):.9f}"
            ),
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "independent_diagnostic_effect",
        },
        {
            "metric": "promotion_decision",
            "value": decision,
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "current_default_unchanged",
        },
    ]


def write_markdown(summary: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Refined Feature Allocator Marine Diagnostic",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-profile replay evaluates the opt-in",
        f"`--abundance-feature-allocator-switch {SWITCH}` on CAMI II marine",
        "exact-split profiles that were not used to choose the refined guard.",
        "The truth namespace is conservative GTDB taxid transfer, so this is",
        "diagnostic evidence and not release-grade evidence.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "diagnostic_scope",
        "baseline_internal_validation",
        "refined_allocator_marine_effect",
        "promotion_decision",
    ]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Method | Samples | Mean F1 | Mean L1 | Mean Pearson | Mean L1 delta | Worsened L1 samples | Switched |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary.to_dict("records"):
        lines.append(
            "| {method} | {samples} | {f1:.6f} | {l1:.6f} | {pearson:.6f} | {delta:.9f} | {worse} | {switched} |".format(
                method=row["method"],
                samples=row["samples"],
                f1=float(row["mean_F1"]),
                l1=float(row["mean_L1"]),
                pearson=float(row["mean_Pearson"]),
                delta=float(row["mean_delta_L1"]),
                worse=int(row["worsened_L1_samples"]),
                switched=int(row["switched_samples"]),
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Keep the current default unchanged.",
            "- Treat this as independent diagnostic stress support only.",
            "- Release-grade promotion still needs clean GTDB holdout validation.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SCORES.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_VALIDATION.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    scores, validation, quality = evaluate()
    summary = summarize(scores, validation, quality)
    audit = audit_rows(summary, validation)

    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    validation.to_csv(OUT_VALIDATION, sep="\t", index=False)
    pd.DataFrame(audit).to_csv(OUT_AUDIT, sep="\t", index=False)
    write_markdown(summary, audit)
    print(pd.DataFrame(audit).to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
