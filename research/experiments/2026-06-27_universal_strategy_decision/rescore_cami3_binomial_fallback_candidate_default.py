#!/usr/bin/env python3
"""Rescore selected-default CAMI3 profiles against binomial-fallback truth.

This is the release-evidence companion to
rescore_cami3_source_readmap_binomial_fallback.py.  It uses the same adjusted
truth table but scores the selected candidate-preset MinCO profiles, so the
result reflects the current no-expertise default instead of an older cached
CAMI3 profile set.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import pandas as pd

import rescore_cami3_source_readmap_binomial_fallback as fallback
import score_cami3_gtdb_source_readmap as source_score
import score_cami3_gtdb_taxid_transfer as taxid_score


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
CANDIDATE_ROOT = Path("/tmp/minco_candidate_preset_replay_20260629")
CANDIDATE_PROFILES = {
    0: CANDIDATE_ROOT / "cami3_sample0.candidate_preset.tsv",
    1: CANDIDATE_ROOT / "cami3_sample1.candidate_preset.tsv",
    2: CANDIDATE_ROOT / "cami3_sample2.candidate_preset.tsv",
}

OUT_SCORES = RESULTS / "cami3_binomial_fallback_candidate_default_scores.tsv"
OUT_SUMMARY = RESULTS / "cami3_binomial_fallback_candidate_default_summary.tsv"
OUT_DELTA = RESULTS / "cami3_binomial_fallback_candidate_default_delta.tsv"
OUT_AUDIT = RESULTS / "cami3_binomial_fallback_candidate_default_audit.tsv"
OUT_MD = EXP / "CAMI3_BINOMIAL_FALLBACK_CANDIDATE_DEFAULT.md"


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
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


def load_candidate_predictions():
    _wgs_to_species, taxid_to_species, name_to_species, _map_diag = source_score.build_transfer_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(source_score.TAXMAP)
    minco_preds: dict[int, pd.DataFrame] = {}
    sylph_preds: dict[int, pd.DataFrame] = {}
    minco_extra: dict[int, dict[str, object]] = {}
    sylph_extra: dict[int, dict[str, object]] = {}

    for sample_id, profile in CANDIDATE_PROFILES.items():
        if not profile.exists():
            raise SystemExit(f"missing candidate-preset CAMI3 profile: {profile}")
        best_ref_species, _best_ref_diag = source_score.best_raw_ref_species_by_taxid(
            source_score.RAW_TABLES[sample_id],
            taxmap,
            by_accession,
            by_core,
        )
        minco_pred, minco_meta = source_score.load_minco_predictions(
            profile,
            taxid_to_species,
            name_to_species,
            best_ref_species,
            by_accession,
            by_core,
        )
        sylph_pred, sylph_meta = taxid_score.load_sylph_predictions(
            source_score.PROFILE_PATHS[sample_id]["sylph"],
            by_accession,
            by_core,
        )
        minco_preds[sample_id] = minco_pred
        sylph_preds[sample_id] = sylph_pred
        minco_extra[sample_id] = {
            **minco_meta,
            "profile": str(profile),
            "selected_default_profile_preset": "candidate",
        }
        sylph_extra[sample_id] = sylph_meta
    return minco_preds, sylph_preds, minco_extra, sylph_extra


def score_candidate_default(adjusted_truth: pd.DataFrame) -> pd.DataFrame:
    minco_preds, sylph_preds, minco_extra, sylph_extra = load_candidate_predictions()
    score_rows = []
    for sample_id in sorted(CANDIDATE_PROFILES):
        truth_df = adjusted_truth.loc[adjusted_truth["sample"].astype(str).eq(str(sample_id))].copy()
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "minco_candidate_default_binomial_fallback",
                minco_preds[sample_id],
                truth_df,
                minco_extra[sample_id],
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_binomial_fallback",
                sylph_preds[sample_id],
                truth_df,
                sylph_extra[sample_id],
            )
        )
    return pd.DataFrame(score_rows)


def build_delta(scores: pd.DataFrame) -> pd.DataFrame:
    previous = pd.read_csv(fallback.OUT_SCORES, sep="\t")
    method_map = {
        "minco_candidate_default_binomial_fallback": "minco_source_readmap_binomial_fallback_rescore",
        "sylph_binomial_fallback": "sylph_source_readmap_binomial_fallback_rescore",
    }
    rows = []
    metrics = ["TP", "FP", "FN", "F1", "L1_union_pp", "Pearson_union"]
    for row in scores.itertuples(index=False):
        method = str(getattr(row, "method"))
        sample = str(getattr(row, "sample"))
        baseline_method = method_map[method]
        base = previous.loc[
            previous["sample"].astype(str).eq(sample) & previous["method"].eq(baseline_method)
        ]
        if base.empty:
            continue
        base_one = base.iloc[0]
        out = {
            "sample": sample,
            "method": method,
            "baseline_method": baseline_method,
        }
        for metric in metrics:
            adjusted = finite(getattr(row, metric))
            old = finite(base_one.get(metric))
            out[f"baseline_{metric}"] = old
            out[f"selected_default_{metric}"] = adjusted
            out[f"delta_{metric}"] = adjusted - old
        rows.append(out)
    return pd.DataFrame(rows)


def build_audit(summary: pd.DataFrame, delta: pd.DataFrame) -> list[dict[str, object]]:
    minco = summary.loc[summary["method"].eq("minco_candidate_default_binomial_fallback")]
    sylph = summary.loc[summary["method"].eq("sylph_binomial_fallback")]
    minco_f1 = finite(minco.iloc[0]["mean_F1"]) if not minco.empty else 0.0
    sylph_f1 = finite(sylph.iloc[0]["mean_F1"]) if not sylph.empty else 0.0
    minco_l1 = finite(minco.iloc[0]["mean_L1_union_pp"]) if not minco.empty else 0.0
    sylph_l1 = finite(sylph.iloc[0]["mean_L1_union_pp"]) if not sylph.empty else 0.0
    minco_pearson = finite(minco.iloc[0]["mean_Pearson_union"]) if not minco.empty else 0.0
    sylph_pearson = finite(sylph.iloc[0]["mean_Pearson_union"]) if not sylph.empty else 0.0
    minco_delta = delta.loc[delta["method"].eq("minco_candidate_default_binomial_fallback")]
    return [
        {
            "metric": "candidate_profiles_scored",
            "value": len(CANDIDATE_PROFILES),
            "evidence": str(CANDIDATE_ROOT),
            "decision": "selected_default_profile_set",
        },
        {
            "metric": "minco_vs_sylph_mean_F1",
            "value": f"{minco_f1:.6f} vs {sylph_f1:.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "minco_higher_F1" if minco_f1 > sylph_f1 else "sylph_higher_F1",
        },
        {
            "metric": "minco_vs_sylph_mean_L1",
            "value": f"{minco_l1:.6f} vs {sylph_l1:.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "minco_lower_L1" if minco_l1 < sylph_l1 else "sylph_lower_L1",
        },
        {
            "metric": "minco_vs_sylph_mean_Pearson",
            "value": f"{minco_pearson:.6f} vs {sylph_pearson:.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "minco_higher_Pearson" if minco_pearson > sylph_pearson else "sylph_higher_Pearson",
        },
        {
            "metric": "candidate_delta_vs_prior_minco_rescore",
            "value": (
                f"mean_F1_delta={minco_delta['delta_F1'].mean():.6f};"
                f"mean_L1_delta={minco_delta['delta_L1_union_pp'].mean():.6f};"
                f"mean_Pearson_delta={minco_delta['delta_Pearson_union'].mean():.6f}"
            ),
            "evidence": str(OUT_DELTA.relative_to(EXP)),
            "decision": "candidate_default_rescore_replaces_prior_minco_rescore",
        },
        {
            "metric": "promotion_decision",
            "value": "release_evidence_candidate_if_truth_policy_accepts_binomial_fallback",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "truth_policy_is_remaining_gate",
        },
    ]


def write_markdown(scores: pd.DataFrame, summary: pd.DataFrame, delta: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    lines = [
        "# CAMI3 Binomial Fallback Candidate Default",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-only audit scores the selected `candidate` default MinCO",
        "profiles against the CAMI3 exact-binomial fallback truth table. It",
        "uses cached profile outputs and does not change MinCO calls.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Summary Scores",
            "",
            "| Method | Mean F1 | Pooled F1 | Mean L1 pp | Mean Pearson |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary.itertuples(index=False):
        lines.append(
            "| {method} | {mean_f1:.6f} | {pooled_f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                method=getattr(row, "method"),
                mean_f1=float(getattr(row, "mean_F1")),
                pooled_f1=float(getattr(row, "pooled_F1")),
                l1=float(getattr(row, "mean_L1_union_pp")),
                pearson=float(getattr(row, "mean_Pearson_union")),
            )
        )
    lines.extend(
        [
            "",
            "## Delta Versus Prior MinCO Rescore",
            "",
            "| Sample | Method | Delta F1 | Delta L1 pp | Delta Pearson |",
            "|---:|---|---:|---:|---:|",
        ]
    )
    for row in delta.itertuples(index=False):
        lines.append(
            "| {sample} | {method} | {f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                sample=getattr(row, "sample"),
                method=getattr(row, "method"),
                f1=float(getattr(row, "delta_F1")),
                l1=float(getattr(row, "delta_L1_union_pp")),
                pearson=float(getattr(row, "delta_Pearson_union")),
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- This is the CAMI3 binomial-fallback score that matches the selected",
            "  no-expertise MinCO default.",
            "- It can replace the older MinCO rescore only if the exact-binomial",
            "  fallback truth policy is accepted.",
            "- It supports MinCO on F1 for this panel, but abundance remains weaker",
            "  than Sylph.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_SCORES.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_DELTA.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    adjusted_truth = pd.read_csv(fallback.OUT_TRUTH, sep="\t")
    scores = score_candidate_default(adjusted_truth)
    summary = taxid_score.summarize(scores)
    delta = build_delta(scores)
    audit = build_audit(summary, delta)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    delta.to_csv(OUT_DELTA, sep="\t", index=False)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(scores, summary, delta, audit)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
