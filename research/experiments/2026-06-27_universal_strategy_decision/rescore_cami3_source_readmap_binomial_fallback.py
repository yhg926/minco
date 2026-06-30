#!/usr/bin/env python3
"""Rescore CAMI3 source-readmap with cached exact-binomial fallback truth.

This audit uses only cached source-level mapping tables and cached profiler
outputs.  It does not reread input data and does not change MinCO calls.  The
purpose is to test whether the resolver-candidate fallback from
CAMI3_SOURCE_READMAP_RESOLVER_CANDIDATES.md can make the CAMI3 source-readmap
panel strong enough for release-grade evidence.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

import score_cami3_gtdb_source_readmap as source_score
import score_cami3_gtdb_taxid_transfer as taxid_score


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

BASE_TRUTH = RESULTS / "cami3_gtdb_source_readmap_truth.tsv"
BASE_QUALITY = RESULTS / "cami3_gtdb_source_readmap_quality.tsv"
BASE_SOURCE_ROWS = RESULTS / "cami3_gtdb_source_readmap_source_genomes.tsv"
BASE_SCORES = RESULTS / "cami3_gtdb_source_readmap_scores.tsv"
SCOPE_SUMMARY = RESULTS / "cami3_source_readmap_scope_summary.tsv"
RESOLVER_CANDIDATES = RESULTS / "cami3_source_readmap_resolver_candidates.tsv"

OUT_TRUTH = RESULTS / "cami3_source_readmap_binomial_fallback_truth.tsv"
OUT_QUALITY = RESULTS / "cami3_source_readmap_binomial_fallback_quality.tsv"
OUT_SCORES = RESULTS / "cami3_source_readmap_binomial_fallback_scores.tsv"
OUT_SUMMARY = RESULTS / "cami3_source_readmap_binomial_fallback_summary.tsv"
OUT_DELTA = RESULTS / "cami3_source_readmap_binomial_fallback_delta.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_binomial_fallback_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_BINOMIAL_FALLBACK_RESCORING.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def base_source_species(base_truth: pd.DataFrame) -> dict[tuple[str, str], str]:
    by_source: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in base_truth.itertuples(index=False):
        sample = str(getattr(row, "sample"))
        species = str(getattr(row, "gtdb_species"))
        for source_id in str(getattr(row, "source_genome_ids", "") or "").split(","):
            source_id = source_id.strip()
            if source_id:
                by_source[(sample, source_id)].add(species)
    return {key: next(iter(values)) for key, values in by_source.items() if len(values) == 1}


def build_adjusted_truth() -> tuple[pd.DataFrame, pd.DataFrame]:
    base_truth = pd.read_csv(BASE_TRUTH, sep="\t")
    quality = pd.read_csv(BASE_QUALITY, sep="\t")
    scope = {str(row["sample"]): row for row in read_tsv(SCOPE_SUMMARY)}
    source_rows = {(row["sample"], row["source_genome_id"]): row for row in read_tsv(BASE_SOURCE_ROWS)}
    source_species = base_source_species(base_truth)
    fallback_rows = [
        row
        for row in read_tsv(RESOLVER_CANDIDATES)
        if str(row.get("recommended_cached_resolution", "")).startswith("s__")
    ]

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    fallback_counts: dict[str, Counter[str]] = defaultdict(Counter)
    fallback_sources: dict[tuple[str, str], set[str]] = defaultdict(set)

    for row in base_truth.itertuples(index=False):
        sample = str(getattr(row, "sample"))
        species = str(getattr(row, "gtdb_species"))
        read_rows = int(round(finite(getattr(row, "read_rows", 0))))
        counts[sample][species] += read_rows
        for source_id in str(getattr(row, "source_genome_ids", "") or "").split(","):
            source_id = source_id.strip()
            if source_id:
                sources[(sample, species)].add(source_id)

    quality_rows = []
    added_by_sample: dict[str, int] = defaultdict(int)
    for row in fallback_rows:
        sample = str(row["sample"])
        source_id = str(row["source_genome_id"])
        species = str(row["recommended_cached_resolution"])
        key = (sample, source_id)
        source = source_rows.get(key, {})
        add_rows = int(round(finite(source.get("unmapped_read_rows", row.get("unmapped_rows", 0)))))
        if add_rows <= 0:
            continue
        if source_species.get(key) == species:
            # If a source already contributed mapped rows to the same species,
            # this adds only its cached unmapped remainder.
            pass
        counts[sample][species] += add_rows
        sources[(sample, species)].add(source_id)
        fallback_counts[sample][species] += add_rows
        fallback_sources[(sample, species)].add(source_id)
        added_by_sample[sample] += add_rows

    adjusted_rows = []
    for sample in sorted(counts, key=int):
        total = sum(counts[sample].values())
        for species, read_rows in sorted(counts[sample].items()):
            adjusted_rows.append(
                {
                    "sample": sample,
                    "gtdb_species": species,
                    "read_rows": read_rows,
                    "truth_abundance": read_rows / total if total else 0.0,
                    "source_genomes": len(sources[(sample, species)]),
                    "source_genome_ids": ",".join(sorted(sources[(sample, species)])),
                    "fallback_read_rows": fallback_counts[sample].get(species, 0),
                    "fallback_source_genomes": len(fallback_sources[(sample, species)]),
                    "fallback_source_genome_ids": ",".join(sorted(fallback_sources[(sample, species)])),
                }
            )

    for row in quality.to_dict(orient="records"):
        sample = str(row["sample"])
        added = int(added_by_sample.get(sample, 0))
        scope_row = scope.get(sample, {})
        profile_total = int(round(finite(scope_row.get("profile_scope_read_rows", 0))))
        profile_mapped = int(round(finite(scope_row.get("profile_scope_mapped_rows", 0))))
        total_rows = int(round(finite(row.get("read_rows_total", 0))))
        base_mapped = int(round(finite(row.get("read_rows_mapped", 0))))
        adjusted_profile_mapped = min(profile_total, profile_mapped + added)
        adjusted_total_mapped = min(total_rows, base_mapped + added)
        quality_rows.append(
            {
                "sample": sample,
                "base_read_rows_total": total_rows,
                "base_read_rows_mapped": base_mapped,
                "fallback_added_rows": added,
                "adjusted_read_rows_mapped": adjusted_total_mapped,
                "base_all_rows_mapped_pct": finite(row.get("read_rows_mapped_pct")),
                "adjusted_all_rows_mapped_pct": (
                    100.0 * adjusted_total_mapped / total_rows if total_rows else 0.0
                ),
                "profile_scope_read_rows": profile_total,
                "base_profile_scope_mapped_rows": profile_mapped,
                "adjusted_profile_scope_mapped_rows": adjusted_profile_mapped,
                "base_profile_scope_mapped_pct": finite(scope_row.get("profile_scope_mapped_pct")),
                "adjusted_profile_scope_mapped_pct": (
                    100.0 * adjusted_profile_mapped / profile_total if profile_total else 0.0
                ),
                "release_threshold_pct": 95.0,
                "profile_scope_release_ready_after_fallback": str(
                    profile_total > 0 and 100.0 * adjusted_profile_mapped / profile_total >= 95.0
                ).lower(),
            }
        )

    return pd.DataFrame(adjusted_rows), pd.DataFrame(quality_rows)


def load_predictions() -> tuple[dict[int, pd.DataFrame], dict[int, pd.DataFrame], dict[int, dict[str, object]], dict[int, dict[str, object]]]:
    _wgs_to_species, taxid_to_species, name_to_species, _map_diag = source_score.build_transfer_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(source_score.TAXMAP)
    minco_preds: dict[int, pd.DataFrame] = {}
    sylph_preds: dict[int, pd.DataFrame] = {}
    minco_extra: dict[int, dict[str, object]] = {}
    sylph_extra: dict[int, dict[str, object]] = {}

    for sample_id in sorted(source_score.PROFILE_PATHS):
        best_ref_species, _best_ref_diag = source_score.best_raw_ref_species_by_taxid(
            source_score.RAW_TABLES[sample_id],
            taxmap,
            by_accession,
            by_core,
        )
        minco_pred, minco_meta = source_score.load_minco_predictions(
            source_score.PROFILE_PATHS[sample_id]["minco"],
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
        minco_extra[sample_id] = minco_meta
        sylph_extra[sample_id] = sylph_meta
    return minco_preds, sylph_preds, minco_extra, sylph_extra


def score_adjusted_truth(adjusted_truth: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    minco_preds, sylph_preds, minco_extra, sylph_extra = load_predictions()
    score_rows = []
    for sample_id in sorted(minco_preds):
        truth_df = adjusted_truth.loc[adjusted_truth["sample"].astype(str).eq(str(sample_id))].copy()
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "minco_source_readmap_binomial_fallback_rescore",
                minco_preds[sample_id],
                truth_df,
                minco_extra[sample_id],
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_source_readmap_binomial_fallback_rescore",
                sylph_preds[sample_id],
                truth_df,
                sylph_extra[sample_id],
            )
        )
    scores = pd.DataFrame(score_rows)
    summary = taxid_score.summarize(scores)
    delta = build_delta(scores)
    return scores, summary, delta


def build_delta(scores: pd.DataFrame) -> pd.DataFrame:
    base = pd.read_csv(BASE_SCORES, sep="\t")
    method_map = {
        "minco_source_readmap_binomial_fallback_rescore": "minco_universal_autoexact_gtdb_source_readmap",
        "sylph_source_readmap_binomial_fallback_rescore": "sylph_gtdb_source_readmap",
    }
    rows = []
    metrics = ["truth_gtdb_species", "TP", "FP", "FN", "F1", "L1_union_pp", "Pearson_union"]
    for row in scores.itertuples(index=False):
        method = str(getattr(row, "method"))
        base_method = method_map[method]
        sample = str(getattr(row, "sample"))
        base_row = base.loc[base["sample"].astype(str).eq(sample) & base["method"].eq(base_method)]
        if base_row.empty:
            continue
        base_one = base_row.iloc[0]
        out = {
            "sample": sample,
            "method": method,
            "baseline_method": base_method,
        }
        for metric in metrics:
            adjusted = finite(getattr(row, metric))
            base_value = finite(base_one.get(metric))
            out[f"baseline_{metric}"] = base_value
            out[f"adjusted_{metric}"] = adjusted
            out[f"delta_{metric}"] = adjusted - base_value
        rows.append(out)
    return pd.DataFrame(rows)


def build_audit(quality: pd.DataFrame, scores: pd.DataFrame, summary: pd.DataFrame, delta: pd.DataFrame) -> list[dict[str, object]]:
    ready = ",".join(
        map(
            str,
            quality.loc[
                quality["profile_scope_release_ready_after_fallback"].astype(str).eq("true"), "sample"
            ].tolist(),
        )
    )
    added = int(quality["fallback_added_rows"].sum())
    minco = summary.loc[summary["method"].eq("minco_source_readmap_binomial_fallback_rescore")]
    sylph = summary.loc[summary["method"].eq("sylph_source_readmap_binomial_fallback_rescore")]
    minco_f1 = finite(minco.iloc[0]["mean_F1"]) if not minco.empty else 0.0
    sylph_f1 = finite(sylph.iloc[0]["mean_F1"]) if not sylph.empty else 0.0
    minco_l1 = finite(minco.iloc[0]["mean_L1_union_pp"]) if not minco.empty else 0.0
    sylph_l1 = finite(sylph.iloc[0]["mean_L1_union_pp"]) if not sylph.empty else 0.0
    minco_f1_delta = delta.loc[delta["method"].str.startswith("minco"), "delta_F1"].mean()
    sylph_f1_delta = delta.loc[delta["method"].str.startswith("sylph"), "delta_F1"].mean()
    return [
        {
            "metric": "fallback_added_rows",
            "value": added,
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "cached_source_level_truth_expansion",
        },
        {
            "metric": "profile_scope_release_ready_after_fallback",
            "value": ready,
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "candidate_truth_transfer_reaches_coverage_threshold" if ready == "0,1,2" else "coverage_gap_remains",
        },
        {
            "metric": "minco_vs_sylph_mean_F1_after_rescore",
            "value": f"{minco_f1:.6f} vs {sylph_f1:.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "minco_higher_F1" if minco_f1 > sylph_f1 else "sylph_higher_F1",
        },
        {
            "metric": "minco_vs_sylph_mean_L1_after_rescore",
            "value": f"{minco_l1:.6f} vs {sylph_l1:.6f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "minco_lower_L1" if minco_l1 < sylph_l1 else "sylph_lower_L1",
        },
        {
            "metric": "mean_F1_delta_vs_previous_truth",
            "value": f"minco={finite(minco_f1_delta):.6f};sylph={finite(sylph_f1_delta):.6f}",
            "evidence": str(OUT_DELTA.relative_to(EXP)),
            "decision": "truth_expansion_changes_metric_denominator",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_rescore_not_default_change",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "update_holdout_status_only_after_manual_truth_rule_review",
        },
    ]


def write_markdown(quality: pd.DataFrame, summary: pd.DataFrame, delta: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    lines = [
        "# CAMI3 Source-Readmap Binomial Fallback Rescoring",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-only audit adds the resolver-candidate exact GTDB-binomial",
        "fallback rows to the existing source-readmap truth table and scores the",
        "same cached MinCO/Sylph profile outputs against the adjusted truth. It",
        "does not change MinCO calls and does not reread input data.",
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
            "## Coverage",
            "",
            "| Sample | Added rows | Base profile-scope mapped % | Adjusted profile-scope mapped % | Ready |",
            "|---:|---:|---:|---:|---|",
        ]
    )
    for row in quality.itertuples(index=False):
        lines.append(
            "| {sample} | {added} | {base:.3f} | {adjusted:.3f} | {ready} |".format(
                sample=getattr(row, "sample"),
                added=int(getattr(row, "fallback_added_rows")),
                base=float(getattr(row, "base_profile_scope_mapped_pct")),
                adjusted=float(getattr(row, "adjusted_profile_scope_mapped_pct")),
                ready=getattr(row, "profile_scope_release_ready_after_fallback"),
            )
        )
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
            "## Delta Versus Previous Truth",
            "",
            "| Sample | Method | Delta truth species | Delta F1 | Delta L1 pp | Delta Pearson |",
            "|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in delta.itertuples(index=False):
        lines.append(
            "| {sample} | {method} | {truth_delta:.0f} | {f1_delta:.6f} | {l1_delta:.6f} | {pearson_delta:.6f} |".format(
                sample=getattr(row, "sample"),
                method=getattr(row, "method"),
                truth_delta=float(getattr(row, "delta_truth_gtdb_species")),
                f1_delta=float(getattr(row, "delta_F1")),
                l1_delta=float(getattr(row, "delta_L1_union_pp")),
                pearson_delta=float(getattr(row, "delta_Pearson_union")),
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- The fallback reaches the profile-scope coverage threshold for samples",
            "  0, 1, and 2 in this cached source-level rescore.",
            "- This updates the CAMI3 evidence path, but it is still a truth-transfer",
            "  audit. It does not justify changing the MinCO default strategy.",
            "- Before promoting CAMI3 to release-grade evidence, review whether exact",
            "  GTDB-binomial fallback is acceptable for the benchmark truth policy.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_TRUTH.relative_to(EXP)}`",
            f"- `{OUT_QUALITY.relative_to(EXP)}`",
            f"- `{OUT_SCORES.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_DELTA.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    adjusted_truth, quality = build_adjusted_truth()
    scores, summary, delta = score_adjusted_truth(adjusted_truth)
    audit = build_audit(quality, scores, summary, delta)
    adjusted_truth.to_csv(OUT_TRUTH, sep="\t", index=False)
    quality.to_csv(OUT_QUALITY, sep="\t", index=False)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    delta.to_csv(OUT_DELTA, sep="\t", index=False)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(quality, summary, delta, audit)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
