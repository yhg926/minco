#!/usr/bin/env python3
"""Audit whether CAMI3 exact-binomial fallback can be release truth.

This policy audit is intentionally separate from scoring.  It checks whether
the exact GTDB-binomial fallback is deterministic and conservative enough to
use as clean GTDB-species truth for the CAMI3 source-readmap holdout.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

RESOLVER = RESULTS / "cami3_source_readmap_resolver_candidates.tsv"
SOURCE_ROWS = RESULTS / "cami3_gtdb_source_readmap_source_genomes.tsv"
BASE_TRUTH = RESULTS / "cami3_gtdb_source_readmap_truth.tsv"
QUALITY = RESULTS / "cami3_source_readmap_binomial_fallback_quality.tsv"
CANDIDATE_AUDIT = RESULTS / "cami3_binomial_fallback_candidate_default_audit.tsv"

OUT_DETAIL = RESULTS / "cami3_binomial_fallback_truth_policy_detail.tsv"
OUT_AUDIT = RESULTS / "cami3_binomial_fallback_truth_policy_audit.tsv"
OUT_MD = EXP / "CAMI3_BINOMIAL_FALLBACK_TRUTH_POLICY.md"


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


def table_by(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_text(value: bool) -> str:
    return str(bool(value)).lower()


def build_detail() -> list[dict[str, object]]:
    source_by_key = {(row["sample"], row["source_genome_id"]): row for row in read_tsv(SOURCE_ROWS)}
    prior_species_by_key: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in read_tsv(BASE_TRUTH):
        for source_id in str(row.get("source_genome_ids", "")).split(","):
            source_id = source_id.strip()
            if source_id:
                prior_species_by_key[(row["sample"], source_id)].add(row["gtdb_species"])
    detail = []
    for row in read_tsv(RESOLVER):
        recommended = str(row.get("recommended_cached_resolution", ""))
        if not recommended.startswith("s__"):
            continue
        source = source_by_key.get((row["sample"], row["source_genome_id"]), {})
        prior_mapped = int(finite(source.get("wgs_prefix_read_rows"))) + int(
            finite(source.get("unique_taxid_read_rows"))
        )
        prior_species = sorted(prior_species_by_key.get((row["sample"], row["source_genome_id"]), set()))
        added = int(finite(source.get("unmapped_read_rows", row.get("unmapped_rows", 0))))
        same_species_remainder = prior_mapped > 0 and prior_species == [recommended]
        priority_ok = (prior_mapped == 0 and added > 0) or same_species_remainder
        detail.append(
            {
                "sample": row["sample"],
                "source_genome_id": row["source_genome_id"],
                "species_label": row["species_label"],
                "recommended_gtdb_species": recommended,
                "recommended_resolution_rule": row.get("recommended_resolution_rule", ""),
                "gtdb_binomial_species_count": int(finite(row.get("gtdb_binomial_species_count"))),
                "taxid_gtdb_species_count": int(finite(row.get("taxid_gtdb_species_count"))),
                "ncbi_name_gtdb_species_count": int(finite(row.get("ncbi_name_gtdb_species_count"))),
                "prior_wgs_or_unique_taxid_rows": prior_mapped,
                "prior_mapped_gtdb_species": ",".join(prior_species),
                "fallback_added_rows": added,
                "source_priority_case": "unmapped_source" if prior_mapped == 0 else "same_species_remainder"
                if same_species_remainder
                else "priority_conflict",
                "source_priority_ok": bool_text(priority_ok),
            }
        )
    return detail


def build_audit(detail: list[dict[str, object]]) -> list[dict[str, object]]:
    quality = read_tsv(QUALITY)
    candidate = table_by(CANDIDATE_AUDIT, "metric")
    fallback_rows = len(detail)
    added_rows = sum(int(row["fallback_added_rows"]) for row in detail)
    unique_rule_rows = sum(
        1
        for row in detail
        if row["recommended_resolution_rule"] == "unique_gtdb_binomial"
        and int(row["gtdb_binomial_species_count"]) == 1
    )
    priority_conflicts = sum(1 for row in detail if str(row["source_priority_ok"]) != "true")
    same_species_remainders = sum(1 for row in detail if row["source_priority_case"] == "same_species_remainder")
    ready_samples = [
        row["sample"]
        for row in quality
        if str(row.get("profile_scope_release_ready_after_fallback", "")).lower() == "true"
    ]
    min_ready_pct = min(finite(row.get("adjusted_profile_scope_mapped_pct")) for row in quality)
    candidate_profile_decision = candidate.get("promotion_decision", {}).get("decision", "")
    selected_default_scored = (
        candidate.get("candidate_profiles_scored", {}).get("value") == "3"
        and candidate_profile_decision == "truth_policy_is_remaining_gate"
    )
    deterministic = fallback_rows > 0 and unique_rule_rows == fallback_rows
    conservative_priority = priority_conflicts == 0
    coverage_ready = len(ready_samples) == 3 and min_ready_pct >= 95.0
    accepted = deterministic and conservative_priority and coverage_ready and selected_default_scored
    rows = [
        {
            "metric": "fallback_rows_reviewed",
            "value": fallback_rows,
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "nonzero_fallback_set" if fallback_rows else "no_fallback_rows",
        },
        {
            "metric": "fallback_added_rows",
            "value": added_rows,
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "truth_rows_added_by_policy",
        },
        {
            "metric": "deterministic_unique_binomial",
            "value": f"{unique_rule_rows}/{fallback_rows}",
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "pass" if deterministic else "fail",
        },
        {
            "metric": "source_mapping_priority",
            "value": (
                f"priority_conflicts={priority_conflicts};"
                f"same_species_remainders={same_species_remainders}"
            ),
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "pass" if conservative_priority else "fail",
        },
        {
            "metric": "profile_scope_coverage",
            "value": f"ready_samples={','.join(ready_samples)};min_adjusted_pct={min_ready_pct:.3f}",
            "evidence": str(QUALITY.relative_to(EXP)),
            "decision": "pass" if coverage_ready else "fail",
        },
        {
            "metric": "selected_default_scored",
            "value": candidate.get("candidate_profiles_scored", {}).get("value", "NA"),
            "evidence": str(CANDIDATE_AUDIT.relative_to(EXP)),
            "decision": "pass" if selected_default_scored else "fail",
        },
        {
            "metric": "truth_policy_decision",
            "value": "accept_exact_binomial_fallback_for_release_truth" if accepted else "do_not_accept",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "release_truth_policy_pass" if accepted else "release_truth_policy_fail",
        },
    ]
    return rows


def write_markdown(detail: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# CAMI3 Binomial Fallback Truth Policy",
        "",
        "Date: 2026-06-29",
        "",
        "This cached audit decides whether the exact GTDB-binomial fallback is",
        "acceptable as CAMI3 source-readmap GTDB-species truth. It does not rerun",
        "profiling and does not change MinCO calls.",
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
            "## Detail Summary",
            "",
            "| Fallback rows | Added rows | Unique-binomial rows | Prior-mapped violations |",
            "|---:|---:|---:|---:|",
            "| {rows} | {added} | {unique} | {violations} |".format(
                rows=audit_by_metric["fallback_rows_reviewed"]["value"],
                added=audit_by_metric["fallback_added_rows"]["value"],
                unique=str(audit_by_metric["deterministic_unique_binomial"]["value"]).split("/")[0],
                violations=str(audit_by_metric["source_mapping_priority"]["value"]).split(";")[0].split("=")[-1],
            ),
            "",
            "## Decision",
            "",
            "- Accept the exact GTDB-binomial fallback as clean CAMI3",
            "  source-readmap GTDB-species truth if all audit checks pass.",
            "- This policy fills rows left unmapped by the stronger WGS-prefix",
            "  and unique-taxid transfer routes. If a source already mapped by",
            "  WGS/unique-taxid, the fallback is allowed only when it assigns the",
            "  unmapped remainder to that same GTDB species.",
            "- The policy does not imply MinCO abundance beats Sylph; it only decides",
            "  whether the CAMI3 panel is clean enough to use as release evidence.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_DETAIL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    detail = build_detail()
    audit = build_audit(detail)
    write_tsv(
        OUT_DETAIL,
        detail,
        [
            "sample",
            "source_genome_id",
            "species_label",
            "recommended_gtdb_species",
            "recommended_resolution_rule",
            "gtdb_binomial_species_count",
            "taxid_gtdb_species_count",
            "ncbi_name_gtdb_species_count",
            "prior_wgs_or_unique_taxid_rows",
            "prior_mapped_gtdb_species",
            "fallback_added_rows",
            "source_priority_case",
            "source_priority_ok",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(detail, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
