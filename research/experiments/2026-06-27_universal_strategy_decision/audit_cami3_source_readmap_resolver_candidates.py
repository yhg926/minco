#!/usr/bin/env python3
"""Audit cached CAMI3 source-readmap mapping fallbacks.

The source-readmap truth transfer currently uses source contig WGS prefixes,
then unique NCBI taxid mappings.  This cached-only audit checks whether the
remaining in-scope source rows can be resolved by deterministic metadata
fallbacks without touching raw reads or rerunning profilers.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import score_cami3_gtdb_source_readmap as source_score


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
TOP_IN_SCOPE = RESULTS / "cami3_source_readmap_top_unmapped_in_scope.tsv"
SCOPE_SUMMARY = RESULTS / "cami3_source_readmap_scope_summary.tsv"
OUT_DETAIL = RESULTS / "cami3_source_readmap_resolver_candidates.tsv"
OUT_SUMMARY = RESULTS / "cami3_source_readmap_resolver_candidate_summary.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_resolver_candidate_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_RESOLVER_CANDIDATES.md"

TOP_N_PER_SAMPLE = 20


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def build_candidate_maps() -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]]]:
    taxid_to_species: dict[str, set[str]] = defaultdict(set)
    ncbi_name_to_species: dict[str, set[str]] = defaultdict(set)
    gtdb_binomial_to_species: dict[str, set[str]] = defaultdict(set)

    for path in source_score.truth.GTDB_METADATA:
        with source_score.truth.open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                species = source_score.gtdb_species(row.get("gtdb_taxonomy", ""))
                if not species:
                    continue
                for key in ["ncbi_species_taxid", "ncbi_taxid"]:
                    taxid = str(row.get(key, "")).strip()
                    if taxid:
                        taxid_to_species[taxid].add(species)
                for name in [
                    source_score.ncbi_species_from_taxonomy(row.get("ncbi_taxonomy", "")),
                    row.get("ncbi_organism_name", ""),
                ]:
                    normalized = source_score.normalize_name(name)
                    if normalized:
                        ncbi_name_to_species[normalized].add(species)
                gtdb_label = species[3:] if species.startswith("s__") else species
                gtdb_parts = gtdb_label.split()
                if len(gtdb_parts) >= 2:
                    gtdb_binomial_to_species[source_score.normalize_name(" ".join(gtdb_parts[:2]))].add(species)

    return taxid_to_species, ncbi_name_to_species, gtdb_binomial_to_species


def sample_top_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_sample[str(row.get("sample", ""))].append(row)
    out: list[dict[str, str]] = []
    for sample in sorted(by_sample, key=lambda value: int(value)):
        sorted_rows = sorted(
            by_sample[sample],
            key=lambda row: -int(float(row.get("unmapped_rows", "0") or 0)),
        )
        out.extend(sorted_rows[:TOP_N_PER_SAMPLE])
    return out


def one_or_empty(values: set[str]) -> str:
    return next(iter(values)) if len(values) == 1 else ""


def summarize(detail_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_sample: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    scope_by_sample = {str(row["sample"]): row for row in read_tsv(SCOPE_SUMMARY)}
    for row in detail_rows:
        sample = str(row["sample"])
        rows = int(row["unmapped_rows"])
        by_sample[sample]["reviewed_sources"] += 1
        by_sample[sample]["reviewed_unmapped_rows"] += rows
        if row["unique_taxid_species"]:
            by_sample[sample]["unique_taxid_resolvable_rows"] += rows
        if row["unique_ncbi_name_species"]:
            by_sample[sample]["unique_ncbi_name_resolvable_rows"] += rows
        if row["unique_gtdb_binomial_species"]:
            by_sample[sample]["unique_gtdb_binomial_resolvable_rows"] += rows
        if row["recommended_cached_resolution"]:
            by_sample[sample]["recommended_resolvable_rows"] += rows

    summary_rows: list[dict[str, object]] = []
    for sample in sorted(by_sample, key=int):
        data = by_sample[sample]
        scope = scope_by_sample.get(sample, {})
        scope_total = int(float(scope.get("profile_scope_read_rows", "0") or 0))
        scope_mapped = int(float(scope.get("profile_scope_mapped_rows", "0") or 0))
        current_pct = 100.0 * scope_mapped / scope_total if scope_total else 0.0
        projected_mapped = min(scope_total, scope_mapped + data["recommended_resolvable_rows"])
        projected_pct = 100.0 * projected_mapped / scope_total if scope_total else 0.0
        total = data["reviewed_unmapped_rows"]
        summary_rows.append(
            {
                "sample": sample,
                "current_profile_scope_mapped_pct": current_pct,
                "additional_rows_needed_for_95pct": int(
                    float(scope.get("additional_profile_scope_rows_needed_for_95pct", "0") or 0)
                ),
                "reviewed_sources": data["reviewed_sources"],
                "reviewed_unmapped_rows": total,
                "unique_taxid_resolvable_rows": data["unique_taxid_resolvable_rows"],
                "unique_ncbi_name_resolvable_rows": data["unique_ncbi_name_resolvable_rows"],
                "unique_gtdb_binomial_resolvable_rows": data["unique_gtdb_binomial_resolvable_rows"],
                "recommended_resolvable_rows": data["recommended_resolvable_rows"],
                "recommended_resolvable_pct_of_reviewed": (
                    100.0 * data["recommended_resolvable_rows"] / total if total else 0.0
                ),
                "projected_profile_scope_mapped_pct_after_recommended": projected_pct,
                "projected_profile_scope_release_ready": str(projected_pct >= 95.0).lower(),
            }
        )

    total_reviewed = sum(int(row["reviewed_unmapped_rows"]) for row in summary_rows)
    total_recommended = sum(int(row["recommended_resolvable_rows"]) for row in summary_rows)
    audit_rows = [
        {
            "metric": "reviewed_top_sources_per_sample",
            "value": TOP_N_PER_SAMPLE,
            "evidence": str(TOP_IN_SCOPE.relative_to(EXP)),
            "decision": "cached_only_no_profile_rerun",
        },
        {
            "metric": "reviewed_unmapped_rows",
            "value": total_reviewed,
            "evidence": str(OUT_DETAIL.relative_to(EXP)),
            "decision": "top_in_scope_source_rows",
        },
        {
            "metric": "recommended_resolvable_rows",
            "value": total_recommended,
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "candidate_truth_transfer_fallback",
        },
        {
            "metric": "recommended_resolvable_pct_of_reviewed",
            "value": f"{(100.0 * total_recommended / total_reviewed if total_reviewed else 0.0):.3f}",
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "large_enough_to_test_release_upgrade",
        },
        {
            "metric": "projected_samples_release_ready_after_candidate",
            "value": ",".join(
                str(row["sample"])
                for row in summary_rows
                if str(row["projected_profile_scope_release_ready"]).lower() == "true"
            ),
            "evidence": str(OUT_SUMMARY.relative_to(EXP)),
            "decision": "rescored_truth_test_is_worth_running",
        },
        {
            "metric": "promotion_decision",
            "value": "do_not_promote_without_rescored_truth_audit",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "exact_binomial_fallback_is_candidate_not_current_truth",
        },
    ]
    return summary_rows, audit_rows


def write_markdown(
    detail_rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
    audit_rows: list[dict[str, object]],
) -> None:
    lines = [
        "# CAMI3 Source-Readmap Resolver Candidates",
        "",
        "Date: 2026-06-29",
        "",
        "This cached-only audit checks deterministic metadata fallbacks for the",
        "top in-scope CAMI3 source-readmap rows that remain unmapped after the",
        "current WGS-prefix and unique-taxid transfer. It does not use raw",
        "input data and does not rerun MinCO or Sylph.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit_rows:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Sample | Reviewed rows | Unique taxid rows | Unique NCBI-name rows | Unique GTDB-binomial rows | Recommended rows |",
            "|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in summary_rows:
        lines.append(
            "| {sample} | {reviewed} | {taxid} | {name} | {binomial} | {recommended} |".format(
                sample=row["sample"],
                reviewed=row["reviewed_unmapped_rows"],
                taxid=row["unique_taxid_resolvable_rows"],
                name=row["unique_ncbi_name_resolvable_rows"],
                binomial=row["unique_gtdb_binomial_resolvable_rows"],
                recommended=row["recommended_resolvable_rows"],
            )
        )

    lines.extend(
        [
            "",
            "## Release-Coverage Projection",
            "",
            "| Sample | Current mapped % | Rows needed for 95% | Recommended rows | Projected mapped % | Projected ready |",
            "|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in summary_rows:
        lines.append(
            "| {sample} | {current:.3f} | {needed} | {recommended} | {projected:.3f} | {ready} |".format(
                sample=row["sample"],
                current=float(row["current_profile_scope_mapped_pct"]),
                needed=row["additional_rows_needed_for_95pct"],
                recommended=row["recommended_resolvable_rows"],
                projected=float(row["projected_profile_scope_mapped_pct_after_recommended"]),
                ready=row["projected_profile_scope_release_ready"],
            )
        )

    lines.extend(
        [
            "",
            "## Top Candidate Resolutions",
            "",
            "| Sample | Source | Species label | Unmapped rows | Recommended GTDB species | Rule | Ambiguity note |",
            "|---:|---|---|---:|---|---|---|",
        ]
    )
    for row in detail_rows[:30]:
        lines.append(
            "| {sample} | {source} | {label} | {rows} | {species} | {rule} | {note} |".format(
                sample=row["sample"],
                source=row["source_genome_id"],
                label=str(row["species_label"]).replace("|", "/"),
                rows=row["unmapped_rows"],
                species=row["recommended_cached_resolution"],
                rule=row["recommended_resolution_rule"],
                note=row["ambiguity_note"],
            )
        )

    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Do not promote CAMI3 source-readmap to release-grade evidence yet.",
            "- The exact GTDB-binomial fallback is promising enough to test in the",
            "  scorer because it can resolve many top in-scope rows, but it needs",
            "  a rescored truth audit before it becomes evidence.",
            "- This is a metadata/truth-transfer task, not a MinCO threshold task.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_DETAIL.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    taxid_map, ncbi_name_map, gtdb_binomial_map = build_candidate_maps()
    top_rows = sample_top_rows(read_tsv(TOP_IN_SCOPE))
    detail_rows: list[dict[str, object]] = []

    for row in top_rows:
        taxid = str(row.get("taxid", "")).strip()
        label = str(row.get("species_label", "")).strip()
        normalized_label = source_score.normalize_name(label)
        taxid_candidates = taxid_map.get(taxid, set())
        ncbi_name_candidates = ncbi_name_map.get(normalized_label, set())
        gtdb_binomial_candidates = gtdb_binomial_map.get(normalized_label, set())
        unique_taxid = one_or_empty(taxid_candidates)
        unique_ncbi_name = one_or_empty(ncbi_name_candidates)
        unique_gtdb_binomial = one_or_empty(gtdb_binomial_candidates)
        recommended = unique_taxid or unique_ncbi_name or unique_gtdb_binomial
        if unique_taxid:
            rule = "unique_taxid"
        elif unique_ncbi_name:
            rule = "unique_ncbi_name"
        elif unique_gtdb_binomial:
            rule = "unique_gtdb_binomial"
        else:
            rule = ""
        ambiguity_parts = []
        if len(taxid_candidates) > 1:
            ambiguity_parts.append(f"taxid_candidates={len(taxid_candidates)}")
        if len(ncbi_name_candidates) > 1:
            ambiguity_parts.append(f"ncbi_name_candidates={len(ncbi_name_candidates)}")
        if len(gtdb_binomial_candidates) > 1:
            ambiguity_parts.append(f"gtdb_binomial_candidates={len(gtdb_binomial_candidates)}")
        detail_rows.append(
            {
                "sample": row.get("sample", ""),
                "source_genome_id": row.get("source_genome_id", ""),
                "taxid": taxid,
                "species_label": label,
                "unmapped_rows": int(float(row.get("unmapped_rows", "0") or 0)),
                "taxid_gtdb_species_count": len(taxid_candidates),
                "ncbi_name_gtdb_species_count": len(ncbi_name_candidates),
                "gtdb_binomial_species_count": len(gtdb_binomial_candidates),
                "unique_taxid_species": unique_taxid,
                "unique_ncbi_name_species": unique_ncbi_name,
                "unique_gtdb_binomial_species": unique_gtdb_binomial,
                "recommended_cached_resolution": recommended,
                "recommended_resolution_rule": rule,
                "ambiguity_note": ";".join(ambiguity_parts),
            }
        )

    summary_rows, audit_rows = summarize(detail_rows)
    detail_fields = [
        "sample",
        "source_genome_id",
        "taxid",
        "species_label",
        "unmapped_rows",
        "taxid_gtdb_species_count",
        "ncbi_name_gtdb_species_count",
        "gtdb_binomial_species_count",
        "unique_taxid_species",
        "unique_ncbi_name_species",
        "unique_gtdb_binomial_species",
        "recommended_cached_resolution",
        "recommended_resolution_rule",
        "ambiguity_note",
    ]
    summary_fields = [
        "sample",
        "current_profile_scope_mapped_pct",
        "additional_rows_needed_for_95pct",
        "reviewed_sources",
        "reviewed_unmapped_rows",
        "unique_taxid_resolvable_rows",
        "unique_ncbi_name_resolvable_rows",
        "unique_gtdb_binomial_resolvable_rows",
        "recommended_resolvable_rows",
        "recommended_resolvable_pct_of_reviewed",
        "projected_profile_scope_mapped_pct_after_recommended",
        "projected_profile_scope_release_ready",
    ]
    audit_fields = ["metric", "value", "evidence", "decision"]
    write_tsv(OUT_DETAIL, detail_rows, detail_fields)
    write_tsv(OUT_SUMMARY, summary_rows, summary_fields)
    write_tsv(OUT_AUDIT, audit_rows, audit_fields)
    write_markdown(detail_rows, summary_rows, audit_rows)
    for row in summary_rows:
        print(
            "{sample}\t{recommended_resolvable_rows}\t{recommended_resolvable_pct_of_reviewed:.3f}".format(
                **row
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
