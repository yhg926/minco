#!/usr/bin/env python3
"""Audit exact-binomial GTDB transfer for the CAMI II marine panel.

The existing marine scorer uses only unique NCBI-taxid transfer and leaves a
large ambiguous truth fraction.  This cached-only audit tests a stricter
fallback: use unique taxid first, then exact GTDB binomial when the gold-profile
species label maps to exactly one GTDB species.  It does not rerun profilers.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd

import score_cami3_gtdb_taxid_transfer as taxid_score
import score_marine_gtdb_taxid_transfer as marine


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

OUT_TRUTH = RESULTS / "marine_binomial_transfer_truth.tsv"
OUT_QUALITY = RESULTS / "marine_binomial_transfer_quality.tsv"
OUT_DETAIL = RESULTS / "marine_binomial_transfer_detail.tsv"
OUT_AUDIT = RESULTS / "marine_binomial_transfer_audit.tsv"
OUT_MD = EXP / "MARINE_BINOMIAL_TRANSFER_AUDIT.md"


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


def normalize_name(value: object) -> str:
    return " ".join(str(value or "").strip().lower().split())


def build_gtdb_binomial_map(by_accession: dict[str, dict[str, str]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for rec in by_accession.values():
        species = str(rec.get("gtdb_species", ""))
        if not species.startswith("s__"):
            continue
        label = species[3:]
        parts = label.split()
        if len(parts) >= 2:
            out[normalize_name(" ".join(parts[:2]))].add(species)
    return out


def transfer_sample(
    sample_id: int,
    taxid_to_gtdb: dict[str, str],
    ambiguous_taxids: dict[str, list[str]],
    gtdb_binomial: dict[str, set[str]],
) -> tuple[list[dict[str, object]], dict[str, object], list[dict[str, object]]]:
    ba_truth, scope = marine.read_marine_ba_truth(sample_id)
    counts: dict[str, float] = defaultdict(float)
    source_taxids: dict[str, set[str]] = defaultdict(set)
    detail_rows: list[dict[str, object]] = []
    method_mass = defaultdict(float)
    method_taxids = defaultdict(int)

    for row in ba_truth.itertuples(index=False):
        taxid = str(row.ncbi_species_taxid)
        name = str(row.truth_name)
        abundance = float(row.abundance_pct_all)
        method = ""
        species = taxid_to_gtdb.get(taxid, "")
        if species:
            method = "unique_taxid"
        else:
            candidates = gtdb_binomial.get(normalize_name(name), set())
            if len(candidates) == 1:
                species = next(iter(candidates))
                method = "unique_gtdb_binomial"
            elif taxid in ambiguous_taxids:
                method = "ambiguous_taxid"
            else:
                method = "unmapped"
        method_mass[method] += abundance
        method_taxids[method] += 1
        detail_rows.append(
            {
                "sample": sample_id,
                "ncbi_species_taxid": taxid,
                "truth_name": name,
                "abundance_pct_all": abundance,
                "transfer_method": method,
                "gtdb_species": species,
                "taxid_gtdb_species_count": len(ambiguous_taxids.get(taxid, [])) if taxid in ambiguous_taxids else (1 if species and method == "unique_taxid" else 0),
                "binomial_gtdb_species_count": len(gtdb_binomial.get(normalize_name(name), set())),
            }
        )
        if species:
            counts[species] += abundance
            source_taxids[species].add(taxid)

    total_mapped = sum(counts.values())
    truth_rows = [
        {
            "sample": sample_id,
            "gtdb_species": species,
            "abundance_pct_all": abundance,
            "truth_abundance": abundance / total_mapped if total_mapped else 0.0,
            "source_taxids": ",".join(sorted(source_taxids[species])),
        }
        for species, abundance in sorted(counts.items())
    ]
    ba_mass = float(scope["gold_bacteria_archaea_species_mass_pct_all"])
    quality = {
        "sample": sample_id,
        **scope,
        "truth_taxids_total": int(len(ba_truth)),
        "truth_taxids_unique_taxid": int(method_taxids["unique_taxid"]),
        "truth_taxids_binomial_fallback": int(method_taxids["unique_gtdb_binomial"]),
        "truth_taxids_ambiguous_after_fallback": int(method_taxids["ambiguous_taxid"]),
        "truth_taxids_unmapped_after_fallback": int(method_taxids["unmapped"]),
        "truth_mass_unique_taxid_pct_all": method_mass["unique_taxid"],
        "truth_mass_binomial_fallback_pct_all": method_mass["unique_gtdb_binomial"],
        "truth_mass_ambiguous_after_fallback_pct_all": method_mass["ambiguous_taxid"],
        "truth_mass_unmapped_after_fallback_pct_all": method_mass["unmapped"],
        "truth_mass_mapped_pct_all": total_mapped,
        "truth_mass_mapped_pct_bacteria_archaea": 100.0 * total_mapped / ba_mass if ba_mass else 0.0,
        "release_grade_threshold_pct_bacteria_archaea": 95.0,
        "release_grade_ready_bacteria_archaea": str((100.0 * total_mapped / ba_mass if ba_mass else 0.0) >= 95.0).lower(),
        "scored_with_profiles": str(sample_id in marine.SAMPLE_PATHS).lower(),
    }
    return truth_rows, quality, detail_rows


def build_audit(quality_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    scored = [row for row in quality_rows if str(row["scored_with_profiles"]).lower() == "true"]
    all_ready = all(str(row["release_grade_ready_bacteria_archaea"]).lower() == "true" for row in quality_rows)
    scored_ready = all(str(row["release_grade_ready_bacteria_archaea"]).lower() == "true" for row in scored)
    min_all = min(finite(row["truth_mass_mapped_pct_bacteria_archaea"]) for row in quality_rows)
    min_scored = min(finite(row["truth_mass_mapped_pct_bacteria_archaea"]) for row in scored) if scored else 0.0
    total_fallback_mass = sum(finite(row["truth_mass_binomial_fallback_pct_all"]) for row in quality_rows)
    return [
        {
            "metric": "samples_audited",
            "value": len(quality_rows),
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "all_gold_samples",
        },
        {
            "metric": "profile_scored_samples",
            "value": ",".join(str(row["sample"]) for row in scored),
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "cached_profile_subset",
        },
        {
            "metric": "fallback_mapped_mass_pct_all_sum",
            "value": f"{total_fallback_mass:.6f}",
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "binomial_fallback_rescues_truth_mass",
        },
        {
            "metric": "all_samples_min_mapped_pct_bacteria_archaea",
            "value": f"{min_all:.6f}",
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "pass" if all_ready else "below_release_threshold",
        },
        {
            "metric": "scored_samples_min_mapped_pct_bacteria_archaea",
            "value": f"{min_scored:.6f}",
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "pass" if scored_ready else "below_release_threshold",
        },
        {
            "metric": "promotion_decision",
            "value": "truth_transfer_ready_profiles_not_selected_default"
            if scored_ready
            else "do_not_promote",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "need_selected_default_profiles_before_release"
            if scored_ready
            else "truth_transfer_still_partial",
        },
    ]


def write_markdown(quality_rows: list[dict[str, object]], audit: list[dict[str, object]]) -> None:
    lines = [
        "# Marine Binomial Transfer Audit",
        "",
        "Date: 2026-06-29",
        "",
        "This cached audit tests whether exact GTDB-binomial fallback can upgrade",
        "CAMI II marine GTDB truth transfer. It uses gold-profile species labels",
        "and cached GTDB metadata only; no profilers are rerun.",
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
            "## Quality",
            "",
            "| Sample | Mapped B/A % | Unique-taxid mass | Binomial fallback mass | Remaining ambiguous mass | Ready | Scored |",
            "|---:|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in quality_rows:
        lines.append(
            "| {sample} | {mapped:.3f} | {unique:.4f} | {fallback:.4f} | {ambig:.4f} | {ready} | {scored} |".format(
                sample=row["sample"],
                mapped=float(row["truth_mass_mapped_pct_bacteria_archaea"]),
                unique=float(row["truth_mass_unique_taxid_pct_all"]),
                fallback=float(row["truth_mass_binomial_fallback_pct_all"]),
                ambig=float(row["truth_mass_ambiguous_after_fallback_pct_all"]),
                ready=row["release_grade_ready_bacteria_archaea"],
                scored=row["scored_with_profiles"],
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Exact-binomial fallback is useful for marine truth transfer.",
            "- Do not promote the marine panel unless both truth coverage is high",
            "  and selected-default MinCO/Sylph profiles are available in the same",
            "  benchmark namespace.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_TRUTH.relative_to(EXP)}`",
            f"- `{OUT_QUALITY.relative_to(EXP)}`",
            f"- `{OUT_DETAIL.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    by_accession, _by_core = marine.truth.load_gtdb_metadata(marine.truth.GTDB_METADATA)
    taxid_to_gtdb, ambiguous_taxids = taxid_score.build_taxid_transfer(by_accession)
    gtdb_binomial = build_gtdb_binomial_map(by_accession)

    truth_rows: list[dict[str, object]] = []
    quality_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    for sample_id in marine.all_gold_sample_ids(marine.GOLD):
        sample_truth, quality, detail = transfer_sample(
            sample_id,
            taxid_to_gtdb,
            ambiguous_taxids,
            gtdb_binomial,
        )
        truth_rows.extend(sample_truth)
        quality_rows.append(quality)
        detail_rows.extend(detail)
    audit = build_audit(quality_rows)
    write_tsv(
        OUT_TRUTH,
        truth_rows,
        ["sample", "gtdb_species", "abundance_pct_all", "truth_abundance", "source_taxids"],
    )
    write_tsv(
        OUT_QUALITY,
        quality_rows,
        [
            "sample",
            "gold_species_mass_pct_all",
            "gold_bacteria_archaea_species_mass_pct_all",
            "truth_taxids_total",
            "truth_taxids_unique_taxid",
            "truth_taxids_binomial_fallback",
            "truth_taxids_ambiguous_after_fallback",
            "truth_taxids_unmapped_after_fallback",
            "truth_mass_unique_taxid_pct_all",
            "truth_mass_binomial_fallback_pct_all",
            "truth_mass_ambiguous_after_fallback_pct_all",
            "truth_mass_unmapped_after_fallback_pct_all",
            "truth_mass_mapped_pct_all",
            "truth_mass_mapped_pct_bacteria_archaea",
            "release_grade_threshold_pct_bacteria_archaea",
            "release_grade_ready_bacteria_archaea",
            "scored_with_profiles",
        ],
    )
    write_tsv(
        OUT_DETAIL,
        detail_rows,
        [
            "sample",
            "ncbi_species_taxid",
            "truth_name",
            "abundance_pct_all",
            "transfer_method",
            "gtdb_species",
            "taxid_gtdb_species_count",
            "binomial_gtdb_species_count",
        ],
    )
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(quality_rows, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
