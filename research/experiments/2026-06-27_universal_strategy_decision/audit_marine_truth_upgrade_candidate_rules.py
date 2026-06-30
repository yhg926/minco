#!/usr/bin/env python3
"""Audit candidate marine truth-upgrade rules without changing release truth.

The accepted marine transfer already uses unique NCBI taxid, then exact unique
GTDB-binomial fallback. This script tests whether a weaker ambiguous-taxid
epithet rule could clear the 95% mapped-truth threshold. It is diagnostic only:
an NCBI species split across multiple GTDB species is still not release-grade
unless the GTDB species assignment is deterministic.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

import score_cami3_gtdb_taxid_transfer as taxid_score
import score_marine_gtdb_taxid_transfer as marine


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

BASE_QUALITY = RESULTS / "marine_binomial_transfer_quality.tsv"
BASE_DETAIL = RESULTS / "marine_binomial_transfer_detail.tsv"

OUT_RULES = RESULTS / "marine_truth_upgrade_candidate_rules.tsv"
OUT_RESCUES = RESULTS / "marine_truth_upgrade_candidate_rescues.tsv"
OUT_BLOCKERS = RESULTS / "marine_truth_upgrade_remaining_blockers.tsv"
OUT_AUDIT = RESULTS / "marine_truth_upgrade_candidate_audit.tsv"
OUT_MD = EXP / "MARINE_TRUTH_UPGRADE_CANDIDATE_RULES.md"


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
    text = str(value or "").replace("[", "").replace("]", "").lower()
    return " ".join(text.split())


def truth_epithet(name: object) -> str:
    parts = [part for part in normalize_name(name).split() if part != "candidatus"]
    if len(parts) < 2 or parts[1] in {"sp", "sp."}:
        return ""
    return parts[1]


def gtdb_epithet(species: object) -> str:
    label = str(species or "")
    if label.startswith("s__"):
        label = label[3:]
    parts = label.split()
    if len(parts) < 2:
        return ""
    epithet = re.sub(r"_[A-Z]+$", "", parts[1]).lower()
    return "" if epithet.startswith("sp") else epithet


def build_unique_epithet_rescues() -> list[dict[str, object]]:
    by_accession, _by_core = marine.truth.load_gtdb_metadata(marine.truth.GTDB_METADATA)
    _taxid_to_gtdb, ambiguous_taxids = taxid_score.build_taxid_transfer(by_accession)
    detail = pd.read_csv(BASE_DETAIL, sep="\t")
    rescues: list[dict[str, object]] = []
    ambiguous = detail[detail["transfer_method"] == "ambiguous_taxid"]
    for row in ambiguous.itertuples(index=False):
        taxid = str(row.ncbi_species_taxid)
        epithet = truth_epithet(row.truth_name)
        candidates = sorted(set(ambiguous_taxids.get(taxid, [])))
        matches = [species for species in candidates if epithet and gtdb_epithet(species) == epithet]
        matches = sorted(set(matches))
        if len(matches) != 1:
            continue
        rescues.append(
            {
                "sample": int(row.sample),
                "ncbi_species_taxid": taxid,
                "truth_name": row.truth_name,
                "abundance_pct_all": finite(row.abundance_pct_all),
                "candidate_gtdb_species": matches[0],
                "candidate_count_for_taxid": len(candidates),
                "matched_epithet": epithet,
                "rule_id": "diagnostic_unique_ambiguous_epithet",
            }
        )
    return rescues


def summarize_rules(rescues: list[dict[str, object]]) -> list[dict[str, object]]:
    quality = pd.read_csv(BASE_QUALITY, sep="\t")
    detail = pd.read_csv(BASE_DETAIL, sep="\t")
    rescue_df = pd.DataFrame(rescues)
    rescue_mass = (
        rescue_df.groupby("sample")["abundance_pct_all"].sum().to_dict()
        if not rescue_df.empty
        else {}
    )
    ambiguous_mass = (
        detail[detail["transfer_method"] == "ambiguous_taxid"]
        .groupby("sample")["abundance_pct_all"]
        .sum()
        .to_dict()
    )
    rows: list[dict[str, object]] = []
    for row in quality.itertuples(index=False):
        sample = int(row.sample)
        ba_mass = finite(row.gold_bacteria_archaea_species_mass_pct_all)
        base_mapped = finite(row.truth_mass_mapped_pct_all)
        base_pct = finite(row.truth_mass_mapped_pct_bacteria_archaea)
        diag_extra = finite(rescue_mass.get(sample, 0.0))
        unsafe_extra = finite(ambiguous_mass.get(sample, 0.0))
        for rule_id, extra, release_allowed, description in [
            (
                "accepted_exact_binomial",
                0.0,
                True,
                "unique NCBI taxid plus exact unique GTDB-binomial fallback",
            ),
            (
                "diagnostic_unique_ambiguous_epithet",
                diag_extra,
                False,
                "within ambiguous taxid, one GTDB candidate has matching species epithet",
            ),
            (
                "unsafe_all_ambiguous_assigned",
                unsafe_extra,
                False,
                "upper bound if every ambiguous taxid could be assigned",
            ),
        ]:
            mapped_pct_all = base_mapped + extra
            mapped_pct_ba = 100.0 * mapped_pct_all / ba_mass if ba_mass else 0.0
            rows.append(
                {
                    "sample": sample,
                    "rule_id": rule_id,
                    "release_grade_allowed": str(release_allowed).lower(),
                    "description": description,
                    "extra_mapped_pct_all": extra,
                    "mapped_pct_all": mapped_pct_all,
                    "mapped_pct_bacteria_archaea": mapped_pct_ba,
                    "ready_at_95pct_bacteria_archaea": str(mapped_pct_ba >= 95.0).lower(),
                    "base_mapped_pct_bacteria_archaea": base_pct,
                }
            )
    return rows


def build_blockers() -> list[dict[str, object]]:
    detail = pd.read_csv(BASE_DETAIL, sep="\t")
    rem = detail[detail["transfer_method"].isin(["ambiguous_taxid", "unmapped"])].copy()
    grouped = (
        rem.groupby(["transfer_method", "ncbi_species_taxid", "truth_name"], as_index=False)
        .agg(
            total_abundance_pct_all=("abundance_pct_all", "sum"),
            sample_count=("sample", "nunique"),
            max_taxid_gtdb_species_count=("taxid_gtdb_species_count", "max"),
            max_binomial_gtdb_species_count=("binomial_gtdb_species_count", "max"),
        )
        .sort_values("total_abundance_pct_all", ascending=False)
    )
    return grouped.head(40).to_dict(orient="records")


def build_audit(rule_rows: list[dict[str, object]], rescues: list[dict[str, object]]) -> list[dict[str, object]]:
    df = pd.DataFrame(rule_rows)
    audit_rows: list[dict[str, object]] = [
        {
            "metric": "candidate_rules_tested",
            "value": ",".join(df["rule_id"].drop_duplicates().tolist()),
            "evidence": str(OUT_RULES.relative_to(EXP)),
            "decision": "diagnostic_rules_only",
        }
    ]
    for rule_id, group in df.groupby("rule_id", sort=False):
        min_mapped = float(group["mapped_pct_bacteria_archaea"].min())
        ready_n = int((group["mapped_pct_bacteria_archaea"] >= 95.0).sum())
        allowed = str(group["release_grade_allowed"].iloc[0]).lower() == "true"
        audit_rows.append(
            {
                "metric": f"{rule_id}_summary",
                "value": (
                    f"min_mapped_pct_bacteria_archaea={min_mapped:.6f};"
                    f"ready_samples={ready_n}/{len(group)}"
                ),
                "evidence": str(OUT_RULES.relative_to(EXP)),
                "decision": (
                    "release_rule_fails_threshold"
                    if allowed and ready_n < len(group)
                    else "diagnostic_not_release_rule"
                    if not allowed
                    else "release_rule_passes_threshold"
                ),
            }
        )
    diag = df[df["rule_id"] == "diagnostic_unique_ambiguous_epithet"]
    diag_min = float(diag["mapped_pct_bacteria_archaea"].min())
    diag_ready = int((diag["mapped_pct_bacteria_archaea"] >= 95.0).sum())
    decision = (
        "candidate_rules_do_not_clear_release_threshold"
        if diag_ready < len(diag)
        else "diagnostic_rule_clears_threshold_but_not_release_grade"
    )
    audit_rows.append(
        {
            "metric": "diagnostic_unique_epithet_rescues",
            "value": (
                f"rescued_rows={len(rescues)};"
                f"rescued_mass_pct_all={sum(finite(row['abundance_pct_all']) for row in rescues):.6f};"
                f"min_mapped_pct_bacteria_archaea={diag_min:.6f};"
                f"ready_samples={diag_ready}/{len(diag)}"
            ),
            "evidence": str(OUT_RESCUES.relative_to(EXP)),
            "decision": decision,
        }
    )
    audit_rows.append(
        {
            "metric": "promotion_decision",
            "value": "do_not_promote_marine_truth_upgrade",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "truth_ambiguity_not_resolved_by_current_candidate_rules",
        }
    )
    return audit_rows


def write_markdown(audit_rows: list[dict[str, object]], blockers: list[dict[str, object]]) -> None:
    lines = [
        "# Marine Truth Upgrade Candidate Rules",
        "",
        "Date: 2026-06-29",
        "",
        "This cached audit tests whether deterministic-looking fallback rules can",
        "raise marine GTDB truth mapping to the 95% release threshold. It does not",
        "rerun profilers and does not change the accepted release truth policy.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---|---|",
    ]
    for row in audit_rows:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Main Remaining Blockers",
            "",
            "| Method | Taxid | Name | Total mass | Samples | Candidate count |",
            "|---|---:|---|---:|---:|---:|",
        ]
    )
    for row in blockers[:20]:
        lines.append(
            "| {transfer_method} | {ncbi_species_taxid} | {truth_name} | {mass:.4f} | {samples} | {candidates} |".format(
                transfer_method=row["transfer_method"],
                ncbi_species_taxid=row["ncbi_species_taxid"],
                truth_name=row["truth_name"],
                mass=float(row["total_abundance_pct_all"]),
                samples=int(row["sample_count"]),
                candidates=int(row["max_taxid_gtdb_species_count"]),
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "The current marine route should not be promoted from cached truth alone.",
            "The relaxed epithet fallback is diagnostic only and still leaves most",
            "samples below 95% mapped Bacteria/Archaea truth mass.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_RULES.relative_to(EXP)}`",
            f"- `{OUT_RESCUES.relative_to(EXP)}`",
            f"- `{OUT_BLOCKERS.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rescues = build_unique_epithet_rescues()
    rule_rows = summarize_rules(rescues)
    blockers = build_blockers()
    audit_rows = build_audit(rule_rows, rescues)
    write_tsv(
        OUT_RULES,
        rule_rows,
        [
            "sample",
            "rule_id",
            "release_grade_allowed",
            "description",
            "extra_mapped_pct_all",
            "mapped_pct_all",
            "mapped_pct_bacteria_archaea",
            "ready_at_95pct_bacteria_archaea",
            "base_mapped_pct_bacteria_archaea",
        ],
    )
    write_tsv(
        OUT_RESCUES,
        rescues,
        [
            "sample",
            "ncbi_species_taxid",
            "truth_name",
            "abundance_pct_all",
            "candidate_gtdb_species",
            "candidate_count_for_taxid",
            "matched_epithet",
            "rule_id",
        ],
    )
    write_tsv(
        OUT_BLOCKERS,
        blockers,
        [
            "transfer_method",
            "ncbi_species_taxid",
            "truth_name",
            "total_abundance_pct_all",
            "sample_count",
            "max_taxid_gtdb_species_count",
            "max_binomial_gtdb_species_count",
        ],
    )
    write_tsv(OUT_AUDIT, audit_rows, ["metric", "value", "evidence", "decision"])
    write_markdown(audit_rows, blockers)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
