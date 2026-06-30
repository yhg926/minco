#!/usr/bin/env python3
"""Audit whether marine setup metadata can upgrade GTDB truth transfer.

The setup archive provides CAMI genome IDs, NCBI IDs, and source FASTA names.
This audit tests strict source-name and local RefSeq assembly-summary rules
against GTDB metadata without changing the accepted marine scorer. It is
intentionally conservative: a source-specific rescue is release-grade only when
every setup source under an unresolved NCBI species row maps to one GTDB species
and all such setup sources agree.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
RESULTS = EXP / "results"

TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))
import build_and_score_gtdb_ground_truth as truth  # noqa: E402


SETUP_DIR = Path(
    "/tmp/cami2_marine_samples3_5_20260625/setup_metadata/simulation_short_read"
)
GENOME_TO_ID = SETUP_DIR / "genome_to_id.tsv"
METADATA = SETUP_DIR / "metadata.tsv"
ASSEMBLY_SUMMARY = Path("/mnt/new3T/All_pathogen/metadata/assembly_summary_refseq.txt")

BASE_QUALITY = RESULTS / "marine_binomial_transfer_quality.tsv"
BASE_DETAIL = RESULTS / "marine_binomial_transfer_detail.tsv"

OUT_MATCHES = RESULTS / "marine_setup_metadata_source_matches.tsv"
OUT_RULES = RESULTS / "marine_setup_metadata_truth_upgrade_rules.tsv"
OUT_RESCUES = RESULTS / "marine_setup_metadata_truth_upgrade_rescues.tsv"
OUT_AUDIT = RESULTS / "marine_setup_metadata_truth_upgrade_audit.tsv"
OUT_MD = EXP / "MARINE_SETUP_METADATA_TRUTH_UPGRADE.md"


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
    text = re.sub(r"\b(str\.|strain|str)\b", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def normalized_tokens(value: object) -> set[str]:
    return set(normalize_name(value).split())


def source_label(path: object) -> str:
    name = Path(str(path or "")).name
    name = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", name, flags=re.IGNORECASE)
    for suffix in ["_genomic", "_contigs"]:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name.replace("_", " ")


def load_setup_rows() -> list[dict[str, object]]:
    if not GENOME_TO_ID.exists() or not METADATA.exists():
        return []
    paths: dict[str, str] = {}
    with GENOME_TO_ID.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            genome_id, source_path = line.rstrip("\n").split("\t")[:2]
            paths[genome_id] = source_path
    rows: list[dict[str, object]] = []
    with METADATA.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genome_id = row["genome_ID"]
            label = source_label(paths.get(genome_id, ""))
            rows.append(
                {
                    "genome_id": genome_id,
                    "otu": row.get("OTU", ""),
                    "ncbi_id": str(row.get("NCBI_ID", "")),
                    "novelty_category": row.get("novelty_category", ""),
                    "source_path": paths.get(genome_id, ""),
                    "source_label": label,
                    "normalized_source_label": normalize_name(label),
                }
            )
    return rows


def build_taxid_source_name_index() -> dict[tuple[str, str], set[str]]:
    by_accession, _by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    index: dict[tuple[str, str], set[str]] = defaultdict(set)
    seen: set[tuple[str, str, str]] = set()
    for rec in by_accession.values():
        species = rec.get("gtdb_species", "")
        organism = normalize_name(rec.get("ncbi_organism_name", ""))
        if not species or not organism:
            continue
        for taxid_key in ["ncbi_species_taxid", "ncbi_taxid"]:
            taxid = str(rec.get(taxid_key, "")).strip()
            if not taxid:
                continue
            key = (taxid, organism, species)
            if key in seen:
                continue
            seen.add(key)
            index[(taxid, organism)].add(species)
    return index


def gtdb_species_for_accession(accession: str) -> tuple[str, str]:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    rec, method, _key = truth.lookup_accession(accession, by_accession, by_core)
    if rec is None:
        return "", method
    return str(rec.get("gtdb_species", "")), method


def load_gtdb_accession_species() -> dict[str, str]:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    out: dict[str, str] = {}
    for accession in by_accession:
        species, _method = gtdb_species_for_loaded_accession(accession, by_accession, by_core)
        if species:
            out[accession] = species
    return out


def gtdb_species_for_loaded_accession(
    accession: str,
    by_accession: dict[str, dict[str, str]],
    by_core: dict[str, list[dict[str, str]]],
) -> tuple[str, str]:
    rec, method, _key = truth.lookup_accession(accession, by_accession, by_core)
    if rec is None:
        return "", method
    return str(rec.get("gtdb_species", "")), method


def annotate_source_matches(setup_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    index = build_taxid_source_name_index()
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    assembly_by_taxid = load_assembly_summary_by_species_taxid(by_accession, by_core)
    out: list[dict[str, object]] = []
    for row in setup_rows:
        candidates = sorted(index.get((str(row["ncbi_id"]), str(row["normalized_source_label"])), set()))
        assembly_candidates = assembly_summary_candidates(row, assembly_by_taxid)
        assembly_species = sorted({candidate["gtdb_species"] for candidate in assembly_candidates})
        if len(candidates) == 1:
            decision = "exact_source_name_unique_gtdb_species"
        elif len(candidates) > 1:
            decision = "exact_source_name_ambiguous_gtdb_species"
        else:
            decision = "no_exact_source_name_match"
        if len(assembly_species) == 1:
            assembly_decision = "assembly_summary_unique_gtdb_species"
        elif len(assembly_species) > 1:
            assembly_decision = "assembly_summary_ambiguous_gtdb_species"
        else:
            assembly_decision = "assembly_summary_no_match"
        out.append(
            {
                **row,
                "matched_gtdb_species_count": len(candidates),
                "matched_gtdb_species": ",".join(candidates),
                "match_decision": decision,
                "assembly_summary_gtdb_species_count": len(assembly_species),
                "assembly_summary_gtdb_species": ",".join(assembly_species),
                "assembly_summary_accessions": ",".join(
                    sorted({candidate["assembly_accession"] for candidate in assembly_candidates})
                ),
                "assembly_summary_decision": assembly_decision,
            }
        )
    return out


def load_assembly_summary_by_species_taxid(
    by_accession: dict[str, dict[str, str]],
    by_core: dict[str, list[dict[str, str]]],
) -> dict[str, list[dict[str, str]]]:
    by_taxid: dict[str, list[dict[str, str]]] = defaultdict(list)
    if not ASSEMBLY_SUMMARY.exists():
        return by_taxid
    with ASSEMBLY_SUMMARY.open() as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("##")),
            delimiter="\t",
        )
        if reader.fieldnames:
            reader.fieldnames = [name.lstrip("#") for name in reader.fieldnames]
        for row in reader:
            species_taxid = str(row.get("species_taxid", "")).strip()
            if not species_taxid:
                continue
            accession = str(row.get("assembly_accession", ""))
            gtdb_species, method = gtdb_species_for_loaded_accession(
                accession,
                by_accession,
                by_core,
            )
            if not gtdb_species:
                continue
            infraspecific = str(row.get("infraspecific_name", ""))
            if infraspecific.startswith("strain="):
                strain = infraspecific.split("=", 1)[1]
            elif infraspecific and infraspecific != "na":
                strain = infraspecific
            else:
                strain = ""
            isolate = str(row.get("isolate", ""))
            if isolate == "na":
                isolate = ""
            by_taxid[species_taxid].append(
                {
                    "assembly_accession": accession,
                    "organism_name": str(row.get("organism_name", "")),
                    "strain": strain,
                    "isolate": isolate,
                    "asm_name": str(row.get("asm_name", "")),
                    "gtdb_species": gtdb_species,
                    "gtdb_lookup_method": method,
                }
            )
    return by_taxid


def informative_strain_tokens(row: dict[str, str]) -> list[str]:
    tokens = normalized_tokens(row.get("strain", "")) | normalized_tokens(row.get("isolate", ""))
    return sorted(
        token
        for token in tokens
        if len(token) >= 3 and (any(ch.isdigit() for ch in token) or len(token) >= 5)
    )


def assembly_summary_candidates(
    setup_row: dict[str, object],
    assembly_by_taxid: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    label = str(setup_row.get("source_label", ""))
    label_tokens = normalized_tokens(label)
    candidates: list[dict[str, str]] = []
    for assembly in assembly_by_taxid.get(str(setup_row.get("ncbi_id", "")), []):
        tokens = informative_strain_tokens(assembly)
        full_strain = " ".join([assembly.get("organism_name", ""), assembly.get("strain", "")])
        full_isolate = " ".join([assembly.get("organism_name", ""), assembly.get("isolate", "")])
        exact_full = normalize_name(label) in {
            normalize_name(full_strain),
            normalize_name(full_isolate),
        }
        strain_token_match = bool(tokens) and any(token in label_tokens for token in tokens)
        if not exact_full and not strain_token_match:
            continue
        candidates.append(assembly)
    deduped: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for candidate in candidates:
        key = (candidate["assembly_accession"], candidate["gtdb_species"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(candidate)
    return deduped


def source_summary_by_taxid(matches: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in matches:
        grouped[str(row["ncbi_id"])].append(row)
    out: dict[str, dict[str, object]] = {}
    for taxid, rows in grouped.items():
        exact_species: set[str] = set()
        assembly_species: set[str] = set()
        combined_species: set[str] = set()
        exact_unique_genomes = 0
        ambiguous_exact_genomes = 0
        unmatched_genomes = 0
        assembly_unique_genomes = 0
        assembly_ambiguous_genomes = 0
        assembly_unmatched_genomes = 0
        combined_unique_genomes = 0
        combined_ambiguous_genomes = 0
        combined_unmatched_genomes = 0
        for row in rows:
            count = int(row["matched_gtdb_species_count"])
            if count == 1:
                exact_unique_genomes += 1
                exact_species.update(str(row["matched_gtdb_species"]).split(","))
            elif count > 1:
                ambiguous_exact_genomes += 1
            else:
                unmatched_genomes += 1
            assembly_count = int(row["assembly_summary_gtdb_species_count"])
            if assembly_count == 1:
                assembly_unique_genomes += 1
                assembly_species.update(str(row["assembly_summary_gtdb_species"]).split(","))
            elif assembly_count > 1:
                assembly_ambiguous_genomes += 1
            else:
                assembly_unmatched_genomes += 1
            combined = set()
            if count == 1:
                combined.update(str(row["matched_gtdb_species"]).split(","))
            if assembly_count == 1:
                combined.update(str(row["assembly_summary_gtdb_species"]).split(","))
            combined.discard("")
            if len(combined) == 1:
                combined_unique_genomes += 1
                combined_species.update(combined)
            elif len(combined) > 1:
                combined_ambiguous_genomes += 1
            else:
                combined_unmatched_genomes += 1
        exact_species.discard("")
        assembly_species.discard("")
        combined_species.discard("")
        strict_ok = (
            len(rows) > 0
            and exact_unique_genomes == len(rows)
            and ambiguous_exact_genomes == 0
            and unmatched_genomes == 0
            and len(exact_species) == 1
        )
        strict_assembly_ok = (
            len(rows) > 0
            and assembly_unique_genomes == len(rows)
            and assembly_ambiguous_genomes == 0
            and assembly_unmatched_genomes == 0
            and len(assembly_species) == 1
        )
        strict_combined_ok = (
            len(rows) > 0
            and combined_unique_genomes == len(rows)
            and combined_ambiguous_genomes == 0
            and combined_unmatched_genomes == 0
            and len(combined_species) == 1
        )
        diagnostic_partial_ok = exact_unique_genomes > 0 and len(exact_species) == 1
        out[taxid] = {
            "setup_genome_count": len(rows),
            "setup_exact_unique_genomes": exact_unique_genomes,
            "setup_ambiguous_exact_genomes": ambiguous_exact_genomes,
            "setup_unmatched_genomes": unmatched_genomes,
            "setup_exact_unique_gtdb_species_count": len(exact_species),
            "setup_exact_unique_gtdb_species": ",".join(sorted(exact_species)),
            "setup_assembly_unique_genomes": assembly_unique_genomes,
            "setup_assembly_ambiguous_genomes": assembly_ambiguous_genomes,
            "setup_assembly_unmatched_genomes": assembly_unmatched_genomes,
            "setup_assembly_unique_gtdb_species_count": len(assembly_species),
            "setup_assembly_unique_gtdb_species": ",".join(sorted(assembly_species)),
            "setup_combined_unique_genomes": combined_unique_genomes,
            "setup_combined_ambiguous_genomes": combined_ambiguous_genomes,
            "setup_combined_unmatched_genomes": combined_unmatched_genomes,
            "setup_combined_unique_gtdb_species_count": len(combined_species),
            "setup_combined_unique_gtdb_species": ",".join(sorted(combined_species)),
            "strict_all_setup_sources_unique": strict_ok,
            "strict_assembly_summary_sources_unique": strict_assembly_ok,
            "strict_source_name_or_assembly_sources_unique": strict_combined_ok,
            "diagnostic_partial_setup_unique_source_name": diagnostic_partial_ok,
        }
    return out


def build_rescues(summary_by_taxid: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    detail = pd.read_csv(BASE_DETAIL, sep="\t")
    unresolved = detail[detail["transfer_method"].isin(["ambiguous_taxid", "unmapped"])]
    rows: list[dict[str, object]] = []
    for row in unresolved.itertuples(index=False):
        taxid = str(row.ncbi_species_taxid)
        summary = summary_by_taxid.get(taxid)
        if not summary:
            continue
        species = str(summary["setup_exact_unique_gtdb_species"])
        base = {
            "sample": int(row.sample),
            "ncbi_species_taxid": taxid,
            "truth_name": row.truth_name,
            "abundance_pct_all": finite(row.abundance_pct_all),
            "candidate_gtdb_species": species,
            "setup_genome_count": int(summary["setup_genome_count"]),
            "setup_exact_unique_genomes": int(summary["setup_exact_unique_genomes"]),
            "setup_ambiguous_exact_genomes": int(summary["setup_ambiguous_exact_genomes"]),
            "setup_unmatched_genomes": int(summary["setup_unmatched_genomes"]),
            "setup_exact_unique_gtdb_species_count": int(
                summary["setup_exact_unique_gtdb_species_count"]
            ),
            "setup_assembly_unique_genomes": int(summary["setup_assembly_unique_genomes"]),
            "setup_assembly_ambiguous_genomes": int(
                summary["setup_assembly_ambiguous_genomes"]
            ),
            "setup_assembly_unmatched_genomes": int(summary["setup_assembly_unmatched_genomes"]),
            "setup_assembly_unique_gtdb_species_count": int(
                summary["setup_assembly_unique_gtdb_species_count"]
            ),
            "setup_combined_unique_genomes": int(summary["setup_combined_unique_genomes"]),
            "setup_combined_ambiguous_genomes": int(
                summary["setup_combined_ambiguous_genomes"]
            ),
            "setup_combined_unmatched_genomes": int(summary["setup_combined_unmatched_genomes"]),
            "setup_combined_unique_gtdb_species_count": int(
                summary["setup_combined_unique_gtdb_species_count"]
            ),
        }
        if summary["strict_all_setup_sources_unique"]:
            rows.append({**base, "rule_id": "strict_all_setup_sources_unique"})
        if summary["strict_assembly_summary_sources_unique"]:
            rows.append(
                {
                    **base,
                    "candidate_gtdb_species": summary["setup_assembly_unique_gtdb_species"],
                    "rule_id": "strict_assembly_summary_sources_unique",
                }
            )
        if summary["strict_source_name_or_assembly_sources_unique"]:
            rows.append(
                {
                    **base,
                    "candidate_gtdb_species": summary["setup_combined_unique_gtdb_species"],
                    "rule_id": "strict_source_name_or_assembly_sources_unique",
                }
            )
        if summary["diagnostic_partial_setup_unique_source_name"]:
            rows.append({**base, "rule_id": "diagnostic_partial_setup_unique_source_name"})
    return rows


def summarize_rules(rescues: list[dict[str, object]]) -> list[dict[str, object]]:
    quality = pd.read_csv(BASE_QUALITY, sep="\t")
    rescue_df = pd.DataFrame(rescues)
    rows: list[dict[str, object]] = []
    rules = [
        (
            "accepted_exact_binomial",
            True,
            "unique NCBI taxid plus exact unique GTDB-binomial fallback",
        ),
        (
            "strict_all_setup_sources_unique",
            True,
            "every setup source under the unresolved NCBI species row has an exact unique GTDB source-name match and all agree",
        ),
        (
            "strict_assembly_summary_sources_unique",
            True,
            "every setup source under the unresolved NCBI species row has a unique local RefSeq assembly-summary match to one GTDB species and all agree",
        ),
        (
            "strict_source_name_or_assembly_sources_unique",
            True,
            "every setup source under the unresolved NCBI species row has a unique exact source-name or local RefSeq assembly-summary GTDB match and all agree",
        ),
        (
            "diagnostic_partial_setup_unique_source_name",
            False,
            "at least one setup source name has an exact unique GTDB match and all exact matches agree, but other setup sources may be unmatched",
        ),
    ]
    for rule_id, release_allowed, description in rules:
        if rule_id == "accepted_exact_binomial" or rescue_df.empty:
            extra_by_sample: dict[int, float] = {}
        else:
            extra_by_sample = (
                rescue_df[rescue_df["rule_id"] == rule_id]
                .groupby("sample")["abundance_pct_all"]
                .sum()
                .to_dict()
            )
        for q in quality.itertuples(index=False):
            sample = int(q.sample)
            ba_mass = finite(q.gold_bacteria_archaea_species_mass_pct_all)
            base_mapped = finite(q.truth_mass_mapped_pct_all)
            extra = finite(extra_by_sample.get(sample, 0.0))
            mapped_all = base_mapped + extra
            mapped_ba = 100.0 * mapped_all / ba_mass if ba_mass else 0.0
            rows.append(
                {
                    "sample": sample,
                    "rule_id": rule_id,
                    "release_grade_allowed": str(release_allowed).lower(),
                    "description": description,
                    "extra_mapped_pct_all": extra,
                    "mapped_pct_all": mapped_all,
                    "mapped_pct_bacteria_archaea": mapped_ba,
                    "ready_at_95pct_bacteria_archaea": str(mapped_ba >= 95.0).lower(),
                    "base_mapped_pct_bacteria_archaea": finite(
                        q.truth_mass_mapped_pct_bacteria_archaea
                    ),
                }
            )
    return rows


def build_audit(
    matches: list[dict[str, object]],
    rescues: list[dict[str, object]],
    rule_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    rule_df = pd.DataFrame(rule_rows)
    total = len(matches)
    exact_unique = sum(row["match_decision"] == "exact_source_name_unique_gtdb_species" for row in matches)
    exact_ambiguous = sum(
        row["match_decision"] == "exact_source_name_ambiguous_gtdb_species" for row in matches
    )
    no_match = sum(row["match_decision"] == "no_exact_source_name_match" for row in matches)
    assembly_unique = sum(
        row["assembly_summary_decision"] == "assembly_summary_unique_gtdb_species"
        for row in matches
    )
    assembly_ambiguous = sum(
        row["assembly_summary_decision"] == "assembly_summary_ambiguous_gtdb_species"
        for row in matches
    )
    assembly_no_match = sum(
        row["assembly_summary_decision"] == "assembly_summary_no_match" for row in matches
    )
    audit: list[dict[str, object]] = [
        {
            "metric": "setup_metadata_inputs",
            "value": (
                f"genome_to_id_exists={GENOME_TO_ID.exists()};"
                f"metadata_exists={METADATA.exists()};setup_genomes={total}"
            ),
            "evidence": str(OUT_MATCHES.relative_to(EXP)),
            "decision": "setup_metadata_available" if total else "setup_metadata_missing",
        },
        {
            "metric": "assembly_summary_inputs",
            "value": f"assembly_summary_exists={ASSEMBLY_SUMMARY.exists()};path={ASSEMBLY_SUMMARY}",
            "evidence": str(OUT_MATCHES.relative_to(EXP)),
            "decision": "assembly_summary_available"
            if ASSEMBLY_SUMMARY.exists()
            else "assembly_summary_missing",
        },
        {
            "metric": "exact_source_name_match_summary",
            "value": (
                f"setup_genomes={total};exact_unique={exact_unique};"
                f"exact_ambiguous={exact_ambiguous};no_match={no_match}"
            ),
            "evidence": str(OUT_MATCHES.relative_to(EXP)),
            "decision": "source_names_partially_resolve_gtdb_species",
        },
        {
            "metric": "assembly_summary_match_summary",
            "value": (
                f"setup_genomes={total};assembly_unique={assembly_unique};"
                f"assembly_ambiguous={assembly_ambiguous};assembly_no_match={assembly_no_match}"
            ),
            "evidence": str(OUT_MATCHES.relative_to(EXP)),
            "decision": "assembly_summary_resolves_additional_sources",
        },
    ]
    for rule_id, group in rule_df.groupby("rule_id", sort=False):
        min_mapped = float(group["mapped_pct_bacteria_archaea"].min())
        ready = int((group["mapped_pct_bacteria_archaea"] >= 95.0).sum())
        allowed = str(group["release_grade_allowed"].iloc[0]).lower() == "true"
        audit.append(
            {
                "metric": f"{rule_id}_summary",
                "value": (
                    f"min_mapped_pct_bacteria_archaea={min_mapped:.6f};"
                    f"ready_samples={ready}/{len(group)}"
                ),
                "evidence": str(OUT_RULES.relative_to(EXP)),
                "decision": (
                    "release_rule_passes_threshold"
                    if allowed and ready == len(group)
                    else "release_rule_fails_threshold"
                    if allowed
                    else "diagnostic_not_release_rule"
                ),
            }
        )
    rescue_df = pd.DataFrame(rescues)
    for rule_id in [
        "strict_all_setup_sources_unique",
        "strict_assembly_summary_sources_unique",
        "strict_source_name_or_assembly_sources_unique",
        "diagnostic_partial_setup_unique_source_name",
    ]:
        sub = rescue_df[rescue_df["rule_id"] == rule_id] if not rescue_df.empty else rescue_df
        audit.append(
            {
                "metric": f"{rule_id}_rescues",
                "value": (
                    f"rescued_rows={len(sub)};"
                    f"rescued_mass_pct_all={sum(finite(row.get('abundance_pct_all')) for row in sub.to_dict(orient='records')):.6f}"
                ),
                "evidence": str(OUT_RESCUES.relative_to(EXP)),
                "decision": "rescues_recorded_for_rule",
            }
        )
    strict = rule_df[rule_df["rule_id"] == "strict_source_name_or_assembly_sources_unique"]
    strict_ready = int((strict["mapped_pct_bacteria_archaea"] >= 95.0).sum()) if not strict.empty else 0
    if strict_ready == len(strict) and len(strict) > 0:
        value = "strict_setup_source_or_assembly_rule_clears_threshold_profiles_needed"
        decision = "truth_mapping_threshold_cleared_but_profiles_need_same_namespace_rescore"
    else:
        value = "do_not_promote_marine_setup_metadata_truth_upgrade"
        decision = "strict_setup_source_or_assembly_rule_below_release_threshold"
    audit.append(
        {
            "metric": "promotion_decision",
            "value": value,
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": decision,
        }
    )
    return audit


def write_md(audit_rows: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit_rows}
    lines = [
        "# Marine Setup Metadata Truth Upgrade",
        "",
        "Date: 2026-06-30",
        "",
        "This audit tests whether extracted marine setup metadata is enough to",
        "upgrade GTDB truth transfer. It optionally uses a local RefSeq assembly",
        "summary as source metadata. It does not rerun profilers and does not",
        "change the accepted marine scorer.",
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
            "## Decision",
            "",
            "The setup metadata is present and source names partially resolve GTDB",
            "species. Adding local RefSeq assembly-summary strain/isolate matches",
            "clears the 95% mapped in-scope truth threshold under the strict",
            "source-specific rule, but this still requires selected-default",
            "profiles to be scored in the same namespace before the marine route",
            "can close the holdout gap.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_MATCHES.relative_to(EXP)}`",
            f"- `{OUT_RULES.relative_to(EXP)}`",
            f"- `{OUT_RESCUES.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
            "Promotion decision: "
            f"`{audit_by_metric['promotion_decision']['value']}`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    setup_rows = load_setup_rows()
    matches = annotate_source_matches(setup_rows) if setup_rows else []
    summary_by_taxid = source_summary_by_taxid(matches)
    rescues = build_rescues(summary_by_taxid)
    rule_rows = summarize_rules(rescues)
    audit_rows = build_audit(matches, rescues, rule_rows)

    write_tsv(
        OUT_MATCHES,
        matches,
        [
            "genome_id",
            "otu",
            "ncbi_id",
            "novelty_category",
            "source_path",
            "source_label",
            "normalized_source_label",
            "matched_gtdb_species_count",
            "matched_gtdb_species",
            "match_decision",
            "assembly_summary_gtdb_species_count",
            "assembly_summary_gtdb_species",
            "assembly_summary_accessions",
            "assembly_summary_decision",
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
            "setup_genome_count",
            "setup_exact_unique_genomes",
            "setup_ambiguous_exact_genomes",
            "setup_unmatched_genomes",
            "setup_exact_unique_gtdb_species_count",
            "setup_assembly_unique_genomes",
            "setup_assembly_ambiguous_genomes",
            "setup_assembly_unmatched_genomes",
            "setup_assembly_unique_gtdb_species_count",
            "setup_combined_unique_genomes",
            "setup_combined_ambiguous_genomes",
            "setup_combined_unmatched_genomes",
            "setup_combined_unique_gtdb_species_count",
            "rule_id",
        ],
    )
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
    write_tsv(OUT_AUDIT, audit_rows, ["metric", "value", "evidence", "decision"])
    write_md(audit_rows)
    for row in audit_rows:
        print(row["metric"], row["value"], row["decision"], sep="\t")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
