#!/usr/bin/env python3
"""Audit whether exact split work can be restricted to candidate references.

This is a runtime-design audit, not a new profiler.  It asks two separate
questions:

1. If the full exact table already exists, how much of that table belongs to
   final called taxids?  This estimates post-hoc output reduction only.
2. Do the existing aggregate tables prove that a future exact rerun against
   only candidate refs would produce the same exact evidence?  They do not,
   because best-diff assignment depends on per-context competitors that are not
   present in the aggregate per-reference outputs.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")

CASES = [
    {
        "dataset": "hmp_gastrooral_r232",
        "sample": "0",
        "class": "exact_needed_sidecar_matches_rerun",
        "scope": "bacteria",
        "taxmap": "/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv",
        "profile": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/minco_sample0_sidecar.tsv",
        "unique": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work/minco.best_diff_unique.unfiltered.tsv",
        "block": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work/minco.best_diff_split.unfiltered.tsv",
        "exact": "/tmp/minco_hmp_gastrooral_sidecar_speed_20260627/sample0_work/minco.best_diff_split.exact.unfiltered.tsv",
    },
    {
        "dataset": "cami2_toy_mouse",
        "sample": "6",
        "class": "exact_sidecar_available_but_current_default_skips",
        "scope": "bacteria",
        "taxmap": "/tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv",
        "profile": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16.tsv",
        "unique": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work/minco.best_diff_unique.unfiltered.tsv",
        "block": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work/minco.best_diff_split.unfiltered.tsv",
        "exact": "/tmp/minco_exact_split_sidecar_wrapper_20260627/sample6_sidecar_p16_work/minco.best_diff_split.exact.unfiltered.tsv",
    },
]


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def extract_accession(text: str) -> str:
    match = ACC_RE.search(text or "")
    return match.group(1) if match else ""


def as_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return 0.0
    return out if math.isfinite(out) else 0.0


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def collect_accessions(paths: Iterable[Path]) -> set[str]:
    accessions: set[str] = set()
    for path in paths:
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                acc = extract_accession(row.get("Ref", ""))
                if acc:
                    accessions.add(acc)
    return accessions


def scope_from_taxpath(taxpath: str, taxpathsn: str) -> str:
    combined = f"{taxpath}|{taxpathsn}"
    parts = [part.strip() for part in combined.split("|") if part.strip()]
    if any(part in {"d__Bacteria", "Bacteria"} for part in parts):
        return "bacteria"
    if any(part in {"d__Archaea", "Archaea"} for part in parts):
        return "archaea"
    if any(part in {"d__Viruses", "Viruses", "Virus"} for part in parts):
        return "virus"
    return "other"


def keep_scope(taxid: str, wanted: str, taxid_scope: dict[str, str]) -> bool:
    if wanted == "all":
        return True
    observed = taxid_scope.get(taxid, "other")
    if wanted == "prokaryote":
        return observed in {"bacteria", "archaea"}
    return observed == wanted


def read_taxmap_subset(path: Path, needed_accessions: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    ref_to_taxid: dict[str, str] = {}
    taxid_scope: dict[str, str] = {}
    if not path.exists():
        return ref_to_taxid, taxid_scope
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < 5:
                continue
            key, taxid, rank, taxpath, taxpathsn = fields[:5]
            if key.lower() in {"key", "ref_key", "accession"} or rank != "species":
                continue
            acc = extract_accession(key)
            if acc and acc in needed_accessions:
                ref_to_taxid.setdefault(acc, taxid)
                taxid_scope.setdefault(taxid, scope_from_taxpath(taxpath, taxpathsn))
    return ref_to_taxid, taxid_scope


def read_profile(path: Path) -> tuple[int, set[str]]:
    if not path.exists():
        return 0, set()
    called: set[str] = set()
    rows = 0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            taxid = str(row.get("taxid", "")).strip()
            if taxid and truthy(row.get("calibrated_call", "")):
                called.add(taxid)
    return rows, called


def summarize_raw_table(
    path: Path,
    ref_to_taxid: dict[str, str],
    taxid_scope: dict[str, str],
    scope: str,
    called_taxids: set[str],
) -> dict[str, object]:
    out: dict[str, object] = {
        "exists": str(path.exists()).lower(),
        "rows": 0,
        "mapped_rows": 0,
        "scoped_rows": 0,
        "scoped_taxids": 0,
        "strict_candidate_rows": 0,
        "strict_candidate_taxids": 0,
        "called_taxid_rows": 0,
        "called_taxids_seen": 0,
        "called_taxids_missing": len(called_taxids),
    }
    if not path.exists():
        return out
    scoped_taxids: set[str] = set()
    strict_taxids: set[str] = set()
    called_seen: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out["rows"] = int(out["rows"]) + 1
            acc = extract_accession(row.get("Ref", ""))
            taxid = ref_to_taxid.get(acc, "")
            if not taxid:
                continue
            out["mapped_rows"] = int(out["mapped_rows"]) + 1
            if not keep_scope(taxid, scope, taxid_scope):
                continue
            out["scoped_rows"] = int(out["scoped_rows"]) + 1
            scoped_taxids.add(taxid)
            if taxid in called_taxids:
                out["called_taxid_rows"] = int(out["called_taxid_rows"]) + 1
                called_seen.add(taxid)
            strict = (
                as_float(row.get("XnY_ctx", 0.0)) >= 10.0
                and as_float(row.get("ANI", 0.0)) >= 0.95
                and as_float(row.get("Real_min_align_fraction", 0.0)) >= 0.05
            )
            if strict:
                out["strict_candidate_rows"] = int(out["strict_candidate_rows"]) + 1
                strict_taxids.add(taxid)
    out["scoped_taxids"] = len(scoped_taxids)
    out["strict_candidate_taxids"] = len(strict_taxids)
    out["called_taxids_seen"] = len(called_seen)
    out["called_taxids_missing"] = len(called_taxids - called_seen)
    return out


def ratio(num: float, den: float) -> str:
    if den <= 0:
        return "NA"
    return f"{num / den:.8f}"


def case_summary(spec: dict[str, str]) -> dict[str, object]:
    profile = Path(spec["profile"])
    unique = Path(spec["unique"])
    block = Path(spec["block"])
    exact = Path(spec["exact"])
    taxmap = Path(spec["taxmap"])
    profile_rows, called_taxids = read_profile(profile)
    needed_accs = collect_accessions([unique, block, exact])
    ref_to_taxid, taxid_scope = read_taxmap_subset(taxmap, needed_accs)
    unique_sum = summarize_raw_table(unique, ref_to_taxid, taxid_scope, spec["scope"], called_taxids)
    block_sum = summarize_raw_table(block, ref_to_taxid, taxid_scope, spec["scope"], called_taxids)
    exact_sum = summarize_raw_table(exact, ref_to_taxid, taxid_scope, spec["scope"], called_taxids)

    exact_scoped = int(exact_sum["scoped_rows"])
    exact_called_rows = int(exact_sum["called_taxid_rows"])
    exact_eliminable = max(0, exact_scoped - exact_called_rows)
    strict_exact_rows = int(exact_sum["strict_candidate_rows"])
    strict_block_rows = int(block_sum["strict_candidate_rows"])
    return {
        "dataset": spec["dataset"],
        "sample": spec["sample"],
        "class": spec["class"],
        "scope": spec["scope"],
        "profile_exists": str(profile.exists()).lower(),
        "profile_rows": profile_rows,
        "final_called_taxids": len(called_taxids),
        "needed_ref_accessions": len(needed_accs),
        "taxmap_accessions_mapped": len(ref_to_taxid),
        "unique_rows": unique_sum["rows"],
        "unique_scoped_rows": unique_sum["scoped_rows"],
        "unique_strict_candidate_taxids": unique_sum["strict_candidate_taxids"],
        "block_rows": block_sum["rows"],
        "block_scoped_rows": block_sum["scoped_rows"],
        "block_strict_candidate_rows": strict_block_rows,
        "block_strict_candidate_taxids": block_sum["strict_candidate_taxids"],
        "exact_rows": exact_sum["rows"],
        "exact_scoped_rows": exact_scoped,
        "exact_strict_candidate_rows": strict_exact_rows,
        "exact_strict_candidate_taxids": exact_sum["strict_candidate_taxids"],
        "exact_rows_with_final_called_taxids": exact_called_rows,
        "exact_called_taxid_keep_fraction": ratio(exact_called_rows, exact_scoped),
        "exact_rows_eliminable_posthoc": exact_eliminable,
        "exact_rows_eliminable_posthoc_fraction": ratio(exact_eliminable, exact_scoped),
        "final_called_taxids_seen_in_exact": exact_sum["called_taxids_seen"],
        "final_called_taxids_missing_in_exact": exact_sum["called_taxids_missing"],
        "strict_exact_rows_vs_block_delta": strict_exact_rows - strict_block_rows,
        "posthoc_filter_saves_output_rows": str(exact_eliminable > 0).lower(),
        "true_subset_exact_rerun_safety_from_aggregate_outputs": "not_proven",
        "decision": "posthoc_output_filter_only_not_safe_scan_restriction",
        "reason": "aggregate ref rows omit per-context competitor closure needed for best-diff exact assignment",
        "profile": str(profile),
        "unique": str(unique),
        "block": str(block),
        "exact": str(exact),
        "taxmap": str(taxmap),
    }


def main() -> int:
    summaries = [case_summary(spec) for spec in CASES]
    evaluated = [row for row in summaries if row["profile_exists"] == "true" and int(row["exact_rows"]) > 0]
    if evaluated:
        min_keep = min(
            float(row["exact_called_taxid_keep_fraction"])
            for row in evaluated
            if row["exact_called_taxid_keep_fraction"] != "NA"
        )
        max_elim_frac = max(
            float(row["exact_rows_eliminable_posthoc_fraction"])
            for row in evaluated
            if row["exact_rows_eliminable_posthoc_fraction"] != "NA"
        )
        total_missing_called = sum(int(row["final_called_taxids_missing_in_exact"]) for row in evaluated)
    else:
        min_keep = math.nan
        max_elim_frac = math.nan
        total_missing_called = 0

    audit_rows = [
        {
            "metric": "evaluated_cases",
            "value": len(evaluated),
            "evidence": "candidate_restricted_exact_feasibility.tsv",
            "decision": "cached_real_exact_tables",
        },
        {
            "metric": "min_exact_called_taxid_keep_fraction",
            "value": f"{min_keep:.8f}" if math.isfinite(min_keep) else "NA",
            "evidence": "candidate_restricted_exact_feasibility.tsv",
            "decision": "posthoc_filter_can_remove_many_ref_rows",
        },
        {
            "metric": "max_exact_rows_eliminable_posthoc_fraction",
            "value": f"{max_elim_frac:.8f}" if math.isfinite(max_elim_frac) else "NA",
            "evidence": "candidate_restricted_exact_feasibility.tsv",
            "decision": "output_size_reduction_not_scan_cost_reduction",
        },
        {
            "metric": "final_called_taxids_missing_in_exact_total",
            "value": total_missing_called,
            "evidence": "candidate_restricted_exact_feasibility.tsv",
            "decision": "sanity_check_final_calls_have_exact_or_unique_support",
        },
        {
            "metric": "candidate_restricted_exact_decision",
            "value": "posthoc_filter_only_not_safe_as_candidate_subset_rerun",
            "evidence": "candidate_restricted_exact_feasibility.tsv; command_ani.c best-diff candidate competition",
            "decision": "do_not_promote_candidate_only_exact_rerun",
        },
        {
            "metric": "required_evidence_for_safe_subset_exact",
            "value": "context_or_read_level_competitor_closure",
            "evidence": "aggregate raw MinCO tables lack per-context alternative refs",
            "decision": "needs_C_level_ambiguity_groups_or_full_candidate_scan_sidecar",
        },
    ]
    fields = [
        "dataset",
        "sample",
        "class",
        "scope",
        "profile_exists",
        "profile_rows",
        "final_called_taxids",
        "needed_ref_accessions",
        "taxmap_accessions_mapped",
        "unique_rows",
        "unique_scoped_rows",
        "unique_strict_candidate_taxids",
        "block_rows",
        "block_scoped_rows",
        "block_strict_candidate_rows",
        "block_strict_candidate_taxids",
        "exact_rows",
        "exact_scoped_rows",
        "exact_strict_candidate_rows",
        "exact_strict_candidate_taxids",
        "exact_rows_with_final_called_taxids",
        "exact_called_taxid_keep_fraction",
        "exact_rows_eliminable_posthoc",
        "exact_rows_eliminable_posthoc_fraction",
        "final_called_taxids_seen_in_exact",
        "final_called_taxids_missing_in_exact",
        "strict_exact_rows_vs_block_delta",
        "posthoc_filter_saves_output_rows",
        "true_subset_exact_rerun_safety_from_aggregate_outputs",
        "decision",
        "reason",
        "profile",
        "unique",
        "block",
        "exact",
        "taxmap",
    ]
    write_tsv(RESULTS / "candidate_restricted_exact_feasibility.tsv", summaries, fields)
    write_tsv(
        RESULTS / "candidate_restricted_exact_feasibility_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
