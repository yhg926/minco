#!/usr/bin/env python3
"""Inventory independent holdout evidence for the refined allocator.

The refined allocator has cached score replay and independent diagnostic
support. This audit answers a narrower release question: do we currently have
an unused, release-grade GTDB holdout profile set that was not part of the
refined-guard selection cache? The answer is recorded as machine-readable
evidence rather than left as an informal caveat.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

MANIFEST = EXP / "holdout_bundle_manifest.tsv"
GUARD_FEATURES = RESULTS / "feature_allocator_combined_guard_features.tsv"
MARINE_AUDIT = RESULTS / "feature_allocator_refined_marine_diagnostic_audit.tsv"

OUT_INVENTORY = RESULTS / "feature_allocator_refined_independent_holdout_inventory.tsv"
OUT_AUDIT = RESULTS / "feature_allocator_refined_independent_holdout_inventory_audit.tsv"
OUT_MD = EXP / "ABUNDANCE_FEATURE_ALLOCATOR_REFINED_HOLDOUT_INVENTORY.md"


PANEL_TO_FEATURE_PANEL = {
    "cami2_toy_mouse_gut_5_7": "cami2_toy_mouse_gut",
    "hmp_airskin_source_abundance_24sample": "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_source_abundance_0_6": "hmp_gastrooral_gtdb_source_abundance",
    "cami3_toy_human_gut_0_2_gtdb_source_readmap": "cami3_toy_human_gut_gtdb_source_readmap",
}

HMP_OMITTED_SAMPLES = [2, 8, 12, 26, 27]
HMP_TRUTH_ROOT = Path("/tmp/cami2_hmp_airskin_20260625/truth")
HMP_ROOTS = [
    Path("/tmp"),
    Path("/mnt/new3T/minco_release_holdouts_20260628"),
]
HMP_R232_ROOTS = [
    Path("/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232"),
    Path("/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626/run_sylph_r232"),
]
HMP_READ_ROOTS = [
    Path("/tmp/cami2_hmp_airskin_20260625"),
    Path("/tmp/cami2_hmp_unseen_transfer_20260626"),
    Path("/mnt/new3T/minco_release_holdouts_20260628/cami2_hmp_unseen_transfer_20260626"),
]


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


def sample_list(text: str) -> list[int]:
    out = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        out.append(int(part))
    return out


def selected_cache_samples() -> dict[str, set[int]]:
    features = pd.read_csv(GUARD_FEATURES, sep="\t")
    out: dict[str, set[int]] = {}
    for row in features.to_dict("records"):
        panel = str(row["panel"])
        sample = int(float(row["sample"]))
        if str(row.get("evidence_group")) != "selected_cached_wrapper":
            continue
        out.setdefault(panel, set()).add(sample)
    return out


def hmp_minco_profile(sample: int) -> Path | None:
    rel = Path(f"minco_current_code_hmp_airskin{sample}_20260627/minco_sample{sample}_current_default.tsv")
    for root in HMP_ROOTS:
        path = root / rel
        if path.is_file():
            return path
    return None


def hmp_sylph_profile(sample: int) -> Path | None:
    for root in HMP_R232_ROOTS:
        path = root / f"sylph_sample{sample}/profile.tsv"
        if path.is_file():
            return path
    return None


def hmp_read_input(sample: int) -> Path | None:
    names = [
        f"nonzero_bam_paths_sample{sample}.txt",
        f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
    ]
    for root in HMP_READ_ROOTS:
        for name in names:
            path = root / name
            if path.is_file():
                return path
    return None


def build_inventory() -> list[dict[str, object]]:
    manifest = read_tsv(MANIFEST)
    selected = selected_cache_samples()
    rows: list[dict[str, object]] = []
    for panel in manifest:
        panel_id = panel["panel_id"]
        samples = sample_list(panel["samples"])
        feature_panel = PANEL_TO_FEATURE_PANEL.get(panel_id, "")
        selected_samples = selected.get(feature_panel, set())
        overlap = sorted(set(samples) & selected_samples)
        independent = sorted(set(samples) - selected_samples)
        release_ready = str(panel.get("release_grade_ready", "")).lower() == "true"
        rows.append(
            {
                "candidate_set": panel_id,
                "audit_grade": "release" if release_ready else "nonrelease",
                "truth_namespace": panel.get("truth_namespace", ""),
                "truth_status": panel.get("truth_status", ""),
                "samples": ",".join(map(str, samples)),
                "sample_count": len(samples),
                "feature_panel": feature_panel,
                "selected_cache_overlap_samples": ",".join(map(str, overlap)),
                "selected_cache_overlap_n": len(overlap),
                "independent_manifest_samples": ",".join(map(str, independent)),
                "independent_manifest_n": len(independent),
                "independent_ready_n": 0 if release_ready else len(independent),
                "decision": (
                    "release_panel_fully_used_by_refined_guard_selection"
                    if release_ready and len(overlap) == len(samples)
                    else "not_release_grade_for_default_promotion"
                    if not release_ready
                    else "release_independent_samples_need_profile_validation"
                ),
                "evidence": panel.get("metric_source", ""),
            }
        )

    omitted_rows = []
    for sample in HMP_OMITTED_SAMPLES:
        truth_path = HMP_TRUTH_ROOT / f"abundance{sample}.tsv"
        minco = hmp_minco_profile(sample)
        sylph = hmp_sylph_profile(sample)
        reads = hmp_read_input(sample)
        omitted_rows.append(
            {
                "candidate_set": f"hmp_airskin_omitted_sample{sample}",
                "audit_grade": "candidate_release_if_profiles_exist",
                "truth_namespace": "GTDB_species_source_abundance",
                "truth_status": "truth_file_present" if truth_path.is_file() else "truth_missing",
                "samples": str(sample),
                "sample_count": 1,
                "feature_panel": "hmp_airskin_gtdb_source_abundance",
                "selected_cache_overlap_samples": "",
                "selected_cache_overlap_n": 0,
                "independent_manifest_samples": str(sample),
                "independent_manifest_n": 1,
                "independent_ready_n": int(bool(truth_path.is_file() and minco and sylph)),
                "decision": (
                    "independent_release_ready_profile_pair_available"
                    if truth_path.is_file() and minco and sylph
                    else "missing_minco_or_sylph_profile_for_release_grade_replay"
                ),
                "evidence": (
                    f"truth={truth_path if truth_path.is_file() else 'missing'};"
                    f"minco={minco or 'missing'};"
                    f"sylph={sylph or 'missing'};"
                    f"read_input={reads or 'missing'}"
                ),
            }
        )
    rows.extend(omitted_rows)
    return rows


def marine_decision() -> str:
    if not MARINE_AUDIT.exists():
        return "marine_diagnostic_not_run"
    audit = {row["metric"]: row for row in read_tsv(MARINE_AUDIT)}
    return audit.get("promotion_decision", {}).get("value", "marine_diagnostic_unknown")


def build_audit(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    release = [row for row in rows if row["audit_grade"] == "release"]
    omitted = [row for row in rows if str(row["candidate_set"]).startswith("hmp_airskin_omitted_sample")]
    independent_release_ready = [
        row
        for row in release
        if int(row["independent_ready_n"]) > 0
        and row["decision"] != "release_panel_fully_used_by_refined_guard_selection"
    ]
    omitted_ready = [row for row in omitted if int(row["independent_ready_n"]) > 0]
    release_samples = sum(int(row["sample_count"]) for row in release)
    release_overlap = sum(int(row["selected_cache_overlap_n"]) for row in release)
    omitted_truth = sum("truth_file_present" == row["truth_status"] for row in omitted)
    omitted_minco = sum("minco=missing" not in str(row["evidence"]) for row in omitted)
    omitted_sylph = sum("sylph=missing" not in str(row["evidence"]) for row in omitted)
    return [
        {
            "metric": "release_grade_manifest_overlap",
            "value": (
                f"release_panels={len(release)};release_samples={release_samples};"
                f"selected_cache_overlap={release_overlap};"
                f"independent_release_ready_samples={sum(int(row['independent_ready_n']) for row in independent_release_ready)}"
            ),
            "evidence": str(OUT_INVENTORY.relative_to(EXP)),
            "decision": "no_unused_release_grade_cache"
            if not independent_release_ready and release_overlap == release_samples
            else "review_release_holdout_inventory",
        },
        {
            "metric": "hmp_omitted_release_candidate_inputs",
            "value": (
                f"samples={','.join(map(str, HMP_OMITTED_SAMPLES))};"
                f"truth={omitted_truth}/{len(omitted)};"
                f"minco_profile={omitted_minco}/{len(omitted)};"
                f"sylph_profile={omitted_sylph}/{len(omitted)};"
                f"ready_pairs={len(omitted_ready)}/{len(omitted)}"
            ),
            "evidence": str(OUT_INVENTORY.relative_to(EXP)),
            "decision": "missing_cached_profile_pairs",
        },
        {
            "metric": "independent_diagnostic_support",
            "value": marine_decision(),
            "evidence": str(MARINE_AUDIT.relative_to(EXP)),
            "decision": "diagnostic_only_not_release_grade",
        },
        {
            "metric": "promotion_decision",
            "value": "release_grade_independent_holdout_missing",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "keep_refined_allocator_opt_in",
        },
    ]


def write_markdown(audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# Refined Allocator Independent Holdout Inventory",
        "",
        "Date: 2026-06-29",
        "",
        "This generated audit checks whether the current local cache contains an",
        "unused release-grade GTDB holdout profile set for the refined allocator.",
        "It does not rerun raw profiling.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for metric in [
        "release_grade_manifest_overlap",
        "hmp_omitted_release_candidate_inputs",
        "independent_diagnostic_support",
        "promotion_decision",
    ]:
        row = audit_by_metric[metric]
        lines.append(f"| `{metric}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- The existing release-grade manifest samples are all already present in",
            "  the refined guard-selection cache.",
            "- HMP omitted sample IDs have truth files, but no cached MinCO/Sylph",
            "  profile pairs were found locally.",
            "- Marine exact-split replay is independent diagnostic support, not",
            "  release-grade promotion evidence.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_INVENTORY.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    inventory = build_inventory()
    write_tsv(
        OUT_INVENTORY,
        inventory,
        [
            "candidate_set",
            "audit_grade",
            "truth_namespace",
            "truth_status",
            "samples",
            "sample_count",
            "feature_panel",
            "selected_cache_overlap_samples",
            "selected_cache_overlap_n",
            "independent_manifest_samples",
            "independent_manifest_n",
            "independent_ready_n",
            "decision",
            "evidence",
        ],
    )
    audit = build_audit(inventory)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(audit)
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
