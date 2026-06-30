#!/usr/bin/env python3
"""Audit the clean-holdout validation bundle for the universal strategy decision."""

from __future__ import annotations

import csv
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
MANIFEST = NOTE_DIR / "holdout_bundle_manifest.tsv"
SCHEMA = NOTE_DIR / "holdout_metric_schema.tsv"
RESULTS = NOTE_DIR / "results"

REQUIRED_FIELDS = [
    "panel_id",
    "dataset_family",
    "samples",
    "truth_namespace",
    "truth_status",
    "minco_result",
    "sylph_result",
    "metric_source",
    "minco_method",
    "sylph_method",
    "release_grade_ready",
    "diagnostic_grade_ready",
]


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


def repo_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def metric_source_has_rows(metric_path: Path, minco_method: str, sylph_method: str) -> tuple[bool, str]:
    if not metric_path.exists():
        return False, "metric_source_missing"
    rows = read_tsv(metric_path)
    if not rows:
        return False, "metric_source_empty"
    fields = set(rows[0].keys())

    has_minco = any(minco_method in row.values() for row in rows)
    has_sylph = any(sylph_method in row.values() for row in rows)
    has_f1 = bool(fields & {"F1", "mean_F1", "pooled_F1"})
    has_l1 = bool(
        fields
        & {
            "L1",
            "mean_L1",
            "mean_L1_pp",
            "l1_pct_points",
            "L1_union_pp",
            "mean_L1_union_pp",
            "mean_L1_truth_only_pp",
        }
    )
    has_pearson = bool(
        fields
        & {
            "Pearson",
            "mean_Pearson",
            "pearson",
            "Pearson_union",
            "mean_Pearson_union",
            "mean_Pearson_truth_only",
        }
    )
    if not has_minco:
        return False, "minco_method_not_found"
    if not has_sylph:
        return False, "sylph_method_not_found"
    if not has_f1:
        return False, "missing_F1_metric"
    if not has_l1:
        return False, "missing_L1_metric"
    if not has_pearson:
        return False, "missing_Pearson_metric"
    return True, "ok"


def main() -> int:
    manifest = read_tsv(MANIFEST)
    schema_rows = read_tsv(SCHEMA)
    required_from_schema = [row["field"] for row in schema_rows if row.get("required", "").lower() == "true"]
    missing_schema = [field for field in REQUIRED_FIELDS if field not in required_from_schema]
    if missing_schema:
        raise SystemExit(f"schema is missing required fields: {missing_schema}")

    audit_rows: list[dict[str, object]] = []
    for row in manifest:
        missing = [field for field in REQUIRED_FIELDS if not row.get(field, "")]
        metric_path = repo_path(row.get("metric_source", ""))
        metrics_ok, metrics_status = metric_source_has_rows(
            metric_path,
            row.get("minco_method", ""),
            row.get("sylph_method", ""),
        )
        gtdb_clean = row.get("truth_namespace") in {
            "GTDB_species",
            "GTDB_species_source_abundance",
            "GTDB_species_source_readmap_binomial_fallback",
        }
        release_flag = row.get("release_grade_ready", "").lower() == "true"
        diagnostic_flag = row.get("diagnostic_grade_ready", "").lower() == "true"
        release_ready = release_flag and gtdb_clean and metrics_ok and not missing
        diagnostic_ready = diagnostic_flag and metrics_ok and not missing and not release_ready
        audit_rows.append(
            {
                "panel_id": row.get("panel_id", ""),
                "dataset_family": row.get("dataset_family", ""),
                "truth_namespace": row.get("truth_namespace", ""),
                "truth_status": row.get("truth_status", ""),
                "metric_source_exists": metric_path.exists(),
                "metrics_ok": metrics_ok,
                "metrics_status": metrics_status,
                "manifest_release_grade_ready": release_flag,
                "audited_release_grade_ready": release_ready,
                "manifest_diagnostic_grade_ready": diagnostic_flag,
                "audited_diagnostic_grade_ready": diagnostic_ready,
                "audit_grade": "release" if release_ready else "diagnostic" if diagnostic_ready else "nonrelease",
                "blocking_reason": ""
                if release_ready
                else ";".join(
                    reason
                    for reason in [
                        "missing_required_fields" if missing else "",
                        "truth_not_clean_gtdb_species" if not gtdb_clean else "",
                        metrics_status if not metrics_ok else "",
                        "manifest_release_flag_false" if not release_flag else "",
                    ]
                    if reason
                ),
                "notes": row.get("notes", ""),
            }
        )

    write_tsv(
        RESULTS / "holdout_bundle_audit.tsv",
        audit_rows,
        [
            "panel_id",
            "dataset_family",
            "truth_namespace",
            "truth_status",
            "metric_source_exists",
            "metrics_ok",
            "metrics_status",
            "manifest_release_grade_ready",
            "audited_release_grade_ready",
            "manifest_diagnostic_grade_ready",
            "audited_diagnostic_grade_ready",
            "audit_grade",
            "blocking_reason",
            "notes",
        ],
    )

    total = len(audit_rows)
    ready = sum(1 for row in audit_rows if row["audited_release_grade_ready"])
    diagnostic = sum(1 for row in audit_rows if row["audited_diagnostic_grade_ready"])
    write_tsv(
        RESULTS / "holdout_bundle_summary.tsv",
        [
            {
                "metric": "panels_total",
                "value": total,
            },
            {
                "metric": "release_grade_panels",
                "value": ready,
            },
            {
                "metric": "diagnostic_grade_panels",
                "value": diagnostic,
            },
            {
                "metric": "informative_panels",
                "value": ready + diagnostic,
            },
            {
                "metric": "release_grade_all_panels",
                "value": str(ready == total).lower(),
            },
            {
                "metric": "main_blocker",
                "value": "not all panels have clean GTDB_species truth with comparable MinCO/Sylph metrics",
            },
        ],
        ["metric", "value"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
