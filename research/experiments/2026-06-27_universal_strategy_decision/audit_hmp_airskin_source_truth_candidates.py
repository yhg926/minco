#!/usr/bin/env python3
"""Audit local CAMI II HMP airskin samples for release-grade GTDB truth candidates."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

import score_hmp_gtdb_source_abundance as hmp


NOTE_DIR = Path(__file__).resolve().parent
RESULTS = NOTE_DIR / "results"
OUT = RESULTS / "hmp_airskin_source_truth_candidate_audit.tsv"

TRUTH_ROOT = hmp.TRUTH_ROOT
AIRSKIN_ROOT = Path("/tmp/cami2_hmp_airskin_20260625")
UNSEEN_ROOT = Path("/tmp/cami2_hmp_unseen_transfer_20260626")
ALT_ROOT = Path("/mnt/new3T/minco_release_holdouts_20260628")
ALT_UNSEEN_ROOT = ALT_ROOT / "cami2_hmp_unseen_transfer_20260626"
CURRENT_REFRESH = Path("/tmp/minco_current_code_hmp_refresh_20260627")
CURRENT_AIRSKIN28 = Path("/tmp/minco_current_code_hmp_airskin28_20260627")

ABUNDANCE_RE = re.compile(r"abundance([0-9]+)\.tsv$")


def abundance_samples() -> list[int]:
    samples: list[int] = []
    for path in TRUTH_ROOT.glob("abundance*.tsv"):
        match = ABUNDANCE_RE.search(path.name)
        if match:
            samples.append(int(match.group(1)))
    return sorted(set(samples))


def first_existing(paths: list[Path]) -> str:
    for path in paths:
        if path.exists():
            return str(path)
    return ""


def first_existing_min_size(paths: list[Path], min_bytes: int) -> str:
    for path in paths:
        if path.exists() and path.stat().st_size >= min_bytes:
            return str(path)
    return ""


def sample_paths(sample: int) -> dict[str, str]:
    current_dynamic = Path(f"/tmp/minco_current_code_hmp_airskin{sample}_20260627")
    alt_current_dynamic = ALT_ROOT / f"minco_current_code_hmp_airskin{sample}_20260627"
    return {
        "local_fastq": first_existing_min_size(
            [
                AIRSKIN_ROOT / f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
                UNSEEN_ROOT / f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
                ALT_UNSEEN_ROOT / f"reads/airskinurogenital_sample{sample}.nonzero.fastq.gz",
            ],
            min_bytes=1024,
        ),
        "local_bam_list": first_existing(
            [
                AIRSKIN_ROOT / f"nonzero_bam_paths_sample{sample}.txt",
                UNSEEN_ROOT / f"nonzero_bam_paths_sample{sample}.txt",
                ALT_UNSEEN_ROOT / f"nonzero_bam_paths_sample{sample}.txt",
            ]
        ),
        "minco_unique_cached": first_existing(
            [
                current_dynamic / "work/minco.best_diff_unique.unfiltered.tsv",
                alt_current_dynamic / "work/minco.best_diff_unique.unfiltered.tsv",
                UNSEEN_ROOT / f"run/minco_sample{sample}_unique_zip_unfiltered.tsv",
                AIRSKIN_ROOT / f"run/minco_sample{sample}_unique_zip_unfiltered.tsv",
            ]
        ),
        "minco_split_cached": first_existing(
            [
                current_dynamic / "work/minco.best_diff_split.unfiltered.tsv",
                alt_current_dynamic / "work/minco.best_diff_split.unfiltered.tsv",
                UNSEEN_ROOT / f"run/minco_sample{sample}_split_zip_unfiltered.tsv",
                AIRSKIN_ROOT / f"run/minco_sample{sample}_split_zip_unfiltered.tsv",
            ]
        ),
        "minco_current_exact_profile": first_existing(
            [
                current_dynamic / f"minco_sample{sample}_current_default.tsv",
                alt_current_dynamic / f"minco_sample{sample}_current_default.tsv",
                CURRENT_REFRESH / f"minco_sample{sample}_current_exact.tsv",
                CURRENT_AIRSKIN28 / f"minco_sample{sample}_current_default.tsv",
                AIRSKIN_ROOT / f"run/minco_sample{sample}_current_exact.tsv",
            ]
        ),
        "sylph_r226_profile": first_existing(
            [
                UNSEEN_ROOT / f"run/sylph_sample{sample}/profile.tsv",
                ALT_UNSEEN_ROOT / f"run/sylph_sample{sample}/profile.tsv",
                AIRSKIN_ROOT / f"run/sylph_sample{sample}/profile.tsv",
            ]
        ),
        "sylph_r232_profile": first_existing(
            [
                UNSEEN_ROOT / f"run_sylph_r232/sylph_sample{sample}/profile.tsv",
                ALT_UNSEEN_ROOT / f"run_sylph_r232/sylph_sample{sample}/profile.tsv",
                AIRSKIN_ROOT / f"run_sylph_r232/sylph_sample{sample}/profile.tsv",
                AIRSKIN_ROOT / f"run_sylph_r232/sylph_sample{sample}/profile.chunked.tsv",
            ]
        ),
        "sylph_sketch": first_existing(
            [
                UNSEEN_ROOT
                / f"run/sylph_sample{sample}/airskinurogenital_sample{sample}.nonzero.fastq.gz.sylsp",
                ALT_UNSEEN_ROOT
                / f"run_sylph_r232/sylph_sample{sample}/airskinurogenital_sample{sample}.nonzero.fastq.gz.sylsp",
                AIRSKIN_ROOT
                / f"run/sylph_sample{sample}/airskinurogenital_sample{sample}.nonzero.fastq.gz.sylsp",
            ]
        ),
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    by_accession, by_core = hmp.truth.load_gtdb_metadata(hmp.truth.GTDB_METADATA)
    source_accessions = hmp.load_source_accessions(TRUTH_ROOT / "genome_to_id.tsv")

    rows: list[dict[str, object]] = []
    for sample in abundance_samples():
        hmp.SAMPLES.setdefault(sample, {"label": f"airskin_sample{sample}"})
        sample_truth, quality, _source_rows = hmp.build_source_truth(
            sample,
            source_accessions,
            by_accession,
            by_core,
        )
        paths = sample_paths(sample)
        positive_truth_species = int(sample_truth.shape[0])
        has_current_inputs = bool(paths["minco_unique_cached"] and paths["minco_split_cached"])
        has_same_release_outputs = bool(paths["minco_current_exact_profile"] and paths["sylph_r232_profile"])
        can_profile_without_download = bool(paths["local_fastq"] and paths["sylph_sketch"])
        can_restore_from_bam_archive = bool(paths["local_bam_list"])
        rows.append(
            {
                "sample": sample,
                "sample_label": f"airskin_sample{sample}",
                "source_genomes_positive": quality["source_genomes_positive"],
                "source_genomes_mapped": quality["source_genomes_mapped"],
                "truth_gtdb_species": positive_truth_species,
                "truth_mass_total": quality["truth_mass_total"],
                "truth_mass_mapped": quality["truth_mass_mapped"],
                "truth_mass_mapped_pct": quality["truth_mass_mapped_pct"],
                "truth_mass_exact_accession": quality["truth_mass_exact_accession"],
                "truth_mass_unique_assembly_core": quality["truth_mass_unique_assembly_core"],
                "truth_mass_ambiguous_assembly_core": quality["truth_mass_ambiguous_assembly_core"],
                "truth_mass_unmapped": quality["truth_mass_unmapped"],
                "release_truth_quality": quality["truth_mass_mapped_pct"] >= 95.0,
                "has_current_minco_inputs": has_current_inputs,
                "has_same_release_outputs": has_same_release_outputs,
                "can_profile_without_download": can_profile_without_download,
                "can_restore_from_bam_archive": can_restore_from_bam_archive,
                "local_fastq": paths["local_fastq"],
                "local_bam_list": paths["local_bam_list"],
                "minco_unique_cached": paths["minco_unique_cached"],
                "minco_split_cached": paths["minco_split_cached"],
                "minco_current_exact_profile": paths["minco_current_exact_profile"],
                "sylph_r226_profile": paths["sylph_r226_profile"],
                "sylph_r232_profile": paths["sylph_r232_profile"],
                "sylph_sketch": paths["sylph_sketch"],
            }
        )

    out = pd.DataFrame(rows).sort_values(
        [
            "release_truth_quality",
            "has_same_release_outputs",
            "has_current_minco_inputs",
            "can_profile_without_download",
            "can_restore_from_bam_archive",
            "truth_mass_mapped_pct",
        ],
        ascending=[False, False, False, False, False, False],
        kind="mergesort",
    )
    out.to_csv(OUT, sep="\t", index=False)

    summary_rows = [
        {
            "metric": "samples_total",
            "value": int(out.shape[0]),
        },
        {
            "metric": "release_truth_quality_samples",
            "value": int(out["release_truth_quality"].sum()),
        },
        {
            "metric": "same_release_outputs_available",
            "value": int(out["has_same_release_outputs"].sum()),
        },
        {
            "metric": "current_minco_inputs_available",
            "value": int(out["has_current_minco_inputs"].sum()),
        },
        {
            "metric": "can_profile_without_download",
            "value": int(out["can_profile_without_download"].sum()),
        },
        {
            "metric": "can_restore_from_bam_archive",
            "value": int(out["can_restore_from_bam_archive"].sum()),
        },
    ]
    pd.DataFrame(summary_rows).to_csv(
        RESULTS / "hmp_airskin_source_truth_candidate_summary.tsv",
        sep="\t",
        index=False,
    )
    print(out.head(12).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
