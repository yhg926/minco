#!/usr/bin/env python3
"""Regression checks for the CAMI3 post-recovery extension scorer."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(EXP))

import score_cami3_source_readmap_extension_after_recovery as ext  # noqa: E402


def test_missing_readmap_audit_is_clean_blocker() -> None:
    rows = [
        {
            "reads_mapping_present": False,
            "selected_default_profile_present": True,
            "raw_unique_present": True,
            "raw_split_present": True,
            "sylph_profile_present": True,
            "blocking_reason": "missing_reads_mapping",
        }
        for _sample in ext.SAMPLES
    ]

    audit = {row["metric"]: row for row in ext.build_missing_audit(rows)}

    assert audit["post_recovery_input_status"]["value"] == (
        "readmaps=0/3;selected_default=3/3;raw_tables=3/3;sylph=3/3"
    )
    assert audit["post_recovery_input_status"]["decision"] == "blocked_until_readmaps_restored"
    assert audit["blocking_reason"]["value"] == "missing_reads_mapping"
    assert audit["blocking_reason"]["decision"] == "requires_readmap_restore_or_stream_extract"
    assert audit["promotion_decision"]["decision"] == "keep_refined_allocator_opt_in"


def test_exact_binomial_fallback_accepts_unmapped_source(monkeypatch) -> None:
    monkeypatch.setattr(
        ext,
        "source_profile_by_source",
        lambda sample: {
            "ASV1": {
                "source_genome_id": "ASV1",
                "profile_scope": "gtdb_profile_scope",
                "species_taxid": "123",
                "species_label": "Example species",
            }
        },
    )
    sources = pd.DataFrame(
        [
            {
                "source_genome_id": "ASV1",
                "read_rows": 10,
                "wgs_prefix_read_rows": 0,
                "unique_taxid_read_rows": 0,
                "unmapped_read_rows": 10,
            }
        ]
    )

    truth, quality, detail = ext.adjusted_truth_rows(
        3,
        pd.DataFrame(),
        sources,
        {"read_rows_mapped_pct": 0.0},
        {},
        {},
        {"example species": {"s__Example species"}},
    )

    assert quality["truth_policy_used"] == "exact_binomial_fallback"
    assert quality["fallback_policy_accepted"] == "true"
    assert quality["fallback_deterministic_unique_binomial"] == "true"
    assert quality["fallback_conservative_priority"] == "true"
    assert quality["adjusted_profile_scope_mapped_pct"] == 100.0
    assert detail[0]["recommended_resolution_rule"] == "unique_gtdb_binomial"
    assert truth.to_dict("records") == [
        {
            "sample": 3,
            "gtdb_species": "s__Example species",
            "read_rows": 10,
            "truth_abundance": 1.0,
            "source_genomes": 1,
            "source_genome_ids": "ASV1",
        }
    ]


def test_exact_binomial_fallback_rejects_priority_conflict(monkeypatch) -> None:
    monkeypatch.setattr(
        ext,
        "source_profile_by_source",
        lambda sample: {
            "ASV1": {
                "source_genome_id": "ASV1",
                "profile_scope": "gtdb_profile_scope",
                "species_taxid": "123",
                "species_label": "Example species",
            }
        },
    )
    prior_truth = pd.DataFrame(
        [
            {
                "sample": 3,
                "gtdb_species": "s__Prior species",
                "read_rows": 5,
                "truth_abundance": 1.0,
                "source_genomes": 1,
                "source_genome_ids": "ASV1",
            }
        ]
    )
    sources = pd.DataFrame(
        [
            {
                "source_genome_id": "ASV1",
                "read_rows": 10,
                "wgs_prefix_read_rows": 5,
                "unique_taxid_read_rows": 0,
                "unmapped_read_rows": 5,
            }
        ]
    )

    truth, quality, detail = ext.adjusted_truth_rows(
        3,
        prior_truth,
        sources,
        {"read_rows_mapped_pct": 50.0},
        {},
        {},
        {"example species": {"s__Example species"}},
    )

    assert quality["truth_policy_used"] == "direct_source_readmap"
    assert quality["fallback_policy_accepted"] == "false"
    assert quality["fallback_conservative_priority"] == "false"
    assert quality["fallback_priority_conflicts"] == 1
    assert detail[0]["source_priority_case"] == "priority_conflict"
    assert detail[0]["source_priority_ok"] == "false"
    pd.testing.assert_frame_equal(truth.reset_index(drop=True), prior_truth.reset_index(drop=True))
