#!/usr/bin/env python3
"""Score CAMI3 samples3-5 after source-readmap recovery.

This helper is intentionally no-op-safe: while samples3-5 readmap files are
missing it writes an audit that records the missing inputs and exits
successfully.  Once `reads_mapping.tsv.gz` is restored for those samples, the
same command builds GTDB source-readmap truth, applies the exact-binomial
fallback only if the policy checks pass, scores the selected candidate MinCO
default, scores the opt-in refined abundance allocator on the same call set,
and scores Sylph in the same GTDB species namespace.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd

EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
sys.path.insert(0, str(REPO_ROOT))

from scripts import minco_profile_calibrated as wrapper

import audit_cami3_source_readmap_recovery_inputs as recovery
import score_cami3_gtdb_source_readmap as source_score
import score_cami3_gtdb_taxid_transfer as taxid_score
import score_cami3_source_profile_binomial_extension as source_profile_ext


RESULTS = EXP / "results"

SAMPLES = [3, 4, 5]
RELEASE_THRESHOLD_PCT = 95.0
METHOD_MINCO_BASE = "minco_candidate_source_readmap_extension"
METHOD_MINCO_REFINED = "minco_refined_allocator_source_readmap_extension"
METHOD_SYLPH = "sylph_source_readmap_extension"
ALLOCATOR_SWITCH = wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230

OUT_INPUTS = RESULTS / "cami3_source_readmap_extension_after_recovery_inputs.tsv"
OUT_TRUTH = RESULTS / "cami3_source_readmap_extension_after_recovery_truth.tsv"
OUT_SOURCE_ROWS = RESULTS / "cami3_source_readmap_extension_after_recovery_source_genomes.tsv"
OUT_FALLBACK_DETAIL = RESULTS / "cami3_source_readmap_extension_after_recovery_fallback_detail.tsv"
OUT_QUALITY = RESULTS / "cami3_source_readmap_extension_after_recovery_quality.tsv"
OUT_SCORES = RESULTS / "cami3_source_readmap_extension_after_recovery_scores.tsv"
OUT_SUMMARY = RESULTS / "cami3_source_readmap_extension_after_recovery_summary.tsv"
OUT_DELTA = RESULTS / "cami3_source_readmap_extension_after_recovery_delta.tsv"
OUT_AUDIT = RESULTS / "cami3_source_readmap_extension_after_recovery_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_READMAP_EXTENSION_AFTER_RECOVERY.md"


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


def first_existing(paths: Iterable[Path]) -> tuple[bool, str]:
    paths = list(paths)
    for path in paths:
        if path.exists():
            return True, str(path)
    return False, ";".join(str(path) for path in paths)


def input_row(sample: int) -> dict[str, object]:
    readmap_ok, readmap_path = first_existing(recovery.mapping_candidates(sample))
    profile = recovery.candidate_profile(sample)
    raw_unique = recovery.raw_unique(sample)
    raw_split = recovery.raw_split(sample)
    sylph = recovery.sylph_profile(sample)
    blockers = []
    if not readmap_ok:
        blockers.append("missing_reads_mapping")
    if not profile.is_file():
        blockers.append("missing_selected_default_profile")
    if not raw_unique.is_file() or not raw_split.is_file():
        blockers.append("missing_raw_best_ref_tables")
    if not sylph.is_file():
        blockers.append("missing_sylph_profile")
    return {
        "sample": sample,
        "reads_mapping_present": readmap_ok,
        "reads_mapping_path_or_checked": readmap_path,
        "selected_default_profile_present": profile.is_file(),
        "selected_default_profile": profile,
        "raw_unique_present": raw_unique.is_file(),
        "raw_unique": raw_unique,
        "raw_split_present": raw_split.is_file(),
        "raw_split": raw_split,
        "sylph_profile_present": sylph.is_file(),
        "sylph_profile": sylph,
        "ready_to_score": not blockers,
        "blocking_reason": ",".join(blockers),
    }


def raw_tables(sample: int) -> dict[str, Path]:
    return {
        "unique": recovery.raw_unique(sample),
        "split": recovery.raw_split(sample),
    }


def sample_truth_path(row: Mapping[str, object]) -> Path:
    return Path(str(row["reads_mapping_path_or_checked"]).split(";", 1)[0])


def source_profile_by_source(sample: int) -> dict[str, dict[str, object]]:
    path = source_profile_ext.source_profile_path(sample)
    if not path.is_file():
        return {}
    return {
        str(row["source_genome_id"]): row
        for row in source_profile_ext.read_source_profile_rows(path)
    }


def prior_species_by_source(sample_truth: pd.DataFrame) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for row in sample_truth.itertuples(index=False):
        species = str(getattr(row, "gtdb_species", ""))
        for source_id in str(getattr(row, "source_genome_ids", "") or "").split(","):
            source_id = source_id.strip()
            if source_id:
                out.setdefault(source_id, set()).add(species)
    return out


def adjusted_truth_rows(
    sample: int,
    sample_truth: pd.DataFrame,
    sample_sources: pd.DataFrame,
    quality: Mapping[str, object],
    taxid_map: Mapping[str, set[str]],
    ncbi_name_map: Mapping[str, set[str]],
    gtdb_binomial_map: Mapping[str, set[str]],
) -> tuple[pd.DataFrame, dict[str, object], list[dict[str, object]]]:
    counts: dict[str, int] = {}
    sources: dict[str, set[str]] = {}
    for row in sample_truth.itertuples(index=False):
        species = str(getattr(row, "gtdb_species", ""))
        counts[species] = counts.get(species, 0) + int(round(finite(getattr(row, "read_rows", 0))))
        for source_id in str(getattr(row, "source_genome_ids", "") or "").split(","):
            source_id = source_id.strip()
            if source_id:
                sources.setdefault(species, set()).add(source_id)

    profile_rows = source_profile_by_source(sample)
    prior_species = prior_species_by_source(sample_truth)
    fallback_detail: list[dict[str, object]] = []
    fallback_added_rows = 0
    policy_added_rows = 0
    policy_rows = 0
    priority_conflicts = 0
    same_species_remainders = 0
    unique_binomial_rows = 0
    profile_scope_rows = 0
    profile_scope_mapped_rows = 0

    for source in sample_sources.to_dict(orient="records"):
        source_id = str(source.get("source_genome_id", ""))
        meta = profile_rows.get(source_id, {})
        in_scope = str(meta.get("profile_scope", "")) == "gtdb_profile_scope"
        read_rows = int(round(finite(source.get("read_rows"))))
        prior_mapped = int(round(finite(source.get("wgs_prefix_read_rows")))) + int(
            round(finite(source.get("unique_taxid_read_rows")))
        )
        unmapped = int(round(finite(source.get("unmapped_read_rows"))))
        if in_scope:
            profile_scope_rows += read_rows
            profile_scope_mapped_rows += prior_mapped
        if not in_scope or unmapped <= 0:
            continue

        species, rule, taxid_n, name_n, binomial_n = source_profile_ext.choose_resolution(
            meta,
            taxid_map,
            ncbi_name_map,
            gtdb_binomial_map,
        )
        if not species:
            continue
        source_prior_species = sorted(prior_species.get(source_id, set()))
        same_species_remainder = prior_mapped > 0 and source_prior_species == [species]
        priority_ok = (prior_mapped == 0 and unmapped > 0) or same_species_remainder
        if same_species_remainder:
            same_species_remainders += 1
        if not priority_ok:
            priority_conflicts += 1
        if rule == "unique_gtdb_binomial" and binomial_n == 1:
            unique_binomial_rows += 1
        fallback_added_rows += unmapped
        fallback_detail.append(
            {
                "sample": sample,
                "source_genome_id": source_id,
                "species_label": meta.get("species_label", ""),
                "recommended_gtdb_species": species,
                "recommended_resolution_rule": rule,
                "taxid_gtdb_species_count": taxid_n,
                "ncbi_name_gtdb_species_count": name_n,
                "gtdb_binomial_species_count": binomial_n,
                "prior_wgs_or_unique_taxid_rows": prior_mapped,
                "prior_mapped_gtdb_species": ",".join(source_prior_species),
                "fallback_added_rows": unmapped,
                "source_priority_case": "unmapped_source"
                if prior_mapped == 0
                else "same_species_remainder"
                if same_species_remainder
                else "priority_conflict",
                "source_priority_ok": str(priority_ok).lower(),
            }
        )
        if priority_ok:
            policy_rows += 1
            policy_added_rows += unmapped
            counts[species] = counts.get(species, 0) + unmapped
            sources.setdefault(species, set()).add(source_id)

    adjusted_profile_scope_mapped_rows = profile_scope_mapped_rows + policy_added_rows
    adjusted_profile_scope_pct = (
        100.0 * adjusted_profile_scope_mapped_rows / profile_scope_rows
        if profile_scope_rows
        else finite(quality.get("read_rows_mapped_pct"))
    )
    deterministic_unique_binomial = len(fallback_detail) > 0 and unique_binomial_rows == len(fallback_detail)
    conservative_priority = priority_conflicts == 0 and policy_rows == len(fallback_detail)
    coverage_ready = adjusted_profile_scope_pct >= RELEASE_THRESHOLD_PCT
    fallback_policy_accepted = deterministic_unique_binomial and conservative_priority and coverage_ready
    truth_policy_used = "exact_binomial_fallback" if fallback_policy_accepted else "direct_source_readmap"

    adjusted_rows = []
    total = sum(counts.values())
    for species, read_rows in sorted(counts.items()):
        adjusted_rows.append(
            {
                "sample": sample,
                "gtdb_species": species,
                "read_rows": read_rows,
                "truth_abundance": read_rows / total if total else 0.0,
                "source_genomes": len(sources.get(species, set())),
                "source_genome_ids": ",".join(sorted(sources.get(species, set()))),
            }
        )
    quality_update = {
        "fallback_rows_reviewed": len(fallback_detail),
        "fallback_added_rows": fallback_added_rows,
        "fallback_policy_added_rows": policy_added_rows,
        "fallback_unique_binomial_rows": unique_binomial_rows,
        "fallback_priority_conflicts": priority_conflicts,
        "fallback_same_species_remainders": same_species_remainders,
        "profile_scope_read_rows": profile_scope_rows,
        "base_profile_scope_mapped_rows": profile_scope_mapped_rows,
        "adjusted_profile_scope_mapped_rows": adjusted_profile_scope_mapped_rows,
        "adjusted_profile_scope_mapped_pct": adjusted_profile_scope_pct,
        "fallback_deterministic_unique_binomial": str(deterministic_unique_binomial).lower(),
        "fallback_conservative_priority": str(conservative_priority).lower(),
        "fallback_policy_accepted": str(fallback_policy_accepted).lower(),
        "truth_policy_used": truth_policy_used,
    }
    truth_for_scoring = pd.DataFrame(adjusted_rows) if fallback_policy_accepted else sample_truth
    return truth_for_scoring, quality_update, fallback_detail


def score_sample(
    sample: int,
    truth_one: pd.DataFrame,
    profile: pd.DataFrame,
    best_ref_species: Mapping[str, str],
    best_ref_diag: Mapping[str, object],
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    allocator_taxmap: Mapping[str, Mapping[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    call = profile["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()
    raw = pd.to_numeric(
        profile.get("calibrated_abundance_raw", profile.get("calibrated_abundance", 0.0)),
        errors="coerce",
    ).fillna(0.0).to_numpy(dtype=float)
    adjusted, allocator_details = wrapper.guarded_feature_allocator_abundance_raw(
        profile,
        call,
        raw,
        ALLOCATOR_SWITCH,
        allocator_taxmap,
    )

    base_pred, base_meta = source_profile_ext.minco_prediction_from_profile(
        profile,
        None,
        taxid_to_species,
        name_to_species,
        best_ref_species,
        by_accession,
        by_core,
    )
    refined_pred, refined_meta = source_profile_ext.minco_prediction_from_profile(
        profile,
        np.asarray(adjusted, dtype=float),
        taxid_to_species,
        name_to_species,
        best_ref_species,
        by_accession,
        by_core,
    )
    sylph_pred, sylph_meta = taxid_score.load_sylph_predictions(
        recovery.sylph_profile(sample),
        by_accession,
        by_core,
    )

    base_score = taxid_score.score_prediction(
        sample,
        METHOD_MINCO_BASE,
        base_pred,
        truth_one,
        {**best_ref_diag, **base_meta},
    )
    refined_score = taxid_score.score_prediction(
        sample,
        METHOD_MINCO_REFINED,
        refined_pred,
        truth_one,
        {**best_ref_diag, **refined_meta, **allocator_details},
    )
    sylph_score = taxid_score.score_prediction(sample, METHOD_SYLPH, sylph_pred, truth_one, sylph_meta)
    delta = {
        "sample": sample,
        "delta_refined_minus_base_F1": finite(refined_score.get("F1")) - finite(base_score.get("F1")),
        "delta_refined_minus_base_L1_union_pp": finite(refined_score.get("L1_union_pp"))
        - finite(base_score.get("L1_union_pp")),
        "delta_refined_minus_base_Pearson_union": finite(refined_score.get("Pearson_union"))
        - finite(base_score.get("Pearson_union")),
        "allocator_applied": allocator_details.get("abundance_feature_allocator_applied", False),
        "allocator_guard_passed": allocator_details.get("abundance_feature_allocator_guard_passed", False),
        "allocator_adjusted_rows_n": allocator_details.get("abundance_feature_allocator_adjusted_rows_n", 0),
        "allocator_s_xny_median": allocator_details.get("abundance_feature_allocator_s_xny_median", 0.0),
    }
    return [base_score, refined_score, sylph_score], delta


def build_missing_audit(inputs: list[dict[str, object]]) -> list[dict[str, object]]:
    readmaps = sum(bool(row["reads_mapping_present"]) for row in inputs)
    profiles = sum(bool(row["selected_default_profile_present"]) for row in inputs)
    raw = sum(bool(row["raw_unique_present"]) and bool(row["raw_split_present"]) for row in inputs)
    sylph = sum(bool(row["sylph_profile_present"]) for row in inputs)
    blockers = sorted(
        {
            blocker
            for row in inputs
            for blocker in str(row.get("blocking_reason", "")).split(",")
            if blocker
        }
    )
    return [
        {
            "metric": "post_recovery_input_status",
            "value": f"readmaps={readmaps}/3;selected_default={profiles}/3;raw_tables={raw}/3;sylph={sylph}/3",
            "evidence": str(OUT_INPUTS.relative_to(EXP)),
            "decision": "blocked_until_readmaps_restored" if blockers else "ready_to_score",
        },
        {
            "metric": "blocking_reason",
            "value": ",".join(blockers),
            "evidence": str(OUT_INPUTS.relative_to(EXP)),
            "decision": "requires_readmap_restore_or_stream_extract" if blockers else "none",
        },
        {
            "metric": "promotion_decision",
            "value": "not_scored_missing_required_inputs" if blockers else "score_ready",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "keep_refined_allocator_opt_in",
        },
    ]


def score_if_ready(inputs: list[dict[str, object]]) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list[dict[str, object]],
]:
    wgs_to_species, taxid_to_species, name_to_species, map_diag = source_score.build_transfer_maps()
    taxid_map, ncbi_name_map, gtdb_binomial_map = source_profile_ext.resolver.build_candidate_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(source_score.TAXMAP)
    allocator_taxmap = wrapper.parse_species_taxmap(source_score.TAXMAP)
    truth_rows: list[dict[str, object]] = []
    source_rows: list[dict[str, object]] = []
    fallback_rows: list[dict[str, object]] = []
    quality_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    delta_rows: list[dict[str, object]] = []

    for row in inputs:
        sample = int(row["sample"])
        sample_truth, quality, sample_sources = source_score.build_source_truth(
            sample,
            sample_truth_path(row),
            wgs_to_species,
            taxid_to_species,
        )
        truth_rows.extend(sample_truth.to_dict(orient="records"))
        source_rows.extend(sample_sources.to_dict(orient="records"))
        scoring_truth, fallback_quality, sample_fallback = adjusted_truth_rows(
            sample,
            sample_truth,
            sample_sources,
            quality,
            taxid_map,
            ncbi_name_map,
            gtdb_binomial_map,
        )
        fallback_rows.extend(sample_fallback)
        quality_rows.append({**map_diag, **quality, **fallback_quality})
        best_ref_species, best_ref_diag = source_score.best_raw_ref_species_by_taxid(
            raw_tables(sample),
            taxmap,
            by_accession,
            by_core,
        )
        quality_rows[-1].update(best_ref_diag)
        scores, delta = score_sample(
            sample,
            scoring_truth,
            pd.read_csv(recovery.candidate_profile(sample), sep="\t"),
            best_ref_species,
            best_ref_diag,
            taxid_to_species,
            name_to_species,
            by_accession,
            by_core,
            allocator_taxmap,
        )
        score_rows.extend(scores)
        delta_rows.append(delta)

    scores = pd.DataFrame(score_rows)
    summary = taxid_score.summarize(scores)
    quality = pd.DataFrame(quality_rows)
    audit = build_scored_audit(quality, scores, summary, pd.DataFrame(delta_rows))
    return (
        pd.DataFrame(truth_rows),
        pd.DataFrame(source_rows),
        pd.DataFrame(fallback_rows),
        quality,
        scores,
        summary,
        pd.DataFrame(delta_rows),
        audit,
    )


def build_scored_audit(
    quality: pd.DataFrame,
    scores: pd.DataFrame,
    summary: pd.DataFrame,
    delta: pd.DataFrame,
) -> list[dict[str, object]]:
    min_mapped = float(quality["read_rows_mapped_pct"].min()) if not quality.empty else 0.0
    fallback_accepted = int(quality["fallback_policy_accepted"].astype(str).str.lower().eq("true").sum())
    min_adjusted_profile_scope = (
        float(quality["adjusted_profile_scope_mapped_pct"].min()) if not quality.empty else 0.0
    )
    ready_samples = ",".join(
        map(
            str,
            quality.loc[quality["read_rows_mapped_pct"] >= RELEASE_THRESHOLD_PCT, "sample"].astype(int).tolist(),
        )
    )
    refined_improved = int((delta["delta_refined_minus_base_L1_union_pp"] < -1e-12).sum())
    refined_worse = int((delta["delta_refined_minus_base_L1_union_pp"] > 1e-12).sum())
    mean_delta_l1 = float(delta["delta_refined_minus_base_L1_union_pp"].mean()) if not delta.empty else 0.0
    base = scores.loc[scores["method"].eq(METHOD_MINCO_BASE)].set_index("sample")
    sylph = scores.loc[scores["method"].eq(METHOD_SYLPH)].set_index("sample")
    f1_wins = sum(
        1 for sample in SAMPLES if finite(base.loc[sample, "F1"]) > finite(sylph.loc[sample, "F1"])
    )
    l1_wins = sum(
        1
        for sample in SAMPLES
        if finite(base.loc[sample, "L1_union_pp"]) < finite(sylph.loc[sample, "L1_union_pp"])
    )
    return [
        {
            "metric": "post_recovery_input_status",
            "value": "readmaps=3/3;selected_default=3/3;raw_tables=3/3;sylph=3/3",
            "evidence": str(OUT_INPUTS.relative_to(EXP)),
            "decision": "ready_to_score",
        },
        {
            "metric": "direct_readmap_truth_quality",
            "value": f"min_mapped_read_rows_pct={min_mapped:.6f};ready_samples={ready_samples}",
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "direct_readmap_truth_reaches_threshold"
            if len(ready_samples.split(",")) == len(SAMPLES)
            else "direct_readmap_truth_needs_fallback_review",
        },
        {
            "metric": "exact_binomial_fallback_policy",
            "value": (
                f"accepted_samples={fallback_accepted}/{len(SAMPLES)};"
                f"min_adjusted_profile_scope_pct={min_adjusted_profile_scope:.6f}"
            ),
            "evidence": f"{OUT_QUALITY.relative_to(EXP)};{OUT_FALLBACK_DETAIL.relative_to(EXP)}",
            "decision": "fallback_truth_policy_pass"
            if fallback_accepted == len(SAMPLES)
            else "fallback_truth_policy_not_accepted_for_all_samples",
        },
        {
            "metric": "candidate_minco_vs_sylph",
            "value": f"minco_F1_wins={f1_wins}/3;minco_L1_wins={l1_wins}/3",
            "evidence": str(OUT_SCORES.relative_to(EXP)),
            "decision": "independent_holdout_comparison",
        },
        {
            "metric": "refined_allocator_effect",
            "value": (
                f"improved={refined_improved};worsened={refined_worse};"
                f"mean_L1_delta_pp={mean_delta_l1:.9f}"
            ),
            "evidence": str(OUT_DELTA.relative_to(EXP)),
            "decision": "candidate_holdout_support"
            if refined_worse == 0 and refined_improved > 0
            else "do_not_promote_from_extension_result",
        },
        {
            "metric": "promotion_decision",
            "value": "post_recovery_extension_scored",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "review_holdout_result_before_default_change",
        },
    ]


def write_markdown(
    inputs: list[dict[str, object]],
    quality: pd.DataFrame,
    summary: pd.DataFrame,
    delta: pd.DataFrame,
    audit: list[dict[str, object]],
) -> None:
    lines = [
        "# CAMI3 Source-Readmap Extension After Recovery",
        "",
        "Date: 2026-06-29",
        "",
        "This helper records whether CAMI3 samples3-5 are ready for the",
        "source-readmap release-extension score. It does not download data. When",
        "the required `reads_mapping.tsv.gz` files are present, it scores the",
        "selected candidate MinCO default, the opt-in refined abundance allocator,",
        "and Sylph in the same GTDB species namespace. The exact-binomial",
        "fallback truth adjustment is used only if all policy checks pass.",
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
            "## Inputs",
            "",
            "| Sample | Readmap | Selected Default | Raw Tables | Sylph | Ready | Blocker |",
            "|---:|---|---|---|---|---|---|",
        ]
    )
    for row in inputs:
        raw_ok = bool(row["raw_unique_present"]) and bool(row["raw_split_present"])
        lines.append(
            "| {sample} | {readmap} | {profile} | {raw} | {sylph} | {ready} | {blocker} |".format(
                sample=row["sample"],
                readmap=row["reads_mapping_present"],
                profile=row["selected_default_profile_present"],
                raw=raw_ok,
                sylph=row["sylph_profile_present"],
                ready=row["ready_to_score"],
                blocker=row["blocking_reason"],
            )
        )
    if not quality.empty:
        lines.extend(
            [
                "",
                "## Truth Quality",
                "",
                "| Sample | Read rows | Mapped % | GTDB species |",
                "|---:|---:|---:|---:|",
            ]
        )
        for row in quality.itertuples(index=False):
            lines.append(
                "| {sample} | {rows} | {mapped:.6f} | {species} |".format(
                    sample=int(getattr(row, "sample")),
                    rows=int(getattr(row, "read_rows_total")),
                    mapped=float(getattr(row, "read_rows_mapped_pct")),
                    species=int(getattr(row, "truth_gtdb_species")),
                )
            )
    if not summary.empty:
        lines.extend(
            [
                "",
                "## Score Summary",
                "",
                "| Method | Mean F1 | Pooled F1 | Mean L1 pp | Mean Pearson |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for row in summary.itertuples(index=False):
            lines.append(
                "| {method} | {mean_f1:.6f} | {pooled_f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                    method=getattr(row, "method"),
                    mean_f1=float(getattr(row, "mean_F1")),
                    pooled_f1=float(getattr(row, "pooled_F1")),
                    l1=float(getattr(row, "mean_L1_union_pp")),
                    pearson=float(getattr(row, "mean_Pearson_union")),
                )
            )
    if not delta.empty:
        lines.extend(
            [
                "",
                "## Refined Allocator Delta",
                "",
                "| Sample | Delta F1 | Delta L1 pp | Delta Pearson | Guard passed |",
                "|---:|---:|---:|---:|---|",
            ]
        )
        for row in delta.itertuples(index=False):
            lines.append(
                "| {sample} | {f1:.6f} | {l1:.6f} | {pearson:.6f} | {guard} |".format(
                    sample=int(getattr(row, "sample")),
                    f1=float(getattr(row, "delta_refined_minus_base_F1")),
                    l1=float(getattr(row, "delta_refined_minus_base_L1_union_pp")),
                    pearson=float(getattr(row, "delta_refined_minus_base_Pearson_union")),
                    guard=getattr(row, "allocator_guard_passed"),
                )
            )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `{OUT_INPUTS.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
        ]
    )
    if not summary.empty:
        lines.extend(
            [
                f"- `{OUT_TRUTH.relative_to(EXP)}`",
            f"- `{OUT_SOURCE_ROWS.relative_to(EXP)}`",
            f"- `{OUT_FALLBACK_DETAIL.relative_to(EXP)}`",
                f"- `{OUT_QUALITY.relative_to(EXP)}`",
                f"- `{OUT_SCORES.relative_to(EXP)}`",
                f"- `{OUT_SUMMARY.relative_to(EXP)}`",
                f"- `{OUT_DELTA.relative_to(EXP)}`",
            ]
        )
    else:
        lines.append("- Score/truth/fallback TSVs are written only after all readmaps are present.")
    lines.append("")
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    inputs = [input_row(sample) for sample in SAMPLES]
    input_fields = [
        "sample",
        "reads_mapping_present",
        "reads_mapping_path_or_checked",
        "selected_default_profile_present",
        "selected_default_profile",
        "raw_unique_present",
        "raw_unique",
        "raw_split_present",
        "raw_split",
        "sylph_profile_present",
        "sylph_profile",
        "ready_to_score",
        "blocking_reason",
    ]
    write_tsv(OUT_INPUTS, inputs, input_fields)
    if not all(bool(row["ready_to_score"]) for row in inputs):
        audit = build_missing_audit(inputs)
        write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
        write_markdown(inputs, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), audit)
        print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
        return 0

    truth_df, source_df, fallback_df, quality, scores, summary, delta, audit = score_if_ready(inputs)
    truth_df.to_csv(OUT_TRUTH, sep="\t", index=False)
    source_df.to_csv(OUT_SOURCE_ROWS, sep="\t", index=False)
    fallback_df.to_csv(OUT_FALLBACK_DETAIL, sep="\t", index=False)
    quality.to_csv(OUT_QUALITY, sep="\t", index=False)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    delta.to_csv(OUT_DELTA, sep="\t", index=False)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(inputs, quality, summary, delta, audit)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
