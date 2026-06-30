#!/usr/bin/env python3
"""Write a compact decision matrix for abundance strategy candidates."""

from __future__ import annotations

import csv
from pathlib import Path


EXP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXP_DIR / "results"

SOURCE_FILES = {
    "candidate_default": RESULTS_DIR / "candidate_default_decision_audit.tsv",
    "candidate_wrapper": RESULTS_DIR / "cross_panel_candidate_abundance_wrapper_audit.tsv",
    "oracle": RESULTS_DIR / "abundance_oracle_bounds_audit.tsv",
    "fixed_sweep": RESULTS_DIR / "abundance_variant_safety_audit.tsv",
    "adaptive_switch": RESULTS_DIR / "adaptive_abundance_switch_audit.tsv",
    "supervised_lopo": RESULTS_DIR / "supervised_abundance_calibrator_audit.tsv",
    "supervised_external": RESULTS_DIR / "supervised_abundance_external_exactsplit_audit.tsv",
    "edge_em": RESULTS_DIR / "edge_em_cross_domain_policy_summary.tsv",
    "blend": RESULTS_DIR / "candidate_preset_genus_xny_blend_audit.tsv",
    "blend_guard": RESULTS_DIR / "candidate_preset_genus_xny_blend_guard_audit.tsv",
    "alpha_sweep": RESULTS_DIR / "candidate_preset_genus_xny_alpha_sweep_audit.tsv",
    "alpha_overall": RESULTS_DIR / "candidate_preset_genus_xny_alpha_sweep_overall.tsv",
    "candidate_callset_oracle": RESULTS_DIR / "candidate_callset_oracle_feasibility_audit.tsv",
    "candidate_callset_summary": RESULTS_DIR / "candidate_callset_oracle_feasibility_summary.tsv",
    "selected_mass_transform": RESULTS_DIR / "selected_call_mass_transform_audit.tsv",
    "selected_mass_transform_guard": RESULTS_DIR / "selected_call_mass_transform_guard_audit.tsv",
    "selected_feature_allocator": RESULTS_DIR / "selected_call_feature_allocator_audit.tsv",
    "selected_feature_allocator_guard": RESULTS_DIR / "selected_call_feature_allocator_guard_audit.tsv",
    "feature_allocator_wrapper_parity": RESULTS_DIR / "feature_allocator_wrapper_parity_audit.tsv",
    "feature_allocator_external": RESULTS_DIR / "feature_allocator_external_exactsplit_audit.tsv",
    "feature_allocator_release_candidate": RESULTS_DIR / "feature_allocator_release_candidate_audit.tsv",
    "feature_allocator_refined_guard": RESULTS_DIR / "feature_allocator_combined_guard_audit.tsv",
    "feature_allocator_refined_parity": RESULTS_DIR / "feature_allocator_refined_guard_parity_audit.tsv",
    "feature_allocator_refined_score": RESULTS_DIR / "feature_allocator_refined_score_replay_audit.tsv",
    "feature_allocator_refined_marine": RESULTS_DIR / "feature_allocator_refined_marine_diagnostic_audit.tsv",
    "feature_allocator_refined_holdout_inventory": RESULTS_DIR
    / "feature_allocator_refined_independent_holdout_inventory_audit.tsv",
    "feature_allocator_refined_recovery": RESULTS_DIR
    / "feature_allocator_refined_holdout_recovery_audit.tsv",
    "cami3_source_profile_extension": RESULTS_DIR
    / "cami3_source_profile_binomial_extension_audit.tsv",
    "cami3_source_readmap_recovery_inputs": RESULTS_DIR
    / "cami3_source_readmap_recovery_inputs_audit.tsv",
}

OUT_TSV = RESULTS_DIR / "abundance_strategy_decision_matrix.tsv"
OUT_MD = EXP_DIR / "ABUNDANCE_STRATEGY_DECISION.md"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_metric_table(path: Path) -> dict[str, dict[str, str]]:
    return {row["metric"]: row for row in read_tsv(path)}


def metric_value(metrics: dict[str, dict[str, str]], metric: str, default: str = "NA") -> str:
    row = metrics.get(metric)
    if not row:
        return default
    return row.get("value", default)


def metric_field(
    metrics: dict[str, dict[str, str]],
    metric: str,
    field: str,
    default: str = "NA",
) -> str:
    row = metrics.get(metric)
    if not row:
        return default
    return row.get(field, default)


def rel(path: Path) -> str:
    return str(path.relative_to(EXP_DIR))


def edge_row(rows: list[dict[str, str]], comparison: str) -> dict[str, str]:
    for row in rows:
        if row.get("comparison") == comparison:
            return row
    return {}


def write_tsv(rows: list[dict[str, str]]) -> None:
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "strategy_family",
        "scope",
        "best_signal",
        "failure_mode",
        "decision",
        "source",
    ]
    with OUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(rows: list[dict[str, str]]) -> None:
    lines = [
        "# Abundance Strategy Decision Matrix",
        "",
        "Date: 2026-06-29",
        "",
        "This tracked note is generated from cached summary TSVs. It records which",
        "abundance-related strategies are part of the selected MinCO default and",
        "which ones remain diagnostic only. It does not rerun raw profiling jobs.",
        "",
        "## Source Tables",
        "",
    ]
    for key in sorted(SOURCE_FILES):
        lines.append(f"- `{rel(SOURCE_FILES[key])}`")

    lines.extend(
        [
            "",
            "## Decision Matrix",
            "",
            "| Strategy family | Scope | Best signal | Failure mode | Decision |",
            "|---|---|---|---|---|",
        ]
    )
    for row in rows:
        lines.append(
            "| {strategy_family} | {scope} | {best_signal} | {failure_mode} | {decision} |".format(
                **row
            )
        )

    lines.extend(
        [
            "",
            "## Current Default",
            "",
            "- Keep `scripts/minco_profile_default.py` default preset as `candidate`.",
            "- Candidate preset components remain: `universal-auto-exact`,",
            "  `emitted-ani90-xny100-br01-af70` rescue,",
            "  `accession-ani90-xny100-br01-af70` surface, and",
            "  `normalized-depth-alpha2` candidate abundance.",
            "- Do not promote any additional abundance allocator yet. Every tested",
            "  nonzero allocator family either has sample regressions, weak external",
            "  validation, or a call-level F1 cost.",
            "- This is a best-current-MinCO default, not a claim that MinCO broadly",
            "  beats Sylph on abundance.",
            "",
            "## Next Abundance Work",
            "",
            "- Prioritize matched-call mass allocation, because that is the largest",
            "  known gap on most cached panels and the selected-call oracle can",
            "  close the selected-candidate L1 gap on 31/32 cached samples.",
            "- Pair allocation work with targeted call recovery for the remaining",
            "  sample where the selected call set is still insufficient.",
            "- Avoid promoting panel-mean-only improvements unless they are also",
            "  sample-safe on independent holdouts.",
            "- The selected-call mass-transform sweep confirms that small",
            "  unguarded base-row transforms are not sample-safe, even when panel",
            "  means improve.",
            "- The selected-call guard audit found in-panel sample-safe rules, but",
            "  leave-one-panel-out still regressed held-out samples, so guarded",
            "  transforms also remain diagnostic.",
            "- The selected-call feature allocator sweep found useful",
            "  panel-level signals from existing depth features, but no strict",
            "  sample-safe variant, so it is not a default change.",
            "- The corrected feature-allocator guard audit excludes cached",
            "  F1/L1/Pearson score columns and is the first allocator branch to",
            "  pass leave-one-panel-out without held-out sample regressions.",
            "- The implemented guarded feature allocator switch reproduces the",
            "  sample-safe cached signal when accession taxmap labels are used",
            "  for genus grouping. The stricter release-candidate audit improves",
            "  34/38 combined samples with unchanged F1, but two external",
            "  diagnostic samples still have small L1 regressions.",
            "- The refined `guarded-genus-hit-breadth-a002-xny230` switch adds an",
            "  output-only split-support guard that removes those cached external",
            "  regressions while preserving 97.7% of the selected-panel gain. Its",
            "  wrapper guard parity and fixed-call score replay are validated, but",
            "  it still needs release-grade independent holdout evidence before",
            "  default promotion. A marine exact-split diagnostic replay is",
            "  supportive, but its GTDB transfer truth is not release-grade.",
            "- The independent-holdout inventory confirms the current cache has",
            "  no unused release-grade profile set for this allocator: all",
            "  release-grade manifest samples are already in the guard-selection",
            "  cache, and omitted HMP IDs lack cached MinCO/Sylph profile pairs.",
            "- The recovery plan ranks CAMI3 source-readmap samples3-5 as the",
            "  best local route: selected-default replay profiles, raw tables,",
            "  and baseline profiles are present for 3/3 samples; the remaining",
            "  blocker is source-readmap truth cache.",
            "- A source-profile fallback using local CAMI taxonomic profiles was",
            "  tested for samples3-5 and rejected as a release substitute because",
            "  in-scope GTDB truth mapping remains below threshold.",
            "- The exact source-readmap recovery audit shows remote archive URLs",
            "  are present in the local manifest and rescoring artifacts are ready,",
            "  but the extracted source-readmap truth files for samples3-5 are",
            "  missing. Anonymous reads are optional for this rescore because the",
            "  MinCO and Sylph profiling artifacts already exist.",
            "- A post-recovery CAMI3 extension scorer is now wired into the",
            "  evidence path. It is currently blocked only by `reads_mapping.tsv.gz`",
            "  files and will score selected-default MinCO, the refined allocator,",
            "  and Sylph once those files are restored.",
            "",
            f"Machine-readable matrix: `{rel(OUT_TSV)}`",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def build_rows() -> list[dict[str, str]]:
    candidate_default = read_metric_table(SOURCE_FILES["candidate_default"])
    candidate_wrapper = read_metric_table(SOURCE_FILES["candidate_wrapper"])
    oracle = read_metric_table(SOURCE_FILES["oracle"])
    fixed_sweep = read_metric_table(SOURCE_FILES["fixed_sweep"])
    adaptive_switch = read_metric_table(SOURCE_FILES["adaptive_switch"])
    supervised_lopo = read_metric_table(SOURCE_FILES["supervised_lopo"])
    supervised_external = read_metric_table(SOURCE_FILES["supervised_external"])
    blend = read_metric_table(SOURCE_FILES["blend"])
    blend_guard = read_metric_table(SOURCE_FILES["blend_guard"])
    alpha_sweep = read_metric_table(SOURCE_FILES["alpha_sweep"])
    candidate_oracle = read_metric_table(SOURCE_FILES["candidate_callset_oracle"])
    selected_mass_transform = read_metric_table(SOURCE_FILES["selected_mass_transform"])
    selected_mass_transform_guard = read_metric_table(
        SOURCE_FILES["selected_mass_transform_guard"]
    )
    selected_feature_allocator = read_metric_table(
        SOURCE_FILES["selected_feature_allocator"]
    )
    selected_feature_allocator_guard = read_metric_table(
        SOURCE_FILES["selected_feature_allocator_guard"]
    )
    feature_allocator_wrapper_parity = read_metric_table(
        SOURCE_FILES["feature_allocator_wrapper_parity"]
    )
    feature_allocator_external = read_metric_table(
        SOURCE_FILES["feature_allocator_external"]
    )
    feature_allocator_release_candidate = read_metric_table(
        SOURCE_FILES["feature_allocator_release_candidate"]
    )
    feature_allocator_refined_guard = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_guard"]
    )
    feature_allocator_refined_parity = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_parity"]
    )
    feature_allocator_refined_score = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_score"]
    )
    feature_allocator_refined_marine = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_marine"]
    )
    feature_allocator_refined_holdout_inventory = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_holdout_inventory"]
    )
    feature_allocator_refined_recovery = read_metric_table(
        SOURCE_FILES["feature_allocator_refined_recovery"]
    )
    cami3_source_profile_extension = read_metric_table(
        SOURCE_FILES["cami3_source_profile_extension"]
    )
    cami3_source_readmap_recovery_inputs = read_metric_table(
        SOURCE_FILES["cami3_source_readmap_recovery_inputs"]
    )
    edge_rows = read_tsv(SOURCE_FILES["edge_em"])
    edge_vs_current = edge_row(edge_rows, "adaptive_edge_em_vs_current_best_minco")
    edge_vs_sylph = edge_row(edge_rows, "adaptive_edge_em_vs_sylph")

    rows = [
        {
            "strategy_family": "selected candidate preset",
            "scope": "default wrapper",
            "best_signal": metric_value(candidate_default, "candidate_vs_current_mean_delta"),
            "failure_mode": metric_value(candidate_default, "candidate_vs_sylph_panel_wins"),
            "decision": metric_value(candidate_default, "default_candidate_decision"),
            "source": rel(SOURCE_FILES["candidate_default"]),
        },
        {
            "strategy_family": "candidate row abundance",
            "scope": "candidate-only rows",
            "best_signal": metric_value(candidate_wrapper, "delta_vs_wrapper_zero_mass"),
            "failure_mode": metric_value(candidate_wrapper, "max_delta_vs_fixed_call_policy_sweep"),
            "decision": "keep_as_candidate_preset_component_not_general_allocator",
            "source": rel(SOURCE_FILES["candidate_wrapper"]),
        },
        {
            "strategy_family": "truth-aware allocation upper bound",
            "scope": "current calls only",
            "best_signal": "recoverable_L1_pp="
            + metric_value(oracle, "mean_current_minus_oracle_truth_renorm_L1_pp"),
            "failure_mode": "oracle_vs_sylph_panels="
            + metric_value(oracle, "oracle_truth_renorm_beats_sylph_panels"),
            "decision": metric_value(oracle, "promotion_decision"),
            "source": rel(SOURCE_FILES["oracle"]),
        },
        {
            "strategy_family": "fixed-call formula sweep",
            "scope": "32 cached scored profiles",
            "best_signal": metric_value(fixed_sweep, "best_sample_mean_variant"),
            "failure_mode": "sample_safe="
            + metric_value(fixed_sweep, "strict_sample_safe_variants"),
            "decision": metric_value(fixed_sweep, "promotion_decision"),
            "source": rel(SOURCE_FILES["fixed_sweep"]),
        },
        {
            "strategy_family": "adaptive output switch",
            "scope": "output-derived thresholds",
            "best_signal": metric_value(adaptive_switch, "best_panel_mean_switch"),
            "failure_mode": metric_field(
                adaptive_switch, "lopo_holdout_panels_improved", "evidence"
            ),
            "decision": metric_value(adaptive_switch, "promotion_decision"),
            "source": rel(SOURCE_FILES["adaptive_switch"]),
        },
        {
            "strategy_family": "supervised row calibrator",
            "scope": "leave-one-panel-out",
            "best_signal": metric_value(supervised_lopo, "best_model"),
            "failure_mode": metric_value(supervised_lopo, "promotion_decision"),
            "decision": metric_value(supervised_lopo, "promotion_decision"),
            "source": rel(SOURCE_FILES["supervised_lopo"]),
        },
        {
            "strategy_family": "supervised external stress test",
            "scope": "cached external exact-split profiles",
            "best_signal": metric_value(supervised_external, "best_external_mean_delta"),
            "failure_mode": "external_safe_blends="
            + metric_value(supervised_external, "strict_external_safe_blends"),
            "decision": metric_value(supervised_external, "promotion_decision"),
            "source": rel(SOURCE_FILES["supervised_external"]),
        },
        {
            "strategy_family": "edge-level mass redistribution",
            "scope": "three spot panels",
            "best_signal": "vs_current_L1_delta="
            + edge_vs_current.get("mean_delta_L1_pp", "NA"),
            "failure_mode": "vs_current_F1_delta="
            + edge_vs_current.get("mean_delta_F1", "NA")
            + "; vs_sylph_F1_delta="
            + edge_vs_sylph.get("mean_delta_F1", "NA"),
            "decision": edge_vs_current.get("decision", "NA"),
            "source": rel(SOURCE_FILES["edge_em"]),
        },
        {
            "strategy_family": "selected-default genus-XnY blend",
            "scope": "posthoc selected-default profiles",
            "best_signal": "mean_L1_delta_pp="
            + metric_value(blend, "mean_L1_delta_vs_candidate_preset_pp"),
            "failure_mode": "sample_direction="
            + metric_value(blend, "sample_L1_direction")
            + "; max_worse_pp="
            + metric_value(blend, "max_worse_L1_delta_vs_candidate_preset_pp"),
            "decision": metric_value(blend, "promotion_decision"),
            "source": rel(SOURCE_FILES["blend"]),
        },
        {
            "strategy_family": "selected-default blend guard",
            "scope": "output-derived guard rules",
            "best_signal": metric_value(blend_guard, "best_sample_safe_guard"),
            "failure_mode": metric_value(blend_guard, "lopo_guard_result"),
            "decision": metric_value(blend_guard, "promotion_decision"),
            "source": rel(SOURCE_FILES["blend_guard"]),
        },
        {
            "strategy_family": "selected-default alpha/cap sweep",
            "scope": "conservative nonzero blends",
            "best_signal": metric_value(alpha_sweep, "best_ranked_variant"),
            "failure_mode": "sample_safe="
            + metric_value(alpha_sweep, "sample_safe_variants"),
            "decision": metric_value(alpha_sweep, "promotion_decision"),
            "source": rel(SOURCE_FILES["alpha_sweep"]),
        },
        {
            "strategy_family": "selected call-set oracle feasibility",
            "scope": "32 selected-candidate cached samples",
            "best_signal": "allocation_only_can_close_samples="
            + metric_value(candidate_oracle, "allocation_only_can_close_samples")
            + "/"
            + metric_value(candidate_oracle, "evaluated_samples"),
            "failure_mode": "call_recovery_required_samples="
            + metric_value(candidate_oracle, "call_recovery_required_samples"),
            "decision": metric_field(
                candidate_oracle, "promotion_decision", "decision"
            ),
            "source": rel(SOURCE_FILES["candidate_callset_oracle"]),
        },
        {
            "strategy_family": "selected-call simple mass transforms",
            "scope": "base called rows; candidate-row mass preserved",
            "best_signal": metric_value(selected_mass_transform, "best_ranked_variant"),
            "failure_mode": "sample_safe="
            + metric_value(selected_mass_transform, "sample_safe_variants"),
            "decision": metric_value(selected_mass_transform, "promotion_decision"),
            "source": rel(SOURCE_FILES["selected_mass_transform"]),
        },
        {
            "strategy_family": "selected-call guarded mass transforms",
            "scope": "output-derived guards over selected-call profiles",
            "best_signal": metric_value(
                selected_mass_transform_guard, "best_sample_safe_guard"
            ),
            "failure_mode": metric_value(
                selected_mass_transform_guard, "lopo_guard_result"
            ),
            "decision": metric_value(
                selected_mass_transform_guard, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["selected_mass_transform_guard"]),
        },
        {
            "strategy_family": "selected-call feature allocators",
            "scope": "fixed calls; output depth/quality target features",
            "best_signal": metric_value(
                selected_feature_allocator, "best_ranked_variant"
            ),
            "failure_mode": "sample_safe="
            + metric_value(selected_feature_allocator, "sample_safe_variants"),
            "decision": metric_value(
                selected_feature_allocator, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["selected_feature_allocator"]),
        },
        {
            "strategy_family": "selected-call guarded feature allocators",
            "scope": "output-derived guards over feature allocators",
            "best_signal": metric_value(
                selected_feature_allocator_guard, "best_sample_safe_guard"
            ),
            "failure_mode": metric_value(
                selected_feature_allocator_guard, "lopo_guard_result"
            ),
            "decision": metric_value(
                selected_feature_allocator_guard, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["selected_feature_allocator_guard"]),
        },
        {
            "strategy_family": "implemented guarded feature allocator",
            "scope": "wrapper function parity on cached selected profiles",
            "best_signal": "mean_L1_delta_pp="
            + metric_value(feature_allocator_wrapper_parity, "mean_L1_delta_pp"),
            "failure_mode": "worsened_samples="
            + metric_value(feature_allocator_wrapper_parity, "worsened_samples")
            + "; needs independent holdout",
            "decision": metric_value(
                feature_allocator_wrapper_parity, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["feature_allocator_wrapper_parity"]),
        },
        {
            "strategy_family": "implemented guarded feature allocator external",
            "scope": "plant/strain cached exact-split diagnostic profiles",
            "best_signal": metric_value(feature_allocator_external, "candidate_effect"),
            "failure_mode": "diagnostic_nonrelease; "
            + metric_value(feature_allocator_external, "baseline_validation"),
            "decision": metric_value(feature_allocator_external, "promotion_decision"),
            "source": rel(SOURCE_FILES["feature_allocator_external"]),
        },
        {
            "strategy_family": "guarded feature allocator release candidate",
            "scope": "wrapper parity plus external exact-split",
            "best_signal": metric_value(
                feature_allocator_release_candidate, "combined_candidate_effect"
            ),
            "failure_mode": metric_value(
                feature_allocator_release_candidate, "external_regression_details"
            ),
            "decision": metric_value(
                feature_allocator_release_candidate, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["feature_allocator_release_candidate"]),
        },
        {
            "strategy_family": "refined guarded feature allocator",
            "scope": "combined output-only guard, wrapper parity, score replay",
            "best_signal": metric_value(
                feature_allocator_refined_score, "score_replay_effect"
            ),
            "failure_mode": "guard_match="
            + metric_value(feature_allocator_refined_score, "score_replay_matches_guard_estimate")
            + "; independent_holdout_missing",
            "decision": metric_value(
                feature_allocator_refined_score, "promotion_decision"
            ),
            "source": (
                rel(SOURCE_FILES["feature_allocator_refined_guard"])
                + ";"
                + rel(SOURCE_FILES["feature_allocator_refined_parity"])
                + ";"
                + rel(SOURCE_FILES["feature_allocator_refined_score"])
            ),
        },
        {
            "strategy_family": "refined guarded feature allocator marine diagnostic",
            "scope": "independent cached marine exact-split profiles; nonrelease GTDB transfer truth",
            "best_signal": metric_value(
                feature_allocator_refined_marine, "refined_allocator_marine_effect"
            ),
            "failure_mode": metric_value(
                feature_allocator_refined_marine, "diagnostic_scope"
            ),
            "decision": metric_value(
                feature_allocator_refined_marine, "promotion_decision"
            ),
            "source": rel(SOURCE_FILES["feature_allocator_refined_marine"]),
        },
        {
            "strategy_family": "refined guarded feature allocator holdout inventory",
            "scope": "release-grade manifest overlap and omitted HMP input check",
            "best_signal": metric_value(
                feature_allocator_refined_holdout_inventory,
                "release_grade_manifest_overlap",
            ),
            "failure_mode": metric_value(
                feature_allocator_refined_holdout_inventory,
                "hmp_omitted_release_candidate_inputs",
            ),
            "decision": metric_value(
                feature_allocator_refined_holdout_inventory,
                "promotion_decision",
            ),
            "source": rel(SOURCE_FILES["feature_allocator_refined_holdout_inventory"]),
        },
        {
            "strategy_family": "refined guarded feature allocator recovery plan",
            "scope": "ranked independent holdout recovery routes",
            "best_signal": metric_value(
                feature_allocator_refined_recovery,
                "completed_local_step",
            ),
            "failure_mode": metric_value(
                feature_allocator_refined_recovery,
                "remaining_blocker",
            ),
            "decision": metric_field(
                feature_allocator_refined_recovery,
                "default_decision",
                "decision",
            ),
            "source": rel(SOURCE_FILES["feature_allocator_refined_recovery"]),
        },
        {
            "strategy_family": "CAMI3 source-profile extension diagnostic",
            "scope": "samples3-5 local taxonomic-profile source rows",
            "best_signal": metric_value(
                cami3_source_profile_extension,
                "source_profile_truth_scope",
            ),
            "failure_mode": metric_value(
                cami3_source_profile_extension,
                "baseline_minco_vs_sylph",
            ),
            "decision": metric_field(
                cami3_source_profile_extension,
                "promotion_decision",
                "decision",
            ),
            "source": rel(SOURCE_FILES["cami3_source_profile_extension"]),
        },
        {
            "strategy_family": "CAMI3 source-readmap recovery inputs",
            "scope": "samples3-5 exact per-read truth recovery",
            "best_signal": metric_value(
                cami3_source_readmap_recovery_inputs,
                "rescoring_artifacts",
            ),
            "failure_mode": metric_value(
                cami3_source_readmap_recovery_inputs,
                "recovery_blockers",
            ),
            "decision": metric_field(
                cami3_source_readmap_recovery_inputs,
                "promotion_decision",
                "decision",
            ),
            "source": rel(SOURCE_FILES["cami3_source_readmap_recovery_inputs"]),
        },
    ]
    return rows


def main() -> None:
    rows = build_rows()
    write_tsv(rows)
    write_markdown(rows)


if __name__ == "__main__":
    main()
