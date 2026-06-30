#!/usr/bin/env python3
"""Audit an output-derived max-called-species guard for candidate-surface rows.

This is a cached replay over existing wrapper score tables. It asks whether
using the combined candidate surface only below a pre-surface call-count limit
would preserve cross-panel F1 gains while avoiding the observed high-complexity
regression samples. It does not change MinCO defaults.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import sweep_raw_candidate_retention_cross_panel as cross


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
SCORES = RESULTS / "raw_retention_surface_wrapper_scores.tsv"
SMOKE = RESULTS / "raw_retention_surface_maxcalls_guard_smoke.tsv"
SMOKE_TIMING = RESULTS / "raw_retention_surface_maxcalls_guard_smoke_timing.tsv"
SWEEP_OUT = RESULTS / "raw_retention_surface_maxcalls_guard_sweep.tsv"
AUDIT_OUT = RESULTS / "raw_retention_surface_maxcalls_guard_audit.tsv"
WRAPPER_METHOD = "wrapper_accession-current-or-ani93-xny650-br20_normalized-depth-alpha2"
THRESHOLDS = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 999999]


def write_sweep() -> pd.DataFrame:
    scores = pd.read_csv(SCORES, sep="\t")
    current = scores.loc[scores["method"] == "current_default"].copy()
    wrapper = scores.loc[scores["method"] == WRAPPER_METHOD].copy()
    current_i = current.set_index(["panel", "sample"])
    wrapper_i = wrapper.set_index(["panel", "sample"])
    rows: list[dict[str, object]] = []
    for threshold in THRESHOLDS:
        hybrid_rows: list[dict[str, object]] = []
        for index, current_row in current_i.iterrows():
            wrapper_row = wrapper_i.loc[index]
            use_wrapper = int(current_row["pred_species"]) <= threshold
            row = (wrapper_row if use_wrapper else current_row).copy()
            row["method"] = (
                f"max_base_calls_le_{threshold}" if threshold < 999999 else "combined_no_guard"
            )
            row["rule"] = f"use combined wrapper if current pred_species <= {threshold}"
            row["panel"] = index[0]
            row["sample"] = index[1]
            hybrid_rows.append(row.to_dict())
        hybrid = pd.DataFrame(hybrid_rows)
        all_scores = pd.concat([current, hybrid], ignore_index=True)
        deltas, panels = cross.summarize(all_scores)
        overall = cross.overall(deltas, panels).iloc[0].to_dict()
        overall["thr"] = threshold
        rows.append(overall)
    sweep = pd.DataFrame(rows)
    sweep.to_csv(SWEEP_OUT, sep="\t", index=False)
    return sweep


def write_audit(sweep: pd.DataFrame) -> pd.DataFrame:
    best = sweep.loc[sweep["thr"] == 250].iloc[0]
    smoke = pd.read_csv(SMOKE, sep="\t") if SMOKE.exists() else pd.DataFrame()
    timing = pd.read_csv(SMOKE_TIMING, sep="\t") if SMOKE_TIMING.exists() else pd.DataFrame()
    if not timing.empty:
        max_rss = int(pd.to_numeric(timing["max_rss_kb"], errors="coerce").max())
        median_wall = str(timing["wall_clock"].tolist()[len(timing) // 2])
        exit_statuses = ",".join(sorted(timing["exit_status"].astype(str).unique()))
    else:
        max_rss = 0
        median_wall = ""
        exit_statuses = ""
    if not smoke.empty:
        smoke_value = (
            f"samples={','.join(map(str, smoke['sample'].tolist()))};"
            f"surface_added_sum={int(smoke['surface_added'].sum())};"
            f"guard_blocked={','.join(sorted(smoke['guard_blocked'].astype(str).unique()))}"
        )
    else:
        smoke_value = "not_run"
    rows = [
        {
            "metric": "scope",
            "value": "43 samples / 7 panels offline guard replay; marine samples3-5 implementation smoke",
            "evidence": f"{SWEEP_OUT.name}; {SMOKE.name}",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "guard_rule",
            "value": "candidate surface active only when pre-surface called species count <= 250",
            "evidence": "--candidate-surface-max-called-species 250",
            "decision": "opt_in_guard_not_default",
        },
        {
            "metric": "offline_best_thr_250",
            "value": (
                f"mean_delta_F1={best['mean_delta_F1']};"
                f"min_delta_F1={best['min_delta_F1']};"
                f"sample_worsen_n={int(best['sample_worsen_n'])};"
                f"panel_worsen_n={int(best['panel_worsen_n'])};"
                f"mean_delta_L1={best['mean_delta_L1_union_pp']};"
                f"TP_delta={int(best['total_TP_delta'])};"
                f"FP_delta={int(best['total_FP_delta'])}"
            ),
            "evidence": "hybrid replay from current and combined wrapper score tables",
            "decision": "strict_pass_diagnostic_candidate",
        },
        {
            "metric": "smoke_profiles",
            "value": smoke_value,
            "evidence": "/tmp/minco_raw_retention_surface_maxcalls_guard_20260630.marine*.tsv",
            "decision": "implementation_blocks_regression_samples"
            if smoke_value != "not_run"
            else "smoke_not_run",
        },
        {
            "metric": "smoke_runtime",
            "value": f"median_wall={median_wall};max_rss_kb_max={max_rss};exit_statuses={exit_statuses}",
            "evidence": "/tmp/minco_raw_retention_surface_maxcalls_guard_20260630.marine*.time.log",
            "decision": "runtime_smoke_complete" if not timing.empty else "runtime_smoke_missing",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "same cached-panel replay plus three implementation smokes; no new independent holdout",
            "decision": "keep_opt_in_until_full_wrapper_replay_or_independent_validation",
        },
    ]
    audit = pd.DataFrame(rows)
    audit.to_csv(AUDIT_OUT, sep="\t", index=False)
    return audit


def main() -> int:
    sweep = write_sweep()
    audit = write_audit(sweep)
    print(sweep.to_string(index=False))
    print("\nAUDIT")
    print(audit.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
