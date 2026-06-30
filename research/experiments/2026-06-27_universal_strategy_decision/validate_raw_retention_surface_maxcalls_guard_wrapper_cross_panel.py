#!/usr/bin/env python3
"""Replay the guarded candidate surface across cached panels.

This validates the implementation path for the combined candidate surface plus
``--candidate-surface-max-called-species 250``. It is a full wrapper-output
replay over the same cached 43-sample / 7-panel panel used by the unguarded
candidate-surface replay. The script is validation only and does not change
MinCO defaults.
"""

from __future__ import annotations

import csv
import subprocess
import time
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

import validate_raw_retention_surface_wrapper_cross_panel as base


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_raw_retention_surface_maxcalls_guard_wrapper_cross_panel_20260630")
SURFACE_SWITCH = "accession-current-or-ani93-xny650-br20"
RESCUE_SWITCH = "emitted-ani90-xny100-br01-af70"
POLICIES = ["zero", "normalized-depth-alpha2"]
MAX_CALLED_SPECIES = 250
PREFIX = "raw_retention_surface_maxcalls_guard_wrapper"


def log(message: str) -> None:
    print(message, flush=True)


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def clean_panel(panel: str) -> str:
    return panel.replace("/", "_").replace(" ", "_")


def out_path(config: Mapping[str, object], policy: str) -> Path:
    panel = clean_panel(str(config["panel"]))
    sample = int(config["sample"])
    clean_policy = policy.replace("-", "_")
    return OUT_DIR / clean_policy / (
        f"{panel}.sample{sample}.{SURFACE_SWITCH}.maxcalls{MAX_CALLED_SPECIES}.tsv"
    )


def run_wrapper(config: Mapping[str, object], policy: str) -> dict[str, object]:
    panel = str(config["panel"])
    sample = int(config["sample"])
    tables = base.table_paths(config)
    unique = Path(tables["unique"])
    split = Path(tables["split"])
    out = out_path(config, policy)
    time_log = out.with_suffix(out.suffix + ".time.log")
    base_row = {
        "panel": panel,
        "sample": sample,
        "policy": policy,
        "profile": str(out),
        "unique": str(unique),
        "split": str(split),
        "split_source": tables["split_source"],
        "has_exact": tables["has_exact"],
        "time_log": str(time_log),
    }
    if not unique.exists() or not split.exists():
        return {
            **base_row,
            "status": "missing_input",
            "command": "",
            "elapsed_s": "",
            "wall_clock": "",
            "max_rss_kb": "",
            "exit_status": "",
        }
    if out.exists():
        timing = base.parse_time_log(time_log)
        return {
            **base_row,
            **timing,
            "status": "reused",
            "command": "reuse existing guarded raw-retention surface wrapper profile",
            "elapsed_s": "0.000",
        }

    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_log),
        "python3",
        "-B",
        "scripts/minco_profile_calibrated.py",
        "--unique-table",
        str(unique),
        "--split-table",
        str(split),
        "--taxmap",
        str(base.normal_taxmap(config)),
        "--candidate-surface-taxmap",
        str(base.surface_taxmap(config)),
        *base.model_args(),
        "--scope",
        "bacteria",
        "--train-pool",
        "train12",
        "--strategy",
        "universal-auto-exact",
        "--candidate-rescue-switch",
        RESCUE_SWITCH,
        "--candidate-surface-switch",
        SURFACE_SWITCH,
        "--candidate-surface-max-called-species",
        str(MAX_CALLED_SPECIES),
        "--candidate-abundance-policy",
        policy,
        "-o",
        str(out),
    ]
    start = time.monotonic()
    subprocess.run(cmd, cwd=base.REPO, check=True)
    elapsed = time.monotonic() - start
    timing = base.parse_time_log(time_log)
    return {
        **base_row,
        **timing,
        "status": "ran",
        "command": " ".join(cmd),
        "elapsed_s": f"{elapsed:.3f}",
    }


def metadata_for_profile(panel: str, sample: int, policy: str, path: Path) -> dict[str, object]:
    row = {
        "panel": panel,
        "sample": sample,
        "policy": policy,
        "profile": str(path),
        "rows": 0,
        "called": 0,
        "candidate_rescue_added": 0,
        "candidate_surface_added": 0,
        "candidate_added": 0,
        "candidate_added_called": 0,
        "candidate_added_raw_mass": 0.0,
        "candidate_added_abundance": 0.0,
        "surface_switch_values": "",
        "abundance_policy_values": "",
        "guard_blocked_values": "",
        "guard_blocked_rows": 0,
        "surface_called_species_n_max": "",
        "surface_max_called_species_values": "",
    }
    if not path.exists():
        return row
    cols = {
        "calibrated_call",
        "calibrated_abundance",
        "calibrated_abundance_raw",
        "candidate_rescue_added",
        "candidate_surface_added",
        "candidate_surface_switch",
        "candidate_abundance_policy",
        "candidate_surface_guard_blocked",
        "candidate_surface_called_species_n",
        "candidate_surface_max_called_species",
    }
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in cols, low_memory=False)
    called = base.bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
    rescue = base.bool_series(df.get("candidate_rescue_added", pd.Series(False, index=df.index)))
    surface = base.bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    blocked = base.bool_series(
        df.get("candidate_surface_guard_blocked", pd.Series(False, index=df.index))
    )
    added = rescue | surface
    raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    called_species_n = pd.to_numeric(
        df.get("candidate_surface_called_species_n", pd.Series(dtype=float)),
        errors="coerce",
    )
    row.update(
        {
            "rows": int(len(df)),
            "called": int(called.sum()),
            "candidate_rescue_added": int(rescue.sum()),
            "candidate_surface_added": int(surface.sum()),
            "candidate_added": int(added.sum()),
            "candidate_added_called": int((added & called).sum()),
            "candidate_added_raw_mass": float(raw.loc[added].sum()),
            "candidate_added_abundance": float(abundance.loc[added].sum()),
            "surface_switch_values": ",".join(
                sorted(set(df.get("candidate_surface_switch", pd.Series(dtype=str)).dropna().astype(str)))
            ),
            "abundance_policy_values": ",".join(
                sorted(set(df.get("candidate_abundance_policy", pd.Series(dtype=str)).dropna().astype(str)))
            ),
            "guard_blocked_values": ",".join(sorted(set(blocked.astype(str)))),
            "guard_blocked_rows": int(blocked.sum()),
            "surface_called_species_n_max": (
                int(called_species_n.max()) if called_species_n.notna().any() else ""
            ),
            "surface_max_called_species_values": ",".join(
                sorted(
                    set(
                        df.get("candidate_surface_max_called_species", pd.Series(dtype=str))
                        .dropna()
                        .astype(str)
                    )
                )
            ),
        }
    )
    return row


def score_profiles(
    configs: list[dict[str, object]],
    mapper: base.cross.failure.GtdbMapper,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for index, config in enumerate(configs, start=1):
        panel = str(config["panel"])
        sample = int(config["sample"])
        log(f"[score {index}/{len(configs)}] {panel} sample{sample}")
        truth = base.cross.truth_for_config(config)
        current_pred = base.cross.current_pred_for_config(config, mapper)
        rows.append(
            base.cross.score_sample(
                panel,
                sample,
                "current_default",
                truth,
                current_pred,
                set(current_pred),
                set(),
                "current selected default profile",
            )
        )
        for policy in POLICIES:
            path = out_path(config, policy)
            if not path.exists():
                continue
            scored_config = dict(config)
            scored_config["profile"] = path
            pred = base.cross.current_pred_for_config(scored_config, mapper)
            rows.append(
                base.cross.score_sample(
                    panel,
                    sample,
                    f"wrapper_{SURFACE_SWITCH}_maxcalls{MAX_CALLED_SPECIES}_{policy}",
                    truth,
                    pred,
                    set(pred),
                    set(),
                    (
                        f"wrapper emitted {SURFACE_SWITCH}; "
                        f"max_called_species={MAX_CALLED_SPECIES}; policy={policy}"
                    ),
                )
            )
    return pd.DataFrame(rows)


def audit_rows(
    overall: pd.DataFrame,
    metadata: pd.DataFrame,
    timing: pd.DataFrame,
    sample_count: int,
    panel_count: int,
) -> list[dict[str, object]]:
    strict = overall.loc[overall["strict_pass"].astype(bool)] if not overall.empty else pd.DataFrame()
    top = overall.iloc[0].to_dict() if not overall.empty else {}
    added = metadata.loc[metadata["policy"].isin(POLICIES)] if not metadata.empty else pd.DataFrame()
    blocked_profiles = int((metadata["guard_blocked_rows"] > 0).sum()) if not metadata.empty else 0
    return [
        {
            "metric": "scope",
            "value": (
                f"samples={sample_count};panels={panel_count};"
                f"policies={','.join(POLICIES)};max_called_species={MAX_CALLED_SPECIES}"
            ),
            "evidence": "cached table-mode guarded wrapper replay",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "completed_profiles",
            "value": int(timing["status"].isin({"ran", "reused"}).sum()) if not timing.empty else 0,
            "evidence": "wrapper timing table",
            "decision": (
                "complete"
                if not timing.empty and timing["status"].isin({"ran", "reused"}).all()
                else "incomplete"
            ),
        },
        {
            "metric": "strict_pass_methods",
            "value": int(len(strict)),
            "evidence": "requires no sample or panel F1 regression versus current default",
            "decision": "candidate_found" if not strict.empty else "no_safe_wrapper_method",
        },
        {
            "metric": "best_method",
            "value": top.get("method", ""),
            "evidence": (
                f"mean_delta_F1={top.get('mean_delta_F1', '')};"
                f"min_delta_F1={top.get('min_delta_F1', '')};"
                f"sample_worsen_n={top.get('sample_worsen_n', '')};"
                f"panel_worsen_n={top.get('panel_worsen_n', '')};"
                f"mean_delta_L1={top.get('mean_delta_L1_union_pp', '')}"
            ),
            "decision": "review_candidate" if bool(top.get("strict_pass", False)) else "diagnostic_only",
        },
        {
            "metric": "candidate_added_rows",
            "value": int(added["candidate_added"].sum()) if not added.empty else 0,
            "evidence": "sum over wrapper metadata rows",
            "decision": "metadata",
        },
        {
            "metric": "guard_blocked_profiles",
            "value": blocked_profiles,
            "evidence": "profiles with candidate_surface_guard_blocked=True",
            "decision": "metadata",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "same cached-panel replay; no independent same-namespace holdout added here",
            "decision": "keep_opt_in_until_independent_validation",
        },
    ]


def main() -> int:
    configs = base.cross.configs()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    timing_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    for index, config in enumerate(configs, start=1):
        panel = str(config["panel"])
        sample = int(config["sample"])
        log(f"[{index}/{len(configs)}] {panel} sample{sample}")
        for policy in POLICIES:
            timing = run_wrapper(config, policy)
            timing_rows.append(timing)
            metadata_rows.append(
                metadata_for_profile(
                    panel,
                    sample,
                    policy,
                    Path(str(timing["profile"])),
                )
            )
    timing = pd.DataFrame(timing_rows)
    metadata = pd.DataFrame(metadata_rows)
    timing.to_csv(RESULTS / f"{PREFIX}_timing.tsv", sep="\t", index=False)
    metadata.to_csv(RESULTS / f"{PREFIX}_metadata.tsv", sep="\t", index=False)
    log(f"wrote timing and metadata for {len(metadata)} profiles")
    mapper = base.cross.failure.GtdbMapper()
    scores = score_profiles(configs, mapper)
    deltas, panels = base.cross.summarize(scores)
    overall = base.cross.overall(deltas, panels)
    audit = audit_rows(
        overall,
        metadata,
        timing,
        len(configs),
        len({str(config["panel"]) for config in configs}),
    )

    scores.to_csv(RESULTS / f"{PREFIX}_scores.tsv", sep="\t", index=False)
    deltas.to_csv(RESULTS / f"{PREFIX}_deltas.tsv", sep="\t", index=False)
    panels.to_csv(RESULTS / f"{PREFIX}_panel_summary.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / f"{PREFIX}_overall.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / f"{PREFIX}_audit.tsv",
        audit,
        ["metric", "value", "evidence", "decision"],
    )

    print(overall.to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
