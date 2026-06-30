#!/usr/bin/env python3
"""Replay the opt-in raw-retention candidate surface across cached panels.

The offline raw-candidate replay added strict raw-visible species on top of the
current selected default. This wrapper replay validates the implementation path:
candidate rescue plus the combined accession-level surface
``accession-current-or-ani93-xny650-br20``. It runs cached table-mode inputs,
then scores emitted profiles with the same GTDB panel scorers used by the
decision note.

This script is validation only. It does not change MinCO defaults.
"""

from __future__ import annotations

import csv
import math
import subprocess
import time
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

import audit_hmp_missed_truth_raw_tables as hmp_raw_audit
import sweep_raw_candidate_retention_cross_panel as cross
import validate_cross_panel_candidate_abundance_wrapper as wrapper_validation
import validate_integrated_candidate_surface_hmp as hmp_surface


EXP = Path(__file__).resolve().parent
REPO = EXP.parents[2]
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_raw_retention_surface_wrapper_cross_panel_20260630")
SURFACE_SWITCH = "accession-current-or-ani93-xny650-br20"
RESCUE_SWITCH = "emitted-ani90-xny100-br01-af70"
POLICIES = ["zero", "normalized-depth-alpha2"]
GTDB_TAXMAP = wrapper_validation.GTDB_TAXMAP
TOY_TAXMAP = wrapper_validation.TOY_TAXMAP
SURFACE_TAXMAP = hmp_surface.SURFACE_TAXMAP
MODEL_CACHE = wrapper_validation.MODEL_CACHE
TRAIN_FEATURES = wrapper_validation.TRAIN_FEATURES


def log(message: str) -> None:
    print(message, flush=True)


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
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


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def clean_panel(panel: str) -> str:
    return panel.replace("/", "_").replace(" ", "_")


def parse_time_log(path: Path) -> dict[str, object]:
    text = path.read_text() if path.exists() else ""
    out = {"wall_clock": "", "max_rss_kb": "", "exit_status": ""}
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            out["wall_clock"] = line.rsplit("):", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            out["max_rss_kb"] = line.rsplit(":", 1)[-1].strip()
        elif "Exit status" in line:
            out["exit_status"] = line.rsplit(":", 1)[-1].strip()
    return out


def model_args() -> list[str]:
    if MODEL_CACHE.exists():
        return ["--model-cache", str(MODEL_CACHE)]
    return ["--train-features", str(TRAIN_FEATURES)]


def normal_taxmap(config: Mapping[str, object]) -> Path:
    return TOY_TAXMAP if str(config["kind"]) == "toy" else GTDB_TAXMAP


def surface_taxmap(config: Mapping[str, object]) -> Path:
    return SURFACE_TAXMAP if SURFACE_TAXMAP.exists() else normal_taxmap(config)


def hmp_panel_for_paths(kind: str) -> str:
    if kind == "hmp_gastro":
        return "hmp_gastrooral_gtdb_source_abundance"
    return "hmp_airskin_gtdb_source_abundance"


def optional_path(value: object) -> Path | None:
    if value is None:
        return None
    text = str(value)
    if not text or text.lower() == "nan":
        return None
    return Path(text)


def table_paths(config: Mapping[str, object]) -> dict[str, object]:
    kind = str(config["kind"])
    sample = int(config["sample"])
    profile = Path(str(config["profile"]))
    if kind in {"hmp_airskin", "hmp_gastro"}:
        table_info = hmp_surface.discover_tables(hmp_panel_for_paths(kind), sample, profile)
        paths = {
            key: value
            for key, value in {
                "unique": table_info.get("unique"),
                "split": table_info.get("split"),
            }.items()
            if value is not None
        }
        if not paths.get("unique") or not paths.get("split"):
            paths.update(hmp_raw_audit.raw_paths_from_profile(hmp_panel_for_paths(kind), sample, profile))
    else:
        paths = {key: Path(value) for key, value in dict(config["raw_paths"]).items()}  # type: ignore[arg-type]
    unique = optional_path(paths.get("unique"))
    split_exact = optional_path(paths.get("split_exact"))
    split_raw = optional_path(paths.get("split"))
    split = split_exact if split_exact is not None and split_exact.exists() else split_raw
    return {
        "unique": unique or Path("/nonexistent/minco_missing_unique.tsv"),
        "split": split or Path("/nonexistent/minco_missing_split.tsv"),
        "split_source": "split_exact" if split_exact is not None and split == split_exact and split_exact.exists() else "split",
        "has_exact": bool(split_exact is not None and split_exact.exists()),
    }


def out_path(config: Mapping[str, object], policy: str) -> Path:
    panel = clean_panel(str(config["panel"]))
    sample = int(config["sample"])
    clean_policy = policy.replace("-", "_")
    return OUT_DIR / clean_policy / f"{panel}.sample{sample}.{SURFACE_SWITCH}.tsv"


def run_wrapper(config: Mapping[str, object], policy: str) -> dict[str, object]:
    panel = str(config["panel"])
    sample = int(config["sample"])
    tables = table_paths(config)
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
        timing = parse_time_log(time_log)
        return {
            **base_row,
            **timing,
            "status": "reused",
            "command": "reuse existing raw-retention surface wrapper profile",
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
        str(normal_taxmap(config)),
        "--candidate-surface-taxmap",
        str(surface_taxmap(config)),
        *model_args(),
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
        "--candidate-abundance-policy",
        policy,
        "-o",
        str(out),
    ]
    start = time.monotonic()
    subprocess.run(cmd, cwd=REPO, check=True)
    elapsed = time.monotonic() - start
    timing = parse_time_log(time_log)
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
    }
    df = pd.read_csv(path, sep="\t", usecols=lambda col: col in cols, low_memory=False)
    called = bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
    rescue = bool_series(df.get("candidate_rescue_added", pd.Series(False, index=df.index)))
    surface = bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    added = rescue | surface
    raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
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
        }
    )
    return row


def score_profiles(
    configs: list[dict[str, object]],
    timing: pd.DataFrame,
    mapper: cross.failure.GtdbMapper,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for config in configs:
        panel = str(config["panel"])
        sample = int(config["sample"])
        truth = cross.truth_for_config(config)
        current_pred = cross.current_pred_for_config(config, mapper)
        rows.append(
            cross.score_sample(
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
            pred = cross.current_pred_for_config(scored_config, mapper)
            rows.append(
                cross.score_sample(
                    panel,
                    sample,
                    f"wrapper_{SURFACE_SWITCH}_{policy}",
                    truth,
                    pred,
                    set(pred),
                    set(),
                    f"wrapper emitted {SURFACE_SWITCH}; policy={policy}",
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
    return [
        {
            "metric": "scope",
            "value": f"samples={sample_count};panels={panel_count};policies={','.join(POLICIES)}",
            "evidence": "cached table-mode wrapper replay",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "completed_profiles",
            "value": int(timing["status"].isin({"ran", "reused"}).sum()) if not timing.empty else 0,
            "evidence": "wrapper timing table",
            "decision": "complete" if not timing.empty and timing["status"].isin({"ran", "reused"}).all() else "incomplete",
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
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "same cached-panel replay; no independent same-namespace holdout added here",
            "decision": "keep_opt_in_until_independent_validation",
        },
    ]


def main() -> int:
    configs = cross.configs()
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
    mapper = cross.failure.GtdbMapper()
    scores = score_profiles(configs, timing, mapper)
    deltas, panels = cross.summarize(scores)
    overall = cross.overall(deltas, panels)
    audit = audit_rows(overall, metadata, timing, len(configs), len({str(c["panel"]) for c in configs}))

    timing.to_csv(RESULTS / "raw_retention_surface_wrapper_timing.tsv", sep="\t", index=False)
    metadata.to_csv(RESULTS / "raw_retention_surface_wrapper_metadata.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "raw_retention_surface_wrapper_scores.tsv", sep="\t", index=False)
    deltas.to_csv(RESULTS / "raw_retention_surface_wrapper_deltas.tsv", sep="\t", index=False)
    panels.to_csv(RESULTS / "raw_retention_surface_wrapper_panel_summary.tsv", sep="\t", index=False)
    overall.to_csv(RESULTS / "raw_retention_surface_wrapper_overall.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "raw_retention_surface_wrapper_audit.tsv",
        audit,
        ["metric", "value", "evidence", "decision"],
    )

    print(overall.to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
