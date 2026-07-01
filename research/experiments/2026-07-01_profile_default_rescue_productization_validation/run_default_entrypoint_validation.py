#!/usr/bin/env python3
"""Validate the productized default profile entry point.

This script runs ``scripts/minco_profile`` from saved unique/split profile
tables and compares:

* default candidate preset
* default candidate preset with ``--no-profile-rescue``

Large profile outputs stay under /tmp; repository results are compact TSVs.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
SWITCH_EXP = ROOT / "research/experiments/2026-07-01_profile_rescue_switch_implementation_validation"
sys.path.insert(0, str(SWITCH_EXP))

import run_switch_validation as switch  # noqa: E402


DEFAULT_RUN_ROOT = Path("/tmp/minco_profile_default_rescue_productization_validation_20260701")
METHODS = [
    ("default_candidate", []),
    ("default_candidate_no_profile_rescue", ["--no-profile-rescue"]),
]


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def profile_path(run_root: Path, sample: switch.Sample, method: str) -> Path:
    return run_root / "profiles" / f"{sample.label}.{method}.tsv"


def time_path(run_root: Path, sample: switch.Sample, method: str) -> Path:
    return run_root / "logs" / f"{sample.label}.{method}.time.log"


def log_path(run_root: Path, sample: switch.Sample, method: str) -> Path:
    return run_root / "logs" / f"{sample.label}.{method}.run.log"


def command(args: argparse.Namespace, sample: switch.Sample, method_args: list[str], out: Path) -> list[str]:
    return [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_path(args.run_root, sample, args.current_method)),
        "python3",
        "-B",
        str(ROOT / "scripts/minco_profile"),
        "--unique-table",
        str(sample.unique_table),
        "--split-table",
        str(sample.split_table),
        "--taxmap",
        str(args.taxmap),
        "--model-cache",
        str(args.model_cache),
        "--scope",
        "bacteria",
        "--report-all",
        *method_args,
        "-o",
        str(out),
    ]


def run_profile(args: argparse.Namespace, sample: switch.Sample, method: str, method_args: list[str]) -> dict[str, object]:
    out = profile_path(args.run_root, sample, method)
    args.current_method = method
    cmd = command(args, sample, method_args, out)
    row = {
        "panel": sample.panel,
        "sample": sample.sample,
        "label": sample.label,
        "method": method,
        "profile": str(out),
        "command": shlex.join(cmd),
    }
    if out.exists() and out.stat().st_size > 0 and not args.force:
        row["run_status"] = "reused_existing_profile"
        row.update(switch.parse_time_log(time_path(args.run_root, sample, method)))
        return row
    if args.skip_run:
        row["run_status"] = "skipped_by_flag"
        row.update(switch.parse_time_log(time_path(args.run_root, sample, method)))
        return row

    out.parent.mkdir(parents=True, exist_ok=True)
    log_path(args.run_root, sample, method).parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    with log_path(args.run_root, sample, method).open("w") as handle:
        handle.write(f"# {shlex.join(cmd)}\n")
        handle.flush()
        try:
            subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=handle, env=env, check=True)
            row["run_status"] = "completed"
        except subprocess.CalledProcessError as exc:
            row["run_status"] = f"failed_exit_{exc.returncode}"
    row.update(switch.parse_time_log(time_path(args.run_root, sample, method)))
    return row


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--taxmap", type=Path, default=switch.TAXMAP)
    parser.add_argument("--model-cache", type=Path, default=switch.MODEL_CACHE)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.set_defaults(current_method="")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True)
    command_rows = []
    for sample in switch.SAMPLES:
        for method, method_args in METHODS:
            print(f"[default-validation] {sample.label} {method}", flush=True)
            command_rows.append(run_profile(args, sample, method, method_args))
            write_tsv(pd.DataFrame(command_rows), RESULTS / "default_profile_runtime.tsv")

    mapping_resources = switch.fresh.build_mapping_resources()
    score_rows = []
    addition_rows = []
    surface_rows = []
    for sample in switch.SAMPLES:
        for method, _method_args in METHODS:
            path = profile_path(args.run_root, sample, method)
            if not path.exists() or path.stat().st_size == 0:
                continue
            score, additions, surface = switch.score_profile(path, sample, method, mapping_resources)
            score_rows.append(score)
            addition_rows.extend(additions)
            surface_rows.append(surface)

    scores = pd.DataFrame(score_rows)
    additions = pd.DataFrame(addition_rows)
    write_tsv(scores, RESULTS / "default_profile_scores.tsv")
    write_tsv(additions, RESULTS / "default_profile_added.tsv")
    write_tsv(pd.DataFrame(surface_rows), RESULTS / "default_profile_surface.tsv")
    write_tsv(switch.summarize_additions(additions), RESULTS / "default_profile_added_summary.tsv")
    if not scores.empty:
        panel_summary, overall = switch.sweep.summarize(scores)
        write_tsv(panel_summary, RESULTS / "default_profile_panel_summary.tsv")
        write_tsv(overall, RESULTS / "default_profile_overall_summary.tsv")
        write_tsv(overall, EXP / "summary.tsv")
    print(f"[default-validation] wrote {RESULTS}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
