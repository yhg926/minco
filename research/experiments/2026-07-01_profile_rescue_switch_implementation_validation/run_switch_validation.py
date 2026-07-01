#!/usr/bin/env python3
"""Validate the implemented high-specificity profile-rescue switch.

This reruns the user-facing calibrated wrapper from saved unique/split tables
created by the fresh raw-run validation, comparing the current no-rescue output
against the implemented explicit rescue switch:

    split-p002-x300-ani095-af06-b025-d1-top1

Large profile tables are written to /tmp; this experiment stores only compact
commands and score summaries in the repository.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
BASE_EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
FRESH_EXP = ROOT / "research/experiments/2026-07-01_profile_rescue_fresh_raw_validation"
sys.path.insert(0, str(BASE_EXP))
sys.path.insert(0, str(FRESH_EXP))

import decompose_abundance_errors as decomp  # noqa: E402
import run_fresh_profile_rescue_validation as fresh  # noqa: E402
import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402
import sweep_cross_panel_abundance_variants as sweep  # noqa: E402


REF = Path(
    "/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/"
    "sketch_T_S1000_anno"
)
TAXMAP = REF / "species_taxmap.tsv"
MODEL_CACHE = REF / "minco_profile_rf_hgb.train12.unfiltered.joblib"
SOURCE_RUN = Path("/tmp/minco_profile_rescue_fresh_raw_validation_20260701")
DEFAULT_RUN_ROOT = Path("/tmp/minco_profile_rescue_switch_implementation_validation_20260701")
STRICT_SWITCH = "split-p002-x300-ani095-af06-b025-d1-top1"


@dataclass(frozen=True)
class Sample:
    panel: str
    sample: int
    label: str
    workdir: Path

    @property
    def unique_table(self) -> Path:
        return self.workdir / "minco.best_diff_unique.unfiltered.tsv"

    @property
    def split_table(self) -> Path:
        return self.workdir / "minco.best_diff_split.unfiltered.tsv"


SAMPLES = [
    Sample("cami3_toy_human_gut_gtdb_source_readmap", 0, "cami3_sample0", SOURCE_RUN / "work/cami3_sample0"),
    Sample("cami3_toy_human_gut_gtdb_source_readmap", 1, "cami3_sample1", SOURCE_RUN / "work/cami3_sample1"),
    Sample("cami3_toy_human_gut_gtdb_source_readmap", 2, "cami3_sample2", SOURCE_RUN / "work/cami3_sample2"),
    Sample("hmp_airskin_gtdb_source_abundance", 13, "hmp_airskin13", SOURCE_RUN / "work/hmp_airskin13"),
]

METHODS = [
    ("current_calibrated_abundance", "off"),
    ("profile_rescue_p002_x300_ani095_af06_b025_d1_top1", STRICT_SWITCH),
]


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_elapsed_seconds(line: str) -> float:
    text = line.split("):", 1)[1].strip() if "):" in line else line.rsplit(":", 1)[1].strip()
    parts = text.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(text)


def parse_time_log(path: Path) -> dict[str, object]:
    row: dict[str, object] = {
        "elapsed_seconds": np.nan,
        "max_rss_kb": np.nan,
        "user_seconds": np.nan,
        "system_seconds": np.nan,
        "exit_status": "",
    }
    if not path.exists():
        return row
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("\tElapsed (wall clock) time"):
            row["elapsed_seconds"] = parse_elapsed_seconds(line)
        elif line.startswith("\tMaximum resident set size"):
            row["max_rss_kb"] = float(line.split(":", 1)[1].strip())
        elif line.startswith("\tUser time"):
            row["user_seconds"] = float(line.split(":", 1)[1].strip())
        elif line.startswith("\tSystem time"):
            row["system_seconds"] = float(line.split(":", 1)[1].strip())
        elif line.startswith("\tExit status"):
            row["exit_status"] = line.split(":", 1)[1].strip()
    return row


def truth_df(panel: str, sample: int) -> pd.DataFrame:
    cfg = decomp.PANELS[panel]
    values = decomp.load_truth(Path(cfg["truth"](sample)), int(sample), str(cfg["truth_abundance_col"]))
    return pd.DataFrame({"gtdb_species": list(values), "truth_abundance": list(values.values())})


def collapse_rule(panel: str) -> str:
    return str(decomp.PANELS[panel]["minco_collapse"])


def profile_path(run_root: Path, sample: Sample, method: str) -> Path:
    return run_root / "profiles" / f"{sample.label}.{method}.tsv"


def time_path(run_root: Path, sample: Sample, method: str) -> Path:
    return run_root / "logs" / f"{sample.label}.{method}.time.log"


def log_path(run_root: Path, sample: Sample, method: str) -> Path:
    return run_root / "logs" / f"{sample.label}.{method}.run.log"


def command(args: argparse.Namespace, sample: Sample, method: str, switch: str) -> list[str]:
    return [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_path(args.run_root, sample, method)),
        "python3",
        "-B",
        str(ROOT / "scripts/minco_profile_calibrated.py"),
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
        "--strategy",
        "universal-auto-exact",
        "--abundance-ani-floor",
        "0.90",
        "--abundance-sparse-depth-cap",
        "poisson-breadth",
        "--abundance-sparse-breadth-max",
        "0.15",
        "--abundance-sparse-depth-ratio-min",
        "200",
        "--candidate-rescue-switch",
        switch,
        "--report-all",
        "-o",
        str(profile_path(args.run_root, sample, method)),
    ]


def run_profile(args: argparse.Namespace, sample: Sample, method: str, switch: str) -> dict[str, object]:
    out = profile_path(args.run_root, sample, method)
    cmd = command(args, sample, method, switch)
    row = {
        "panel": sample.panel,
        "sample": sample.sample,
        "label": sample.label,
        "method": method,
        "candidate_rescue_switch": switch,
        "unique_table": str(sample.unique_table),
        "split_table": str(sample.split_table),
        "profile": str(out),
        "command": shlex.join(cmd),
    }
    if out.exists() and out.stat().st_size > 0 and not args.force:
        row["run_status"] = "reused_existing_profile"
        row.update(parse_time_log(time_path(args.run_root, sample, method)))
        return row
    if args.skip_run:
        row["run_status"] = "skipped_by_flag"
        row.update(parse_time_log(time_path(args.run_root, sample, method)))
        return row

    for path in [out.parent, log_path(args.run_root, sample, method).parent]:
        path.mkdir(parents=True, exist_ok=True)
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
    row.update(parse_time_log(time_path(args.run_root, sample, method)))
    return row


def score_profile(profile: Path, sample: Sample, method: str, mapping_resources) -> tuple[dict[str, object], list[dict[str, object]], dict[str, object]]:
    by_accession, by_core, taxid_to_species, name_to_species = mapping_resources
    raw = pd.read_csv(profile, sep="\t", low_memory=False)
    mapped = fresh.add_gtdb_species(raw, by_accession, by_core, taxid_to_species, name_to_species)
    called = mapped["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = mapped.loc[called].copy()
    raw_abundance = pd.to_numeric(selected.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    pred = sweep.collapse_prediction(selected, raw_abundance.to_numpy(dtype=float), collapse_rule(sample.panel))
    row = taxid_score.score_prediction(sample.sample, method, pred, truth_df(sample.panel, sample.sample), {})
    row["panel"] = sample.panel
    row["label"] = sample.label
    row["collapse_rule"] = collapse_rule(sample.panel)
    row["profile_rows"] = int(len(mapped))
    row["called_rows"] = int(called.sum())
    row["called_species"] = int(selected["gtdb_species"].nunique())
    row["candidate_rescue_added_rows"] = int(
        selected.get("candidate_rescue_added", pd.Series(False, index=selected.index))
        .astype(str)
        .str.lower()
        .isin({"true", "1", "yes"})
        .sum()
    )
    truth = truth_df(sample.panel, sample.sample)
    truth_map = dict(zip(truth["gtdb_species"].astype(str), truth["truth_abundance"].astype(float)))
    additions: list[dict[str, object]] = []
    if "candidate_rescue_added" in selected.columns:
        added_mask = selected["candidate_rescue_added"].astype(str).str.lower().isin({"true", "1", "yes"})
        for rec in selected.loc[added_mask].itertuples(index=False):
            species = str(getattr(rec, "gtdb_species", ""))
            additions.append(
                {
                    "panel": sample.panel,
                    "sample": sample.sample,
                    "label": sample.label,
                    "method": method,
                    "gtdb_species": species,
                    "truth_abundance": truth_map.get(species, 0.0),
                    "is_truth_species": species in truth_map,
                    "calibrated_probability": fresh.finite_float(getattr(rec, "calibrated_probability", 0.0)),
                    "calibrated_abundance_raw": fresh.finite_float(getattr(rec, "calibrated_abundance_raw", 0.0)),
                    "s_XnY_ctx_max": fresh.finite_float(getattr(rec, "s_XnY_ctx_max", 0.0)),
                    "s_ANI_max": fresh.finite_float(getattr(rec, "s_ANI_max", 0.0)),
                    "s_Real_min_align_fraction_max": fresh.finite_float(getattr(rec, "s_Real_min_align_fraction_max", 0.0)),
                    "s_Ref_breadth_max": fresh.finite_float(getattr(rec, "s_Ref_breadth_max", 0.0)),
                    "s_Ref_mean_depth_max": fresh.finite_float(getattr(rec, "s_Ref_mean_depth_max", 0.0)),
                }
            )
    surface = {
        "panel": sample.panel,
        "sample": sample.sample,
        "label": sample.label,
        "method": method,
        "profile_rows": int(len(mapped)),
        "called_rows": int(called.sum()),
        "called_species": int(selected["gtdb_species"].nunique()),
        "rescue_added_rows": int(len(additions)),
        "rescue_added_truth_rows": int(sum(1 for item in additions if item["is_truth_species"])),
        "rescue_added_fp_rows": int(sum(1 for item in additions if not item["is_truth_species"])),
    }
    return row, additions, surface


def summarize_additions(additions: pd.DataFrame) -> pd.DataFrame:
    if additions.empty:
        return pd.DataFrame(columns=["panel", "method", "added_rows", "added_truth_rows", "added_fp_rows", "added_truth_abundance_mass_pct"])
    work = additions.copy()
    work["truth_abundance"] = pd.to_numeric(work["truth_abundance"], errors="coerce").fillna(0.0)
    out = work.groupby(["panel", "method"], as_index=False).agg(
        added_rows=("gtdb_species", "size"),
        added_truth_rows=("is_truth_species", "sum"),
        added_truth_abundance_mass_pct=("truth_abundance", lambda x: float(x.sum() * 100.0)),
    )
    out["added_fp_rows"] = out["added_rows"] - out["added_truth_rows"]
    return out[["panel", "method", "added_rows", "added_truth_rows", "added_fp_rows", "added_truth_abundance_mass_pct"]]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--taxmap", type=Path, default=TAXMAP)
    parser.add_argument("--model-cache", type=Path, default=MODEL_CACHE)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True)
    command_rows = []
    for sample in SAMPLES:
        for method, switch in METHODS:
            if not sample.unique_table.exists() or not sample.split_table.exists():
                raise SystemExit(f"missing saved tables for {sample.label}: {sample.workdir}")
            print(f"[switch-validation] {sample.label} {method}", flush=True)
            command_rows.append(run_profile(args, sample, method, switch))
            write_tsv(pd.DataFrame(command_rows), RESULTS / "switch_profile_runtime.tsv")

    mapping_resources = fresh.build_mapping_resources()
    score_rows: list[dict[str, object]] = []
    addition_rows: list[dict[str, object]] = []
    surface_rows: list[dict[str, object]] = []
    for sample in SAMPLES:
        for method, _switch in METHODS:
            out = profile_path(args.run_root, sample, method)
            if not out.exists() or out.stat().st_size == 0:
                continue
            score, additions, surface = score_profile(out, sample, method, mapping_resources)
            score_rows.append(score)
            addition_rows.extend(additions)
            surface_rows.append(surface)

    scores = pd.DataFrame(score_rows)
    additions = pd.DataFrame(addition_rows)
    write_tsv(scores, RESULTS / "switch_profile_scores.tsv")
    write_tsv(additions, RESULTS / "switch_profile_added.tsv")
    write_tsv(pd.DataFrame(surface_rows), RESULTS / "switch_profile_surface.tsv")
    write_tsv(summarize_additions(additions), RESULTS / "switch_profile_added_summary.tsv")
    if not scores.empty:
        panel_summary, overall = sweep.summarize(scores)
        write_tsv(panel_summary, RESULTS / "switch_profile_panel_summary.tsv")
        write_tsv(overall, RESULTS / "switch_profile_overall_summary.tsv")
        write_tsv(overall, EXP / "summary.tsv")
    print(f"[switch-validation] wrote {RESULTS}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
