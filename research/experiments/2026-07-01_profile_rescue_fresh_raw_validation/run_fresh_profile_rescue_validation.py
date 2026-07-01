#!/usr/bin/env python3
"""Fresh raw-read validation for the high-specificity profile rescue rule.

The cached cross-panel replay showed that a strict uncalled-candidate rescue can
recover high-abundance CAMI3 false negatives without cached false additions. This
script reruns MinCO from raw reads with ``--report-all`` for local samples, maps
profile rows to GTDB species by best reference accession, then scores:

* current calibrated abundance
* current abundance guard
* current abundance guard plus high-specificity profile rescue

Large profile tables are written under /tmp. Only compact summaries should be
kept in the repository.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
BASE_EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
CROSS_EXP = ROOT / "research/experiments/2026-07-01_abundance_guard_cross_panel_validation"

sys.path.insert(0, str(BASE_EXP))
sys.path.insert(0, str(CROSS_EXP))

import decompose_abundance_errors as decomp  # noqa: E402
import diagnose_cami3_fn_candidate_recovery as rescue  # noqa: E402
import score_abundance_guards_cached as guard  # noqa: E402,F401
import score_cami3_gtdb_source_readmap as cami3  # noqa: E402
import sweep_cross_panel_abundance_variants as sweep  # noqa: E402
import validate_profile_rescue_cross_panel as cached_rescue  # noqa: E402


REF = Path(
    "/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/"
    "sketch_T_S1000_anno"
)
TAXMAP = REF / "species_taxmap.tsv"
MODEL_CACHE = REF / "minco_profile_rf_hgb.train12.unfiltered.joblib"
DEFAULT_RUN_ROOT = Path("/tmp/minco_profile_rescue_fresh_raw_validation_20260701")

METHODS = [
    "current_calibrated_abundance",
    "abundance_ani_floor090_sparse_cap",
    "profile_rescue_p002_x300_ani095_af06_b025_d1_top1",
]

ACC_RE = re.compile(r"(GC[AF]_[0-9]+(?:\.[0-9]+)?)")


@dataclass(frozen=True)
class Sample:
    panel: str
    sample: int
    label: str
    reads: Path

    def out_path(self, run_root: Path) -> Path:
        return run_root / "profiles" / f"{self.label}.profile.tsv"

    def workdir(self, run_root: Path) -> Path:
        return run_root / "work" / self.label

    def log_path(self, run_root: Path) -> Path:
        return run_root / "logs" / f"{self.label}.run.log"

    def time_path(self, run_root: Path) -> Path:
        return run_root / "logs" / f"{self.label}.time.log"


SAMPLES = [
    Sample(
        "cami3_toy_human_gut_gtdb_source_readmap",
        0,
        "cami3_sample0",
        Path("/tmp/minco_exact_trigger_validation_20260630/reads/cami3_sample0/sample_0_reads/anonymous_reads.fq.gz"),
    ),
    Sample(
        "cami3_toy_human_gut_gtdb_source_readmap",
        1,
        "cami3_sample1",
        Path("/tmp/minco_exact_trigger_validation_20260630/reads/cami3_sample1/sample_1_reads/anonymous_reads.fq.gz"),
    ),
    Sample(
        "cami3_toy_human_gut_gtdb_source_readmap",
        2,
        "cami3_sample2",
        Path("/tmp/minco_exact_trigger_validation_20260630/reads/cami3_sample2/sample_2_reads/anonymous_reads.fq.gz"),
    ),
    Sample(
        "hmp_airskin_gtdb_source_abundance",
        13,
        "hmp_airskin13",
        Path("/tmp/minco_exact_trigger_validation_20260630/reads/airskinurogenital_sample13.nonzero.fastq.gz"),
    ),
]

MISSING_AUDIT = [
    Sample("cami2_toy_mouse_gut", 5, "toy_mouse5", Path("/tmp/cami2_toymouse_samples5_7_20260625/reads/sample5.nonzero.fastq.gz")),
    Sample("cami2_toy_mouse_gut", 6, "toy_mouse6", Path("/tmp/cami2_toymouse_samples5_7_20260625/reads/sample6.nonzero.fastq.gz")),
    Sample("cami2_toy_mouse_gut", 7, "toy_mouse7", Path("/tmp/cami2_toymouse_samples5_7_20260625/reads/sample7.nonzero.fastq.gz")),
    Sample("hmp_airskin_gtdb_source_abundance", 6, "hmp_airskin6", Path("/tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample6.nonzero.fastq.gz")),
    Sample("hmp_airskin_gtdb_source_abundance", 11, "hmp_airskin11", Path("/tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample11.nonzero.fastq.gz")),
    Sample("hmp_airskin_gtdb_source_abundance", 28, "hmp_airskin28", Path("/tmp/cami2_hmp_airskin_20260625/reads/airskinurogenital_sample28.nonzero.fastq.gz")),
]


def finite_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def extract_accession(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def clean_taxid(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    if re.fullmatch(r"[0-9]+(?:\\.0+)?", text):
        return text.split(".", 1)[0]
    return text


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def parse_elapsed_seconds(value: str) -> float:
    value = value.strip()
    parts = value.split(":")
    try:
        if len(parts) == 3:
            hours, minutes, seconds = parts
            return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
        if len(parts) == 2:
            minutes, seconds = parts
            return int(minutes) * 60 + float(seconds)
        return float(value)
    except ValueError:
        return float("nan")


def parse_time_log(path: Path) -> dict[str, object]:
    out: dict[str, object] = {
        "time_log": str(path),
        "elapsed_seconds": float("nan"),
        "max_rss_kb": float("nan"),
        "user_seconds": float("nan"),
        "system_seconds": float("nan"),
        "exit_status": "",
    }
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("\tElapsed (wall clock) time"):
            elapsed = line.split("):", 1)[1].strip() if "):" in line else line.rsplit(":", 1)[1].strip()
            out["elapsed_seconds"] = parse_elapsed_seconds(elapsed)
        elif line.startswith("\tMaximum resident set size"):
            out["max_rss_kb"] = finite_float(line.split(":", 1)[1].strip(), float("nan"))
        elif line.startswith("\tUser time"):
            out["user_seconds"] = finite_float(line.split(":", 1)[1].strip(), float("nan"))
        elif line.startswith("\tSystem time"):
            out["system_seconds"] = finite_float(line.split(":", 1)[1].strip(), float("nan"))
        elif line.startswith("\tExit status"):
            out["exit_status"] = line.split(":", 1)[1].strip()
    return out


def profile_command(sample: Sample, args: argparse.Namespace) -> list[str]:
    return [
        "/usr/bin/time",
        "-v",
        "-o",
        str(sample.time_path(args.run_root)),
        "python3",
        "-B",
        str(ROOT / "scripts/minco_profile_calibrated.py"),
        "-r",
        str(args.ref),
        "--reads",
        str(sample.reads),
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
        "--workdir",
        str(sample.workdir(args.run_root)),
        "--minco",
        str(args.minco),
        "-p",
        str(args.threads),
        "--report-all",
        "-o",
        str(sample.out_path(args.run_root)),
    ]


def run_profile(sample: Sample, args: argparse.Namespace) -> dict[str, object]:
    out = sample.out_path(args.run_root)
    log = sample.log_path(args.run_root)
    cmd = profile_command(sample, args)
    row = {
        "panel": sample.panel,
        "sample": sample.sample,
        "label": sample.label,
        "reads": str(sample.reads),
        "profile": str(out),
        "log": str(log),
        "command": shlex.join(cmd),
    }
    if out.exists() and out.stat().st_size > 0 and not args.force:
        row["run_status"] = "reused_existing_profile"
        row.update(parse_time_log(sample.time_path(args.run_root)))
        return row
    if args.skip_run:
        row["run_status"] = "skipped_by_flag"
        row.update(parse_time_log(sample.time_path(args.run_root)))
        return row

    out.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    sample.workdir(args.run_root).mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    with log.open("w") as handle:
        handle.write(f"# {shlex.join(cmd)}\n")
        handle.flush()
        try:
            subprocess.run(cmd, cwd=ROOT, stdout=handle, stderr=handle, env=env, check=True)
            row["run_status"] = "completed"
        except subprocess.CalledProcessError as exc:
            row["run_status"] = f"failed_exit_{exc.returncode}"
            row.update(parse_time_log(sample.time_path(args.run_root)))
            return row
    row.update(parse_time_log(sample.time_path(args.run_root)))
    return row


def truth_df(panel: str, sample: int) -> pd.DataFrame:
    cfg = decomp.PANELS[panel]
    values = decomp.load_truth(Path(cfg["truth"](sample)), int(sample), str(cfg["truth_abundance_col"]))
    return pd.DataFrame({"gtdb_species": list(values), "truth_abundance": list(values.values())})


def build_mapping_resources() -> tuple[
    Mapping[str, Mapping[str, str]],
    Mapping[str, list[Mapping[str, str]]],
    Mapping[str, str],
    Mapping[str, str],
]:
    toy_mod = decomp.load_toy_scorer()
    by_accession, by_core = toy_mod.score.truth.load_gtdb_metadata(toy_mod.score.truth.GTDB_METADATA)
    _wgs_to_species, taxid_to_species, name_to_species, _diag = cami3.build_transfer_maps()
    return by_accession, by_core, taxid_to_species, name_to_species


def species_from_accession(
    accession: str,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    cache: dict[str, str],
) -> str:
    if not accession:
        return ""
    if accession not in cache:
        cache[accession] = cami3.gtdb_from_accession(accession, by_accession, by_core)
    return cache[accession]


def add_gtdb_species(
    profile: pd.DataFrame,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
) -> pd.DataFrame:
    acc_cache: dict[str, str] = {}
    species_out: list[str] = []
    method_out: list[str] = []

    for row in profile.itertuples(index=False):
        species = ""
        method = ""
        for col in ("s_best_accession", "u_best_accession"):
            acc = extract_accession(getattr(row, col, ""))
            species = species_from_accession(acc, by_accession, by_core, acc_cache)
            if species:
                method = col
                break
        if not species:
            name = str(getattr(row, "species_name", "") or "").strip()
            if name.startswith("s__"):
                species = name
                method = "species_name_gtdb"
        if not species:
            taxid = clean_taxid(getattr(row, "taxid", ""))
            species = taxid_to_species.get(taxid, "")
            if species:
                method = "taxid_unique_gtdb"
        if not species:
            name = str(getattr(row, "species_name", "") or "").strip()
            normalized = cami3.normalize_name(name)
            species = name_to_species.get(normalized, "")
            if species:
                method = "name_unique_gtdb"
        species_out.append(species)
        method_out.append(method)

    out = profile.copy()
    out["gtdb_species"] = species_out
    out["gtdb_mapping_method"] = method_out
    out["called"] = out["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    return out.loc[out["gtdb_species"].astype(str).astype(bool)].copy()


def load_mapped_profile(
    sample: Sample,
    args: argparse.Namespace,
    by_accession,
    by_core,
    taxid_to_species,
    name_to_species,
) -> tuple[pd.DataFrame, dict[str, object]]:
    path = sample.out_path(args.run_root)
    profile = pd.read_csv(path, sep="\t", low_memory=False)
    mapped = add_gtdb_species(profile, by_accession, by_core, taxid_to_species, name_to_species)
    surface = {
        "panel": sample.panel,
        "sample": sample.sample,
        "label": sample.label,
        "profile": str(path),
        "profile_rows": int(len(profile)),
        "mapped_rows": int(len(mapped)),
        "mapped_species": int(mapped["gtdb_species"].nunique()),
        "called_rows": int(mapped["called"].sum()),
        "called_species": int(mapped.loc[mapped["called"], "gtdb_species"].nunique()),
        "uncalled_rows": int((~mapped["called"]).sum()),
        "uncalled_species": int(mapped.loc[~mapped["called"], "gtdb_species"].nunique()),
    }
    return mapped, surface


def score_profile(sample: Sample, candidates: pd.DataFrame) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    truth = truth_df(sample.panel, sample.sample)
    score_rows: list[dict[str, object]] = []
    addition_rows: list[dict[str, object]] = []
    for method in METHODS:
        selected, added = cached_rescue.select_rows(candidates, method)
        row = cached_rescue.score_selected(sample.panel, sample.sample, method, selected, added, truth)
        row["label"] = sample.label
        if method.startswith("profile_rescue"):
            added_records = cached_rescue.summarize_additions(sample.panel, sample.sample, candidates, method, added)
            addition_rows.extend(added_records)
            row["added_truth_rows"] = int(sum(1 for item in added_records if item["is_truth_species"]))
            row["added_fp_rows"] = int(sum(1 for item in added_records if not item["is_truth_species"]))
        score_rows.append(row)
    return score_rows, addition_rows


def score_completed_profiles(samples: list[Sample], args: argparse.Namespace) -> None:
    completed = [
        sample
        for sample in samples
        if sample.out_path(args.run_root).exists() and sample.out_path(args.run_root).stat().st_size > 0
    ]
    if not completed:
        for path in [
            RESULTS / "fresh_profile_scores.tsv",
            RESULTS / "fresh_profile_added.tsv",
            RESULTS / "fresh_profile_surface.tsv",
            RESULTS / "fresh_profile_added_summary.tsv",
            RESULTS / "fresh_profile_panel_summary.tsv",
            RESULTS / "fresh_profile_overall_summary.tsv",
            EXP / "summary.tsv",
        ]:
            write_tsv(pd.DataFrame(), path)
        return

    by_accession, by_core, taxid_to_species, name_to_species = build_mapping_resources()
    score_rows: list[dict[str, object]] = []
    addition_rows: list[dict[str, object]] = []
    surface_rows: list[dict[str, object]] = []

    for sample in completed:
        candidates, surface = load_mapped_profile(
            sample,
            args,
            by_accession,
            by_core,
            taxid_to_species,
            name_to_species,
        )
        surface_rows.append(surface)
        rows, additions = score_profile(sample, candidates)
        score_rows.extend(rows)
        addition_rows.extend(additions)

    scores = pd.DataFrame(score_rows)
    additions = pd.DataFrame(addition_rows)
    surfaces = pd.DataFrame(surface_rows)
    write_tsv(scores, RESULTS / "fresh_profile_scores.tsv")
    write_tsv(additions, RESULTS / "fresh_profile_added.tsv")
    write_tsv(surfaces, RESULTS / "fresh_profile_surface.tsv")
    write_tsv(cached_rescue.addition_summary(additions), RESULTS / "fresh_profile_added_summary.tsv")
    if not scores.empty:
        panel_summary, overall = sweep.summarize(scores)
        write_tsv(panel_summary, RESULTS / "fresh_profile_panel_summary.tsv")
        write_tsv(overall, RESULTS / "fresh_profile_overall_summary.tsv")
        write_tsv(overall, EXP / "summary.tsv")


def input_audit(args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    for sample in [*SAMPLES, *MISSING_AUDIT]:
        rows.append(
            {
                "panel": sample.panel,
                "sample": sample.sample,
                "label": sample.label,
                "reads": str(sample.reads),
                "reads_exists": sample.reads.exists(),
                "profile": str(sample.out_path(args.run_root)),
                "profile_exists": sample.out_path(args.run_root).exists(),
                "planned_for_fresh_run": sample in SAMPLES,
            }
        )
    return pd.DataFrame(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--ref", type=Path, default=REF)
    parser.add_argument("--taxmap", type=Path, default=TAXMAP)
    parser.add_argument("--model-cache", type=Path, default=MODEL_CACHE)
    parser.add_argument("--minco", type=Path, default=ROOT / "bin/minco")
    parser.add_argument("-p", "--threads", type=int, default=8)
    parser.add_argument("--sample", action="append", default=[], help="Run only these labels; may be repeated.")
    parser.add_argument("--force", action="store_true", help="Rerun profiles even if output TSV already exists.")
    parser.add_argument("--skip-run", action="store_true", help="Score existing profiles without launching MinCO.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    RESULTS.mkdir(parents=True, exist_ok=True)
    for required in [args.ref, args.taxmap, args.model_cache, args.minco]:
        if not required.exists():
            raise SystemExit(f"missing required input: {required}")

    requested = set(args.sample)
    samples = [sample for sample in SAMPLES if not requested or sample.label in requested]
    audit = input_audit(args)
    write_tsv(audit, RESULTS / "fresh_profile_input_audit.tsv")
    available = [sample for sample in samples if sample.reads.exists()]
    if not available:
        raise SystemExit("no requested fresh-run samples have local reads")

    command_rows = []
    for sample in available:
        print(f"[fresh-validation] profiling {sample.label}", flush=True)
        command_rows.append(run_profile(sample, args))
        write_tsv(pd.DataFrame(command_rows), RESULTS / "fresh_profile_runtime.tsv")

    score_completed_profiles(available, args)
    print(f"[fresh-validation] wrote {RESULTS}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
