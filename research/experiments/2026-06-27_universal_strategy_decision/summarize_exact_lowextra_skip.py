#!/usr/bin/env python3
"""Summarize the low-extra exact-skip speed experiment."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pandas as pd


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-27_universal_strategy_decision"
RESULTS = EXP / "results"
SCORER = (
    ROOT
    / "research/experiments/2026-06-26_cami2_toymouse_current_default/"
    "score_sample7_current_default.py"
)
RUN = Path("/tmp/minco_exact_split_sidecar_wrapper_20260627")
READS = Path(
    "/tmp/cami2_toymouse_samples5_7_20260625/data/sample_6/"
    "2017.12.29_11.37.26_sample_6/reads/anonymous_reads.fq.gz"
)
REF = Path(
    "/mnt/new3T/gtdbr220/gtdb232/"
    "GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno"
)
SYLPH_PROFILE = Path("/tmp/cami2_toymouse_samples5_7_20260625/run/sylph_sample6/profile.tsv")
SYLPH_RUNTIME_TABLE = RESULTS / "runtime_memory_minco_vs_sylph.tsv"


MINCO_PROFILES = {
    "minco_lowextra_skip_autoexact_sample6": RUN / "sample6_lowextra_skip_p16.tsv",
    "minco_sidecar_autoexact_sample6": RUN / "sample6_sidecar_p16.tsv",
}
MINCO_RUNTIME = {
    "wrapper_lowextra_skip_default": (
        RUN / "sample6_lowextra_skip_p16.time.log",
        MINCO_PROFILES["minco_lowextra_skip_autoexact_sample6"],
        "skip",
    ),
    "wrapper_same_stream_exact_sidecar": (
        RUN / "sample6_sidecar_p16.time.log",
        MINCO_PROFILES["minco_sidecar_autoexact_sample6"],
        "allow",
    ),
    "wrapper_one_stream_unique_block_exact": (
        RUN / "sample6_onestream_p16.time.log",
        RUN / "sample6_onestream_p16.tsv",
        "allow",
    ),
}


def parse_elapsed(text: str) -> tuple[str, float]:
    match = re.search(r"Elapsed .*: ([0-9:.]+)", text)
    if not match:
        return "", 0.0
    elapsed = match.group(1)
    parts = elapsed.split(":")
    if len(parts) == 3:
        seconds = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        seconds = int(parts[0]) * 60 + float(parts[1])
    else:
        seconds = float(elapsed)
    return elapsed, seconds


def parse_time_log(path: Path) -> tuple[str, float, int, float]:
    text = path.read_text()
    elapsed, seconds = parse_elapsed(text)
    rss_match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    rss = int(rss_match.group(1)) if rss_match else 0
    return elapsed, seconds, rss, rss / 1048576.0 if rss else 0.0


def load_metadata(path: Path) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t", nrows=1)
    if df.empty:
        return {
            "exact_split_requested": "",
            "exact_split_used": "",
            "check": "empty profile",
        }
    row = df.iloc[0]
    requested = row.get("auto_exact_split_requested", "")
    used = row.get("auto_exact_split_used", "")
    reason = str(row.get("auto_exact_split_unavailable_reason", "") or "")
    source = str(row.get("auto_exact_split_source", "") or "")
    if str(used).lower() == "true":
        check = f"auto_exact_split_used=true; source={source}"
    elif reason:
        check = f"auto_exact_split_unavailable_reason={reason}"
    else:
        check = "auto_exact_split_used=false"
    return {
        "exact_split_requested": str(requested).lower(),
        "exact_split_used": str(used).lower(),
        "check": check,
    }


def load_sylph_runtime() -> dict[str, object]:
    if not SYLPH_RUNTIME_TABLE.exists():
        return {"wall_time": "", "seconds": "", "peak_rss_gib": ""}
    table = pd.read_csv(SYLPH_RUNTIME_TABLE, sep="\t")
    row = table.loc[
        table["dataset"].eq("CAMI2 toy mouse")
        & table["sample"].eq("sample6")
        & table["method"].eq("Sylph sketch+profile")
    ]
    if row.empty:
        return {"wall_time": "", "seconds": "", "peak_rss_gib": ""}
    rec = row.iloc[0]
    return {
        "wall_time": rec["wall_time"],
        "seconds": rec["seconds"],
        "peak_rss_gib": rec["peak_rss_gib"],
    }


def score_profiles() -> None:
    spec = importlib.util.spec_from_file_location("toy_scorer", SCORER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import scorer {SCORER}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    score = mod.score
    by_accession, by_core = score.truth.load_gtdb_metadata(score.truth.GTDB_METADATA)
    genome_taxid, genome_location = mod.load_meta()
    score.EXP_DIR = RESULTS
    score.BASE = mod.DATA_BASE
    score.distribution_path = mod.distribution_path
    sample = 6
    source_rows = score.build_gtdb_truth(
        sample, by_accession, by_core, genome_taxid, genome_location
    )
    profile = score.write_truth_files(sample, source_rows)
    gold = {str(row["gtdb_species"]) for row in profile}

    rows: list[dict[str, object]] = []
    methods = [(name, path, "minco") for name, path in MINCO_PROFILES.items()]
    methods.append(("sylph_gtdb_profile_sample6", SYLPH_PROFILE, "sylph"))
    for method, path, kind in methods:
        calls = (
            score.load_sylph_rows(path, by_accession, by_core)
            if kind == "sylph"
            else mod.score_calibrated_minco(path)
        )
        selected = calls.loc[calls["active_gate_pass"]].copy()
        tp, fp, fn, precision, recall, f1 = score.score_sets(selected["gtdb_species"], gold)
        renorm = [
            rec
            for rec in score.abundance_metrics(selected, profile)
            if str(rec.get("renorm_pred")).lower() == "true"
        ][0]
        rows.append(
            {
                "sample": sample,
                "method": method,
                "path": str(path),
                "gold_taxa": len(gold),
                "pred_taxa": len(tp | fp),
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "precision": precision,
                "recall": recall,
                "F1": f1,
                "l1_pct_points": renorm["l1_pct_points"],
                "pearson": renorm["pearson"],
            }
        )
    pd.DataFrame(rows).to_csv(
        RESULTS / "exact_split_lowextra_skip_sample6_score.tsv",
        sep="\t",
        index=False,
    )


def write_runtime() -> None:
    rows: list[dict[str, object]] = []
    for mode, (log_path, output, low_extra_mode) in MINCO_RUNTIME.items():
        wall, seconds, rss_kb, rss_gib = parse_time_log(log_path)
        meta = load_metadata(output)
        if mode == "wrapper_lowextra_skip_default":
            interpretation = (
                "skips exact after block low-extra split rescue; faster than "
                "Sylph sample6 110.98s with much lower RSS, but Sylph remains "
                "more accurate on this sample"
            )
        elif mode == "wrapper_same_stream_exact_sidecar":
            interpretation = "previous exact-sidecar path; same calls as lowextra skip on sample6"
        else:
            interpretation = (
                "forcing one process to write unique, block split, and exact split "
                "sidecars loses p16 pass-level concurrency and is slower"
            )
        rows.append(
            {
                "dataset": "CAMI2 toy mouse",
                "sample": "sample6 wrapper full",
                "mode": mode,
                "threads": 16,
                "reads": str(READS),
                "ref": str(REF),
                "strategy": "universal-auto-exact",
                "exact_low_extra_mode": low_extra_mode,
                "exact_split_requested": meta["exact_split_requested"],
                "exact_split_used": meta["exact_split_used"],
                "wall_time": wall,
                "seconds": f"{seconds:.2f}",
                "peak_rss_kb": rss_kb,
                "peak_rss_gib": f"{rss_gib:.4f}",
                "output": str(output),
                "check": meta["check"],
                "interpretation": interpretation,
            }
        )
    sylph = load_sylph_runtime()
    rows.append(
        {
            "dataset": "CAMI2 toy mouse",
            "sample": "sample6 external baseline",
            "mode": "Sylph sketch+profile",
            "threads": "NA",
            "reads": str(READS),
            "ref": "Sylph GTDB r226 c200 DB",
            "strategy": "Sylph",
            "exact_low_extra_mode": "NA",
            "exact_split_requested": "NA",
            "exact_split_used": "NA",
            "wall_time": sylph.get("wall_time", ""),
            "seconds": sylph.get("seconds", ""),
            "peak_rss_kb": "",
            "peak_rss_gib": sylph.get("peak_rss_gib", ""),
            "output": str(SYLPH_PROFILE),
            "check": "from runtime_memory_minco_vs_sylph.tsv",
            "interpretation": "external speed and accuracy baseline",
        }
    )
    pd.DataFrame(rows).to_csv(
        RESULTS / "exact_split_lowextra_skip_runtime.tsv",
        sep="\t",
        index=False,
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    score_profiles()
    write_runtime()
    print(RESULTS / "exact_split_lowextra_skip_runtime.tsv")
    print(RESULTS / "exact_split_lowextra_skip_sample6_score.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
