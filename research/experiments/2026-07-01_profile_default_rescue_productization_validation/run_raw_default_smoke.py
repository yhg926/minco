#!/usr/bin/env python3
"""Run a tiny raw-input smoke test through scripts/minco_profile.

The fixture is generated under /tmp and is intentionally small. It verifies
that the user-facing default launcher can build readwise tables from raw input,
load a fitted model cache, and emit a profile using the productized defaults.
"""

from __future__ import annotations

import os
import random
import shutil
import subprocess
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
ROOT = EXP.parents[2]
RESULTS = EXP / "results"
RUN_ROOT = Path(os.environ.get("MINCO_RAW_SMOKE_ROOT", "/tmp/minco_default_raw_smoke_20260701"))
MODEL_CACHE = Path(
    os.environ.get(
        "MINCO_RAW_SMOKE_MODEL_CACHE",
        "/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/"
        "sketch_T_S1000_anno/minco_profile_rf_hgb.train12.unfiltered.joblib",
    )
)


def write_wrapped_fasta(path: Path, header: str, seq: str) -> None:
    path.write_text(
        f">{header}\n" + "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n"
    )


def parse_time_log(path: Path) -> dict[str, object]:
    out: dict[str, object] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("User time"):
            out["user_seconds"] = float(line.split(":", 1)[1])
        elif line.startswith("System time"):
            out["system_seconds"] = float(line.split(":", 1)[1])
        elif line.startswith("Elapsed"):
            out["elapsed"] = line.split("):", 1)[1].strip()
        elif line.startswith("Maximum resident set size"):
            out["max_rss_kb"] = int(line.split(":", 1)[1])
        elif line.startswith("Exit status"):
            out["exit_status"] = int(line.split(":", 1)[1])
    return out


def run(label: str, cmd: list[str], cwd: Path, log_dir: Path) -> None:
    print(f"[raw-smoke] {label}: {' '.join(cmd)}", flush=True)
    with (log_dir / f"{label}.log").open("w") as handle:
        subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=True)


def main() -> int:
    if not MODEL_CACHE.is_file():
        raise SystemExit(f"model cache missing: {MODEL_CACHE}")
    if RUN_ROOT.exists():
        shutil.rmtree(RUN_ROOT)
    refs = RUN_ROOT / "refs"
    logs = RUN_ROOT / "logs"
    work = RUN_ROOT / "profile_work"
    refs.mkdir(parents=True)
    logs.mkdir()
    work.mkdir()

    random.seed(1701)
    alphabet = "ACGT"
    seq1 = "".join(random.choice(alphabet) for _ in range(6000))
    seq2 = "".join(random.choice(alphabet) for _ in range(6000))
    ref1 = refs / "GCF_000000001.1.fa"
    ref2 = refs / "GCF_000000002.1.fa"
    write_wrapped_fasta(ref1, "GCF_000000001.1 fixture_1", seq1)
    write_wrapped_fasta(ref2, "GCF_000000002.1 fixture_2", seq2)
    (RUN_ROOT / "ref.list").write_text(f"{ref1}\n{ref2}\n")

    read_rows = []
    for i, start in enumerate(range(0, 2500, 125), 1):
        read = seq1[start : start + 150]
        read_rows.append(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n")
    (RUN_ROOT / "reads.fq").write_text("".join(read_rows))
    (RUN_ROOT / "taxmap.tsv").write_text(
        "ref_key\tTAXID\tRANK\tTAXPATH\tTAXPATHSN\n"
        "GCF_000000001.1\t1001\tspecies\t"
        "d__Bacteria|p__Smoke|c__Smoke|o__Smoke|f__Smoke|g__Smoke|s__fixture_1\t"
        "d__Bacteria|p__Smoke|c__Smoke|o__Smoke|f__Smoke|g__Smoke|s__fixture_1\n"
        "GCF_000000002.1\t1002\tspecies\t"
        "d__Bacteria|p__Smoke|c__Smoke|o__Smoke|f__Smoke|g__Smoke|s__fixture_2\t"
        "d__Bacteria|p__Smoke|c__Smoke|o__Smoke|f__Smoke|g__Smoke|s__fixture_2\n"
    )

    run(
        "sketch",
        [
            str(ROOT / "bin/minco"),
            "sketch",
            "-S",
            "100",
            "-l",
            str(RUN_ROOT / "ref.list"),
            "-o",
            str(RUN_ROOT / "ref.minco"),
            "--anno",
        ],
        ROOT,
        logs,
    )
    run("index", [str(ROOT / "bin/minco"), "sketch", "-i", str(RUN_ROOT / "ref.minco")], ROOT, logs)

    profile = RUN_ROOT / "profile.tsv"
    profile_cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(logs / "default_profile.time.log"),
        "python3",
        "-B",
        str(ROOT / "scripts/minco_profile"),
        "-r",
        str(RUN_ROOT / "ref.minco"),
        "--reads",
        str(RUN_ROOT / "reads.fq"),
        "--taxmap",
        str(RUN_ROOT / "taxmap.tsv"),
        "--model-cache",
        str(MODEL_CACHE),
        "--scope",
        "bacteria",
        "--report-all",
        "--workdir",
        str(work),
        "-p",
        "1",
        "-o",
        str(profile),
    ]
    run("default_profile", profile_cmd, ROOT, logs)

    table = pd.read_csv(profile, sep="\t")
    first = table.iloc[0].to_dict() if not table.empty else {}
    summary = {
        "run_root": str(RUN_ROOT),
        "model_cache": str(MODEL_CACHE),
        "profile_rows": len(table),
        "called_rows": int(table["calibrated_call"].astype(bool).sum()) if "calibrated_call" in table else 0,
        "candidate_rescue_switch": first.get("candidate_rescue_switch", ""),
        "candidate_surface_switch": first.get("candidate_surface_switch", ""),
        "candidate_abundance_policy": first.get("candidate_abundance_policy", ""),
        "abundance_sparse_depth_cap_switch": first.get("abundance_sparse_depth_cap_switch", ""),
        "auto_exact_split_used": first.get("auto_exact_split_used", ""),
        "auto_exact_split_source": first.get("auto_exact_split_source", ""),
        "abundance_exact_hit_applied": first.get("abundance_exact_hit_applied", ""),
    }
    summary.update(parse_time_log(logs / "default_profile.time.log"))
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(RESULTS / "raw_default_smoke_summary.tsv", sep="\t", index=False)
    print(f"[raw-smoke] wrote {RESULTS / 'raw_default_smoke_summary.tsv'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
