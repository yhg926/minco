#!/usr/bin/env python3
"""Whole-genome check for crispatus trace contexts contributed by helveticus reads."""

from __future__ import annotations

import gzip
import math
import shutil
import subprocess
from pathlib import Path

import pandas as pd

import analyze_ctx_trace_sources as trace_sources
from check_origin_marker_membership import SketchCtxLookup


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-23_low_abundance_ctx_contamination"
SOURCE_REP = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_source_rep_ani.tsv"
WORK = Path("/tmp/crispatus_helv_wholegenome_hashbottomk_20260623")
FASTA_DIR = WORK / "fasta"
DELTA_DIR = WORK / "delta"
MINCO = ROOT / "minco_core/bin_hashbottomk/minco"

HEL_SOURCE_ACC = "GCF_001702095.1"
HEL_REP_ACC = "GCF_000160855.1"
CRISPATUS_REP_ACC = "GCF_018987235.1"


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def copy_fasta(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size > 0:
        return
    if src.suffix == ".gz":
        with gzip.open(src, "rt") as inp, dst.open("w") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(src, dst)


def fasta_length(path: Path) -> int:
    total = 0
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def parse_filtered_delta(path: Path) -> tuple[int, int]:
    aln_len = 0
    sim_errors = 0
    with path.open() as fh:
        for raw in fh:
            fields = raw.strip().split()
            if not fields or fields[0] == "NUCMER" or fields[0].startswith(">"):
                continue
            if len(fields) == 7:
                aln_len += abs(int(fields[1]) - int(fields[0])) + 1
                sim_errors += int(fields[4])
    return aln_len, sim_errors


def anim_pair(ref_path: Path, qry_path: Path, label: str) -> dict[str, object]:
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    DELTA_DIR.mkdir(parents=True, exist_ok=True)
    ref_fa = FASTA_DIR / f"{label}.ref.fa"
    qry_fa = FASTA_DIR / f"{label}.qry.fa"
    copy_fasta(ref_path, ref_fa)
    copy_fasta(qry_path, qry_fa)
    prefix = DELTA_DIR / label
    delta = Path(str(prefix) + ".delta")
    filt = Path(str(prefix) + ".filter")
    if not filt.exists() or filt.stat().st_size == 0:
        subprocess.run(
            ["nucmer", "--mum", "-p", str(prefix), str(ref_fa), str(qry_fa)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        with filt.open("w") as out:
            subprocess.run(["delta-filter", "-1", str(delta)], check=True, stdout=out)
    aln_len, sim_errors = parse_filtered_delta(filt)
    ref_len = fasta_length(ref_fa)
    qry_len = fasta_length(qry_fa)
    return {
        "pair": label,
        "ANIm_ANI": 1.0 - sim_errors / aln_len if aln_len else math.nan,
        "aligned_bp": aln_len,
        "similarity_errors": sim_errors,
        "ref_aligned_fraction": aln_len / ref_len if ref_len else math.nan,
        "qry_aligned_fraction": aln_len / qry_len if qry_len else math.nan,
        "ref_length": ref_len,
        "qry_length": qry_len,
    }


def load_paths() -> dict[str, Path]:
    df = pd.read_csv(SOURCE_REP, sep="\t", dtype=str)
    hel = df[df["source_accession"].eq(HEL_SOURCE_ACC)].iloc[0]
    crisp = df[df["gtdb_representative_accession"].eq(CRISPATUS_REP_ACC)].iloc[0]
    return {
        "helveticus_source": Path(hel["source_path"]),
        "helveticus_rep": Path(hel["representative_path"]),
        "crispatus_source": Path(crisp["source_path"]),
        "crispatus_rep": Path(crisp["representative_path"]),
    }


def crispatus_helveticus_events() -> pd.DataFrame:
    cached = EXP / "crispatus_helveticus_wholegenome_membership_events.tsv"
    if cached.exists() and cached.stat().st_size > 0:
        df = pd.read_csv(cached, sep="\t", dtype={"read_id": str})
        df = df[["read_id", "qctx", "diff"]].copy()
        df["qctx"] = pd.to_numeric(df["qctx"], errors="coerce").astype("int64")
        df["diff"] = pd.to_numeric(df["diff"], errors="coerce")
        return df

    truth = pd.read_csv(trace_sources.SOURCE_TRUTH, sep="\t", dtype=str)
    target = next(t for t in trace_sources.TARGETS if t["target_label"] == "crispatus_low_bad")
    needed = trace_sources.collect_needed_read_ids()
    mapping = trace_sources.load_read_mapping(needed)
    trace = pd.read_csv(trace_sources.trace_path_for(target), sep="\t", dtype={"read_id": str})
    classified = trace_sources.classify_trace(
        trace,
        mapping,
        truth,
        target["gtdb_species"],
        set(truth[truth["gtdb_species"].eq(target["gtdb_species"])]["genome_id"].astype(str)),
    )
    out = classified[classified["source_accession"].eq(HEL_SOURCE_ACC)].copy()
    out["qctx"] = pd.to_numeric(out["qctx"], errors="coerce").astype("int64")
    out["diff"] = pd.to_numeric(out["diff"], errors="coerce")
    return out


def build_sketch(paths: dict[str, Path], outdir: Path, conflict: bool) -> None:
    if (outdir / "minco.ctxobj64").exists() and (outdir / "minco.stat").exists():
        return
    if outdir.exists():
        shutil.rmtree(outdir)
    WORK.mkdir(parents=True, exist_ok=True)
    list_path = WORK / ("whole_conflict.list" if conflict else "whole_default.list")
    list_path.write_text(
        "\n".join(str(paths[k]) for k in ["helveticus_source", "helveticus_rep", "crispatus_rep"])
        + "\n"
    )
    cmd = [
        str(MINCO),
        "sketch",
        "-p",
        "4",
        "-S",
        "10000000",
        "--ctxmeta",
        "none",
        "-l",
        str(list_path),
        "-o",
        str(outdir),
    ]
    if conflict:
        cmd.insert(2, "--conflict")
    subprocess.run(cmd, check=True)


def summarize_membership(events: pd.DataFrame, default_sketch: Path, conflict_sketch: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    default = SketchCtxLookup(default_sketch)
    conflict = SketchCtxLookup(conflict_sketch)
    accs = {
        "helveticus_source": HEL_SOURCE_ACC,
        "helveticus_rep": HEL_REP_ACC,
        "crispatus_rep": CRISPATUS_REP_ACC,
    }
    rows = []
    unique_rows = []
    for _, event in events.iterrows():
        row = {
            "read_id": event["read_id"],
            "qctx": int(event["qctx"]),
            "diff": event["diff"],
        }
        for label, acc in accs.items():
            row[f"{label}_default_present"] = default.contains(acc, int(event["qctx"]))
            row[f"{label}_conflict_present"] = conflict.contains(acc, int(event["qctx"]))
            row[f"{label}_default_size"] = default.size(acc)
            row[f"{label}_conflict_size"] = conflict.size(acc)
        rows.append(row)
    detail = pd.DataFrame(rows)
    for qctx, g in detail.groupby("qctx", sort=True):
        rec = {"qctx": qctx, "event_n": len(g), "mean_diff": g["diff"].mean()}
        for label in accs:
            rec[f"{label}_default_present"] = bool(g[f"{label}_default_present"].iloc[0])
            rec[f"{label}_conflict_present"] = bool(g[f"{label}_conflict_present"].iloc[0])
        unique_rows.append(rec)
    return detail, pd.DataFrame(unique_rows)


def fraction_summary(detail: pd.DataFrame, unique: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label in ["helveticus_source", "helveticus_rep", "crispatus_rep"]:
        for mode in ["default", "conflict"]:
            col = f"{label}_{mode}_present"
            rows.append(
                {
                    "genome": label,
                    "mode": mode,
                    "event_present_n": int(detail[col].sum()),
                    "event_total": len(detail),
                    "event_present_frac": float(detail[col].mean()) if len(detail) else math.nan,
                    "unique_qctx_present_n": int(unique[col].sum()),
                    "unique_qctx_total": len(unique),
                    "unique_qctx_present_frac": float(unique[col].mean()) if len(unique) else math.nan,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    paths = load_paths()
    events = crispatus_helveticus_events()
    default_sketch = WORK / "helv_crisp_whole_default"
    conflict_sketch = WORK / "helv_crisp_whole_conflict"
    build_sketch(paths, default_sketch, conflict=False)
    build_sketch(paths, conflict_sketch, conflict=True)
    detail, unique = summarize_membership(events, default_sketch, conflict_sketch)
    summary = fraction_summary(detail, unique)
    ani = pd.DataFrame(
        [
            anim_pair(paths["helveticus_source"], paths["helveticus_rep"], "helveticus_source_vs_helveticus_rep"),
            anim_pair(paths["helveticus_source"], paths["crispatus_rep"], "helveticus_source_vs_crispatus_rep"),
            anim_pair(paths["helveticus_source"], paths["crispatus_source"], "helveticus_source_vs_crispatus_source"),
            anim_pair(paths["helveticus_rep"], paths["crispatus_rep"], "helveticus_rep_vs_crispatus_rep"),
        ]
    )
    detail.to_csv(EXP / "crispatus_helveticus_wholegenome_membership_events.tsv", sep="\t", index=False)
    unique.to_csv(EXP / "crispatus_helveticus_wholegenome_membership_unique_qctx.tsv", sep="\t", index=False)
    summary.to_csv(EXP / "crispatus_helveticus_wholegenome_membership_summary.tsv", sep="\t", index=False)
    ani.to_csv(EXP / "crispatus_helveticus_anim.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
