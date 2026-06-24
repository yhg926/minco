#!/usr/bin/env python3
"""Assembly source-to-GTDB-representative MinCO naive ANI checks."""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani"
WORK = EXP / "assembly_minco_naive_work"

TRUTH_C11 = EXP / "toymouse_sample0_source_rep_ani.tsv"
TRUTH_C15 = EXP / "toymouse_sample0_source_rep_ani_coden15.tsv"

C11_BIN = ROOT / "minco_core/bin/minco"
C15_BIN = ROOT / "minco_core/bin_coden15/minco"


def run(cmd: list[str], stdout_path: Path | None = None, stderr_path: Path | None = None) -> None:
    stdout = stdout_path.open("w") if stdout_path else subprocess.DEVNULL
    stderr = stderr_path.open("w") if stderr_path else subprocess.DEVNULL
    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
    finally:
        if stdout_path:
            stdout.close()
        if stderr_path:
            stderr.close()


def write_lists(df: pd.DataFrame) -> tuple[Path, Path]:
    WORK.mkdir(parents=True, exist_ok=True)
    src_list = WORK / "source_genomes_75.list"
    rep_list = WORK / "representative_genomes_unique.list"
    src_list.write_text("".join(f"{p}\n" for p in df["source_path"]))
    reps = list(dict.fromkeys(str(p) for p in df["representative_path"]))
    rep_list.write_text("".join(f"{p}\n" for p in reps))
    return src_list, rep_list


def build_and_compare(label: str, binary: Path, src_list: Path, rep_list: Path) -> Path:
    src_sketch = WORK / f"sources_{label}_S2000"
    rep_sketch = WORK / f"representatives_{label}_S2000"
    out = WORK / f"source_vs_rep_all_{label}.tsv"
    logs = WORK / "logs"
    logs.mkdir(exist_ok=True)

    if not (src_sketch / "minco.stat").exists():
        run(
            [str(binary), "sketch", "-p8", "-S", "2000", "-l", str(src_list), "-o", str(src_sketch)],
            stderr_path=logs / f"sketch_sources_{label}.stderr.log",
        )
    if not (rep_sketch / "minco.stat").exists():
        run(
            [str(binary), "sketch", "-p8", "-S", "2000", "-l", str(rep_list), "-o", str(rep_sketch)],
            stderr_path=logs / f"sketch_reps_{label}.stderr.log",
        )
    index_name = "minco.refindex.ctx64gid32obj32" if label == "coden15" else "minco.refindex.ctxgid64obj32"
    if not (rep_sketch / index_name).exists():
        run(
            [str(binary), "sketch", "-i", str(rep_sketch)],
            stderr_path=logs / f"index_reps_{label}.stderr.log",
        )
    if not out.exists() or out.stat().st_size == 0:
        run(
            [
                str(binary),
                "ani",
                "-r",
                str(rep_sketch),
                "-q",
                str(src_sketch),
                "-v",
                "-s4",
                "-m0",
                "-f0",
                "-n0",
                "-t0",
                "-o",
                str(out),
            ],
            stderr_path=logs / f"ani_{label}.stderr.log",
        )
    return out


def load_ani_by_pair(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    df = pd.read_csv(path, sep="\t")
    out: dict[tuple[str, str], dict[str, float]] = {}
    for row in df.itertuples(index=False):
        key = (str(row.Qry), str(row.Ref))
        out[key] = {
            "ANI": float(row.ANI),
            "XnY_ctx": float(row.XnY_ctx),
            "N_diff_obj": float(row.N_diff_obj),
            "N_diff_obj_section": float(row.N_diff_obj_section),
            "Ref_align_fraction": float(row.Ref_align_fraction),
            "Qry_align_fraction": float(row.Qry_align_fraction),
        }
    return out


def metrics(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    rows = []
    for col in cols:
        sub = df.loc[df["ANIm_ANI"].notna() & df[col].notna()].copy()
        err = sub[col] - sub["ANIm_ANI"]
        rows.append(
            {
                "estimator": col,
                "n": len(sub),
                "pearson": sub[col].corr(sub["ANIm_ANI"], method="pearson"),
                "spearman": sub[col].corr(sub["ANIm_ANI"], method="spearman"),
                "mae": err.abs().mean(),
                "mean_error": err.mean(),
                "median_abs_error": err.abs().median(),
                "estimator_mean": sub[col].mean(),
                "ANIm_mean": sub["ANIm_ANI"].mean(),
            }
        )
    return pd.DataFrame(rows).sort_values("mae")


def main() -> int:
    c11 = pd.read_csv(TRUTH_C11, sep="\t")
    c15 = pd.read_csv(TRUTH_C15, sep="\t")
    src_list, rep_list = write_lists(c11)

    c11_all = build_and_compare("coden11", C11_BIN, src_list, rep_list)
    c15_all = build_and_compare("coden15", C15_BIN, src_list, rep_list)

    c11_pairs = load_ani_by_pair(c11_all)
    c15_pairs = load_ani_by_pair(c15_all)

    rows = []
    c15_by_key = {
        (r.source_accession, r.gtdb_representative_accession): r
        for r in c15.itertuples(index=False)
    }
    for row in c11.itertuples(index=False):
        key_paths = (str(row.source_path), str(row.representative_path))
        key_ids = (row.source_accession, row.gtdb_representative_accession)
        c11_pair = c11_pairs.get(key_paths, {})
        c15_pair = c15_pairs.get(key_paths, {})
        c15_row = c15_by_key.get(key_ids)
        out = {
            "gtdb_species": row.gtdb_species,
            "source_accession": row.source_accession,
            "gtdb_representative_accession": row.gtdb_representative_accession,
            "source_abundance": row.source_abundance,
            "ANIm_ANI": row.ANIm_ANI,
            "asm_minco_naive_coden11": c11_pair.get("ANI", math.nan),
            "asm_minco_naive_coden15": c15_pair.get("ANI", math.nan),
            "readwise_minco_naive_coden11": row.minco_naive_ANI_calc,
            "readwise_minco_naive_coden15": getattr(c15_row, "minco_naive_ANI_calc", math.nan),
            "sylph_Adjusted_ANI": row.sylph_Adjusted_ANI,
            "sylph_Naive_ANI": row.sylph_Naive_ANI,
            "asm_XnY_coden11": c11_pair.get("XnY_ctx", math.nan),
            "asm_XnY_coden15": c15_pair.get("XnY_ctx", math.nan),
            "asm_N_diff_obj_coden11": c11_pair.get("N_diff_obj", math.nan),
            "asm_N_diff_obj_coden15": c15_pair.get("N_diff_obj", math.nan),
            "asm_ref_af_coden11": c11_pair.get("Ref_align_fraction", math.nan),
            "asm_ref_af_coden15": c15_pair.get("Ref_align_fraction", math.nan),
        }
        rows.append(out)

    out_df = pd.DataFrame(rows).sort_values(["ANIm_ANI", "source_abundance"])
    out_df.to_csv(EXP / "source_ref_assembly_minco_naive.tsv", sep="\t", index=False)

    summary = metrics(
        out_df,
        [
            "asm_minco_naive_coden11",
            "asm_minco_naive_coden15",
            "readwise_minco_naive_coden11",
            "readwise_minco_naive_coden15",
            "sylph_Adjusted_ANI",
            "sylph_Naive_ANI",
        ],
    )
    summary.to_csv(EXP / "source_ref_assembly_minco_naive_summary.tsv", sep="\t", index=False)
    print(summary.to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
