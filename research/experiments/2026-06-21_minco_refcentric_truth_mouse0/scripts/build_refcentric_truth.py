#!/usr/bin/env python3
"""Build ref-centric truth from source-vs-reference assembly ANI hits."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Sequence

import pandas as pd


ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def accession(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def read_list(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def load_source_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = {"genome_id", "abundance", "source_accession", "local_path"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"{path} missing required columns: {', '.join(missing)}")
    return df


def orient_skani(df: pd.DataFrame, source_paths: set[str], ref_paths: set[str]) -> pd.DataFrame:
    cols = set(df.columns)
    if {"Ref_file", "Query_file"}.issubset(cols):
        a_col, b_col = "Ref_file", "Query_file"
    elif {"Reference_file", "Query_file"}.issubset(cols):
        a_col, b_col = "Reference_file", "Query_file"
    else:
        raise SystemExit("skani output must contain Ref_file/Query_file columns")

    a_source = df[a_col].isin(source_paths).sum()
    b_source = df[b_col].isin(source_paths).sum()
    if b_source >= a_source:
        source_col, ref_col = b_col, a_col
    else:
        source_col, ref_col = a_col, b_col

    out = df.copy()
    out["source_path"] = out[source_col]
    out["truth_ref_path"] = out[ref_col]
    out["source_accession_from_path"] = out["source_path"].map(accession)
    out["truth_ref_accession"] = out["truth_ref_path"].map(accession)
    out["ANI"] = pd.to_numeric(out["ANI"], errors="coerce")
    for col in ["Align_fraction_ref", "Align_fraction_query"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
        else:
            out[col] = 0.0
    out["min_af"] = out[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
    out["max_af"] = out[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
    return out


def best_per_source(skani: pd.DataFrame, source_meta: pd.DataFrame) -> pd.DataFrame:
    merged = source_meta.merge(
        skani,
        left_on="local_path",
        right_on="source_path",
        how="left",
        suffixes=("", "_skani"),
    )
    merged["ANI"] = pd.to_numeric(merged["ANI"], errors="coerce").fillna(-1.0)
    merged["min_af"] = pd.to_numeric(merged["min_af"], errors="coerce").fillna(0.0)
    merged["max_af"] = pd.to_numeric(merged["max_af"], errors="coerce").fillna(0.0)
    merged = merged.sort_values(
        ["local_path", "ANI", "min_af", "max_af", "truth_ref_accession"],
        ascending=[True, False, False, False, True],
    )
    return merged.drop_duplicates("local_path", keep="first").reset_index(drop=True)


def weighted_truth(best: pd.DataFrame) -> pd.DataFrame:
    rows = best.loc[best["truth_ref_accession"].astype(bool)].copy()
    rows["abundance"] = pd.to_numeric(rows["abundance"], errors="coerce").fillna(0.0)
    grouped = rows.groupby("truth_ref_accession", as_index=False).agg(
        truth_ref_path=("truth_ref_path", "first"),
        source_count=("genome_id", "count"),
        source_accessions=("source_accession", lambda x: ",".join(sorted(set(map(str, x))))),
        source_genome_ids=("genome_id", lambda x: ",".join(sorted(set(map(str, x))))),
        abundance_sum=("abundance", "sum"),
        best_source_ani=("ANI", "max"),
        worst_source_ani=("ANI", "min"),
        min_source_min_af=("min_af", "min"),
    )
    return grouped.sort_values(["abundance_sum", "truth_ref_accession"], ascending=[False, True])


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skani", type=Path, required=True)
    ap.add_argument("--source-meta", type=Path, required=True)
    ap.add_argument("--source-list", type=Path, required=True)
    ap.add_argument("--ref-list", type=Path, required=True)
    ap.add_argument("--out-source-map", type=Path, required=True)
    ap.add_argument("--out-truth-refs", type=Path, required=True)
    ap.add_argument("--out-summary", type=Path, required=True)
    args = ap.parse_args(argv)

    source_paths = read_list(args.source_list)
    ref_paths = read_list(args.ref_list)
    source_meta = load_source_meta(args.source_meta)
    skani = orient_skani(pd.read_csv(args.skani, sep="\t"), source_paths, ref_paths)
    best = best_per_source(skani, source_meta)
    truth = weighted_truth(best)

    args.out_source_map.parent.mkdir(parents=True, exist_ok=True)
    best.to_csv(args.out_source_map, sep="\t", index=False)
    truth.to_csv(args.out_truth_refs, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "source_genomes": len(source_meta),
                "sources_with_skani_hit": int(best["truth_ref_accession"].astype(bool).sum()),
                "unique_truth_refs": int(len(truth)),
                "min_top_ani": float(best.loc[best["ANI"] >= 0, "ANI"].min()) if (best["ANI"] >= 0).any() else 0.0,
                "median_top_ani": float(best.loc[best["ANI"] >= 0, "ANI"].median()) if (best["ANI"] >= 0).any() else 0.0,
                "min_top_min_af": float(best.loc[best["ANI"] >= 0, "min_af"].min()) if (best["ANI"] >= 0).any() else 0.0,
            }
        ]
    )
    summary.to_csv(args.out_summary, sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
