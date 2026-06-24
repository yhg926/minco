#!/usr/bin/env python3
"""Score minco/Sylph profiles against CAMI marine species gold profiles."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, Mapping, Set

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[4]
PREV_SCRIPT_DIR = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(PREV_SCRIPT_DIR))

from analyze_readwise_corrections import (  # noqa: E402
    load_minco,
    load_sylph,
    parse_gold_profile,
    parse_species_taxmap,
    score_taxids,
)


def _safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _truth_map(gold: pd.DataFrame) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for _, row in gold.iterrows():
        taxid = str(row["taxid"])
        out[taxid] = _safe_float(row["gold_percentage"]) / 100.0
    return out


def _pred_abundance(rows: pd.DataFrame, taxids: Iterable[str], column: str) -> Dict[str, float]:
    pred: Dict[str, float] = {}
    if column not in rows.columns:
        return pred
    subset = rows.loc[rows["taxid"].astype(bool)].copy()
    subset[column] = pd.to_numeric(subset[column], errors="coerce").fillna(0.0)
    for taxid in taxids:
        vals = subset.loc[subset["taxid"].astype(str) == taxid, column]
        pred[taxid] = float(vals.max()) if len(vals) else 0.0
    if column in {"Taxonomic_abundance", "Sequence_abundance"} and max(pred.values(), default=0.0) > 1.0:
        pred = {k: v / 100.0 for k, v in pred.items()}
    total = sum(v for v in pred.values() if math.isfinite(v) and v > 0.0)
    if total > 0.0 and column not in {"Normalized_abundance_depth", "Taxonomic_abundance"}:
        pred = {k: v / total for k, v in pred.items()}
    return pred


def _abundance_stats(gold: pd.DataFrame, rows: pd.DataFrame, pred_taxids: Set[str], column: str) -> Mapping[str, object]:
    gold_ab = _truth_map(gold)
    tp_taxids = sorted(pred_taxids & set(gold_ab))
    pred_ab = _pred_abundance(rows, tp_taxids, column)
    truth = np.array([gold_ab[t] for t in tp_taxids], dtype=float)
    pred = np.array([pred_ab.get(t, 0.0) for t in tp_taxids], dtype=float)
    mask = np.isfinite(truth) & np.isfinite(pred)
    truth = truth[mask]
    pred = pred[mask]
    if len(truth) >= 2 and float(np.std(truth)) > 0.0 and float(np.std(pred)) > 0.0:
        pearson = float(np.corrcoef(truth, pred)[0, 1])
    else:
        pearson = float("nan")
    mae = float(np.mean(np.abs(pred - truth))) if len(truth) else float("nan")
    l1 = float(np.sum(np.abs(pred - truth))) if len(truth) else float("nan")
    return {
        "abundance_column": column,
        "tp_abundance_n": int(len(truth)),
        "tp_abundance_pearson": pearson,
        "tp_abundance_mae": mae,
        "tp_abundance_l1": l1,
    }


def _score_rows(
    label: str,
    sample_id: str,
    gold: pd.DataFrame,
    rows: pd.DataFrame,
    abundance_column: str,
) -> Dict[str, object]:
    gold_taxids = set(gold["taxid"].astype(str))
    selected = rows.loc[rows["taxid"].astype(bool)]
    pred_taxids = set(selected["taxid"].astype(str))
    stats = score_taxids(pred_taxids, gold_taxids)
    out: Dict[str, object] = {
        "sample_id": sample_id,
        "label": label,
        "gold_taxa": len(gold_taxids),
        "pred_rows": int(len(selected)),
        **stats,
    }
    out.update(_abundance_stats(gold, rows, pred_taxids, abundance_column))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--minco", type=Path)
    ap.add_argument("--sylph", type=Path)
    ap.add_argument("--gold", type=Path, default=Path("/tmp/gs_marine_short.profile"))
    ap.add_argument("--taxmap", type=Path, default=Path("/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ctx-k", type=float, default=11.0)
    args = ap.parse_args()

    gold = parse_gold_profile(args.gold, args.sample_id)
    taxmap = parse_species_taxmap(args.taxmap)
    records = []

    if args.minco:
        rows = load_minco(args.minco, taxmap, args.ctx_k)
        records.append(_score_rows("minco", args.sample_id, gold, rows, "Normalized_abundance_depth"))

    if args.sylph:
        rows = load_sylph(args.sylph, taxmap)
        records.append(_score_rows("sylph", args.sample_id, gold, rows, "Taxonomic_abundance"))

    if not records:
        raise SystemExit("provide --minco, --sylph, or both")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(args.out, sep="\t", index=False)
    print(pd.DataFrame(records).to_string(index=False))


if __name__ == "__main__":
    main()
