#!/usr/bin/env python3
"""Refine the blend weight for the best mouse0 filtered edge-EM row."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
SEARCH = EXP / "search_mouse0_filtered_edge_em.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


search = load_module("mouse0_filtered_edge_em_search", SEARCH)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    _, marker_norm, truth, _, _ = search.base.load_marker_values()
    edge = pd.read_csv(RESULTS / "mouse0_marker_candidate_edge_species.tsv", sep="\t")
    edge = search.enrich_edges(edge)
    cfg = next(c for c in search.make_configs() if c.name == "selected_group_species2")
    filtered = search.apply_filter(edge, cfg)
    edge_norm = search.em_counts_filtered(filtered, marker_norm, alpha=0.0, weight_mode="event")

    beta_values = [
        0.0,
        0.0001,
        0.0003,
        0.0005,
        0.0007,
        0.001,
        0.0015,
        0.002,
        0.0025,
        0.003,
        0.0035,
        0.004,
        0.0045,
        0.005,
        0.006,
        0.008,
        0.01,
        0.02,
    ]
    rows = []
    for beta in beta_values:
        pred = search.base.blend(marker_norm, edge_norm, beta)
        score = search.base.score_prediction(
            truth,
            pred,
            f"blend_selected_group_species2_a0_event_b{beta:g}",
        )
        score.update(
            {
                "filter_id": cfg.name,
                "edge_mode": "blend",
                "alpha": 0.0,
                "weight_mode": "event",
                "beta": beta,
                "filtered_rows": len(filtered),
                "filtered_groups": filtered.groupby(["unit_id", "qctx"]).ngroups,
                "filtered_species": int(filtered["target_id"].nunique()),
            }
        )
        rows.append(score)

    out = pd.DataFrame(rows).sort_values(["l1_pct_points", "F1"], ascending=[True, False])
    out.to_csv(RESULTS / "mouse0_filtered_edge_em_beta_refine.tsv", sep="\t", index=False)
    print(out.head(12).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
