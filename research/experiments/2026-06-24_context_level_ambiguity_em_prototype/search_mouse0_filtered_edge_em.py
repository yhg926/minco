#!/usr/bin/env python3
"""Search filters for mouse0 context-level ambiguity EM."""

from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
BASE = EXP / "score_mouse0_edge_em.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


base = load_module("mouse0_edge_em_base", BASE)


@dataclass(frozen=True)
class FilterConfig:
    name: str
    selected_only: bool
    max_candidate_refs: int
    max_selected_refs: int
    max_best_diff: int
    max_diff_delta: int
    max_qctx_depth: int
    max_product: int
    min_group_species: int


def normalize(values: dict[str, float]) -> dict[str, float]:
    total = sum(v for v in values.values() if v > 0.0 and math.isfinite(v))
    if total <= 0.0:
        return {}
    return {k: v / total for k, v in values.items() if v > 0.0 and math.isfinite(v)}


def enrich_edges(edge: pd.DataFrame) -> pd.DataFrame:
    edge = edge.copy()
    group_index = pd.MultiIndex.from_frame(edge[["unit_id", "qctx"]])
    edge["group_id"] = pd.factorize(group_index, sort=False)[0].astype(np.int64)
    qctx_depth = (
        edge[["unit_id", "qctx"]]
        .drop_duplicates()
        .groupby("qctx")
        .size()
        .rename("qctx_depth")
    )
    group_species = edge.groupby("group_id")["target_id"].nunique().rename("group_species_n")
    edge["qctx_depth"] = edge["qctx"].map(qctx_depth).astype(np.int32)
    edge["group_species_n"] = edge["group_id"].map(group_species).astype(np.int16)
    edge["diff_delta"] = (edge["diff"] - edge["best_diff"]).clip(lower=0).astype(np.int16)
    edge["product"] = (edge["best_diff"].astype(np.int64) * edge["qctx_depth"].astype(np.int64)).clip(
        upper=np.iinfo(np.int32).max
    ).astype(np.int32)
    return edge


def make_configs() -> list[FilterConfig]:
    big = 10**9
    configs: list[FilterConfig] = []

    def add(name: str, *, selected_only: bool = True, max_candidate_refs: int = big,
            max_selected_refs: int = big, max_best_diff: int = big,
            max_diff_delta: int = big, max_qctx_depth: int = big,
            max_product: int = big, min_group_species: int = 1) -> None:
        configs.append(
            FilterConfig(
                name=name,
                selected_only=selected_only,
                max_candidate_refs=max_candidate_refs,
                max_selected_refs=max_selected_refs,
                max_best_diff=max_best_diff,
                max_diff_delta=max_diff_delta,
                max_qctx_depth=max_qctx_depth,
                max_product=max_product,
                min_group_species=min_group_species,
            )
        )

    add("selected_raw", selected_only=True)
    add("allcand_raw", selected_only=False)
    add("selected_group_species2", selected_only=True, min_group_species=2)

    for x in [2, 3, 5, 10, 20]:
        add(f"selected_cand_le{x}", selected_only=True, max_candidate_refs=x)
    for x in [1, 2, 5]:
        add(f"selected_selrefs_le{x}", selected_only=True, max_selected_refs=x)
    for x in [0, 1, 2, 4]:
        add(f"selected_bestdiff_le{x}", selected_only=True, max_best_diff=x)
    for x in [1, 2, 5, 10, 20, 50]:
        add(f"selected_qdepth_le{x}", selected_only=True, max_qctx_depth=x)
    for x in [0, 5, 10, 20, 50]:
        add(f"selected_product_le{x}", selected_only=True, max_product=x)

    for cand in [2, 3, 5, 10]:
        for best in [0, 1, 2]:
            for qdepth in [1, 2, 5, 10, 20]:
                add(
                    f"strict_c{cand}_b{best}_qd{qdepth}",
                    selected_only=True,
                    max_candidate_refs=cand,
                    max_selected_refs=1,
                    max_best_diff=best,
                    max_qctx_depth=qdepth,
                )

    for cand in [2, 3, 5, 10]:
        for delta in [0, 1]:
            for qdepth in [1, 2, 5, 10]:
                add(
                    f"allcand_c{cand}_dd{delta}_qd{qdepth}",
                    selected_only=False,
                    max_candidate_refs=cand,
                    max_diff_delta=delta,
                    max_qctx_depth=qdepth,
                )

    # Deduplicate by full config in case labels overlap after future edits.
    out: list[FilterConfig] = []
    seen: set[tuple[object, ...]] = set()
    for cfg in configs:
        key = cfg.__dict__.copy()
        key.pop("name")
        tup = tuple(key.items())
        if tup in seen:
            continue
        seen.add(tup)
        out.append(cfg)
    return out


def apply_filter(edge: pd.DataFrame, cfg: FilterConfig) -> pd.DataFrame:
    mask = (
        (edge["candidate_refs"] <= cfg.max_candidate_refs)
        & (edge["selected_refs"] <= cfg.max_selected_refs)
        & (edge["best_diff"] <= cfg.max_best_diff)
        & (edge["diff_delta"] <= cfg.max_diff_delta)
        & (edge["qctx_depth"] <= cfg.max_qctx_depth)
        & (edge["product"] <= cfg.max_product)
        & (edge["group_species_n"] >= cfg.min_group_species)
    )
    if cfg.selected_only:
        mask &= edge["selected_by_mode"] > 0
    return edge.loc[mask].copy()


def em_counts_filtered(
    rows: pd.DataFrame,
    marker_norm: dict[str, float],
    *,
    alpha: float,
    weight_mode: str,
    iterations: int = 20,
) -> dict[str, float]:
    if rows.empty:
        return {}
    marker_species = sorted(marker_norm)
    target_to_code = {target: i for i, target in enumerate(marker_species)}
    rows = rows.loc[rows["target_id"].isin(target_to_code)].copy()
    if rows.empty:
        return {}

    target_codes = rows["target_id"].map(target_to_code).to_numpy(dtype=np.int64)
    group_codes, _ = pd.factorize(pd.MultiIndex.from_frame(rows[["unit_id", "qctx"]]), sort=False)
    n_groups = int(group_codes.max()) + 1
    diff_delta = rows["diff_delta"].to_numpy(dtype=float)
    likelihood = np.exp(-alpha * diff_delta)

    qdepth = rows["qctx_depth"].to_numpy(dtype=float)
    group_qdepth = np.zeros(n_groups, dtype=float)
    np.maximum.at(group_qdepth, group_codes, qdepth)
    if weight_mode == "event":
        group_weight = np.ones(n_groups, dtype=float)
    elif weight_mode == "qctx_equal":
        group_weight = 1.0 / np.maximum(group_qdepth, 1.0)
    elif weight_mode == "qctx_sqrt":
        group_weight = 1.0 / np.sqrt(np.maximum(group_qdepth, 1.0))
    else:
        raise ValueError(weight_mode)

    prior = np.array([max(marker_norm.get(target, 0.0), 1e-12) for target in marker_species], dtype=float)
    prior /= prior.sum()
    abundance = prior.copy()
    for _ in range(iterations):
        weights = abundance[target_codes] * likelihood
        denom = np.bincount(group_codes, weights=weights, minlength=n_groups)
        valid = denom[group_codes] > 0.0
        post = np.zeros_like(weights)
        post[valid] = weights[valid] / denom[group_codes][valid]
        weighted_post = post * group_weight[group_codes]
        new_abundance = np.bincount(target_codes, weights=weighted_post, minlength=len(marker_species))
        total = new_abundance.sum()
        if total <= 0.0:
            break
        new_abundance /= total
        if np.max(np.abs(new_abundance - abundance)) < 1e-10:
            abundance = new_abundance
            break
        abundance = new_abundance
    return {target: float(abundance[i]) for i, target in enumerate(marker_species) if abundance[i] > 0.0}


def score_edge_filter_search() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    marker_raw, marker_norm, truth, _, _ = base.load_marker_values()
    edge = pd.read_csv(RESULTS / "mouse0_marker_candidate_edge_species.tsv", sep="\t")
    edge = enrich_edges(edge)

    diagnostics = []
    qstats = (
        edge[["qctx", "unit_id", "best_diff"]]
        .drop_duplicates(["qctx", "unit_id"])
        .groupby("qctx")
        .agg(event_depth=("unit_id", "size"), min_best_diff=("best_diff", "min"))
        .reset_index()
    )
    qstats["product"] = qstats["event_depth"] * qstats["min_best_diff"]
    diagnostics.append({
        "metric": "qctx_count",
        "value": len(qstats),
    })
    for col in ["event_depth", "product"]:
        for q in [0.5, 0.9, 0.99, 0.999]:
            diagnostics.append({
                "metric": f"{col}_q{q}",
                "value": float(qstats[col].quantile(q)),
            })
        diagnostics.append({
            "metric": f"{col}_max",
            "value": float(qstats[col].max()),
        })

    rows = []
    baseline = base.score_prediction(truth, marker_norm, "marker_l1_per_read_profile")
    baseline.update({
        "filter_id": "baseline",
        "edge_mode": "none",
        "alpha": math.nan,
        "weight_mode": "none",
        "beta": 0.0,
        "filtered_rows": 0,
        "filtered_groups": 0,
        "filtered_species": 0,
    })
    rows.append(baseline)

    configs = make_configs()
    alpha_values = [0.0, 1.0]
    weight_modes = ["event", "qctx_equal", "qctx_sqrt"]
    beta_values = [0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 1.0]
    for i, cfg in enumerate(configs, start=1):
        filtered = apply_filter(edge, cfg)
        group_n = filtered.groupby(["unit_id", "qctx"]).ngroups if len(filtered) else 0
        species_n = int(filtered["target_id"].nunique()) if len(filtered) else 0
        if group_n < 10 or species_n < 2:
            continue
        for alpha in alpha_values:
            for weight_mode in weight_modes:
                edge_norm = em_counts_filtered(filtered, marker_norm, alpha=alpha, weight_mode=weight_mode)
                edge_score = base.score_prediction(
                    truth,
                    edge_norm,
                    f"edge_only_{cfg.name}_a{alpha:g}_{weight_mode}",
                )
                edge_score.update({
                    "filter_id": cfg.name,
                    "edge_mode": "edge_only",
                    "alpha": alpha,
                    "weight_mode": weight_mode,
                    "beta": 1.0,
                    "filtered_rows": len(filtered),
                    "filtered_groups": group_n,
                    "filtered_species": species_n,
                    **{f"cfg_{k}": v for k, v in cfg.__dict__.items() if k != "name"},
                })
                rows.append(edge_score)
                for beta in beta_values:
                    pred = base.blend(marker_norm, edge_norm, beta)
                    score = base.score_prediction(
                        truth,
                        pred,
                        f"blend_{cfg.name}_a{alpha:g}_{weight_mode}_b{beta:g}",
                    )
                    score.update({
                        "filter_id": cfg.name,
                        "edge_mode": "blend",
                        "alpha": alpha,
                        "weight_mode": weight_mode,
                        "beta": beta,
                        "filtered_rows": len(filtered),
                        "filtered_groups": group_n,
                        "filtered_species": species_n,
                        **{f"cfg_{k}": v for k, v in cfg.__dict__.items() if k != "name"},
                    })
                    rows.append(score)
        if i % 25 == 0:
            print(f"scored filters {i}/{len(configs)}", file=sys.stderr)

    scored = pd.DataFrame(rows).sort_values(["l1_pct_points", "F1"], ascending=[True, False])
    diag = pd.DataFrame(diagnostics)
    return scored, diag, edge


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    scored, diag, edge = score_edge_filter_search()
    scored.to_csv(RESULTS / "mouse0_filtered_edge_em_grid.tsv", sep="\t", index=False)
    diag.to_csv(RESULTS / "mouse0_edge_filter_diagnostics.tsv", sep="\t", index=False)
    scored.head(50).to_csv(RESULTS / "mouse0_filtered_edge_em_top50.tsv", sep="\t", index=False)
    scored.head(20).to_csv(EXP / "filtered_summary.tsv", sep="\t", index=False)

    # Species-level diagnostic for the best edge-only row and marker baseline.
    best = scored.iloc[0]
    pd.DataFrame([best]).to_csv(RESULTS / "mouse0_filtered_edge_em_best.tsv", sep="\t", index=False)
    print(scored.head(30).to_csv(sep="\t", index=False), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
