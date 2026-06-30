#!/usr/bin/env python3
"""Sweep simple raw-side candidate retention rules on negative routes.

The preceding failure-mode audit showed that most missed truth species are
visible in raw best-diff rows before final profile emission. This script tests
whether simple truth-blind raw-row rules can add those species back without
hurting F1 on the completed marine and CAMI3 routes.

This is diagnostic only: it does not change MinCO defaults.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable, Iterable, Mapping

import pandas as pd

import audit_negative_route_failure_modes as failure


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def f1(tp: int, fp: int, fn: int) -> float:
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    return 2.0 * precision * recall / (precision + recall) if (precision + recall) else 0.0


def truth_species(config: Mapping[str, object], sample: int) -> set[str]:
    df = pd.read_csv(Path(str(config["truth"])), sep="\t")
    sub = df.loc[df["sample"].astype(int).eq(sample)]
    return set(sub["gtdb_species"].astype(str))


def called_species(config: Mapping[str, object], sample: int, mapper: failure.GtdbMapper) -> set[str]:
    calls = failure.profile_calls(Path(config["profile"](sample)), mapper)  # type: ignore[operator]
    if calls.empty:
        return set()
    return set(calls["gtdb_species"].astype(str))


def best_raw(config: Mapping[str, object], sample: int, mapper: failure.GtdbMapper) -> pd.DataFrame:
    return failure.best_raw_rows(config["raw_paths"](sample), mapper)  # type: ignore[operator]


def truthy_major(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("major")


def build_rules() -> list[tuple[str, Callable[[pd.DataFrame], pd.Series]]]:
    rules: list[tuple[str, Callable[[pd.DataFrame], pd.Series]]] = [
        ("baseline_no_raw_add", lambda raw: pd.Series(False, index=raw.index)),
        ("raw_default_major", lambda raw: truthy_major(raw["Default_call"])),
    ]
    for ani in [0.93, 0.95, 0.97]:
        for xny in [100, 250, 500, 650]:
            for breadth in [0.05, 0.10, 0.20, 0.30, 0.50]:
                name = f"raw_ani{ani:.2f}_xny{xny}_breadth{breadth:.2f}"
                rules.append(
                    (
                        name,
                        lambda raw, ani=ani, xny=xny, breadth=breadth: (
                            (raw["ANI"] >= ani)
                            & (raw["XnY_ctx"] >= xny)
                            & (raw["Ref_breadth"] >= breadth)
                        ),
                    )
                )
    for ani in [0.93, 0.95]:
        for xny in [500, 650]:
            for breadth in [0.20, 0.30]:
                for depth in [2.0, 5.0]:
                    name = f"raw_ani{ani:.2f}_xny{xny}_breadth{breadth:.2f}_depth{depth:g}"
                    rules.append(
                        (
                            name,
                            lambda raw, ani=ani, xny=xny, breadth=breadth, depth=depth: (
                                (raw["ANI"] >= ani)
                                & (raw["XnY_ctx"] >= xny)
                                & (raw["Ref_breadth"] >= breadth)
                                & (raw["Ref_hit_mean_depth"] >= depth)
                            ),
                        )
                    )
    return rules


def score_rule(
    route: str,
    sample: int,
    truth_set: set[str],
    base_calls: set[str],
    raw: pd.DataFrame,
    rule_name: str,
    rule: Callable[[pd.DataFrame], pd.Series],
) -> dict[str, object]:
    if raw.empty:
        add_set: set[str] = set()
    else:
        mask = rule(raw).fillna(False)
        add_set = set(raw.loc[mask, "gtdb_species"].astype(str)) - base_calls
    calls = base_calls | add_set
    tp = len(calls & truth_set)
    fp = len(calls - truth_set)
    fn = len(truth_set - calls)
    base_tp = len(base_calls & truth_set)
    base_fp = len(base_calls - truth_set)
    base_fn = len(truth_set - base_calls)
    return {
        "route": route,
        "sample": sample,
        "rule": rule_name,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "F1": f1(tp, fp, fn),
        "base_TP": base_tp,
        "base_FP": base_fp,
        "base_FN": base_fn,
        "base_F1": f1(base_tp, base_fp, base_fn),
        "delta_F1": f1(tp, fp, fn) - f1(base_tp, base_fp, base_fn),
        "added_species": len(add_set),
        "added_TP": len(add_set & truth_set),
        "added_FP": len(add_set - truth_set),
    }


def route_summary(sample_scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (route, rule), sub in sample_scores.groupby(["route", "rule"], sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        base_tp = int(sub["base_TP"].sum())
        base_fp = int(sub["base_FP"].sum())
        base_fn = int(sub["base_FN"].sum())
        rows.append(
            {
                "route": route,
                "rule": rule,
                "pooled_F1": f1(tp, fp, fn),
                "base_pooled_F1": f1(base_tp, base_fp, base_fn),
                "delta_pooled_F1": f1(tp, fp, fn) - f1(base_tp, base_fp, base_fn),
                "mean_delta_F1": sub["delta_F1"].mean(),
                "worsened_samples": int((sub["delta_F1"] < -1e-12).sum()),
                "improved_samples": int((sub["delta_F1"] > 1e-12).sum()),
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "added_species": int(sub["added_species"].sum()),
                "added_TP": int(sub["added_TP"].sum()),
                "added_FP": int(sub["added_FP"].sum()),
            }
        )
    return pd.DataFrame(rows)


def overall_summary(route_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for rule, sub in route_df.groupby("rule", sort=True):
        rows.append(
            {
                "rule": rule,
                "mean_route_delta_pooled_F1": sub["delta_pooled_F1"].mean(),
                "min_route_delta_pooled_F1": sub["delta_pooled_F1"].min(),
                "route_worsen_n": int((sub["delta_pooled_F1"] < -1e-12).sum()),
                "route_improve_n": int((sub["delta_pooled_F1"] > 1e-12).sum()),
                "sample_worsen_n": int(sub["worsened_samples"].sum()),
                "sample_improve_n": int(sub["improved_samples"].sum()),
                "added_species": int(sub["added_species"].sum()),
                "added_TP": int(sub["added_TP"].sum()),
                "added_FP": int(sub["added_FP"].sum()),
                "strict_pass": bool(
                    (sub["delta_pooled_F1"] >= -1e-12).all()
                    and (sub["worsened_samples"] == 0).all()
                    and (sub["delta_pooled_F1"] > 1e-12).any()
                ),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        [
            "strict_pass",
            "min_route_delta_pooled_F1",
            "mean_route_delta_pooled_F1",
            "sample_worsen_n",
            "added_FP",
        ],
        ascending=[False, False, False, True, True],
        kind="mergesort",
    )


def build_audit(overall: pd.DataFrame) -> list[dict[str, object]]:
    strict = overall.loc[overall["strict_pass"].astype(bool)] if not overall.empty else pd.DataFrame()
    best = overall.iloc[0].to_dict() if not overall.empty else {}
    return [
        {
            "metric": "tested_rules",
            "value": int(len(overall)),
            "evidence": "raw_default_major plus ANI/XnY/breadth/depth grid",
            "decision": "diagnostic_scope",
        },
        {
            "metric": "strict_pass_rules",
            "value": int(len(strict)),
            "evidence": "requires no route or sample F1 regression on marine and CAMI3 routes",
            "decision": "candidate_found" if not strict.empty else "no_simple_rule_passes",
        },
        {
            "metric": "best_rule",
            "value": best.get("rule", ""),
            "evidence": (
                f"min_route_delta={best.get('min_route_delta_pooled_F1', '')};"
                f"mean_route_delta={best.get('mean_route_delta_pooled_F1', '')};"
                f"sample_worsen_n={best.get('sample_worsen_n', '')};"
                f"added_TP={best.get('added_TP', '')};added_FP={best.get('added_FP', '')}"
            ),
            "decision": "review_candidate" if bool(best.get("strict_pass", False)) else "diagnostic_only",
        },
        {
            "metric": "promotion_decision",
            "value": "diagnostic_only_not_default",
            "evidence": "truth-aware raw-side sweep over two completed negative routes only",
            "decision": "needs_cross_panel_wrapper_implementation_and_independent_validation"
            if not strict.empty
            else "simple_raw_side_rules_not_enough",
        },
    ]


def main() -> int:
    mapper = failure.GtdbMapper()
    rules = build_rules()
    sample_rows: list[dict[str, object]] = []
    cache: dict[tuple[str, int], tuple[set[str], set[str], pd.DataFrame]] = {}
    for config in failure.route_configs():
        route = str(config["route"])
        for sample in config["samples"]:  # type: ignore[index]
            sample = int(sample)
            key = (route, sample)
            cache[key] = (
                truth_species(config, sample),
                called_species(config, sample, mapper),
                best_raw(config, sample, mapper),
            )
            truth_set, base_calls, raw = cache[key]
            for name, rule in rules:
                sample_rows.append(score_rule(route, sample, truth_set, base_calls, raw, name, rule))

    sample_df = pd.DataFrame(sample_rows)
    route_df = route_summary(sample_df)
    overall_df = overall_summary(route_df)
    audit_rows = build_audit(overall_df)

    sample_df.to_csv(RESULTS / "raw_candidate_retention_negative_routes_sample_scores.tsv", sep="\t", index=False)
    route_df.to_csv(RESULTS / "raw_candidate_retention_negative_routes_route_summary.tsv", sep="\t", index=False)
    overall_df.to_csv(RESULTS / "raw_candidate_retention_negative_routes_overall.tsv", sep="\t", index=False)
    write_tsv(
        RESULTS / "raw_candidate_retention_negative_routes_audit.tsv",
        audit_rows,
        ["metric", "value", "evidence", "decision"],
    )
    print(overall_df.head(20).to_string(index=False))
    print("\nAUDIT")
    print(pd.DataFrame(audit_rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
