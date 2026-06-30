#!/usr/bin/env python3
"""Validate nonzero candidate-surface abundance in the calibrated wrapper.

This replays the same cached HMP unique/split tables used by the integrated
candidate-surface validation, but enables
--candidate-abundance-policy normalized-depth-alpha2. The expected comparison
is the fixed-call postprocessor sweep method normalized_depth_alpha2.
"""

from __future__ import annotations

import math
import subprocess
import time
from pathlib import Path

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp
import sweep_hmp_raw_candidate_rescue as raw_rescue
import validate_integrated_candidate_surface_hmp as zero_surface


EXP = Path(__file__).resolve().parent
REPO = EXP.parents[2]
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_candidate_surface_integrated_hmp_abundance_20260629")
METHOD = "minco_integrated_candidate_surface_hmp_normalized_depth_alpha2"
EXPECTED_METHOD = "normalized_depth_alpha2"
POLICY = "normalized-depth-alpha2"
PANELS = zero_surface.PANELS


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def run_wrapper(panel: str, sample: int, unique: Path, split: Path, out: Path) -> dict[str, object]:
    time_log = OUT_DIR / f"{panel}.sample{sample}.abundance.time.log"
    cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_log),
        "python3",
        "-B",
        "scripts/minco_profile_calibrated.py",
        "--unique-table",
        str(unique),
        "--split-table",
        str(split),
        "--taxmap",
        str(zero_surface.NORMAL_TAXMAP),
        "--model-cache",
        str(zero_surface.MODEL_CACHE),
        "--strategy",
        "universal",
        "--scope",
        "bacteria",
        "--train-pool",
        "train12",
        "--candidate-surface-switch",
        "accession-ani90-xny100-br01-af70",
        "--candidate-abundance-policy",
        POLICY,
        "--candidate-surface-taxmap",
        str(zero_surface.SURFACE_TAXMAP),
        "-o",
        str(out),
    ]
    start = time.monotonic()
    subprocess.run(cmd, cwd=REPO, check=True)
    elapsed = time.monotonic() - start
    timing = zero_surface.parse_time_log(time_log)
    timing.update(
        {
            "panel": panel,
            "sample": sample,
            "command": " ".join(cmd),
            "elapsed_s": f"{elapsed:.3f}",
            "time_log": str(time_log),
        }
    )
    return timing


def abundance_policy_audit(profile: Path) -> dict[str, object]:
    df = pd.read_csv(
        profile,
        sep="\t",
        usecols=lambda col: col
        in {
            "candidate_surface_added",
            "calibrated_call",
            "calibrated_abundance",
            "calibrated_abundance_raw",
            "species_name",
            "candidate_surface_accession",
            "candidate_abundance_policy",
            "s_Normalized_abundance_depth_max",
            "u_Normalized_abundance_depth_max",
        },
        low_memory=False,
    )
    added = zero_surface.bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    called = zero_surface.bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
    abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    split_norm = pd.to_numeric(
        df.get("s_Normalized_abundance_depth_max", 0.0),
        errors="coerce",
    ).fillna(0.0)
    unique_norm = pd.to_numeric(
        df.get("u_Normalized_abundance_depth_max", 0.0),
        errors="coerce",
    ).fillna(0.0)
    base_called_raw_mass = float(raw.loc[called & ~added].sum())
    if base_called_raw_mass <= 0.0:
        base_called_raw_mass = 1.0
    expected_raw = (
        2.0 * pd.concat([split_norm, unique_norm], axis=1).max(axis=1) * base_called_raw_mass
    )
    raw_delta = (raw - expected_raw).abs()
    added_policy = (
        ",".join(sorted(set(df.loc[added, "candidate_abundance_policy"].astype(str))))
        if "candidate_abundance_policy" in df.columns
        else ""
    )
    return {
        "added_rows": int(added.sum()),
        "nonzero_added_rows": int((raw.loc[added] > 0.0).sum()),
        "added_raw_mass": float(raw.loc[added].sum()),
        "added_abundance_mass": float(abundance.loc[added].sum()),
        "base_called_raw_mass": base_called_raw_mass,
        "max_added_raw_mass_delta": float(raw_delta.loc[added].max()) if added.any() else 0.0,
        "added_policy_values": added_policy,
        "added_species": ",".join(df.loc[added, "species_name"].astype(str).sort_values().tolist()),
        "added_accessions": ",".join(
            df.loc[added, "candidate_surface_accession"].astype(str).sort_values().tolist()
        ),
    }


def score_profile(panel: str, sample: int, profile: Path, mapper: raw_audit.RawMapper) -> dict[str, object]:
    cfg = decomp.PANELS[panel]
    truth = decomp.load_truth(
        Path(cfg["truth"](sample)),
        sample,
        str(cfg["truth_abundance_col"]),
    )
    pred = raw_rescue.load_current_pred_gtdb(
        panel,
        profile,
        str(cfg["minco_collapse"]),
        mapper,
    )
    return raw_rescue.score(
        panel,
        sample,
        METHOD,
        truth,
        set(pred),
        pred,
        "integrated wrapper candidate-surface normalized-depth-alpha2 abundance policy",
    )


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for panel, sub in scores.groupby("panel", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append(
            {
                "panel": panel,
                "method": METHOD,
                "samples": ",".join(map(str, sorted(sub["sample"].astype(int).unique()))),
                "sample_count": int(len(sub)),
                "mean_F1": sub["F1"].mean(),
                "pooled_F1": pooled_f1,
                "pooled_TP": tp,
                "pooled_FP": fp,
                "pooled_FN": fn,
                "mean_L1_union_pp": sub["L1_union_pp"].mean(),
                "mean_Pearson_union": sub["Pearson_union"].mean(),
                "added_rows": int(sub["added_rows"].sum()),
                "nonzero_added_rows": int(sub["nonzero_added_rows"].sum()),
                "base_called_raw_mass": sub["base_called_raw_mass"].sum(),
                "added_raw_mass": sub["added_raw_mass"].sum(),
                "added_abundance_mass": sub["added_abundance_mass"].sum(),
                "max_added_raw_mass_delta": sub["max_added_raw_mass_delta"].max(),
            }
        )
    return pd.DataFrame(rows)


def compare_to_policy_sweep(scores: pd.DataFrame) -> pd.DataFrame:
    path = RESULTS / "candidate_surface_abundance_policy_scores.tsv"
    if not path.exists():
        return pd.DataFrame(
            [{"panel": "", "sample": "", "status": "missing", "path": str(path)}]
        )
    expected = pd.read_csv(path, sep="\t")
    rows = []
    for row in scores.itertuples(index=False):
        panel = str(row.panel)
        sample = int(row.sample)
        exp = expected.loc[
            expected["panel"].astype(str).eq(panel)
            & expected["sample"].astype(int).eq(sample)
            & expected["method"].astype(str).eq(EXPECTED_METHOD)
        ]
        if exp.empty:
            rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "status": "missing_expected",
                    "path": str(path),
                }
            )
            continue
        exp_row = exp.iloc[0]
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "status": "compared",
                "TP_delta": int(row.TP) - int(exp_row["TP"]),
                "FP_delta": int(row.FP) - int(exp_row["FP"]),
                "FN_delta": int(row.FN) - int(exp_row["FN"]),
                "F1_delta": finite(row.F1) - finite(exp_row["F1"]),
                "L1_delta_pp": finite(row.L1_union_pp) - finite(exp_row["L1_union_pp"]),
                "Pearson_delta": finite(row.Pearson_union) - finite(exp_row["Pearson_union"]),
            }
        )
    return pd.DataFrame(rows)


def audit(summary: pd.DataFrame, vs_policy: pd.DataFrame) -> pd.DataFrame:
    compared = vs_policy.loc[vs_policy["status"].astype(str).eq("compared")].copy()
    if compared.empty:
        comparison_value = "missing"
        comparison_decision = "missing_comparison"
    else:
        max_count_delta = int(compared[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max())
        max_l1_delta = float(compared["L1_delta_pp"].abs().max())
        max_pearson_delta = float(compared["Pearson_delta"].abs().max())
        comparison_value = (
            f"counts={max_count_delta};L1={max_l1_delta:.3g};Pearson={max_pearson_delta:.3g}"
        )
        comparison_decision = (
            "matches"
            if max_count_delta == 0 and max_l1_delta < 1e-9 and max_pearson_delta < 1e-12
            else "differs"
        )
    total_profiles = int(summary["sample_count"].sum()) if not summary.empty else 0
    total_added = int(summary["added_rows"].sum()) if not summary.empty else 0
    total_nonzero = int(summary["nonzero_added_rows"].sum()) if not summary.empty else 0
    max_raw_delta = float(summary["max_added_raw_mass_delta"].max()) if not summary.empty else 0.0
    return pd.DataFrame(
        [
            {
                "metric": "profiles_replayed",
                "value": total_profiles,
                "evidence": str(OUT_DIR),
                "decision": "cached_hmp_integrated_replay",
            },
            {
                "metric": "added_rows",
                "value": total_added,
                "evidence": "candidate_surface_added rows",
                "decision": "candidate_surface_behavior",
            },
            {
                "metric": "nonzero_added_rows",
                "value": total_nonzero,
                "evidence": POLICY,
                "decision": "abundance_policy_behavior",
            },
            {
                "metric": "max_added_raw_mass_delta",
                "value": max_raw_delta,
                "evidence": "calibrated_abundance_raw versus 2*max(s,u Normalized_abundance_depth)*base_called_raw_mass",
                "decision": "pass" if max_raw_delta < 1e-12 else "fail",
            },
            {
                "metric": "max_delta_vs_fixed_call_policy_sweep",
                "value": comparison_value,
                "evidence": "candidate_surface_integrated_hmp_abundance_vs_policy_sweep.tsv",
                "decision": comparison_decision,
            },
            {
                "metric": "promotion_decision",
                "value": "experimental_not_default",
                "evidence": "HMP wrapper validation only; cross-panel wrapper validation still incomplete",
                "decision": "do_not_change_current_default",
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mapper = raw_audit.RawMapper()
    apply_rows = []
    score_rows = []
    timing_rows = []
    missing_rows = []
    for panel in PANELS:
        cfg = decomp.PANELS[panel]
        for sample_raw in cfg["samples"]:
            sample = int(sample_raw)
            profile = Path(cfg["minco"](sample))
            table_info = zero_surface.discover_tables(panel, sample, profile)
            unique = table_info.get("unique")
            split = table_info.get("split")
            out = OUT_DIR / f"{panel}.sample{sample}.integrated_candidate_surface_abundance.tsv"
            apply_row = {
                "panel": panel,
                "sample": sample,
                "profile": str(profile),
                "unique": str(unique or ""),
                "split": str(split or ""),
                "out": str(out),
                "table_source": table_info.get("table_source", ""),
                "auto_exact_split_used": table_info.get("auto_exact_split_used", ""),
                "candidate_abundance_policy": POLICY,
            }
            if not unique or not Path(unique).exists() or not split or not Path(split).exists():
                apply_row.update({"status": "missing_tables", "added_rows": "", "nonzero_added_rows": ""})
                missing_rows.append(apply_row.copy())
                apply_rows.append(apply_row)
                continue
            timing_rows.append(run_wrapper(panel, sample, Path(unique), Path(split), out))
            mass_audit = abundance_policy_audit(out)
            apply_row.update({"status": "replayed", **mass_audit})
            apply_rows.append(apply_row)
            score = score_profile(panel, sample, out, mapper)
            score.update(mass_audit)
            score["profile"] = str(out)
            score_rows.append(score)
    apply_df = pd.DataFrame(apply_rows)
    scores = pd.DataFrame(score_rows)
    summary = summarize(scores) if not scores.empty else pd.DataFrame()
    timing = pd.DataFrame(timing_rows)
    missing = pd.DataFrame(missing_rows)
    vs_policy = compare_to_policy_sweep(scores)
    audit_df = audit(summary, vs_policy)

    apply_df.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_apply.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_summary.tsv", sep="\t", index=False)
    timing.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_runtime.tsv", sep="\t", index=False)
    missing.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_missing.tsv", sep="\t", index=False)
    vs_policy.to_csv(
        RESULTS / "candidate_surface_integrated_hmp_abundance_vs_policy_sweep.tsv",
        sep="\t",
        index=False,
    )
    audit_df.to_csv(RESULTS / "candidate_surface_integrated_hmp_abundance_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    if not missing.empty:
        print("\nMISSING")
        print(missing.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
