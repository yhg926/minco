#!/usr/bin/env python3
"""Validate the integrated accession-level candidate surface on cached HMP panels.

This replays saved unique/split candidate tables through
scripts/minco_profile_calibrated.py with the off-by-default
--candidate-surface-switch. The normal taxmap remains the model feature map; an
accession-level GTDB taxmap is supplied only for candidate-surface labels.
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import decompose_abundance_errors as decomp
import sweep_hmp_raw_candidate_rescue as raw_rescue


EXP = Path(__file__).resolve().parent
REPO = EXP.parents[2]
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_candidate_surface_integrated_hmp_20260629")
NORMAL_TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
SURFACE_TAXMAP = Path("/tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv")
MODEL_CACHE = Path("/tmp/minco_model_cache_bench_20260627/minco_profile_rf_hgb.train12.unfiltered.joblib")
METHOD = "minco_integrated_candidate_surface_hmp"
OFFLINE_METHOD = "cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass"
PANELS = [
    "hmp_airskin_gtdb_source_abundance",
    "hmp_gastrooral_gtdb_source_abundance",
]

# Current HMP airskin samples6/11 were refreshed from saved table-mode inputs
# and only retained the exact split table in their current output directory.
TABLE_OVERRIDES = {
    ("hmp_airskin_gtdb_source_abundance", 6): {
        "unique": Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample6_unique_zip_unfiltered.tsv"),
        "split": Path("/tmp/minco_current_code_hmp_refresh_20260627/minco_sample6_split_exact_unfiltered.tsv"),
    },
    ("hmp_airskin_gtdb_source_abundance", 11): {
        "unique": Path("/tmp/cami2_hmp_unseen_transfer_20260626/run/minco_sample11_unique_zip_unfiltered.tsv"),
        "split": Path("/tmp/minco_current_code_hmp_refresh_20260627/minco_sample11_split_exact_unfiltered.tsv"),
    },
}


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def profile_metadata(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        row = pd.read_csv(path, sep="\t", nrows=1).iloc[0].to_dict()
    except Exception:
        return {}
    return row


def discover_tables(panel: str, sample: int, profile: Path) -> dict[str, object]:
    override = TABLE_OVERRIDES.get((panel, sample))
    if override:
        return {
            "unique": override["unique"],
            "split": override["split"],
            "table_source": "override_saved_table",
        }
    meta = profile_metadata(profile)
    parent = profile.parent
    candidates = []
    for subdir in ["work", f"sample{sample}_work"]:
        base = parent / subdir
        candidates.append(
            {
                "unique": base / "minco.best_diff_unique.unfiltered.tsv",
                "split": base / "minco.best_diff_split.unfiltered.tsv",
                "exact": base / "minco.best_diff_split.exact.unfiltered.tsv",
                "source": subdir,
            }
        )
    unique = None
    split = None
    exact = None
    source = ""
    for cand in candidates:
        if unique is None and cand["unique"].exists():
            unique = cand["unique"]
            source = cand["source"]
        if split is None and cand["split"].exists():
            split = cand["split"]
        if exact is None and cand["exact"].exists():
            exact = cand["exact"]
    exact_path_raw = str(meta.get("auto_exact_split_path", "") or "")
    if exact_path_raw and exact_path_raw.lower() != "nan" and Path(exact_path_raw).exists():
        exact = Path(exact_path_raw)
    exact_used = str(meta.get("auto_exact_split_used", "")).strip().lower() == "true"
    chosen_split = exact if exact_used and exact is not None else split
    return {
        "unique": unique,
        "split": chosen_split,
        "table_source": source,
        "auto_exact_split_used": exact_used,
    }


def parse_time_log(path: Path) -> dict[str, object]:
    text = path.read_text() if path.exists() else ""
    wall = ""
    rss = ""
    status = ""
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            wall = line.rsplit("):", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            rss = line.rsplit(":", 1)[1].strip()
        elif "Exit status" in line:
            status = line.rsplit(":", 1)[1].strip()
    return {"wall_clock": wall, "max_rss_kb": rss, "exit_status": status}


def run_wrapper(panel: str, sample: int, unique: Path, split: Path, out: Path) -> dict[str, object]:
    time_log = OUT_DIR / f"{panel}.sample{sample}.time.log"
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
        str(NORMAL_TAXMAP),
        "--model-cache",
        str(MODEL_CACHE),
        "--strategy",
        "universal",
        "--scope",
        "bacteria",
        "--train-pool",
        "train12",
        "--candidate-surface-switch",
        "accession-ani90-xny100-br01-af70",
        "--candidate-surface-taxmap",
        str(SURFACE_TAXMAP),
        "-o",
        str(out),
    ]
    start = time.monotonic()
    subprocess.run(cmd, cwd=REPO, check=True)
    elapsed = time.monotonic() - start
    timing = parse_time_log(time_log)
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


def zero_mass_audit(profile: Path) -> dict[str, object]:
    df = pd.read_csv(
        profile,
        sep="\t",
        usecols=lambda col: col
        in {
            "candidate_surface_added",
            "calibrated_abundance",
            "calibrated_abundance_raw",
            "species_name",
            "candidate_surface_accession",
        },
        low_memory=False,
    )
    added = bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
    abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
    raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
    return {
        "added_rows": int(added.sum()),
        "zero_mass_violations": int(((abundance.abs() > 1e-15) | (raw.abs() > 1e-15))[added].sum()),
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
        "integrated wrapper accession-level zero-mass candidate surface",
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
                "zero_mass_violations": int(sub["zero_mass_violations"].sum()),
            }
        )
    return pd.DataFrame(rows)


def compare(scores: pd.DataFrame, path: Path, method: str, label: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            [{"panel": "", "sample": "", "comparison": label, "status": "missing", "path": str(path)}]
        )
    other = pd.read_csv(path, sep="\t")
    rows = []
    for row in scores.itertuples(index=False):
        panel = str(row.panel)
        sample = int(row.sample)
        exp = other.loc[
            other["panel"].astype(str).eq(panel)
            & other["sample"].astype(int).eq(sample)
            & other["method"].astype(str).eq(method)
        ]
        if exp.empty:
            rows.append(
                {
                    "panel": panel,
                    "sample": sample,
                    "comparison": label,
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
                "comparison": label,
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


def audit(summary: pd.DataFrame, vs_post: pd.DataFrame, vs_offline: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "metric": "profiles_replayed",
            "value": int(summary["sample_count"].sum()) if not summary.empty else 0,
            "evidence": str(OUT_DIR),
            "decision": "all_cached_hmp_integrated_replay",
        },
        {
            "metric": "added_rows",
            "value": int(summary["added_rows"].sum()) if not summary.empty else 0,
            "evidence": "candidate_surface_added rows",
            "decision": "candidate_surface_behavior",
        },
        {
            "metric": "zero_mass_violations",
            "value": int(summary["zero_mass_violations"].sum()) if not summary.empty else 0,
            "evidence": "calibrated_abundance and calibrated_abundance_raw on added rows",
            "decision": "pass" if int(summary["zero_mass_violations"].sum()) == 0 else "fail",
        },
    ]
    for label, evidence_stem, frame in [
        ("postprocessor", "postprocessor", vs_post),
        ("offline_raw_cache_rule", "offline", vs_offline),
    ]:
        compared = frame.loc[frame["status"].astype(str).eq("compared")].copy()
        if compared.empty:
            value = "missing"
            decision = "missing_comparison"
        else:
            max_count_delta = int(compared[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max())
            max_f1_delta = float(compared["F1_delta"].abs().max())
            value = f"counts={max_count_delta};F1={max_f1_delta:.3g}"
            decision = "matches" if max_count_delta == 0 and max_f1_delta < 1e-12 else "differs"
        rows.append(
            {
                "metric": f"max_delta_vs_{label}",
                "value": value,
                "evidence": f"candidate_surface_integrated_hmp_vs_{evidence_stem}.tsv",
                "decision": decision,
            }
        )
    rows.append(
        {
            "metric": "promotion_decision",
            "value": "experimental_not_default",
            "evidence": "zero-mass F1-only surface; nonzero abundance policy not validated",
            "decision": "do_not_change_current_default",
        }
    )
    return pd.DataFrame(rows)


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
            table_info = discover_tables(panel, sample, profile)
            unique = table_info.get("unique")
            split = table_info.get("split")
            out = OUT_DIR / f"{panel}.sample{sample}.integrated_candidate_surface.tsv"
            apply_row = {
                "panel": panel,
                "sample": sample,
                "profile": str(profile),
                "unique": str(unique or ""),
                "split": str(split or ""),
                "out": str(out),
                "table_source": table_info.get("table_source", ""),
                "auto_exact_split_used": table_info.get("auto_exact_split_used", ""),
            }
            if not unique or not Path(unique).exists() or not split or not Path(split).exists():
                apply_row.update({"status": "missing_tables", "added_rows": "", "zero_mass_violations": ""})
                missing_rows.append(apply_row.copy())
                apply_rows.append(apply_row)
                continue
            timing_rows.append(run_wrapper(panel, sample, Path(unique), Path(split), out))
            zero = zero_mass_audit(out)
            apply_row.update({"status": "replayed", **zero})
            apply_rows.append(apply_row)
            score = score_profile(panel, sample, out, mapper)
            score.update(zero)
            score["profile"] = str(out)
            score_rows.append(score)
    apply_df = pd.DataFrame(apply_rows)
    scores = pd.DataFrame(score_rows)
    summary = summarize(scores) if not scores.empty else pd.DataFrame()
    timing = pd.DataFrame(timing_rows)
    missing = pd.DataFrame(missing_rows)
    vs_post = compare(
        scores,
        RESULTS / "raw_side_candidate_surface_hmp_scores.tsv",
        "raw_side_candidate_surface_hmp",
        "postprocessor",
    )
    vs_offline = compare(
        scores,
        Path("/tmp/minco_cross_panel_candidate_rescue_scores.tsv"),
        OFFLINE_METHOD,
        "offline_raw_cache_rule",
    )
    audit_df = audit(summary, vs_post, vs_offline)

    apply_df.to_csv(RESULTS / "candidate_surface_integrated_hmp_apply.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "candidate_surface_integrated_hmp_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "candidate_surface_integrated_hmp_summary.tsv", sep="\t", index=False)
    timing.to_csv(RESULTS / "candidate_surface_integrated_hmp_runtime.tsv", sep="\t", index=False)
    missing.to_csv(RESULTS / "candidate_surface_integrated_hmp_missing.tsv", sep="\t", index=False)
    vs_post.to_csv(RESULTS / "candidate_surface_integrated_hmp_vs_postprocessor.tsv", sep="\t", index=False)
    vs_offline.to_csv(RESULTS / "candidate_surface_integrated_hmp_vs_offline.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / "candidate_surface_integrated_hmp_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    if not missing.empty:
        print("\nMISSING")
        print(missing.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
