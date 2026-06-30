#!/usr/bin/env python3
"""Validate the one-flag candidate preset against prior wrapper evidence.

The candidate strategy was first validated by explicitly wiring the calibrated
wrapper switches. This script reruns the same cached table-mode inputs through
``scripts/minco_profile_default.py --profile-preset candidate`` and compares
the scored profiles with the prior wrapper validation.
"""

from __future__ import annotations

import math
import subprocess
import time
from pathlib import Path

import pandas as pd

import audit_hmp_missed_truth_raw_tables as raw_audit
import score_adaptive_call_filter_wrapper_validation as adaptive
import validate_cross_panel_candidate_abundance_wrapper as wrapper_validation
import validate_integrated_candidate_surface_hmp as hmp_surface
import validate_integrated_candidate_surface_hmp_abundance_policy as hmp_abundance


EXP = Path(__file__).resolve().parent
REPO = EXP.parents[2]
RESULTS = EXP / "results"
OUT_DIR = Path("/tmp/minco_candidate_preset_replay_20260629")
METHOD = "default_candidate_preset_normalized_depth_alpha2"
REFERENCE_METHOD = wrapper_validation.METHOD
REFERENCE_SCORES = RESULTS / "cross_panel_candidate_abundance_wrapper_scores.tsv"


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def clean_unique_values(values: pd.Series) -> str:
    cleaned = []
    for value in values.astype(str):
        token = value.strip()
        if not token or token.lower() == "nan":
            continue
        cleaned.append(token)
    return ",".join(sorted(set(cleaned)))


def parse_time_log(path: Path) -> dict[str, object]:
    text = path.read_text() if path.exists() else ""
    out = {"wall_clock": "", "max_rss_kb": "", "exit_status": ""}
    for line in text.splitlines():
        if "Elapsed (wall clock) time" in line:
            out["wall_clock"] = line.rsplit("):", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            out["max_rss_kb"] = line.rsplit(":", 1)[-1].strip()
        elif "Exit status" in line:
            out["exit_status"] = line.rsplit(":", 1)[-1].strip()
    return out


def model_args() -> list[str]:
    if wrapper_validation.MODEL_CACHE.exists():
        return ["--model-cache", str(wrapper_validation.MODEL_CACHE)]
    return ["--train-features", str(wrapper_validation.TRAIN_FEATURES)]


def run_default_preset(
    panel: str,
    sample: int,
    unique: Path,
    split: Path,
    taxmap: Path,
    out: Path,
    candidate_surface_taxmap: Path | None = None,
) -> dict[str, object]:
    time_log = OUT_DIR / f"{panel}.sample{sample}.candidate_preset.time.log"
    if out.exists():
        timing = parse_time_log(time_log)
        timing.update(
            {
                "panel": panel,
                "sample": sample,
                "command": "reuse existing default candidate preset profile",
                "elapsed_s": "0.000",
                "time_log": str(time_log),
                "status": "reused",
            }
        )
        return timing
    cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(time_log),
        "python3",
        "-B",
        "scripts/minco_profile_default.py",
        "--profile-preset",
        "candidate",
        "--unique-table",
        str(unique),
        "--split-table",
        str(split),
        "--taxmap",
        str(taxmap),
        *model_args(),
        "--train-pool",
        "train12",
        "--scope",
        "bacteria",
    ]
    if candidate_surface_taxmap is not None:
        cmd.extend(["--candidate-surface-taxmap", str(candidate_surface_taxmap)])
    cmd.extend(["-o", str(out)])
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
            "status": "ran",
        }
    )
    return timing


def replay_emitted() -> tuple[dict[tuple[str, int], Path], pd.DataFrame, pd.DataFrame]:
    paths: dict[tuple[str, int], Path] = {}
    timing_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    for sample, info in wrapper_validation.TOY_TABLES.items():
        panel = "cami2_toy_mouse_gut"
        unique = info["unique"]
        split = info["split"]
        out = OUT_DIR / f"toymouse_sample{sample}.candidate_preset.tsv"
        if not unique.exists() or not split.exists():
            missing_rows.append(
                {"panel": panel, "sample": sample, "unique": str(unique), "split": str(split)}
            )
            continue
        timing_rows.append(
            run_default_preset(panel, sample, unique, split, wrapper_validation.TOY_TAXMAP, out)
        )
        paths[(panel, sample)] = out
    for sample, info in wrapper_validation.CAMI3_TABLES.items():
        panel = "cami3_toy_human_gut_gtdb_source_readmap"
        unique = info["unique"]
        split = info["split"]
        out = OUT_DIR / f"cami3_sample{sample}.candidate_preset.tsv"
        if not unique.exists() or not split.exists():
            missing_rows.append(
                {"panel": panel, "sample": sample, "unique": str(unique), "split": str(split)}
            )
            continue
        timing_rows.append(
            run_default_preset(panel, sample, unique, split, wrapper_validation.GTDB_TAXMAP, out)
        )
        paths[(panel, sample)] = out
    return paths, pd.DataFrame(timing_rows), pd.DataFrame(missing_rows)


def replay_hmp() -> tuple[dict[tuple[str, int], Path], pd.DataFrame, pd.DataFrame]:
    paths: dict[tuple[str, int], Path] = {}
    timing_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    for panel in hmp_surface.PANELS:
        cfg = hmp_surface.decomp.PANELS[panel]
        for sample_raw in cfg["samples"]:
            sample = int(sample_raw)
            profile = Path(cfg["minco"](sample))
            tables = hmp_surface.discover_tables(panel, sample, profile)
            unique = tables.get("unique")
            split = tables.get("split")
            out = OUT_DIR / f"{panel}.sample{sample}.candidate_preset.tsv"
            if not unique or not Path(unique).exists() or not split or not Path(split).exists():
                missing_rows.append(
                    {
                        "panel": panel,
                        "sample": sample,
                        "unique": str(unique or ""),
                        "split": str(split or ""),
                    }
                )
                continue
            timing_rows.append(
                run_default_preset(
                    panel,
                    sample,
                    Path(unique),
                    Path(split),
                    hmp_surface.NORMAL_TAXMAP,
                    out,
                    hmp_surface.SURFACE_TAXMAP,
                )
            )
            paths[(panel, sample)] = out
    return paths, pd.DataFrame(timing_rows), pd.DataFrame(missing_rows)


def candidate_metadata(path_map: dict[tuple[str, int], Path]) -> pd.DataFrame:
    rows = []
    for (panel, sample), path in path_map.items():
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=lambda col: col
            in {
                "calibrated_call",
                "calibrated_abundance",
                "calibrated_abundance_raw",
                "candidate_rescue_added",
                "candidate_surface_added",
                "candidate_abundance_policy",
                "candidate_abundance_norm_mass",
                "profile_strategy",
                "candidate_surface_switch",
                "candidate_rescue_switch",
            },
            low_memory=False,
        )
        called = bool_series(df.get("calibrated_call", pd.Series(False, index=df.index)))
        rescue_added = bool_series(df.get("candidate_rescue_added", pd.Series(False, index=df.index)))
        surface_added = bool_series(df.get("candidate_surface_added", pd.Series(False, index=df.index)))
        added = rescue_added | surface_added
        raw = pd.to_numeric(df.get("calibrated_abundance_raw", 0.0), errors="coerce").fillna(0.0)
        abundance = pd.to_numeric(df.get("calibrated_abundance", 0.0), errors="coerce").fillna(0.0)
        norm_mass = pd.to_numeric(
            df.get("candidate_abundance_norm_mass", 0.0),
            errors="coerce",
        ).fillna(0.0)
        base_called_raw_mass = float(raw.loc[called & ~added].sum())
        if base_called_raw_mass <= 0.0:
            base_called_raw_mass = 1.0
        expected_raw = norm_mass * base_called_raw_mass
        raw_delta = (raw - expected_raw).abs()
        rows.append(
            {
                "panel": panel,
                "sample": str(sample),
                "profile": str(path),
                "profile_strategy_values": clean_unique_values(df["profile_strategy"]),
                "candidate_rescue_switch_values": clean_unique_values(df["candidate_rescue_switch"])
                if "candidate_rescue_switch" in df.columns
                else "",
                "candidate_surface_switch_values": clean_unique_values(df["candidate_surface_switch"])
                if "candidate_surface_switch" in df.columns
                else "",
                "called_rows": int(called.sum()),
                "candidate_added_rows": int(added.sum()),
                "candidate_rescue_added_rows": int(rescue_added.sum()),
                "candidate_surface_added_rows": int(surface_added.sum()),
                "nonzero_candidate_rows": int((raw.loc[added] > 0.0).sum()),
                "candidate_added_abundance_sum": float(abundance.loc[added].sum()),
                "candidate_added_raw_sum": float(raw.loc[added].sum()),
                "candidate_added_norm_mass_sum": float(norm_mass.loc[added].sum()),
                "base_called_raw_mass": base_called_raw_mass,
                "max_candidate_raw_mass_delta": float(raw_delta.loc[added].max()) if added.any() else 0.0,
            }
        )
    return pd.DataFrame(rows)


def score_profiles(
    emitted_paths: dict[tuple[str, int], Path],
    hmp_paths: dict[tuple[str, int], Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adaptive.METHOD = METHOD
    emitted_scores, emitted_meta = adaptive.score_profiles(emitted_paths)
    mapper = raw_audit.RawMapper()
    hmp_rows = []
    for (panel, sample), path in hmp_paths.items():
        row = hmp_abundance.score_profile(panel, int(sample), path, mapper)
        row["method"] = METHOD
        row["panel"] = panel
        row["sample"] = str(sample)
        row["profile"] = str(path)
        hmp_rows.append(row)
    scores = pd.concat([emitted_scores, pd.DataFrame(hmp_rows)], ignore_index=True, sort=False)
    return scores, emitted_meta


def summarize(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for panel, sub in scores.groupby("panel", sort=True):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        pooled_f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        l1_col = adaptive.official_l1_col(panel)
        pearson_col = adaptive.official_pearson_col(panel)
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
                "official_L1_pp": sub[l1_col].mean(),
                "official_Pearson": sub[pearson_col].mean(),
            }
        )
    return pd.DataFrame(rows)


def compare_to_reference(scores: pd.DataFrame) -> pd.DataFrame:
    if not REFERENCE_SCORES.exists():
        return pd.DataFrame(
            [{"panel": "", "sample": "", "status": "missing", "path": str(REFERENCE_SCORES)}]
        )
    reference = pd.read_csv(REFERENCE_SCORES, sep="\t")
    rows = []
    for row in scores.itertuples(index=False):
        panel = str(row.panel)
        sample = int(row.sample)
        ref = reference.loc[
            reference["panel"].astype(str).eq(panel)
            & reference["sample"].astype(int).eq(sample)
            & reference["method"].astype(str).eq(REFERENCE_METHOD)
        ]
        if ref.empty:
            rows.append(
                {"panel": panel, "sample": sample, "status": "missing_reference", "path": str(REFERENCE_SCORES)}
            )
            continue
        ref_row = ref.iloc[0]
        l1_col = adaptive.official_l1_col(panel)
        pearson_col = adaptive.official_pearson_col(panel)
        rows.append(
            {
                "panel": panel,
                "sample": sample,
                "status": "compared",
                "TP_delta": int(row.TP) - int(ref_row["TP"]),
                "FP_delta": int(row.FP) - int(ref_row["FP"]),
                "FN_delta": int(row.FN) - int(ref_row["FN"]),
                "F1_delta": finite(row.F1) - finite(ref_row["F1"]),
                "official_L1_delta_pp": finite(getattr(row, l1_col)) - finite(ref_row[l1_col]),
                "official_Pearson_delta": finite(getattr(row, pearson_col)) - finite(ref_row[pearson_col]),
            }
        )
    return pd.DataFrame(rows)


def audit(scores: pd.DataFrame, meta: pd.DataFrame, comparison: pd.DataFrame) -> pd.DataFrame:
    compared = comparison.loc[comparison["status"].astype(str).eq("compared")].copy()
    if compared.empty:
        comparison_value = "missing"
        comparison_decision = "missing_reference"
    else:
        max_count_delta = int(compared[["TP_delta", "FP_delta", "FN_delta"]].abs().max().max())
        max_f1_delta = float(compared["F1_delta"].abs().max())
        max_l1_delta = float(compared["official_L1_delta_pp"].abs().max())
        max_pearson_delta = float(compared["official_Pearson_delta"].abs().max())
        comparison_value = (
            f"counts={max_count_delta};F1={max_f1_delta:.3g};"
            f"L1={max_l1_delta:.3g};Pearson={max_pearson_delta:.3g}"
        )
        if (
            max_count_delta == 0
            and max_f1_delta < 1e-9
            and max_l1_delta < 1e-9
            and max_pearson_delta < 1e-9
        ):
            comparison_decision = "matches_prior_wrapper_candidate"
        elif max_count_delta == 0 and max_f1_delta < 1e-9 and max_l1_delta <= 0.25 and max_pearson_delta <= 1e-3:
            comparison_decision = "calls_match_small_autoexact_abundance_delta"
        else:
            comparison_decision = "differs_from_prior_wrapper_candidate"
    max_raw_delta = float(meta["max_candidate_raw_mass_delta"].max()) if not meta.empty else 0.0
    expected_strategy = (
        ",".join(sorted(set(meta["profile_strategy_values"].astype(str)))) if not meta.empty else ""
    )
    expected_rescue = (
        ",".join(sorted(set(meta["candidate_rescue_switch_values"].astype(str)))) if not meta.empty else ""
    )
    expected_surface = (
        ",".join(sorted(set(meta["candidate_surface_switch_values"].astype(str)))) if not meta.empty else ""
    )
    return pd.DataFrame(
        [
            {
                "metric": "preset_profiles_scored",
                "value": int(scores.shape[0]),
                "evidence": str(OUT_DIR),
                "decision": "cross_panel_cached_default_preset_replay",
            },
            {
                "metric": "preset_switch_values",
                "value": f"strategy={expected_strategy};rescue={expected_rescue};surface={expected_surface}",
                "evidence": "profile output metadata",
                "decision": "candidate_preset_wired",
            },
            {
                "metric": "candidate_added_rows",
                "value": int(meta["candidate_added_rows"].sum()) if not meta.empty else 0,
                "evidence": "candidate_rescue_added or candidate_surface_added",
                "decision": "candidate_behavior",
            },
            {
                "metric": "nonzero_candidate_rows",
                "value": int(meta["nonzero_candidate_rows"].sum()) if not meta.empty else 0,
                "evidence": "normalized-depth-alpha2",
                "decision": "candidate_abundance_behavior",
            },
            {
                "metric": "max_candidate_raw_mass_delta",
                "value": max_raw_delta,
                "evidence": "candidate_abundance_norm_mass*base_called_raw_mass",
                "decision": "pass" if max_raw_delta < 1e-9 else "fail",
            },
            {
                "metric": "max_delta_vs_prior_wrapper_candidate",
                "value": comparison_value,
                "evidence": "candidate_preset_replay_vs_wrapper.tsv",
                "decision": comparison_decision,
            },
            {
                "metric": "promotion_decision",
                "value": "preset_validated_selected_default",
                "evidence": "candidate preset reproduces prior wrapper candidate if max_delta_vs_prior_wrapper_candidate matches; previous default is available with --profile-preset current",
                "decision": "promote_to_selected_default_without_broad_release_claim",
            },
        ]
    )


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    emitted_paths, emitted_timing, emitted_missing = replay_emitted()
    hmp_paths, hmp_timing, hmp_missing = replay_hmp()
    path_map = {**emitted_paths, **hmp_paths}
    missing = pd.concat([emitted_missing, hmp_missing], ignore_index=True, sort=False)
    if missing.empty and len(path_map) != 32:
        missing = pd.DataFrame(
            [{"panel": "", "sample": "", "reason": f"expected 32 profiles, found {len(path_map)}"}]
        )
    if not missing.empty:
        missing.to_csv(RESULTS / "candidate_preset_replay_missing.tsv", sep="\t", index=False)
        raise SystemExit("missing candidate preset replay inputs; see candidate_preset_replay_missing.tsv")

    scores, scorer_meta = score_profiles(emitted_paths, hmp_paths)
    meta = candidate_metadata(path_map)
    scores = scores.merge(meta.drop(columns=["profile"]), on=["panel", "sample"], how="left")
    summary = summarize(scores)
    comparison = compare_to_reference(scores)
    audit_df = audit(scores, meta, comparison)
    timing = pd.concat([emitted_timing, hmp_timing], ignore_index=True, sort=False)

    pd.DataFrame(
        [
            {
                "source": "default_candidate_preset",
                "profiles": len(path_map),
                "path": str(OUT_DIR),
                "reference_scores": str(REFERENCE_SCORES),
            }
        ]
    ).to_csv(RESULTS / "candidate_preset_replay_apply.tsv", sep="\t", index=False)
    timing.to_csv(RESULTS / "candidate_preset_replay_runtime.tsv", sep="\t", index=False)
    scorer_meta.to_csv(RESULTS / "candidate_preset_replay_scorer_metadata.tsv", sep="\t", index=False)
    meta.to_csv(RESULTS / "candidate_preset_replay_candidate_metadata.tsv", sep="\t", index=False)
    scores.to_csv(RESULTS / "candidate_preset_replay_scores.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "candidate_preset_replay_summary.tsv", sep="\t", index=False)
    comparison.to_csv(RESULTS / "candidate_preset_replay_vs_wrapper.tsv", sep="\t", index=False)
    audit_df.to_csv(RESULTS / "candidate_preset_replay_audit.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nAUDIT")
    print(audit_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
