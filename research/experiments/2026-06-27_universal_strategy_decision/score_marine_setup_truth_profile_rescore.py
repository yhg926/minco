#!/usr/bin/env python3
"""Score marine profiles against the setup+assembly GTDB truth rule.

This scorer consumes the truth-rule outputs from
``audit_marine_setup_metadata_truth_upgrade.py``. It is intentionally separate
from ``score_marine_gtdb_taxid_transfer.py`` because the accepted route now uses
the stricter setup source-name or local assembly-summary rule, not only the
older taxid/binomial transfer.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd


EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
RESULTS = EXP / "results"

TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402
import score_marine_gtdb_taxid_transfer as marine_taxid  # noqa: E402


TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
BASE_DETAIL = RESULTS / "marine_binomial_transfer_detail.tsv"
RULES = RESULTS / "marine_setup_metadata_truth_upgrade_rules.tsv"
RESCUES = RESULTS / "marine_setup_metadata_truth_upgrade_rescues.tsv"

DEFAULT_MINCO = {
    0: Path(
        "/tmp/minco_fair_l1_edge_adaptive_20260624/marine_current_best/"
        "minco_s1000_unique_zipaaf_sample0.tsv"
    ),
    3: Path(
        "/tmp/cami2_marine_samples3_5_20260625/run/"
        "minco_sample3_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
    ),
    4: Path(
        "/tmp/cami2_marine_samples3_5_20260625/run/"
        "minco_sample4_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
    ),
    5: Path(
        "/tmp/cami2_marine_samples3_5_20260625/run/"
        "minco_sample5_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
    ),
}

DEFAULT_SYLPH = {
    0: Path("/tmp/minco_fair_l1_edge_adaptive_20260624/marine_sylph/profile.tsv"),
    3: Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample3/profile.tsv"),
    4: Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample4/profile.tsv"),
    5: Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample5/profile.tsv"),
}

OUT_PREFIX = "marine_setup_truth_profile_rescore"


def parse_sample_path_overrides(values: Iterable[str]) -> dict[int, Path]:
    out: dict[int, Path] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"Expected SAMPLE=PATH override, got: {value}")
        raw_sample, raw_path = value.split("=", 1)
        out[int(raw_sample.strip())] = Path(raw_path.strip())
    return out


def parse_samples(value: str) -> list[int]:
    samples = [int(part.strip()) for part in value.split(",") if part.strip()]
    if not samples:
        raise SystemExit("--samples selected no sample IDs")
    return samples


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def gtdb_from_accession(
    accession: object,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[str, str]:
    rec, method, _key = truth.lookup_accession(str(accession or ""), by_accession, by_core)
    if rec is None:
        return "", method
    return str(rec.get("gtdb_species", "")), method


def build_truth_for_rule(rule_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = pd.read_csv(BASE_DETAIL, sep="\t")
    mapped = detail.loc[detail["gtdb_species"].fillna("").astype(str).ne("")].copy()
    rows: list[dict[str, object]] = []
    for row in mapped.itertuples(index=False):
        rows.append(
            {
                "sample": int(row.sample),
                "ncbi_species_taxid": str(row.ncbi_species_taxid),
                "truth_name": row.truth_name,
                "gtdb_species": str(row.gtdb_species),
                "abundance_pct_all": finite(row.abundance_pct_all),
                "truth_transfer_rule": str(row.transfer_method),
            }
        )

    rescues = pd.read_csv(RESCUES, sep="\t")
    if not rescues.empty:
        sub = rescues.loc[rescues["rule_id"].astype(str).eq(rule_id)].copy()
        for row in sub.itertuples(index=False):
            species = str(row.candidate_gtdb_species)
            if not species:
                continue
            rows.append(
                {
                    "sample": int(row.sample),
                    "ncbi_species_taxid": str(row.ncbi_species_taxid),
                    "truth_name": row.truth_name,
                    "gtdb_species": species,
                    "abundance_pct_all": finite(row.abundance_pct_all),
                    "truth_transfer_rule": rule_id,
                }
            )

    truth_rows = pd.DataFrame(rows)
    if truth_rows.empty:
        return (
            pd.DataFrame(
                columns=[
                    "sample",
                    "gtdb_species",
                    "abundance_pct_all",
                    "truth_abundance",
                    "source_taxids",
                    "truth_names",
                    "truth_transfer_rules",
                ]
            ),
            pd.DataFrame(),
        )

    grouped = (
        truth_rows.groupby(["sample", "gtdb_species"], as_index=False)
        .agg(
            abundance_pct_all=("abundance_pct_all", "sum"),
            source_taxids=("ncbi_species_taxid", lambda x: ",".join(sorted(set(map(str, x))))),
            truth_names=("truth_name", lambda x: ";".join(sorted(set(map(str, x))))),
            truth_transfer_rules=(
                "truth_transfer_rule",
                lambda x: ",".join(sorted(set(map(str, x)))),
            ),
        )
        .sort_values(["sample", "gtdb_species"])
    )
    totals = grouped.groupby("sample")["abundance_pct_all"].transform("sum")
    grouped["truth_abundance"] = grouped["abundance_pct_all"] / totals.where(totals.ne(0), 1.0)

    quality = pd.read_csv(RULES, sep="\t")
    quality = quality.loc[quality["rule_id"].astype(str).eq(rule_id)].copy()
    return grouped, quality


def load_current_minco_predictions(
    path: Path,
    taxid_to_gtdb: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = as_bool(raw["calibrated_call"])
    selected = raw.loc[call].copy()
    rows: list[tuple[str, float, float]] = []
    ref_mapped = 0
    taxid_mapped = 0
    unmapped = 0
    for row in selected.itertuples(index=False):
        species = ""
        for col in ["candidate_surface_accession", "s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species, _method = gtdb_from_accession(getattr(row, col), by_accession, by_core)
                if species:
                    ref_mapped += 1
                    break
        if not species:
            species = taxid_to_gtdb.get(str(getattr(row, "taxid", "")), "")
            if species:
                taxid_mapped += 1
        if not species:
            unmapped += 1
            continue
        abundance = finite(getattr(row, "calibrated_abundance", 0.0))
        ani = finite(getattr(row, "reported_ani", 0.0)) or finite(
            getattr(row, "s_Ref_zip_aaf_ani_max", 0.0)
        )
        rows.append((species, abundance, ani))
    return taxid_score.collapse_prediction(rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_ref_mapped": ref_mapped,
        "pred_rows_taxid_fallback_mapped": taxid_mapped,
        "profile_format": "calibrated_selected_default",
    }


def load_minco_any_profile(
    path: Path,
    taxid_to_gtdb: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    columns = pd.read_csv(path, sep="\t", nrows=0).columns
    if "calibrated_call" in columns:
        return load_current_minco_predictions(path, taxid_to_gtdb, by_accession, by_core)
    taxmap = marine_taxid.parse_species_taxmap(TAXMAP)
    pred, extra = marine_taxid.collapse_minco_prediction(
        path,
        taxmap,
        dict(taxid_to_gtdb),
        by_accession,
        by_core,
    )
    extra["profile_format"] = "legacy_readwise_profile"
    return pred, extra


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def build_audit(
    samples: list[int],
    rule_id: str,
    profile_set_label: str,
    score_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    quality_df: pd.DataFrame,
) -> list[dict[str, object]]:
    minco = summary_df.loc[summary_df["method"].astype(str).str.startswith("minco")]
    sylph = summary_df.loc[summary_df["method"].astype(str).str.startswith("sylph")]
    if not minco.empty and not sylph.empty:
        minco_row = minco.iloc[0]
        sylph_row = sylph.iloc[0]
        comparison = (
            f"minco_mean_F1={finite(minco_row['mean_F1']):.6f};"
            f"sylph_mean_F1={finite(sylph_row['mean_F1']):.6f};"
            f"minco_mean_L1_union_pp={finite(minco_row['mean_L1_union_pp']):.6f};"
            f"sylph_mean_L1_union_pp={finite(sylph_row['mean_L1_union_pp']):.6f};"
            f"minco_mean_Pearson_union={finite(minco_row['mean_Pearson_union']):.6f};"
            f"sylph_mean_Pearson_union={finite(sylph_row['mean_Pearson_union']):.6f}"
        )
        comparison_decision = (
            "minco_higher_F1"
            if finite(minco_row["mean_F1"]) >= finite(sylph_row["mean_F1"])
            else "sylph_higher_F1"
        )
    else:
        comparison = "minco_or_sylph_summary_missing"
        comparison_decision = "comparison_incomplete"

    profile_set_is_selected = profile_set_label == "selected_default_same_namespace"
    all_ready = False
    if not quality_df.empty:
        all_ready = bool(quality_df["ready_at_95pct_bacteria_archaea"].astype(str).str.lower().eq("true").all())
    promotion_decision = (
        "selected_default_profile_scored_review_metrics"
        if profile_set_is_selected and all_ready
        else "diagnostic_cached_profile_score_not_selected_default"
    )
    return [
        {
            "metric": "scoring_scope",
            "value": (
                f"samples={','.join(map(str, samples))};rule_id={rule_id};"
                f"profile_set_label={profile_set_label}"
            ),
            "evidence": f"results/{OUT_PREFIX}_scores.tsv",
            "decision": "same_namespace_truth_rule_applied",
        },
        {
            "metric": "truth_rule_quality",
            "value": (
                f"ready_samples={int(quality_df['ready_at_95pct_bacteria_archaea'].astype(str).str.lower().eq('true').sum())}/{len(quality_df)};"
                f"min_mapped_pct_bacteria_archaea={finite(quality_df['mapped_pct_bacteria_archaea'].min() if not quality_df.empty else 0.0):.6f}"
            ),
            "evidence": str(RULES.relative_to(EXP)),
            "decision": "truth_rule_passes_threshold" if all_ready else "truth_rule_not_all_ready",
        },
        {
            "metric": "minco_vs_sylph_summary",
            "value": comparison,
            "evidence": f"results/{OUT_PREFIX}_summary.tsv",
            "decision": comparison_decision,
        },
        {
            "metric": "promotion_decision",
            "value": promotion_decision,
            "evidence": f"results/{OUT_PREFIX}_audit.tsv",
            "decision": (
                "route_profile_evidence_ready_for_review"
                if profile_set_is_selected and all_ready
                else "do_not_close_route_from_cached_nondefault_profiles"
            ),
        },
    ]


def write_md(audit: list[dict[str, object]]) -> None:
    lines = [
        "# Marine Setup Truth Profile Rescore",
        "",
        "Date: 2026-06-30",
        "",
        "This scorer applies the strict setup source-name or local assembly-summary",
        "GTDB truth rule to marine profile outputs. It does not run profiling.",
        "",
        "## Audit",
        "",
        "| Metric | Value | Decision |",
        "|---|---:|---|",
    ]
    for row in audit:
        lines.append(f"| `{row['metric']}` | {row['value']} | {row['decision']} |")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- `results/{OUT_PREFIX}_truth.tsv`",
            f"- `results/{OUT_PREFIX}_quality.tsv`",
            f"- `results/{OUT_PREFIX}_scores.tsv`",
            f"- `results/{OUT_PREFIX}_summary.tsv`",
            f"- `results/{OUT_PREFIX}_audit.tsv`",
            "",
        ]
    )
    (EXP / "MARINE_SETUP_TRUTH_PROFILE_RESCORE.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Score marine MinCO/Sylph profiles against setup+assembly GTDB truth.",
    )
    ap.add_argument("--samples", default="0,3,4,5")
    ap.add_argument(
        "--rule-id",
        default="strict_source_name_or_assembly_sources_unique",
        help="Truth-upgrade rule from marine_setup_metadata_truth_upgrade_rules.tsv.",
    )
    ap.add_argument(
        "--profile-set-label",
        default="cached_existing",
        help="Use selected_default_same_namespace only for same-release selected-default profile sets.",
    )
    ap.add_argument(
        "--minco-method",
        default="minco_cached_marine_setup_truth",
        help="Method label for MinCO summary rows.",
    )
    ap.add_argument(
        "--sylph-method",
        default="sylph_cached_marine_setup_truth",
        help="Method label for Sylph summary rows.",
    )
    ap.add_argument("--minco-sample", action="append", default=[], metavar="SAMPLE=PATH")
    ap.add_argument("--sylph-sample", action="append", default=[], metavar="SAMPLE=PATH")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    samples = parse_samples(args.samples)
    minco_paths = {**DEFAULT_MINCO, **parse_sample_path_overrides(args.minco_sample)}
    sylph_paths = {**DEFAULT_SYLPH, **parse_sample_path_overrides(args.sylph_sample)}

    required = [BASE_DETAIL, RULES, RESCUES]
    missing: list[Path] = []
    for sample in samples:
        if sample not in minco_paths:
            missing.append(Path(f"<missing-minco-sample-{sample}>"))
        else:
            required.append(minco_paths[sample])
        if sample not in sylph_paths:
            missing.append(Path(f"<missing-sylph-sample-{sample}>"))
        else:
            required.append(sylph_paths[sample])
    missing.extend(path for path in required if not path.exists())
    if missing:
        raise SystemExit(
            "missing required marine setup truth profile-rescore inputs:\n"
            + "\n".join(map(str, missing))
        )

    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    taxid_to_gtdb, _ambiguous_taxids = taxid_score.build_taxid_transfer(by_accession)
    truth_df, quality_df_all = build_truth_for_rule(args.rule_id)
    quality_df = quality_df_all.loc[quality_df_all["sample"].isin(samples)].copy()

    truth_out = truth_df.loc[truth_df["sample"].isin(samples)].copy()
    score_rows: list[dict[str, object]] = []
    for sample in samples:
        sample_truth = truth_out.loc[truth_out["sample"].eq(sample)].copy()
        minco_pred, minco_extra = load_minco_any_profile(
            minco_paths[sample],
            taxid_to_gtdb,
            by_accession,
            by_core,
        )
        minco_extra.update({"profile": str(minco_paths[sample])})
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(
            sylph_paths[sample],
            by_accession,
            by_core,
        )
        sylph_extra.update({"profile": str(sylph_paths[sample])})
        score_rows.append(
            taxid_score.score_prediction(
                sample,
                args.minco_method,
                minco_pred,
                sample_truth,
                minco_extra,
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample,
                args.sylph_method,
                sylph_pred,
                sample_truth,
                sylph_extra,
            )
        )

    score_df = pd.DataFrame(score_rows)
    summary_df = taxid_score.summarize(score_df)
    audit = build_audit(
        samples,
        args.rule_id,
        args.profile_set_label,
        score_df,
        summary_df,
        quality_df,
    )

    RESULTS.mkdir(parents=True, exist_ok=True)
    truth_out.to_csv(RESULTS / f"{OUT_PREFIX}_truth.tsv", sep="\t", index=False)
    quality_df.to_csv(RESULTS / f"{OUT_PREFIX}_quality.tsv", sep="\t", index=False)
    score_df.to_csv(RESULTS / f"{OUT_PREFIX}_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / f"{OUT_PREFIX}_summary.tsv", sep="\t", index=False)
    write_tsv(RESULTS / f"{OUT_PREFIX}_audit.tsv", audit, ["metric", "value", "evidence", "decision"])
    write_md(audit)
    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
