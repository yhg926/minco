#!/usr/bin/env python3
"""Summarize MinCO readwise trace rows by CAMISIM source genome."""

from __future__ import annotations

from pathlib import Path
import math

import pandas as pd


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-23_low_abundance_ctx_contamination"
ANI_EXP = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani"
TRUTH_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"

SOURCE_TRUTH = TRUTH_EXP / "mouse0_gtdb_source_species_ground_truth.tsv"
READWISE_BIAS = ANI_EXP / "readwise_bias_depth_abundance_check.tsv"
SOURCE_REP = ANI_EXP / "toymouse_sample0_source_rep_ani.tsv"
READ_MAPPING = Path("/tmp/sample0_reads_mapping.tsv.gz")
TRACE_DIR = EXP / "traces_perread"

TARGETS = [
    {
        "target_label": "kinnaridis_low_bad",
        "gtdb_species": "s__Limosilactobacillus kinnaridis",
        "representative_accession": "GCA_946487795.1",
    },
    {
        "target_label": "johnsonii_low_bad",
        "gtdb_species": "s__Lactobacillus johnsonii",
        "representative_accession": "GCF_000159355.1",
    },
    {
        "target_label": "crispatus_low_bad",
        "gtdb_species": "s__Lactobacillus crispatus",
        "representative_accession": "GCF_018987235.1",
    },
    {
        "target_label": "rhamnosus_low_multisource",
        "gtdb_species": "s__Lacticaseibacillus rhamnosus",
        "representative_accession": "GCF_900636965.1",
    },
    {
        "target_label": "paracasei_high_control",
        "gtdb_species": "s__Lacticaseibacillus paracasei",
        "representative_accession": "GCF_000829035.1",
    },
]


def safe_num(value) -> float:
    try:
        if value == "":
            return math.nan
        return float(value)
    except Exception:
        return math.nan


def frac(num: float, den: float) -> float:
    return float(num) / float(den) if den else math.nan


def trace_path_for(target: dict[str, str]) -> Path:
    return TRACE_DIR / f"{target['target_label']}.trace.tsv"


def collect_needed_read_ids() -> set[str]:
    needed: set[str] = set()
    for target in TARGETS:
        trace_path = trace_path_for(target)
        if not trace_path.exists():
            raise SystemExit(f"missing trace {trace_path}")
        for chunk in pd.read_csv(trace_path, sep="\t", usecols=["read_id"], chunksize=200_000):
            needed.update(chunk["read_id"].dropna().astype(str))
    return needed


def load_read_mapping(needed_read_ids: set[str]) -> pd.DataFrame:
    if not READ_MAPPING.exists():
        raise SystemExit(f"missing {READ_MAPPING}; extract it with commands.sh first")
    if not needed_read_ids:
        return pd.DataFrame(columns=["read_id", "genome_id", "tax_id"])
    kept = []
    for chunk in pd.read_csv(
        READ_MAPPING,
        sep="\t",
        dtype=str,
        chunksize=1_000_000,
        usecols=["#anonymous_read_id", "genome_id", "tax_id"],
    ):
        chunk = chunk.rename(columns={"#anonymous_read_id": "read_id"})
        hit = chunk[chunk["read_id"].isin(needed_read_ids)]
        if not hit.empty:
            kept.append(hit)
    if not kept:
        return pd.DataFrame(columns=["read_id", "genome_id", "tax_id"])
    return (
        pd.concat(kept, ignore_index=True)
        .drop_duplicates("read_id", keep="first")[["read_id", "genome_id", "tax_id"]]
    )


def classify_trace(
    trace: pd.DataFrame,
    mapping: pd.DataFrame,
    truth: pd.DataFrame,
    target_species: str,
    target_source_ids: set[str],
) -> pd.DataFrame:
    x = trace.merge(mapping, on="read_id", how="left")
    source_cols = [
        "genome_id",
        "source_accession",
        "source_taxid",
        "gtdb_species",
        "abundance",
    ]
    x = x.merge(truth[source_cols], on="genome_id", how="left")
    x["is_target_source"] = x["genome_id"].isin(target_source_ids)
    x["is_same_gtdb_species"] = x["gtdb_species"].eq(target_species)
    x["is_other_source"] = ~x["is_same_gtdb_species"].fillna(False)
    x["category"] = "other_or_unmapped"
    x.loc[x["is_same_gtdb_species"], "category"] = "same_gtdb_species"
    x.loc[x["is_target_source"], "category"] = "target_source"
    x["diff"] = pd.to_numeric(x["diff"], errors="coerce")
    x["best_diff"] = pd.to_numeric(x["best_diff"], errors="coerce")
    x["candidate_refs"] = pd.to_numeric(x["candidate_refs"], errors="coerce")
    x["selected_refs"] = pd.to_numeric(x["selected_refs"], errors="coerce")
    x["cov_inc"] = pd.to_numeric(x["cov_inc"], errors="coerce")
    return x


def add_category_stats(prefix: str, rows: pd.DataFrame, out: dict[str, object]) -> None:
    out[f"{prefix}_event_n"] = len(rows)
    out[f"{prefix}_unique_read_n"] = rows["read_id"].nunique() if len(rows) else 0
    out[f"{prefix}_unique_qctx_n"] = rows["qctx"].nunique() if len(rows) else 0
    out[f"{prefix}_event_frac"] = frac(len(rows), out["trace_event_n"])
    out[f"{prefix}_cov_weight_frac"] = frac(rows["cov_inc"].sum(), out["cov_inc_sum"])
    out[f"{prefix}_diff_mean"] = rows["diff"].mean()
    out[f"{prefix}_best_diff_mean"] = rows["best_diff"].mean()
    out[f"{prefix}_diff_gt0_frac"] = frac((rows["diff"] > 0).sum(), len(rows))
    out[f"{prefix}_best_diff_gt0_frac"] = frac((rows["best_diff"] > 0).sum(), len(rows))


def main() -> None:
    truth = pd.read_csv(SOURCE_TRUTH, sep="\t", dtype=str)
    truth["abundance_num"] = truth["abundance"].map(safe_num)
    bias = pd.read_csv(READWISE_BIAS, sep="\t")
    source_rep = pd.read_csv(SOURCE_REP, sep="\t")
    needed_read_ids = collect_needed_read_ids()
    mapping = load_read_mapping(needed_read_ids)

    summary_rows = []
    top_rows = []
    selected_rows = []

    for target in TARGETS:
        label = target["target_label"]
        species = target["gtdb_species"]
        rep_acc = target["representative_accession"]
        species_truth = truth[truth["gtdb_species"].eq(species)].copy()
        target_source_ids = set(species_truth["genome_id"].astype(str))

        rep_rows = source_rep[
            source_rep["gtdb_representative_accession"].astype(str).eq(rep_acc)
            & source_rep["gtdb_species"].astype(str).eq(species)
        ].copy()
        bias_rows = bias[
            bias["gtdb_representative_accession"].astype(str).eq(rep_acc)
            & bias["gtdb_species"].astype(str).eq(species)
        ].copy()
        trace_path = trace_path_for(target)
        if not trace_path.exists():
            raise SystemExit(f"missing trace {trace_path}")
        trace = pd.read_csv(trace_path, sep="\t", dtype={"read_id": str})
        trace["qctx"] = trace["qctx"].astype(str)
        classified = classify_trace(trace, mapping, truth, species, target_source_ids)

        primary = bias_rows.sort_values("source_abundance", ascending=False).iloc[0]
        rep_primary = rep_rows.sort_values("source_abundance", ascending=False).iloc[0]

        selected_rows.append(
            {
                "target_label": label,
                "gtdb_species": species,
                "representative_accession": rep_acc,
                "trace_ref_match": rep_acc,
                "source_genome_ids_same_gtdb_species": ",".join(
                    species_truth["genome_id"].astype(str)
                ),
                "source_accessions_same_gtdb_species": ",".join(
                    species_truth["source_accession"].astype(str)
                ),
                "raw_source_abundance_sum": species_truth["abundance_num"].sum(),
                "representative_path": rep_primary.get("representative_path", ""),
            }
        )

        out = {
            "target_label": label,
            "gtdb_species": species,
            "representative_accession": rep_acc,
            "same_gtdb_source_genome_n": len(species_truth),
            "same_gtdb_raw_source_abundance_sum": species_truth["abundance_num"].sum(),
            "ANIm_ANI_primary_source": safe_num(primary.get("ANIm_ANI")),
            "readwise_minco_naive_coden11": safe_num(
                primary.get("readwise_minco_naive_coden11")
            ),
            "c11_readwise_err_primary_source": safe_num(primary.get("c11_readwise_err")),
            "c11_XnY_ctx": safe_num(primary.get("c11_XnY_ctx")),
            "c11_Ref_breadth": safe_num(primary.get("c11_Ref_breadth")),
            "c11_Ref_mean_depth": safe_num(primary.get("c11_Ref_mean_depth")),
            "c11_Ref_hit_mean_depth": safe_num(primary.get("c11_Ref_hit_mean_depth")),
            "c11_Reliable_Ref_hit_median_depth": safe_num(
                primary.get("c11_Reliable_Ref_hit_median_depth")
            ),
            "trace_event_n": len(classified),
            "trace_unique_read_n": classified["read_id"].nunique(),
            "trace_unique_qctx_n": classified["qctx"].nunique(),
            "cov_inc_sum": classified["cov_inc"].sum(),
            "candidate_refs_mean": classified["candidate_refs"].mean(),
            "selected_refs_mean": classified["selected_refs"].mean(),
            "diff_mean": classified["diff"].mean(),
            "best_diff_mean": classified["best_diff"].mean(),
            "diff_gt0_frac": frac((classified["diff"] > 0).sum(), len(classified)),
            "best_diff_gt0_frac": frac(
                (classified["best_diff"] > 0).sum(), len(classified)
            ),
        }
        add_category_stats(
            "target_source", classified[classified["is_target_source"]], out
        )
        add_category_stats(
            "same_gtdb_species",
            classified[classified["is_same_gtdb_species"]],
            out,
        )
        add_category_stats(
            "other_or_unmapped_source",
            classified[classified["is_other_source"]],
            out,
        )
        summary_rows.append(out)

        group_cols = [
            "genome_id",
            "source_accession",
            "source_taxid",
            "gtdb_species",
            "category",
        ]
        grouped = (
            classified.groupby(group_cols, dropna=False)
            .agg(
                event_n=("read_id", "size"),
                unique_read_n=("read_id", "nunique"),
                unique_qctx_n=("qctx", "nunique"),
                cov_inc_sum=("cov_inc", "sum"),
                diff_mean=("diff", "mean"),
                best_diff_mean=("best_diff", "mean"),
                diff_gt0_n=("diff", lambda s: int((s > 0).sum())),
                best_diff_gt0_n=("best_diff", lambda s: int((s > 0).sum())),
            )
            .reset_index()
        )
        grouped["target_label"] = label
        grouped["target_species"] = species
        grouped["representative_accession"] = rep_acc
        grouped["event_frac"] = grouped["event_n"] / len(classified) if len(classified) else math.nan
        grouped["cov_weight_frac"] = (
            grouped["cov_inc_sum"] / classified["cov_inc"].sum()
            if classified["cov_inc"].sum()
            else math.nan
        )
        grouped["diff_gt0_frac"] = grouped["diff_gt0_n"] / grouped["event_n"]
        grouped["best_diff_gt0_frac"] = grouped["best_diff_gt0_n"] / grouped["event_n"]
        grouped = grouped.sort_values(
            ["target_label", "event_n", "cov_inc_sum"], ascending=[True, False, False]
        )
        top_rows.extend(grouped.head(12).to_dict("records"))

    pd.DataFrame(selected_rows).to_csv(EXP / "selected_trace_refs.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(EXP / "ctx_trace_source_summary.tsv", sep="\t", index=False)
    pd.DataFrame(top_rows).to_csv(EXP / "ctx_trace_top_sources.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
