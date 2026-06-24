#!/usr/bin/env python3
"""Build GTDB truth and score CAMI II Toy Mouse samples 0-2."""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
OLD_EXP = ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(OLD_EXP))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402


BASE = Path("/mnt/new3T/minco_cami2_toymouse_20260621")
RUN_DIR = Path("/tmp/cami2_toymouse_more_gtdb_20260622")
CTX_MARKER = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker")
CTXOBJ_MARKER = Path(
    "/tmp/gtdb232_s2000_ctxobj_marker_20260622/"
    "sketch_T_S2000_aaf003_dedup_ctxobjmarker"
)

SAMPLES = [0, 1, 2]
EFFECTIVE_CTX_LENGTH = 24
ANI_THRESHOLD = 0.95
ACTIVE_CTX_MIN = 15.0
ACTIVE_RELIABLE_ZTP_AF_FLOOR = 0.40
ACTIVE_MEAN_DEPTH_MIN = 3.0
ACTIVE_VMR_MIN = 50.0
ACTIVE_DELTA_MAX = 0.03


def sample_read_path(sample: int) -> Path:
    return BASE / f"sample_{sample}/2017.12.29_11.37.26_sample_{sample}/reads/anonymous_reads.fq.gz"


def distribution_path(sample: int) -> Path:
    sample_path = BASE / f"sample_{sample}/distributions/distribution_{sample}.txt"
    if sample_path.exists():
        return sample_path
    return BASE / f"setup/distributions/distribution_{sample}.txt"


def minco_output_path(sample: int, marker_kind: str) -> Path:
    if sample == 0 and marker_kind == "ctx_only_current_best":
        return OLD_EXP / "toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv"
    if sample == 0 and marker_kind == "ctxobj_markerdb_product0_active":
        return (
            Path("/tmp/gtdb232_s2000_ctxobj_marker_20260622")
            / "toymouse_sample0_ctxobjmarker_split_naive_product_topfrac_median025.tsv"
        )
    suffix = "ctxmarker" if marker_kind == "ctx_only_current_best" else "ctxobjmarker"
    return RUN_DIR / f"toymouse_sample{sample}_{suffix}_split_naive_product_topfrac_median025.tsv"


def sylph_profile_path(sample: int) -> Path:
    return BASE / f"sylph_sample{sample}/profile.tsv"


def read_tsv_dict(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_meta() -> tuple[dict[str, str], dict[str, str]]:
    meta_rows = read_tsv_dict(BASE / "setup/internal/meta_data.tsv")
    genome_taxid = {row["genome_ID"]: row["NCBI_ID"] for row in meta_rows}
    genome_location = {}
    with (BASE / "setup/internal/genome_locations.tsv").open() as fh:
        for line in fh:
            if not line.strip():
                continue
            genome_id, location = line.rstrip("\n").split("\t")[:2]
            genome_location[genome_id] = location
    return genome_taxid, genome_location


def load_distribution(sample: int) -> list[tuple[str, float]]:
    rows = []
    with distribution_path(sample).open() as fh:
        for line in fh:
            if not line.strip():
                continue
            genome_id, abundance_s = line.rstrip("\n").split("\t")[:2]
            abundance = float(abundance_s)
            if abundance > 0.0:
                rows.append((genome_id, abundance))
    return rows


def build_gtdb_truth(
    sample: int,
    by_accession: Mapping[str, dict[str, str]],
    by_core: Mapping[str, list[dict[str, str]]],
    genome_taxid: Mapping[str, str],
    genome_location: Mapping[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for genome_id, abundance in load_distribution(sample):
        source_location = genome_location.get(genome_id, "")
        source_accession = truth.extract_accession(source_location)
        rec, method, key = truth.lookup_accession(source_accession, by_accession, by_core)
        if rec is None:
            gtdb_fields = {
                "metadata_accession": "",
                "metadata_domain": "",
                "metadata_genbank_accession": "",
                "ncbi_taxid_metadata": "",
                "ncbi_species_taxid_metadata": "",
                "ncbi_species_metadata": "",
                "ncbi_organism_name_metadata": "",
                "gtdb_species": "",
                "gtdb_taxonomy": "",
            }
        else:
            gtdb_fields = {
                "metadata_accession": rec["metadata_accession"],
                "metadata_domain": rec["metadata_domain"],
                "metadata_genbank_accession": rec["genbank_accession"],
                "ncbi_taxid_metadata": rec["ncbi_taxid"],
                "ncbi_species_taxid_metadata": rec["ncbi_species_taxid"],
                "ncbi_species_metadata": rec["ncbi_species"],
                "ncbi_organism_name_metadata": rec["ncbi_organism_name"],
                "gtdb_species": rec["gtdb_species"],
                "gtdb_taxonomy": rec["gtdb_taxonomy"],
            }
        rows.append(
            {
                "genome_id": genome_id,
                "abundance": abundance,
                "source_taxid": genome_taxid.get(genome_id, ""),
                "source_accession": source_accession,
                "source_basename": Path(source_location).name,
                "source_lookup_key": key,
                "mapping_method": method,
                **gtdb_fields,
                "tar_path": f"source_genomes/{Path(source_location).name}" if source_location else "",
                "local_path": "",
            }
        )
    return rows


def write_truth_files(
    sample: int,
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    source_path = EXP_DIR / f"mouse{sample}_gtdb_source_species_ground_truth.tsv"
    profile_path = EXP_DIR / f"mouse{sample}_gtdb_species_profile.tsv"
    summary_path = EXP_DIR / f"mouse{sample}_gtdb_truth_summary.tsv"
    truth.write_source_truth(rows, source_path)
    profile_rows = truth.write_species_profile(rows, profile_path)
    method_counts = Counter(str(row.get("mapping_method", "")) for row in rows)
    with summary_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(
            [
                {"metric": "sample", "value": sample},
                {"metric": "positive_genomes", "value": len(rows)},
                {"metric": "raw_abundance_sum", "value": sum(float(r["abundance"]) for r in rows)},
                {
                    "metric": "mapped_positive_genomes",
                    "value": sum(1 for r in rows if str(r.get("gtdb_species", ""))),
                },
                {
                    "metric": "unmapped_positive_genomes",
                    "value": sum(1 for r in rows if not str(r.get("gtdb_species", ""))),
                },
                {"metric": "gtdb_species_count", "value": len(profile_rows)},
            ]
        )
        for method, count in sorted(method_counts.items()):
            writer.writerow({"metric": f"mapping_method:{method}", "value": count})
    return profile_rows


def numeric(rows: pd.DataFrame, col: str) -> pd.Series:
    if col not in rows.columns:
        return pd.Series([0.0] * len(rows), index=rows.index)
    return pd.to_numeric(rows[col], errors="coerce").fillna(0.0)


def ztp_lambda_from_pos_mean(mean_pos: float) -> float:
    if not math.isfinite(mean_pos) or mean_pos <= 0.0:
        return 0.0
    if mean_pos <= 1.0 + 1e-10:
        return 1e-10
    lo = 1e-10
    hi = max(2.0, mean_pos * 2.0)

    def cond_mean(lam: float) -> float:
        return lam / (1.0 - math.exp(-lam))

    while cond_mean(hi) < mean_pos and hi < 1e6:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if cond_mean(mid) < mean_pos:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ztp_adjusted_af(breadth: float, mean_pos: float) -> float:
    if not math.isfinite(breadth) or breadth <= 0.0:
        return 0.0
    lam = ztp_lambda_from_pos_mean(mean_pos)
    p_nonzero = 1.0 - math.exp(-lam) if lam > 0.0 else 0.0
    if p_nonzero <= 0.0:
        return 1.0
    return min(1.0, breadth / p_nonzero)


def score_sets(predicted: Iterable[str], gold: Iterable[str]):
    pred = {str(x) for x in predicted if str(x)}
    truth_set = {str(x) for x in gold if str(x)}
    tp = pred & truth_set
    fp = pred - truth_set
    fn = truth_set - pred
    precision = len(tp) / len(pred) if pred else 0.0
    recall = len(tp) / len(truth_set) if truth_set else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return tp, fp, fn, precision, recall, f1


def load_minco_active_rows(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = truth.load_minco_with_gtdb_species(path, by_accession, by_core)
    rows = truth.add_naive_ani(rows)
    for col in [
        "XnY_ctx",
        "ANI_naive_calc",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_depth_variance",
        "Normalized_abundance_depth",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
    ]:
        rows[col] = numeric(rows, col)
    rows["Reliable_ztp_af"] = [
        ztp_adjusted_af(float(b), float(m))
        for b, m in zip(rows["Reliable_Ref_breadth"], rows["Reliable_Ref_hit_mean_depth"])
    ]
    rows["ANI_from_Reliable_ztp_af"] = [
        1.0 + math.log(max(float(af), 1e-300)) / EFFECTIVE_CTX_LENGTH
        for af in rows["Reliable_ztp_af"]
    ]
    rows["ANI_AF_delta"] = rows["ANI_naive_calc"] - rows["ANI_from_Reliable_ztp_af"]
    rows["Reliable_depth_vmr"] = [
        (float(var) / float(mean)) if float(mean) > 0.0 else 0.0
        for mean, var in zip(rows["Reliable_Ref_hit_mean_depth"], rows["Reliable_Ref_hit_depth_variance"])
    ]
    rows["active_delta_trigger"] = (
        (rows["Reliable_Ref_hit_mean_depth"] > ACTIVE_MEAN_DEPTH_MIN)
        & (rows["Reliable_depth_vmr"] > ACTIVE_VMR_MIN)
    )
    rows["active_delta_pass"] = (
        ~rows["active_delta_trigger"] | (rows["ANI_AF_delta"] < ACTIVE_DELTA_MAX)
    )
    rows["active_gate_pass"] = (
        (rows["XnY_ctx"] >= ACTIVE_CTX_MIN)
        & (rows["ANI_naive_calc"] > ANI_THRESHOLD)
        & (rows["Reliable_ztp_af"] >= ACTIVE_RELIABLE_ZTP_AF_FLOOR)
        & rows["active_delta_pass"]
        & rows["gtdb_species"].astype(bool)
    )
    return rows


def load_sylph_rows(path: Path, by_accession, by_core) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows["accession"] = rows["Genome_file"].map(truth.extract_accession)
    mapped = [truth.lookup_accession(acc, by_accession, by_core) for acc in rows["accession"]]
    rows["gtdb_species"] = [rec["gtdb_species"] if rec else "" for rec, _, _ in mapped]
    rows["gtdb_taxonomy"] = [rec["gtdb_taxonomy"] if rec else "" for rec, _, _ in mapped]
    rows["Normalized_abundance_depth"] = numeric(rows, "Taxonomic_abundance") / 100.0
    rows["ANI_naive_calc"] = numeric(rows, "Adjusted_ANI") / 100.0
    rows["XnY_ctx"] = 0.0
    rows["Ref_zip_af"] = 0.0
    rows["Reliable_Ref_breadth"] = 0.0
    rows["Reliable_ztp_af"] = 0.0
    rows["Reliable_Ref_hit_mean_depth"] = numeric(rows, "Mean_cov_geq1")
    rows["Reliable_depth_vmr"] = 0.0
    rows["ANI_AF_delta"] = 0.0
    rows["Ref"] = rows["Genome_file"]
    rows["active_gate_pass"] = rows["gtdb_species"].astype(bool)
    return rows


def append_details(
    detail_rows: list[dict[str, object]],
    sample: int,
    method: str,
    kind: str,
    species_set: set[str],
    rows: pd.DataFrame,
    gold_abundance: Mapping[str, float],
) -> None:
    for species in sorted(species_set):
        rec: dict[str, object] = {
            "sample": sample,
            "method": method,
            "kind": kind,
            "gtdb_species": species,
            "gold_abundance": gold_abundance.get(species, ""),
            "accession": "",
            "Ref": "",
            "ANI_naive_calc": "",
            "XnY_ctx": "",
            "Ref_zip_af": "",
            "Reliable_Ref_breadth": "",
            "Reliable_ztp_af": "",
            "Reliable_Ref_hit_mean_depth": "",
            "Reliable_depth_vmr": "",
            "ANI_AF_delta": "",
            "Normalized_abundance_depth": "",
            "active_gate_pass": "",
        }
        examples = rows.loc[rows["gtdb_species"] == species].copy()
        if not examples.empty:
            best = examples.sort_values(
                ["active_gate_pass", "ANI_naive_calc", "Normalized_abundance_depth"],
                ascending=False,
            ).iloc[0]
            for field in [
                "accession",
                "Ref",
                "ANI_naive_calc",
                "XnY_ctx",
                "Ref_zip_af",
                "Reliable_Ref_breadth",
                "Reliable_ztp_af",
                "Reliable_Ref_hit_mean_depth",
                "Reliable_depth_vmr",
                "ANI_AF_delta",
                "Normalized_abundance_depth",
                "active_gate_pass",
            ]:
                rec[field] = best.get(field, "")
        detail_rows.append(rec)


def abundance_metrics(selected: pd.DataFrame, gold_profile: Sequence[Mapping[str, object]]):
    pred = (
        selected.loc[selected["gtdb_species"].astype(bool)]
        .groupby("gtdb_species")["Normalized_abundance_depth"]
        .max()
        .to_dict()
    )
    gold = {
        str(row["gtdb_species"]): float(row["relative_abundance"])
        for row in gold_profile
        if str(row.get("gtdb_species", ""))
    }
    rows = []
    for renorm in [False, True]:
        values = pred.copy()
        if renorm:
            total = sum(values.values())
            if total > 0.0:
                values = {k: v / total for k, v in values.items()}
        y_true = []
        y_pred = []
        for species, abundance in gold.items():
            y_true.append(float(abundance))
            y_pred.append(float(values.get(species, 0.0)))
        true_s = pd.Series(y_true, dtype=float)
        pred_s = pd.Series(y_pred, dtype=float)
        rows.append(
            {
                "renorm_pred": renorm,
                "truth_species": len(y_true),
                "pred_sum_on_truth": sum(y_pred),
                "pred_sum_all": sum(values.values()),
                "pearson": pred_s.corr(true_s, method="pearson"),
                "spearman": pred_s.corr(true_s, method="spearman"),
                "mae_pct_points": (pred_s - true_s).abs().mean() * 100.0,
                "l1_pct_points": (pred_s - true_s).abs().sum() * 100.0,
            }
        )
    return rows


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    genome_taxid, genome_location = load_meta()
    score_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []
    truth_summary_rows: list[dict[str, object]] = []

    profiles: dict[int, list[dict[str, object]]] = {}
    for sample in SAMPLES:
        source_rows = build_gtdb_truth(sample, by_accession, by_core, genome_taxid, genome_location)
        profiles[sample] = write_truth_files(sample, source_rows)
        unmapped = [r for r in source_rows if not str(r.get("gtdb_species", ""))]
        truth_summary_rows.append(
            {
                "sample": sample,
                "positive_genomes": len(source_rows),
                "mapped_positive_genomes": len(source_rows) - len(unmapped),
                "unmapped_positive_genomes": len(unmapped),
                "gtdb_species": len(profiles[sample]),
                "raw_abundance_sum": sum(float(r["abundance"]) for r in source_rows),
            }
        )

    for sample in SAMPLES:
        profile = profiles[sample]
        gold = {str(row["gtdb_species"]) for row in profile}
        gold_abundance = {str(row["gtdb_species"]): float(row["relative_abundance"]) for row in profile}
        methods = [
            ("ctx_only_current_best", minco_output_path(sample, "ctx_only_current_best"), "minco"),
            (
                "ctxobj_markerdb_product0_active",
                minco_output_path(sample, "ctxobj_markerdb_product0_active"),
                "minco",
            ),
            ("sylph_gtdb_profile", sylph_profile_path(sample), "sylph"),
        ]
        for method, path, kind in methods:
            if not path.exists():
                score_rows.append(
                    {
                        "sample": sample,
                        "method": method,
                        "status": "missing",
                        "path": str(path),
                        "gold_taxa": len(gold),
                        "pred_taxa": 0,
                        "TP": 0,
                        "FP": 0,
                        "FN": len(gold),
                        "precision": 0.0,
                        "recall": 0.0,
                        "F1": 0.0,
                        "selected_rows": 0,
                    }
                )
                continue
            rows = (
                load_minco_active_rows(path, by_accession, by_core)
                if kind == "minco"
                else load_sylph_rows(path, by_accession, by_core)
            )
            selected = rows.loc[rows["active_gate_pass"]].copy()
            tp, fp, fn, precision, recall, f1 = score_sets(selected["gtdb_species"], gold)
            score_rows.append(
                {
                    "sample": sample,
                    "method": method,
                    "status": "ok",
                    "path": str(path),
                    "gold_taxa": len(gold),
                    "pred_taxa": len(tp | fp),
                    "TP": len(tp),
                    "FP": len(fp),
                    "FN": len(fn),
                    "precision": precision,
                    "recall": recall,
                    "F1": f1,
                    "selected_rows": len(selected),
                }
            )
            append_details(detail_rows, sample, method, "FP", fp, rows, gold_abundance)
            append_details(detail_rows, sample, method, "FN", fn, rows, gold_abundance)
            for rec in abundance_metrics(selected, profile):
                abundance_rows.append({"sample": sample, "method": method, **rec})

    score_df = pd.DataFrame(score_rows)
    ok_score = score_df.loc[score_df["status"] == "ok"].copy()
    mean_df = (
        ok_score.groupby("method", as_index=False)[["precision", "recall", "F1", "TP", "FP", "FN"]]
        .mean()
        .rename(
            columns={
                "precision": "mean_precision",
                "recall": "mean_recall",
                "F1": "mean_F1",
                "TP": "mean_TP",
                "FP": "mean_FP",
                "FN": "mean_FN",
            }
        )
    )
    total_rows = []
    for method, sub in ok_score.groupby("method"):
        tp = int(sub["TP"].sum())
        fp = int(sub["FP"].sum())
        fn = int(sub["FN"].sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
        total_rows.append(
            {
                "method": method,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": precision,
                "recall": recall,
                "F1": f1,
            }
        )
    score_df.to_csv(EXP_DIR / "scores.tsv", sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(EXP_DIR / "details.tsv", sep="\t", index=False)
    pd.DataFrame(abundance_rows).to_csv(EXP_DIR / "abundance.tsv", sep="\t", index=False)
    pd.DataFrame(truth_summary_rows).to_csv(EXP_DIR / "truth_summary.tsv", sep="\t", index=False)
    mean_df.to_csv(EXP_DIR / "mean_scores.tsv", sep="\t", index=False)
    pd.DataFrame(total_rows).to_csv(EXP_DIR / "total_scores.tsv", sep="\t", index=False)
    print(score_df.to_csv(sep="\t", index=False), end="")
    print(f"wrote {EXP_DIR / 'scores.tsv'}")
    print(f"wrote {EXP_DIR / 'mean_scores.tsv'}")
    print(f"wrote {EXP_DIR / 'total_scores.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
