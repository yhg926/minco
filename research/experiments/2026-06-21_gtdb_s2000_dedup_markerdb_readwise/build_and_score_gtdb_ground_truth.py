#!/usr/bin/env python3
"""Build Toy Mouse sample0 GTDB-species truth and score direct MinCO calls."""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd


EXP_DIR = Path(__file__).resolve().parent
SOURCE_GENOMES = Path(
    "/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/"
    "positive_source_genomes.tsv"
)
GTDB_METADATA = [
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz"),
    Path("/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz"),
]

METHODS = [
    (
        "s1000_unique_direct",
        Path(
            "/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/"
            "toymouse_sample0_s1000_gtdb_unique_zip_unfiltered.tsv"
        ),
    ),
    (
        "s1000_split_direct",
        Path(
            "/mnt/new3T/minco_cami2_toymouse_20260621/minco_outputs/"
            "toymouse_sample0_s1000_gtdb_split_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_unique_direct",
        Path(
            "/tmp/gtdb232_s2000_dedup_marker.qKJofv/"
            "toymouse_sample0_s2000_dedup_ctxmarker_unique_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_split_direct",
        Path(
            "/tmp/gtdb232_s2000_dedup_marker.qKJofv/"
            "toymouse_sample0_s2000_dedup_ctxmarker_split_zip_unfiltered.tsv"
        ),
    ),
    (
        "s2000_marker_split_naive_unfiltered",
        EXP_DIR / "toymouse_sample0_s2000_marker_split_naive_unfiltered.tsv",
    ),
    (
        "s2000_marker_split_zip_poisson_depth_p005",
        EXP_DIR / "toymouse_sample0_s2000_marker_split_zip_poisson_depth_p005.tsv",
    ),
    (
        "s2000_marker_split_naive_poisson_depth_p005",
        EXP_DIR / "toymouse_sample0_s2000_marker_split_naive_poisson_depth_p005.tsv",
    ),
    (
        "s2000_marker_split_naive_poisson_product_p005",
        EXP_DIR / "toymouse_sample0_s2000_marker_split_naive_poisson_product_p005.tsv",
    ),
]

ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")
EPSILON = 1e-8


def extract_accession(text: object) -> str:
    if not isinstance(text, str):
        return ""
    match = ACC_RE.search(text)
    return match.group(1) if match else ""


def accession_core(accession: str) -> str:
    acc = extract_accession(accession)
    return acc.split("_", 1)[1] if acc else ""


def gtdb_species(taxonomy: str) -> str:
    for part in str(taxonomy).split(";"):
        if part.startswith("s__"):
            return part
    return ""


def ncbi_species_from_taxonomy(taxonomy: str) -> str:
    for part in str(taxonomy).split(";"):
        if part.startswith("s__"):
            name = part[3:].strip()
            if name:
                return name
    return ""


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def metadata_record(row: Mapping[str, str], domain: str) -> Dict[str, str]:
    accession = extract_accession(row.get("accession", ""))
    genbank_accession = extract_accession(row.get("ncbi_genbank_assembly_accession", ""))
    taxonomy = row.get("gtdb_taxonomy", "")
    ncbi_species = ncbi_species_from_taxonomy(row.get("ncbi_taxonomy", ""))
    if not ncbi_species:
        ncbi_species = row.get("ncbi_organism_name", "")
    return {
        "metadata_domain": domain,
        "metadata_accession": row.get("accession", ""),
        "accession": accession,
        "genbank_accession": genbank_accession,
        "gtdb_species": gtdb_species(taxonomy),
        "gtdb_taxonomy": taxonomy,
        "gtdb_representative": row.get("gtdb_representative", ""),
        "gtdb_genome_representative": row.get("gtdb_genome_representative", ""),
        "ncbi_taxid": row.get("ncbi_taxid", ""),
        "ncbi_species_taxid": row.get("ncbi_species_taxid", ""),
        "ncbi_organism_name": row.get("ncbi_organism_name", ""),
        "ncbi_species": ncbi_species,
        "ncbi_taxonomy": row.get("ncbi_taxonomy", ""),
    }


def load_gtdb_metadata(
    paths: Sequence[Path],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, List[Dict[str, str]]]]:
    by_accession: Dict[str, Dict[str, str]] = {}
    by_core: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    seen_core_record: set[Tuple[str, int]] = set()

    for path in paths:
        domain = "ar53" if "ar53" in path.name else "bac120"
        with open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rec = metadata_record(row, domain)
                if not rec["gtdb_species"]:
                    continue
                for acc in [rec["accession"], rec["genbank_accession"]]:
                    if not acc:
                        continue
                    by_accession.setdefault(acc, rec)
                    core = accession_core(acc)
                    key = (core, id(rec))
                    if core and key not in seen_core_record:
                        by_core[core].append(rec)
                        seen_core_record.add(key)
    return by_accession, by_core


def unique_records(records: Iterable[Dict[str, str]]) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    seen: set[Tuple[str, str]] = set()
    for rec in records:
        key = (rec.get("metadata_accession", ""), rec.get("gtdb_taxonomy", ""))
        if key in seen:
            continue
        out.append(rec)
        seen.add(key)
    return out


def lookup_accession(
    accession: str,
    by_accession: Mapping[str, Dict[str, str]],
    by_core: Mapping[str, List[Dict[str, str]]],
) -> Tuple[Optional[Dict[str, str]], str, str]:
    acc = extract_accession(accession)
    if not acc:
        return None, "unmapped", "no_accession"
    rec = by_accession.get(acc)
    if rec is not None:
        return rec, "exact_accession", acc
    candidates = unique_records(by_core.get(accession_core(acc), []))
    if len(candidates) == 1:
        return candidates[0], "unique_assembly_core", accession_core(acc)
    if len(candidates) > 1:
        return None, "ambiguous_assembly_core", accession_core(acc)
    return None, "unmapped", acc


def source_truth_rows(
    by_accession: Mapping[str, Dict[str, str]],
    by_core: Mapping[str, List[Dict[str, str]]],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with SOURCE_GENOMES.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for src in reader:
            abundance = float(src["abundance"])
            rec, method, key = lookup_accession(src["source_accession"], by_accession, by_core)
            if rec is None:
                rows.append(
                    {
                        **src,
                        "abundance": abundance,
                        "source_lookup_key": key,
                        "mapping_method": method,
                        "metadata_accession": "",
                        "metadata_domain": "",
                        "metadata_genbank_accession": "",
                        "gtdb_species": "",
                        "gtdb_taxonomy": "",
                        "ncbi_taxid_metadata": "",
                        "ncbi_species_taxid_metadata": "",
                        "ncbi_species_metadata": "",
                        "ncbi_organism_name_metadata": "",
                    }
                )
                continue
            rows.append(
                {
                    **src,
                    "abundance": abundance,
                    "source_lookup_key": key,
                    "mapping_method": method,
                    "metadata_accession": rec["metadata_accession"],
                    "metadata_domain": rec["metadata_domain"],
                    "metadata_genbank_accession": rec["genbank_accession"],
                    "gtdb_species": rec["gtdb_species"],
                    "gtdb_taxonomy": rec["gtdb_taxonomy"],
                    "ncbi_taxid_metadata": rec["ncbi_taxid"],
                    "ncbi_species_taxid_metadata": rec["ncbi_species_taxid"],
                    "ncbi_species_metadata": rec["ncbi_species"],
                    "ncbi_organism_name_metadata": rec["ncbi_organism_name"],
                }
            )
    return rows


def write_source_truth(rows: Sequence[Mapping[str, object]], path: Path) -> None:
    fields = [
        "genome_id",
        "abundance",
        "source_taxid",
        "source_accession",
        "source_basename",
        "source_lookup_key",
        "mapping_method",
        "metadata_accession",
        "metadata_domain",
        "metadata_genbank_accession",
        "ncbi_taxid_metadata",
        "ncbi_species_taxid_metadata",
        "ncbi_species_metadata",
        "ncbi_organism_name_metadata",
        "gtdb_species",
        "gtdb_taxonomy",
        "tar_path",
        "local_path",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_species_profile(rows: Sequence[Mapping[str, object]], path: Path) -> List[Dict[str, object]]:
    grouped: Dict[str, Dict[str, object]] = {}
    total = sum(float(row["abundance"]) for row in rows if row.get("gtdb_species"))
    for row in rows:
        species = str(row.get("gtdb_species", ""))
        if not species:
            continue
        rec = grouped.setdefault(
            species,
            {
                "gtdb_species": species,
                "gtdb_taxonomy": row.get("gtdb_taxonomy", ""),
                "source_genomes": 0,
                "total_abundance": 0.0,
                "relative_abundance": 0.0,
                "source_accessions": [],
                "source_genome_ids": [],
                "source_taxids": [],
                "source_ncbi_species": [],
                "mapping_methods": [],
            },
        )
        rec["source_genomes"] = int(rec["source_genomes"]) + 1
        rec["total_abundance"] = float(rec["total_abundance"]) + float(row["abundance"])
        for field, key in [
            ("source_accessions", "source_accession"),
            ("source_genome_ids", "genome_id"),
            ("source_taxids", "source_taxid"),
            ("source_ncbi_species", "ncbi_species_metadata"),
            ("mapping_methods", "mapping_method"),
        ]:
            value = str(row.get(key, ""))
            if value and value not in rec[field]:
                rec[field].append(value)

    out = list(grouped.values())
    for rec in out:
        rec["relative_abundance"] = float(rec["total_abundance"]) / total if total else 0.0
    out.sort(key=lambda r: (-float(r["total_abundance"]), str(r["gtdb_species"])))

    fields = [
        "gtdb_species",
        "gtdb_taxonomy",
        "source_genomes",
        "total_abundance",
        "relative_abundance",
        "source_accessions",
        "source_genome_ids",
        "source_taxids",
        "source_ncbi_species",
        "mapping_methods",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rec in out:
            row = rec.copy()
            for field in fields:
                if isinstance(row.get(field), list):
                    row[field] = ",".join(row[field])
            writer.writerow(row)
    return out


def to_float_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([0.0] * len(df), index=df.index)
    return pd.to_numeric(df[col], errors="coerce").fillna(0.0)


def minco_naive_ani(xny_ctx: float, n_diff_obj: float, n_diff_obj_section: float) -> float:
    if xny_ctx <= 0.0:
        return 0.0
    ratio = (n_diff_obj_section + EPSILON) / (n_diff_obj + EPSILON)
    dist0 = n_diff_obj / (xny_ctx + n_diff_obj)
    final_dist = 1.0 - math.pow(1.0 - dist0, ratio)
    if final_dist <= 0.0:
        naive_dist = 0.0
    else:
        naive_dist = final_dist * 0.1544286
    return max(0.0, min(1.0, 1.0 - naive_dist))


def add_naive_ani(rows: pd.DataFrame) -> pd.DataFrame:
    xny = to_float_series(rows, "XnY_ctx")
    ndiff = to_float_series(rows, "N_diff_obj")
    nsection = to_float_series(rows, "N_diff_obj_section")
    rows["ANI_naive_calc"] = [
        minco_naive_ani(float(x), float(d), float(s)) for x, d, s in zip(xny, ndiff, nsection)
    ]
    return rows


def direct_mask(
    rows: pd.DataFrame,
    ani_col: str = "ANI",
    ani_threshold: float = 0.94,
    ani_strict: bool = False,
) -> pd.Series:
    ani_values = to_float_series(rows, ani_col)
    ani_pass = ani_values > ani_threshold if ani_strict else ani_values >= ani_threshold
    return (
        (to_float_series(rows, "XnY_ctx") >= 10.0)
        & ani_pass
        & (to_float_series(rows, "Real_min_align_fraction") >= 0.05)
    )


def breadth_depth_reject(
    rows: pd.DataFrame,
    hit_depth_threshold: float = 4.0,
    breadth_floor: float = 0.5,
) -> pd.Series:
    return (
        (to_float_series(rows, "Ref_hit_mean_depth") >= hit_depth_threshold)
        & (to_float_series(rows, "Ref_breadth") < breadth_floor)
    )


def adjusted_breadth_reject(rows: pd.DataFrame, adjusted_breadth_floor: float = 0.5) -> pd.Series:
    return to_float_series(rows, "Ref_zip_af") < adjusted_breadth_floor


def score_set(predicted: Iterable[str], gold: Iterable[str]) -> Dict[str, object]:
    pred = {str(x) for x in predicted if str(x)}
    truth = {str(x) for x in gold if str(x)}
    tp_set = pred & truth
    fp_set = pred - truth
    fn_set = truth - pred
    precision = len(tp_set) / len(pred) if pred else 0.0
    recall = len(tp_set) / len(truth) if truth else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "pred_taxa": len(pred),
        "TP": len(tp_set),
        "FP": len(fp_set),
        "FN": len(fn_set),
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "FP_set": fp_set,
        "FN_set": fn_set,
    }


def load_minco_with_gtdb_species(
    path: Path,
    by_accession: Mapping[str, Dict[str, str]],
    by_core: Mapping[str, List[Dict[str, str]]],
) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    rows["accession"] = rows["Ref"].map(extract_accession)
    if "Ref_annotation" in rows.columns:
        missing = rows["accession"] == ""
        rows.loc[missing, "accession"] = rows.loc[missing, "Ref_annotation"].map(extract_accession)
    mapped = [lookup_accession(acc, by_accession, by_core) for acc in rows["accession"]]
    rows["gtdb_species"] = [rec["gtdb_species"] if rec else "" for rec, _, _ in mapped]
    rows["gtdb_taxonomy"] = [rec["gtdb_taxonomy"] if rec else "" for rec, _, _ in mapped]
    rows["gtdb_mapping_method"] = [method for _, method, _ in mapped]
    rows["gtdb_mapping_key"] = [key for _, _, key in mapped]
    return rows


def write_score_detail(
    path: Path,
    kind: str,
    details: Sequence[Mapping[str, object]],
) -> None:
    fields = ["method", "kind", "gtdb_species", "gold_abundance", "example_ref", "example_ani", "example_xny"]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in details:
            writer.writerow({**row, "kind": kind})


def score_methods(
    profile_rows: Sequence[Mapping[str, object]],
    by_accession: Mapping[str, Dict[str, str]],
    by_core: Mapping[str, List[Dict[str, str]]],
    score_path: Path,
    fp_path: Path,
    fn_path: Path,
    ani_col: str,
    ani_threshold: float = 0.94,
    ani_strict: bool = False,
    extra_reject=None,
    extra_reject_label: str = "",
) -> None:
    gold_species = {str(row["gtdb_species"]) for row in profile_rows}
    abundance_by_species = {
        str(row["gtdb_species"]): float(row["relative_abundance"]) for row in profile_rows
    }

    score_rows: List[Dict[str, object]] = []
    false_positive_rows: List[Dict[str, object]] = []
    false_negative_rows: List[Dict[str, object]] = []

    for method, path in METHODS:
        rows = load_minco_with_gtdb_species(path, by_accession, by_core)
        rows = add_naive_ani(rows)
        mask = direct_mask(rows, ani_col, ani_threshold, ani_strict)
        rejected_by_extra = extra_reject(rows) & mask if extra_reject is not None else pd.Series(
            [False] * len(rows), index=rows.index
        )
        selected_rows = rows.loc[mask & ~rejected_by_extra].copy()
        selected_rows = selected_rows.loc[selected_rows["gtdb_species"].astype(bool)].copy()
        stats = score_set(selected_rows["gtdb_species"], gold_species)
        score_row = {
            "method": method,
            "gold_taxa": len(gold_species),
            "pred_taxa": stats["pred_taxa"],
            "TP": stats["TP"],
            "FP": stats["FP"],
            "FN": stats["FN"],
            "precision": f"{stats['precision']:.12g}",
            "recall": f"{stats['recall']:.12g}",
            "F1": f"{stats['F1']:.12g}",
            "selected_rows": len(selected_rows),
            "unmapped_rows": int((rows["gtdb_species"] == "").sum()),
            "selected_unmapped_rows": int((rows.loc[mask & ~rejected_by_extra, "gtdb_species"] == "").sum()),
        }
        if extra_reject_label:
            score_row["candidate_rows_before_filter"] = int(mask.sum())
            score_row[f"rejected_by_{extra_reject_label}"] = int(rejected_by_extra.sum())
        score_rows.append(score_row)

        for species in sorted(stats["FP_set"]):
            examples = selected_rows.loc[selected_rows["gtdb_species"] == species].copy()
            examples[ani_col] = to_float_series(examples, ani_col)
            examples["XnY_ctx"] = to_float_series(examples, "XnY_ctx")
            best = examples.sort_values([ani_col, "XnY_ctx"], ascending=False).iloc[0]
            false_positive_rows.append(
                {
                    "method": method,
                    "gtdb_species": species,
                    "gold_abundance": "",
                    "example_ref": best.get("Ref", ""),
                    "example_ani": f"{float(best.get(ani_col, 0.0)):.6f}",
                    "example_xny": f"{float(best.get('XnY_ctx', 0.0)):.0f}",
                }
            )

        predicted_species = set(selected_rows["gtdb_species"])
        for species in sorted(stats["FN_set"]):
            false_negative_rows.append(
                {
                    "method": method,
                    "gtdb_species": species,
                    "gold_abundance": f"{abundance_by_species.get(species, 0.0):.12g}",
                    "example_ref": "",
                    "example_ani": "",
                    "example_xny": "",
                }
            )
        assert predicted_species <= {str(row["gtdb_species"]) for row in profile_rows} | set(
            stats["FP_set"]
        )

    with score_path.open("w", newline="") as fh:
        fields = list(score_rows[0])
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(score_rows)
    write_score_detail(fp_path, "FP", false_positive_rows)
    write_score_detail(fn_path, "FN", false_negative_rows)


def write_summary(
    source_rows: Sequence[Mapping[str, object]],
    profile_rows: Sequence[Mapping[str, object]],
    path: Path,
) -> None:
    total_abundance = sum(float(row["abundance"]) for row in source_rows)
    mapped_rows = [row for row in source_rows if row.get("gtdb_species")]
    exact = sum(1 for row in mapped_rows if row.get("mapping_method") == "exact_accession")
    core = sum(1 for row in mapped_rows if row.get("mapping_method") == "unique_assembly_core")
    profile_total = sum(float(row["relative_abundance"]) for row in profile_rows)
    fields = ["metric", "value", "notes"]
    rows = [
        {
            "metric": "source_genomes",
            "value": len(source_rows),
            "notes": str(SOURCE_GENOMES),
        },
        {
            "metric": "mapped_source_genomes",
            "value": len(mapped_rows),
            "notes": f"exact_accession={exact}; unique_assembly_core={core}",
        },
        {
            "metric": "unmapped_source_genomes",
            "value": len(source_rows) - len(mapped_rows),
            "notes": "must be zero before GTDB benchmark scoring",
        },
        {
            "metric": "gtdb_species",
            "value": len(profile_rows),
            "notes": "truth taxa after aggregating source genomes to GTDB species",
        },
        {
            "metric": "raw_abundance_sum",
            "value": f"{total_abundance:.12g}",
            "notes": "CAMISIM source abundance units",
        },
        {
            "metric": "relative_abundance_sum",
            "value": f"{profile_total:.12g}",
            "notes": "normalized GTDB species profile",
        },
        {
            "metric": "gtdb_metadata",
            "value": "r232",
            "notes": ",".join(str(p) for p in GTDB_METADATA),
        },
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    by_accession, by_core = load_gtdb_metadata(GTDB_METADATA)
    source_rows = source_truth_rows(by_accession, by_core)
    if any(not row.get("gtdb_species") for row in source_rows):
        missing = [str(row.get("source_accession", "")) for row in source_rows if not row.get("gtdb_species")]
        raise RuntimeError(f"unmapped source genomes in GTDB metadata: {', '.join(missing)}")

    write_source_truth(source_rows, EXP_DIR / "mouse0_gtdb_source_species_ground_truth.tsv")
    profile_rows = write_species_profile(source_rows, EXP_DIR / "mouse0_gtdb_species_profile.tsv")
    write_summary(source_rows, profile_rows, EXP_DIR / "mouse0_gtdb_truth_summary.tsv")
    score_methods(
        profile_rows,
        by_accession,
        by_core,
        EXP_DIR / "gtdb_direct_scores.tsv",
        EXP_DIR / "gtdb_direct_score_false_positives.tsv",
        EXP_DIR / "gtdb_direct_score_false_negatives.tsv",
        "ANI",
    )
    score_methods(
        profile_rows,
        by_accession,
        by_core,
        EXP_DIR / "gtdb_direct_scores_naive_ani.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_false_positives.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_false_negatives.tsv",
        "ANI_naive_calc",
    )
    score_methods(
        profile_rows,
        by_accession,
        by_core,
        EXP_DIR / "gtdb_direct_scores_naive_ani_gt095.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_false_positives.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_false_negatives.tsv",
        "ANI_naive_calc",
        0.95,
        True,
    )
    score_methods(
        profile_rows,
        by_accession,
        by_core,
        EXP_DIR / "gtdb_direct_scores_naive_ani_gt095_breadth_depth.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_breadth_depth_false_positives.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_breadth_depth_false_negatives.tsv",
        "ANI_naive_calc",
        0.95,
        True,
        lambda rows: breadth_depth_reject(rows, 4.0, 0.5),
        "breadth_depth",
    )
    score_methods(
        profile_rows,
        by_accession,
        by_core,
        EXP_DIR / "gtdb_direct_scores_naive_ani_gt095_zipaf05.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_zipaf05_false_positives.tsv",
        EXP_DIR / "gtdb_direct_score_naive_ani_gt095_zipaf05_false_negatives.tsv",
        "ANI_naive_calc",
        0.95,
        True,
        lambda rows: adjusted_breadth_reject(rows, 0.5),
        "zipaf05",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
