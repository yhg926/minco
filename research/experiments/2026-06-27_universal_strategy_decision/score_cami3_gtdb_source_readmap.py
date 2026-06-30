#!/usr/bin/env python3
"""Score CAMI3 Toy Human Gut samples0-2 using source-readmap GTDB truth.

This scorer is stronger than taxid-only transfer because it maps CAMISIM source
read rows to GTDB species through the source contig WGS prefix when possible.
Remaining rows fall back to unique GTDB mappings of the CAMI/NCBI taxid. Read
rows that cannot be mapped by either route are reported and excluded from the
normalized transferred truth set.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"
TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
sys.path.insert(0, str(TRUTH_HELPER))
import build_and_score_gtdb_ground_truth as truth  # noqa: E402

READWISE_HELPER = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(READWISE_HELPER))
from analyze_readwise_corrections import load_minco, parse_species_taxmap  # noqa: E402

import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402


CAMI3_ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620")
CAMI3_EXTRA = Path("/mnt/new3T/minco_cami3_toygut_extra_20260621")
RUN = Path("/tmp/cami3_toy_human_gut_20260626/run")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")

READ_MAPPINGS = {
    0: CAMI3_ROOT / "sample_0_reads_mapping.tsv.gz",
    1: CAMI3_EXTRA / "sample_1_reads/reads_mapping.tsv.gz",
    2: CAMI3_EXTRA / "sample_2_reads/reads_mapping.tsv.gz",
}

PROFILE_PATHS = {
    0: {
        "minco": RUN / "sample0_universal_autoexact.tsv",
        "sylph": CAMI3_ROOT / "sylph_sample0/profile.tsv",
    },
    1: {
        "minco": RUN / "sample1_universal_autoexact.tsv",
        "sylph": Path("/tmp/cami3_toygut_current_minco_vs_sylph_20260622/sylph_sample1/profile.tsv"),
    },
    2: {
        "minco": RUN / "sample2_universal_autoexact_tail_p025.tsv",
        "sylph": RUN / "sylph_sample2/profile.tsv",
    },
}

RAW_TABLES = {
    0: {
        "unique": RUN / "sample0_autoexact/minco.best_diff_unique.unfiltered.tsv",
        "split": RUN / "sample0_autoexact/minco.best_diff_split.unfiltered.tsv",
    },
    1: {
        "unique": RUN / "sample1_autoexact/minco.best_diff_unique.unfiltered.tsv",
        "split": RUN / "sample1_autoexact/minco.best_diff_split.unfiltered.tsv",
    },
    2: {
        "unique": RUN / "sample2_autoexact/minco.best_diff_unique.unfiltered.tsv",
        "split": RUN / "sample2_autoexact/minco.best_diff_split.unfiltered.tsv",
    },
}

PREFIX_RE = re.compile(r"([A-Z]+)")


def gtdb_species(taxonomy: str) -> str:
    for part in str(taxonomy).split(";"):
        if part.startswith("s__"):
            return part
    return ""


def ncbi_species_from_taxonomy(taxonomy: str) -> str:
    for part in str(taxonomy).split(";"):
        if part.startswith("s__"):
            return part[3:].strip()
    return ""


def normalize_name(value: object) -> str:
    return " ".join(str(value or "").strip().lower().split())


def accession_prefix(value: object) -> str:
    text = str(value or "").replace("NZ_", "")
    if not text or text == "na":
        return ""
    core = text.split(".", 1)[0]
    match = PREFIX_RE.match(core)
    return match.group(1) if match else core


def contig_from_read_id(read_id: str) -> str:
    core = str(read_id).strip().split("/", 1)[0]
    if "-" in core:
        core = core.rsplit("-", 1)[0]
    return core


def unique_map(values: Mapping[str, set[str]]) -> dict[str, str]:
    return {key: next(iter(species)) for key, species in values.items() if len(species) == 1}


def build_transfer_maps() -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, object]]:
    wgs_to_species: dict[str, set[str]] = defaultdict(set)
    taxid_to_species: dict[str, set[str]] = defaultdict(set)
    name_to_species: dict[str, set[str]] = defaultdict(set)
    rows = 0

    for path in truth.GTDB_METADATA:
        with truth.open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                species = gtdb_species(row.get("gtdb_taxonomy", ""))
                if not species:
                    continue
                rows += 1
                wgs_prefix = accession_prefix(row.get("ncbi_wgs_master", ""))
                if wgs_prefix:
                    wgs_to_species[wgs_prefix].add(species)
                for key in ["ncbi_species_taxid", "ncbi_taxid"]:
                    taxid = str(row.get(key, "")).strip()
                    if taxid:
                        taxid_to_species[taxid].add(species)
                for name in [
                    ncbi_species_from_taxonomy(row.get("ncbi_taxonomy", "")),
                    row.get("ncbi_organism_name", ""),
                ]:
                    normalized = normalize_name(name)
                    if normalized:
                        name_to_species[normalized].add(species)

    wgs_unique = unique_map(wgs_to_species)
    taxid_unique = unique_map(taxid_to_species)
    name_unique = unique_map(name_to_species)
    diagnostics = {
        "gtdb_metadata_rows_with_species": rows,
        "wgs_prefixes_total": len(wgs_to_species),
        "wgs_prefixes_unique": len(wgs_unique),
        "taxids_total": len(taxid_to_species),
        "taxids_unique": len(taxid_unique),
        "names_total": len(name_to_species),
        "names_unique": len(name_unique),
    }
    return wgs_unique, taxid_unique, name_unique, diagnostics


def build_source_truth(
    sample_id: int,
    path: Path,
    wgs_to_species: Mapping[str, str],
    taxid_to_species: Mapping[str, str],
) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame]:
    truth_counts: Counter[str] = Counter()
    source_genomes_by_species: dict[str, set[str]] = defaultdict(set)
    source_counts: dict[str, Counter[str]] = defaultdict(Counter)
    genome_taxid: dict[str, str] = {}
    cache: dict[tuple[str, str], tuple[str, str]] = {}
    total_rows = 0
    method_rows = Counter()

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            _read_name, genome_id, taxid, read_id = fields[:4]
            total_rows += 1
            prefix = accession_prefix(contig_from_read_id(read_id))
            cache_key = (taxid, prefix)
            species, method = cache.get(cache_key, ("", ""))
            if not method:
                species = wgs_to_species.get(prefix, "")
                method = "wgs_prefix" if species else ""
                if not species:
                    species = taxid_to_species.get(str(taxid), "")
                    method = "unique_taxid" if species else "unmapped"
                cache[cache_key] = (species, method)

            method_rows[method] += 1
            source_counts[genome_id][method] += 1
            genome_taxid.setdefault(genome_id, str(taxid))
            if species:
                truth_counts[species] += 1
                source_genomes_by_species[species].add(genome_id)

    mapped_rows = sum(truth_counts.values())
    truth_rows = []
    for species, count in sorted(truth_counts.items()):
        truth_rows.append(
            {
                "sample": sample_id,
                "gtdb_species": species,
                "read_rows": count,
                "truth_abundance": count / mapped_rows if mapped_rows else 0.0,
                "source_genomes": len(source_genomes_by_species[species]),
                "source_genome_ids": ",".join(sorted(source_genomes_by_species[species])),
            }
        )

    source_rows = []
    for genome_id, counts in sorted(source_counts.items()):
        total = sum(counts.values())
        method, method_count = counts.most_common(1)[0]
        source_rows.append(
            {
                "sample": sample_id,
                "source_genome_id": genome_id,
                "taxid": genome_taxid.get(genome_id, ""),
                "read_rows": total,
                "primary_mapping_method": method,
                "primary_mapping_read_rows": method_count,
                "wgs_prefix_read_rows": counts.get("wgs_prefix", 0),
                "unique_taxid_read_rows": counts.get("unique_taxid", 0),
                "unmapped_read_rows": counts.get("unmapped", 0),
            }
        )

    quality = {
        "sample": sample_id,
        "read_rows_total": total_rows,
        "read_rows_wgs_prefix": method_rows.get("wgs_prefix", 0),
        "read_rows_unique_taxid": method_rows.get("unique_taxid", 0),
        "read_rows_unmapped": method_rows.get("unmapped", 0),
        "read_rows_mapped": mapped_rows,
        "read_rows_mapped_pct": mapped_rows / total_rows * 100.0 if total_rows else 0.0,
        "truth_gtdb_species": len(truth_rows),
        "source_genomes_total": len(source_counts),
        "source_genomes_primary_wgs_prefix": sum(
            1 for counts in source_counts.values() if counts.most_common(1)[0][0] == "wgs_prefix"
        ),
        "source_genomes_primary_unique_taxid": sum(
            1 for counts in source_counts.values() if counts.most_common(1)[0][0] == "unique_taxid"
        ),
        "source_genomes_primary_unmapped": sum(
            1 for counts in source_counts.values() if counts.most_common(1)[0][0] == "unmapped"
        ),
    }
    return pd.DataFrame(truth_rows), quality, pd.DataFrame(source_rows)


def gtdb_from_accession(
    accession: str,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> str:
    rec, _method, _key = truth.lookup_accession(str(accession), by_accession, by_core)
    if rec is None:
        return ""
    return str(rec.get("gtdb_species", ""))


def best_raw_ref_species_by_taxid(
    raw_paths: Mapping[str, Path],
    taxmap: Mapping[str, Mapping[str, str]],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[dict[str, str], dict[str, object]]:
    frames = []
    for mode in ["split", "unique"]:
        path = raw_paths.get(mode)
        if not path or not path.exists():
            continue
        rows = load_minco(path, taxmap, 11.0)
        if rows.empty:
            continue
        rows = rows.loc[rows["taxid"].astype(bool)].copy()
        rows["raw_mode"] = mode
        frames.append(rows)
    if not frames:
        return {}, {"raw_ref_rows": 0, "raw_ref_taxids": 0, "raw_ref_taxids_gtdb_mapped": 0}

    work = pd.concat(frames, ignore_index=True, sort=False)
    for col in ["XnY_ctx", "Real_min_align_fraction", "ANI", "Ref_breadth", "Ref_mean_depth"]:
        work[col] = pd.to_numeric(work.get(col, 0.0), errors="coerce").fillna(0.0)
    work["raw_mode_rank"] = work["raw_mode"].map({"split": 1, "unique": 0}).fillna(0)
    work = work.sort_values(
        ["taxid", "raw_mode_rank", "XnY_ctx", "Real_min_align_fraction", "ANI", "Ref_breadth", "Ref_mean_depth"],
        ascending=[True, False, False, False, False, False, False],
        kind="mergesort",
    )
    best = work.drop_duplicates("taxid", keep="first")
    out: dict[str, str] = {}
    for row in best.itertuples(index=False):
        species = gtdb_from_accession(getattr(row, "accession", ""), by_accession, by_core)
        if species:
            out[str(getattr(row, "taxid"))] = species
    return out, {
        "raw_ref_rows": int(len(work)),
        "raw_ref_taxids": int(best["taxid"].astype(str).nunique()),
        "raw_ref_taxids_gtdb_mapped": int(len(out)),
    }


def load_minco_predictions(
    path: Path,
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
    best_ref_species_by_taxid: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t")
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[call].copy()
    rows = []
    unmapped = 0
    name_mapped = 0
    taxid_mapped = 0
    ref_mapped = 0
    for row in selected.itertuples(index=False):
        taxid = str(getattr(row, "taxid", ""))
        species_name = normalize_name(getattr(row, "species_name", ""))
        gtdb_name = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                gtdb_name = gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if gtdb_name:
                    ref_mapped += 1
                    break
        if not gtdb_name:
            gtdb_name = best_ref_species_by_taxid.get(taxid, "")
            if gtdb_name:
                ref_mapped += 1
        if not gtdb_name:
            gtdb_name = taxid_to_species.get(taxid, "")
            if gtdb_name:
                taxid_mapped += 1
        if not gtdb_name and species_name:
            gtdb_name = name_to_species.get(species_name, "")
            if gtdb_name:
                name_mapped += 1
        if not gtdb_name:
            unmapped += 1
            continue
        abundance = float(getattr(row, "calibrated_abundance", 0.0) or 0.0)
        ani = float(getattr(row, "s_Ref_zip_aaf_ani_max", 0.0) or 0.0)
        rows.append((gtdb_name, abundance, ani))
    return taxid_score.collapse_prediction(rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_ref_mapped": ref_mapped,
        "pred_rows_taxid_mapped": taxid_mapped,
        "pred_rows_name_mapped": name_mapped,
    }


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    for sample_id, readmap in READ_MAPPINGS.items():
        if not readmap.exists():
            raise SystemExit(f"missing read mapping for sample {sample_id}: {readmap}")
        for label, profile in PROFILE_PATHS[sample_id].items():
            if not profile.exists():
                raise SystemExit(f"missing {label} profile for sample {sample_id}: {profile}")
        for label, raw_table in RAW_TABLES[sample_id].items():
            if not raw_table.exists():
                raise SystemExit(f"missing {label} raw MinCO table for sample {sample_id}: {raw_table}")
    if not TAXMAP.exists():
        raise SystemExit(f"missing profiler taxmap: {TAXMAP}")

    wgs_to_species, taxid_to_species, name_to_species, map_diag = build_transfer_maps()
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    taxmap = parse_species_taxmap(TAXMAP)

    score_rows = []
    quality_rows = []
    truth_rows = []
    source_rows = []
    for sample_id in sorted(READ_MAPPINGS):
        sample_truth, quality, sample_sources = build_source_truth(
            sample_id,
            READ_MAPPINGS[sample_id],
            wgs_to_species,
            taxid_to_species,
        )
        quality_rows.append({**map_diag, **quality})
        truth_rows.extend(sample_truth.to_dict(orient="records"))
        source_rows.extend(sample_sources.to_dict(orient="records"))
        best_ref_species, best_ref_diag = best_raw_ref_species_by_taxid(
            RAW_TABLES[sample_id],
            taxmap,
            by_accession,
            by_core,
        )
        quality_rows[-1].update(best_ref_diag)

        minco_pred, minco_extra = load_minco_predictions(
            PROFILE_PATHS[sample_id]["minco"],
            taxid_to_species,
            name_to_species,
            best_ref_species,
            by_accession,
            by_core,
        )
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(
            PROFILE_PATHS[sample_id]["sylph"],
            by_accession,
            by_core,
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "minco_universal_autoexact_gtdb_source_readmap",
                minco_pred,
                sample_truth,
                minco_extra,
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_gtdb_source_readmap",
                sylph_pred,
                sample_truth,
                sylph_extra,
            )
        )

    score_df = pd.DataFrame(score_rows)
    quality_df = pd.DataFrame(quality_rows)
    truth_df = pd.DataFrame(truth_rows)
    source_df = pd.DataFrame(source_rows)
    summary_df = taxid_score.summarize(score_df)

    truth_df.to_csv(RESULTS / "cami3_gtdb_source_readmap_truth.tsv", sep="\t", index=False)
    quality_df.to_csv(RESULTS / "cami3_gtdb_source_readmap_quality.tsv", sep="\t", index=False)
    source_df.to_csv(RESULTS / "cami3_gtdb_source_readmap_source_genomes.tsv", sep="\t", index=False)
    score_df.to_csv(RESULTS / "cami3_gtdb_source_readmap_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / "cami3_gtdb_source_readmap_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
