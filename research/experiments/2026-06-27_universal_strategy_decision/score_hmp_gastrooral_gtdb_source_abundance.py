#!/usr/bin/env python3
"""Score CAMI II HMP gastrooral samples in GTDB source-abundance space.

The earlier HMP gastrooral pilot used CAMI/NCBI taxid truth. This scorer maps
CAMISIM source genome accessions to GTDB species and compares current MinCO
universal table-mode outputs with same-release chunked Sylph r232 profiles.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping

import pandas as pd


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

TRUTH_HELPER = REPO_ROOT / "research/experiments/2026-06-21_gtdb_s2000_dedup_markerdb_readwise"
READWISE_HELPER = REPO_ROOT / "research/experiments/2026-06-20_minco_readwise_correction_research/scripts"
sys.path.insert(0, str(TRUTH_HELPER))
sys.path.insert(0, str(READWISE_HELPER))

import build_and_score_gtdb_ground_truth as truth  # noqa: E402
from analyze_readwise_corrections import load_minco, parse_species_taxmap  # noqa: E402
import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402


TRUTH_ROOT = Path("/tmp/cami2_hmp_pilot_20260625/truth")
RUN = Path("/tmp/cami2_hmp_pilot_20260625/run")
R232_RUN = Path("/tmp/cami2_hmp_pilot_20260625/run_sylph_r232")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")

SAMPLES = {
    0: {
        "label": "gastrooral0",
        "minco": RUN / "minco_sample0_universal_autoexact_tablemode_check.tsv",
        "unique": RUN / "minco_sample0_unique_zip_unfiltered.tsv",
        "split": RUN / "minco_sample0_split_zip_unfiltered.tsv",
        "sylph_r232": R232_RUN / "sylph_sample0/profile.chunked.tsv",
    },
    6: {
        "label": "gastrooral6",
        "minco": RUN / "minco_sample6_universal_autoexact_tablemode_current.tsv",
        "unique": RUN / "minco_sample6_unique_zip_unfiltered.tsv",
        "split": RUN / "minco_sample6_split_zip_unfiltered.tsv",
        "sylph_r232": R232_RUN / "sylph_sample6/profile.chunked.tsv",
    },
}

ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def extract_accession(value: object) -> str:
    match = ACC_RE.search(str(value or ""))
    return match.group(1) if match else ""


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").replace([math.inf, -math.inf], pd.NA).fillna(0.0)


def load_source_accessions(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                out[fields[0]] = extract_accession(fields[1])
    return out


def gtdb_from_accession(
    accession: str,
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[str, str]:
    rec, method, _key = truth.lookup_accession(str(accession), by_accession, by_core)
    if rec is None:
        return "", method
    return str(rec.get("gtdb_species", "")), method


def build_source_truth(
    sample_id: int,
    sample_label: str,
    source_accessions: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame]:
    abundance_path = TRUTH_ROOT / f"abundance{sample_id}.tsv"
    counts: Counter[str] = Counter()
    source_rows = []
    method_mass: Counter[str] = Counter()
    total_mass = 0.0
    mapped_mass = 0.0
    positive_sources = 0
    with abundance_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            genome_id, raw_value = line.rstrip("\n").split("\t")[:2]
            value = float(raw_value)
            if value <= 0.0:
                continue
            positive_sources += 1
            total_mass += value
            accession = source_accessions.get(genome_id, "")
            species, method = gtdb_from_accession(accession, by_accession, by_core)
            method_mass[method] += value
            if species:
                counts[species] += value
                mapped_mass += value
            source_rows.append(
                {
                    "sample": sample_id,
                    "source_genome_id": genome_id,
                    "source_accession": accession,
                    "raw_abundance": value,
                    "gtdb_species": species,
                    "mapping_method": method,
                }
            )

    truth_rows = [
        {
            "sample": sample_id,
            "gtdb_species": species,
            "truth_raw_abundance": value,
            "truth_abundance": value / mapped_mass if mapped_mass > 0.0 else 0.0,
        }
        for species, value in sorted(counts.items())
    ]
    quality = {
        "sample": sample_id,
        "sample_label": sample_label,
        "source_genomes_positive": positive_sources,
        "source_genomes_mapped": sum(1 for row in source_rows if row["gtdb_species"]),
        "truth_gtdb_species": len(truth_rows),
        "truth_mass_total": total_mass,
        "truth_mass_mapped": mapped_mass,
        "truth_mass_mapped_pct": mapped_mass / total_mass * 100.0 if total_mass else 0.0,
        "truth_mass_exact_accession": method_mass.get("exact_accession", 0.0),
        "truth_mass_unique_assembly_core": method_mass.get("unique_assembly_core", 0.0),
        "truth_mass_ambiguous_assembly_core": method_mass.get("ambiguous_assembly_core", 0.0),
        "truth_mass_unmapped": method_mass.get("unmapped", 0.0),
    }
    return pd.DataFrame(truth_rows), quality, pd.DataFrame(source_rows)


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
        species, _method = gtdb_from_accession(getattr(row, "accession", ""), by_accession, by_core)
        if species:
            out[str(getattr(row, "taxid"))] = species
    return out, {
        "raw_ref_rows": int(len(work)),
        "raw_ref_taxids": int(best["taxid"].astype(str).nunique()),
        "raw_ref_taxids_gtdb_mapped": int(len(out)),
    }


def load_minco_predictions(
    path: Path,
    best_ref_species_by_taxid: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    call = raw["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = raw.loc[call].copy()
    rows = []
    ref_mapped = 0
    fallback_mapped = 0
    unmapped = 0
    for row in selected.itertuples(index=False):
        species = ""
        for col in ["s_best_accession", "u_best_accession"]:
            if hasattr(row, col):
                species, _method = gtdb_from_accession(getattr(row, col, ""), by_accession, by_core)
                if species:
                    ref_mapped += 1
                    break
        if not species:
            species = best_ref_species_by_taxid.get(str(getattr(row, "taxid", "")), "")
            if species:
                fallback_mapped += 1
        if not species:
            unmapped += 1
            continue
        abundance = float(getattr(row, "calibrated_abundance", 0.0) or 0.0)
        ani = float(getattr(row, "reported_ani", 0.0) or getattr(row, "s_Ref_zip_aaf_ani_max", 0.0) or 0.0)
        rows.append((species, abundance, ani))
    return taxid_score.collapse_prediction(rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_ref_mapped": ref_mapped,
        "pred_rows_fallback_mapped": fallback_mapped,
        "minco_profile": str(path),
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--minco-sample0",
        type=Path,
        help="Optional explicit MinCO profile TSV for gastrooral sample0.",
    )
    ap.add_argument(
        "--minco-sample6",
        type=Path,
        help="Optional explicit MinCO profile TSV for gastrooral sample6.",
    )
    ap.add_argument(
        "--minco-method",
        default="minco_current_universal_gastrooral_source_abundance",
        help="Method label for MinCO rows.",
    )
    ap.add_argument(
        "--output-prefix",
        default="hmp_gastrooral_r232_source_abundance",
        help="Prefix for result TSVs under this experiment's results directory.",
    )
    ap.add_argument(
        "--sample",
        action="append",
        default=[],
        metavar="ID,LABEL,MINCO,UNIQUE,SPLIT,SYLPH",
        help=(
            "Optional dynamic sample definition. May be repeated. When supplied, "
            "the built-in sample0/sample6 set is replaced. Paths are TSV profile, "
            "unique raw table, split raw table, and Sylph r232 profile."
        ),
    )
    return ap.parse_args()


def sample_paths_with_overrides(args: argparse.Namespace) -> dict[int, dict[str, object]]:
    if args.sample:
        out: dict[int, dict[str, object]] = {}
        for spec in args.sample:
            fields = spec.split(",", 5)
            if len(fields) != 6:
                raise SystemExit(
                    "--sample must have 6 comma-separated fields: "
                    "ID,LABEL,MINCO,UNIQUE,SPLIT,SYLPH"
                )
            sample_text, label, minco, unique, split, sylph = fields
            sample_id = int(sample_text)
            if sample_id in out:
                raise SystemExit(f"duplicate --sample ID: {sample_id}")
            out[sample_id] = {
                "label": label,
                "minco": Path(minco),
                "unique": Path(unique),
                "split": Path(split),
                "sylph_r232": Path(sylph),
            }
        return out

    out = {sample: dict(paths) for sample, paths in SAMPLES.items()}
    if args.minco_sample0:
        out[0]["minco"] = args.minco_sample0
    if args.minco_sample6:
        out[6]["minco"] = args.minco_sample6
    return out


def main() -> int:
    args = parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    sample_paths = sample_paths_with_overrides(args)
    required = [TRUTH_ROOT / "genome_to_id.tsv", TAXMAP]
    for sample_id, paths in sample_paths.items():
        required.extend(
            [
                TRUTH_ROOT / f"abundance{sample_id}.tsv",
                Path(paths["minco"]),
                Path(paths["unique"]),
                Path(paths["split"]),
                Path(paths["sylph_r232"]),
            ]
        )
    missing = [path for path in required if not Path(path).exists()]
    if missing:
        raise SystemExit("missing required HMP gastrooral source-GTDB inputs:\n" + "\n".join(map(str, missing)))

    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    source_accessions = load_source_accessions(TRUTH_ROOT / "genome_to_id.tsv")
    taxmap = parse_species_taxmap(TAXMAP)

    score_rows = []
    quality_rows = []
    truth_rows = []
    source_rows = []
    for sample_id, paths in sample_paths.items():
        sample_truth, quality, sample_sources = build_source_truth(
            sample_id,
            str(paths["label"]),
            source_accessions,
            by_accession,
            by_core,
        )
        best_ref_species, best_ref_diag = best_raw_ref_species_by_taxid(
            {"unique": paths["unique"], "split": paths["split"]},
            taxmap,
            by_accession,
            by_core,
        )
        quality.update(best_ref_diag)
        quality.update(
            {
                "minco_reference_release": "gtdb_r232",
                "minco_profile": str(paths["minco"]),
                "sylph_reference_release": "gtdb_r232_chunked",
                "sylph_profile": str(paths["sylph_r232"]),
                "same_release_ready": True,
                "release_grade_candidate": True,
            }
        )
        quality_rows.append(quality)
        truth_rows.extend(sample_truth.to_dict(orient="records"))
        source_rows.extend(sample_sources.to_dict(orient="records"))

        minco_pred, minco_extra = load_minco_predictions(
            paths["minco"],
            best_ref_species,
            by_accession,
            by_core,
        )
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(
            paths["sylph_r232"],
            by_accession,
            by_core,
        )
        sylph_extra.update(
            {
                "sylph_reference_release": "gtdb_r232_chunked",
                "release_grade_candidate": True,
            }
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                args.minco_method,
                minco_pred,
                sample_truth,
                minco_extra,
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_gtdb_r232_gastrooral_source_abundance",
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

    prefix = args.output_prefix
    truth_df.to_csv(RESULTS / f"{prefix}_truth.tsv", sep="\t", index=False)
    quality_df.to_csv(RESULTS / f"{prefix}_quality.tsv", sep="\t", index=False)
    source_df.to_csv(RESULTS / f"{prefix}_source_genomes.tsv", sep="\t", index=False)
    score_df.to_csv(RESULTS / f"{prefix}_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / f"{prefix}_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))
    print("\nQUALITY")
    print(quality_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
