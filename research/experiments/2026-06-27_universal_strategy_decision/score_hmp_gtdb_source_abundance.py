#!/usr/bin/env python3
"""Score CAMI II HMP airskin samples in a source-accession GTDB namespace.

The CAMISIM setup provides per-source-genome abundance files and
``genome_to_id.tsv`` with source assembly accessions. This scorer maps those
source accessions to GTDB species through GTDB r232 metadata, then maps MinCO
and Sylph predictions to the same GTDB species labels by best reference
accession. With ``--sylph-source r232`` and explicit refreshed MinCO profiles,
the scorer produces the same-release HMP airskin release-grade comparison used
by the universal-strategy decision note; cached r226 Sylph profiles remain
available for historical diagnostics.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping

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


TRUTH_ROOT = Path("/tmp/cami2_hmp_airskin_20260625/truth")
RUN = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run")
R232_RUN = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232")
AIRSKIN_RUN = Path("/tmp/cami2_hmp_airskin_20260625/run")
AIRSKIN_R232_RUN = Path("/tmp/cami2_hmp_airskin_20260625/run_sylph_r232")
AIRSKIN28_CURRENT = Path("/tmp/minco_current_code_hmp_airskin28_20260627")
HMP_REFRESH_CURRENT = Path("/tmp/minco_current_code_hmp_refresh_20260627")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")

SAMPLES = {
    6: {
        "label": "airskin_urogenital6",
        "minco": HMP_REFRESH_CURRENT / "minco_sample6_current_exact.tsv",
        "unique": RUN / "minco_sample6_unique_zip_unfiltered.tsv",
        "split": RUN / "minco_sample6_split_zip_unfiltered.tsv",
        "sylph_r226": RUN / "sylph_sample6/profile.tsv",
        "sylph_r232": R232_RUN / "sylph_sample6/profile.tsv",
    },
    11: {
        "label": "airskin_airways11",
        "minco": HMP_REFRESH_CURRENT / "minco_sample11_current_exact.tsv",
        "unique": RUN / "minco_sample11_unique_zip_unfiltered.tsv",
        "split": RUN / "minco_sample11_split_zip_unfiltered.tsv",
        "sylph_r226": RUN / "sylph_sample11/profile.tsv",
        "sylph_r232": R232_RUN / "sylph_sample11/profile.tsv",
    },
    28: {
        "label": "airskin_sample28",
        "minco": AIRSKIN28_CURRENT / "minco_sample28_current_default.tsv",
        "unique": AIRSKIN28_CURRENT / "work/minco.best_diff_unique.unfiltered.tsv",
        "split": AIRSKIN28_CURRENT / "work/minco.best_diff_split.unfiltered.tsv",
        "sylph_r226": AIRSKIN_RUN / "sylph_sample28/profile.tsv",
        "sylph_r232": AIRSKIN_R232_RUN / "sylph_sample28/profile.tsv",
    },
}

ACC_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)")


def dynamic_sample_root(sample_id: int) -> Path:
    return Path(f"/tmp/minco_current_code_hmp_airskin{sample_id}_20260627")


def ensure_sample_record(sample_id: int) -> None:
    if sample_id in SAMPLES:
        return
    root = dynamic_sample_root(sample_id)
    SAMPLES[sample_id] = {
        "label": f"airskin_sample{sample_id}",
        "minco": root / f"minco_sample{sample_id}_current_default.tsv",
        "unique": root / "work/minco.best_diff_unique.unfiltered.tsv",
        "split": root / "work/minco.best_diff_split.unfiltered.tsv",
        "sylph_r226": AIRSKIN_RUN / f"sylph_sample{sample_id}/profile.tsv",
        "sylph_r232": R232_RUN / f"sylph_sample{sample_id}/profile.tsv",
    }


def extract_accession(value: object) -> str:
    match = ACC_RE.search(str(value or ""))
    return match.group(1) if match else ""


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").replace([math.inf, -math.inf], pd.NA).fillna(0.0)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Score HMP airskin MinCO/Sylph predictions in GTDB r232 source-abundance truth space.",
    )
    ap.add_argument(
        "--sylph-source",
        choices=["auto", "r226", "r232"],
        default="auto",
        help="Sylph profile source. auto uses r232 profiles when all required files exist, otherwise cached r226.",
    )
    ap.add_argument(
        "--r232-run",
        type=Path,
        default=R232_RUN,
        help="Directory containing sylph_sample6/profile.tsv and sylph_sample11/profile.tsv for same-release r232 scoring.",
    )
    ap.add_argument(
        "--samples",
        default="6,11",
        help="Comma-separated HMP airskin sample IDs to score. Default keeps the original release panel: 6,11.",
    )
    ap.add_argument(
        "--minco-profile-dir",
        type=Path,
        help="Optional directory with minco_sample6.tsv and minco_sample11.tsv refreshed profile outputs.",
    )
    ap.add_argument(
        "--minco-sample6",
        type=Path,
        help="Optional explicit MinCO profile TSV for sample6.",
    )
    ap.add_argument(
        "--minco-sample11",
        type=Path,
        help="Optional explicit MinCO profile TSV for sample11.",
    )
    ap.add_argument(
        "--minco-sample28",
        type=Path,
        help="Optional explicit MinCO profile TSV for sample28.",
    )
    ap.add_argument(
        "--minco-sample",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help="Optional generic MinCO profile override, repeatable for dynamic samples.",
    )
    ap.add_argument(
        "--sylph-sample28",
        type=Path,
        help="Optional explicit Sylph profile TSV for sample28.",
    )
    ap.add_argument(
        "--sylph-sample",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help="Optional generic Sylph profile override, repeatable for dynamic samples.",
    )
    ap.add_argument(
        "--minco-unique-sample",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help="Optional generic MinCO best-diff-unique raw table override.",
    )
    ap.add_argument(
        "--minco-split-sample",
        action="append",
        default=[],
        metavar="SAMPLE=PATH",
        help="Optional generic MinCO best-diff-split raw table override.",
    )
    ap.add_argument(
        "--minco-method",
        default="minco_universal_strategy_gtdb_source_abundance",
        help="Method label for the MinCO rows.",
    )
    ap.add_argument(
        "--output-prefix",
        default="hmp_gtdb_source_abundance",
        help="Prefix for result TSVs under the decision-note results directory.",
    )
    return ap.parse_args()


def selected_sample_ids(args: argparse.Namespace) -> list[int]:
    sample_ids = []
    for raw in str(args.samples).split(","):
        raw = raw.strip()
        if not raw:
            continue
        sample_id = int(raw)
        ensure_sample_record(sample_id)
        sample_ids.append(sample_id)
    if not sample_ids:
        raise SystemExit("--samples selected no sample IDs")
    return sample_ids


def parse_sample_path_overrides(values: Iterable[str]) -> dict[int, Path]:
    overrides: dict[int, Path] = {}
    for value in values:
        if "=" not in value:
            raise SystemExit(f"Expected SAMPLE=PATH override, got: {value}")
        raw_sample, raw_path = value.split("=", 1)
        raw_sample = raw_sample.strip()
        raw_path = raw_path.strip()
        if not raw_sample or not raw_path:
            raise SystemExit(f"Expected non-empty SAMPLE=PATH override, got: {value}")
        overrides[int(raw_sample)] = Path(raw_path)
    return overrides


def select_sylph_profiles(args: argparse.Namespace, sample_ids: Iterable[int]) -> tuple[dict[int, Path], str, bool]:
    r232_profiles = {
        sample_id: args.r232_run / f"sylph_sample{sample_id}/profile.tsv"
        for sample_id in sample_ids
    }
    if args.sylph_sample28:
        r232_profiles[28] = args.sylph_sample28
    r232_profiles.update(parse_sample_path_overrides(args.sylph_sample))
    r232_complete = all(path.exists() for path in r232_profiles.values())
    if args.sylph_source == "r232":
        return r232_profiles, "gtdb_r232", True
    if args.sylph_source == "auto" and r232_complete:
        return r232_profiles, "gtdb_r232", True
    return (
        {sample_id: Path(SAMPLES[sample_id]["sylph_r226"]) for sample_id in sample_ids},
        "gtdb_r226_cached",
        False,
    )


def select_minco_profiles(args: argparse.Namespace, sample_ids: Iterable[int]) -> dict[int, Path]:
    out = {sample_id: Path(SAMPLES[sample_id]["minco"]) for sample_id in sample_ids}
    if args.minco_profile_dir:
        for sample_id in sample_ids:
            out[sample_id] = args.minco_profile_dir / f"minco_sample{sample_id}.tsv"
    if args.minco_sample6:
        out[6] = args.minco_sample6
    if args.minco_sample11:
        out[11] = args.minco_sample11
    if args.minco_sample28:
        out[28] = args.minco_sample28
    out.update(parse_sample_path_overrides(args.minco_sample))
    return out


def select_raw_tables(args: argparse.Namespace, sample_ids: Iterable[int]) -> dict[int, dict[str, Path]]:
    out = {
        sample_id: {
            "unique": Path(SAMPLES[sample_id]["unique"]),
            "split": Path(SAMPLES[sample_id]["split"]),
        }
        for sample_id in sample_ids
    }
    for sample_id, path in parse_sample_path_overrides(args.minco_unique_sample).items():
        ensure_sample_record(sample_id)
        out.setdefault(sample_id, {})["unique"] = path
    for sample_id, path in parse_sample_path_overrides(args.minco_split_sample).items():
        ensure_sample_record(sample_id)
        out.setdefault(sample_id, {})["split"] = path
    return out


def load_source_accessions(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
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
    source_accessions: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame]:
    abundance_path = TRUTH_ROOT / f"abundance{sample_id}.tsv"
    counts: Counter[str] = Counter()
    source_rows = []
    total_mass = 0.0
    mapped_mass = 0.0
    method_mass: Counter[str] = Counter()
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

    truth_rows = []
    for species, value in sorted(counts.items()):
        truth_rows.append(
            {
                "sample": sample_id,
                "gtdb_species": species,
                "truth_raw_abundance": value,
                "truth_abundance": value / mapped_mass if mapped_mass > 0.0 else 0.0,
            }
        )
    quality = {
        "sample": sample_id,
        "sample_label": SAMPLES[sample_id]["label"],
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
    }


def main() -> int:
    args = parse_args()
    sample_ids = selected_sample_ids(args)
    sylph_profiles, sylph_release, same_release_ready = select_sylph_profiles(args, sample_ids)
    minco_profiles = select_minco_profiles(args, sample_ids)
    raw_tables = select_raw_tables(args, sample_ids)
    RESULTS.mkdir(parents=True, exist_ok=True)
    required = [TRUTH_ROOT / "genome_to_id.tsv", TAXMAP]
    for sample_id in sample_ids:
        required.extend(
            [
                TRUTH_ROOT / f"abundance{sample_id}.tsv",
                minco_profiles[sample_id],
                raw_tables[sample_id]["unique"],
                raw_tables[sample_id]["split"],
                sylph_profiles[sample_id],
            ]
        )
    missing = [path for path in required if not Path(path).exists()]
    if missing:
        raise SystemExit("missing required HMP source-GTDB inputs:\n" + "\n".join(map(str, missing)))

    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    source_accessions = load_source_accessions(TRUTH_ROOT / "genome_to_id.tsv")
    taxmap = parse_species_taxmap(TAXMAP)

    score_rows = []
    quality_rows = []
    truth_rows = []
    source_rows = []
    for sample_id in sample_ids:
        sample_truth, quality, sample_sources = build_source_truth(
            sample_id,
            source_accessions,
            by_accession,
            by_core,
        )
        best_ref_species, best_ref_diag = best_raw_ref_species_by_taxid(
            raw_tables[sample_id],
            taxmap,
            by_accession,
            by_core,
        )
        quality.update(best_ref_diag)
        quality.update(
            {
                "minco_reference_release": "gtdb_r232",
                "sylph_reference_release": sylph_release,
                "sylph_profile": str(sylph_profiles[sample_id]),
                "same_release_ready": same_release_ready,
                "release_grade_candidate": same_release_ready,
            }
        )
        quality_rows.append(quality)
        truth_rows.extend(sample_truth.to_dict(orient="records"))
        source_rows.extend(sample_sources.to_dict(orient="records"))

        minco_pred, minco_extra = load_minco_predictions(
            minco_profiles[sample_id],
            best_ref_species,
            by_accession,
            by_core,
        )
        minco_extra["minco_profile"] = str(minco_profiles[sample_id])
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(
            sylph_profiles[sample_id],
            by_accession,
            by_core,
        )
        sylph_extra.update(
            {
                "sylph_reference_release": sylph_release,
                "release_grade_candidate": same_release_ready,
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
                "sylph_gtdb_source_abundance",
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
