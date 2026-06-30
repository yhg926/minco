#!/usr/bin/env python3
"""Score CAMI II marine profiles after conservative GTDB-species transfer.

The CAMI II marine gold profile is NCBI-taxid based and contains large
non-bacterial/non-archaeal components. This scorer transfers only positive
Bacteria/Archaea species rows to GTDB species when the NCBI taxid maps to one
GTDB species in r232 metadata. Excluded and unmapped mass are reported so the
panel is not mistaken for clean release-grade truth.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

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


GOLD = Path("/tmp/gs_marine_short.profile")
TAXMAP = Path("/tmp/cami2_marine_samples3_5_20260625/run/gtdb232_refseqvirus_species_taxmap.tsv")
SAMPLE_PATHS = {
    0: {
        "minco": Path(
            "/tmp/minco_fair_l1_edge_adaptive_20260624/marine_current_best/"
            "minco_s1000_unique_zipaaf_sample0.tsv"
        ),
        "sylph": Path("/tmp/minco_fair_l1_edge_adaptive_20260624/marine_sylph/profile.tsv"),
    },
    3: {
        "minco": Path(
            "/tmp/cami2_marine_samples3_5_20260625/run/"
            "minco_sample3_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
        ),
        "sylph": Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample3/profile.tsv"),
    },
    4: {
        "minco": Path(
            "/tmp/cami2_marine_samples3_5_20260625/run/"
            "minco_sample4_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
        ),
        "sylph": Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample4/profile.tsv"),
    },
    5: {
        "minco": Path(
            "/tmp/cami2_marine_samples3_5_20260625/run/"
            "minco_sample5_s1000_unique_zipaaf_f0.05_n0.94_t10.tsv"
        ),
        "sylph": Path("/tmp/cami2_marine_samples3_5_20260625/run/sylph_sample5/profile.tsv"),
    },
}


DROP_NAMES = {"unidentified", "unidentified plasmid", "unidentified virus"}


def sample_name(sample_id: int) -> str:
    return f"marmgCAMI2_short_read_sample_{sample_id}"


def numeric_value(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return 0.0
    return out if math.isfinite(out) else 0.0


def iter_gold_species_rows(path: Path, sample_id: int):
    target = sample_name(sample_id)
    active = False
    seen = False
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("@SampleID:"):
                sid = line.split(":", 1)[1].strip()
                if seen and sid != target:
                    break
                active = sid == target
                seen = seen or active
                continue
            if not active or not line or line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[0] == "TAXID" or fields[1] != "species":
                continue
            taxid, _rank, taxpath, taxpathsn, pct = fields[:5]
            abundance = numeric_value(pct)
            if abundance <= 0.0:
                continue
            name = taxpathsn.split("|")[-1].strip() if taxpathsn else ""
            yield {
                "taxid": str(taxid),
                "taxpath": taxpath,
                "taxpathsn": taxpathsn,
                "name": name,
                "abundance_pct_all": abundance,
            }


def read_marine_ba_truth(sample_id: int) -> tuple[pd.DataFrame, dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_species = 0.0
    ba_species = 0.0
    dropped_named = 0.0
    virus_species = 0.0
    plasmid_or_unknown = 0.0
    other_species = 0.0

    for row in iter_gold_species_rows(GOLD, sample_id):
        abundance = float(row["abundance_pct_all"])
        total_species += abundance
        taxpath = str(row["taxpath"])
        name = str(row["name"]).lower()
        if name in DROP_NAMES:
            dropped_named += abundance
            if "plasmid" in name or "unidentified" in name:
                plasmid_or_unknown += abundance
            continue
        if taxpath.startswith("2|") or taxpath == "2" or taxpath.startswith("2157|") or taxpath == "2157":
            ba_species += abundance
            rows.append(
                {
                    "ncbi_species_taxid": row["taxid"],
                    "truth_name": row["name"],
                    "abundance_pct_all": abundance,
                }
            )
        elif taxpath.startswith("10239|") or taxpath == "10239":
            virus_species += abundance
        else:
            other_species += abundance

    df = pd.DataFrame(rows)
    if not df.empty:
        df = (
            df.groupby("ncbi_species_taxid", as_index=False)
            .agg(truth_name=("truth_name", "first"), abundance_pct_all=("abundance_pct_all", "sum"))
        )
    else:
        df = pd.DataFrame(columns=["ncbi_species_taxid", "truth_name", "abundance_pct_all"])

    scope = {
        "sample": sample_id,
        "gold_species_mass_pct_all": total_species,
        "gold_bacteria_archaea_species_mass_pct_all": ba_species,
        "gold_dropped_named_species_mass_pct_all": dropped_named,
        "gold_virus_species_mass_pct_all": virus_species,
        "gold_plasmid_or_unknown_species_mass_pct_all": plasmid_or_unknown,
        "gold_other_species_mass_pct_all": other_species,
        "gold_bacteria_archaea_species_taxids": int(len(df)),
    }
    return df, scope


def gtdb_from_accession(accession: str, by_accession, by_core) -> tuple[str, str]:
    rec, method, _key = truth.lookup_accession(accession, by_accession, by_core)
    if rec is None:
        return "", method
    return str(rec.get("gtdb_species", "")), method


def collapse_minco_prediction(
    path: Path,
    taxmap,
    taxid_to_gtdb: dict[str, str],
    by_accession,
    by_core,
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = load_minco(path, taxmap, 11.0)
    rows = []
    accession_mapped = 0
    taxid_mapped = 0
    unmapped = 0
    for row in raw.itertuples(index=False):
        species, method = gtdb_from_accession(str(getattr(row, "accession", "")), by_accession, by_core)
        if species:
            accession_mapped += 1
        else:
            species = taxid_to_gtdb.get(str(getattr(row, "taxid", "")), "")
            if species:
                taxid_mapped += 1
        if not species:
            unmapped += 1
            continue
        abundance = numeric_value(getattr(row, "Normalized_abundance_depth", 0.0))
        ani = numeric_value(getattr(row, "Ref_zip_aaf_ani", 0.0)) or numeric_value(getattr(row, "ANI", 0.0))
        rows.append((species, abundance, ani))
    return taxid_score.collapse_prediction(rows), {
        "pred_rows_called": int(len(raw)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_accession_mapped": accession_mapped,
        "pred_rows_taxid_fallback_mapped": taxid_mapped,
    }


def all_gold_sample_ids(path: Path) -> list[int]:
    out = []
    prefix = "marmgCAMI2_short_read_sample_"
    with path.open() as handle:
        for line in handle:
            if not line.startswith("@SampleID:"):
                continue
            name = line.split(":", 1)[1].strip()
            if name.startswith(prefix):
                out.append(int(name[len(prefix) :]))
    return out


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    required = [GOLD, TAXMAP]
    for paths in SAMPLE_PATHS.values():
        required.extend([paths["minco"], paths["sylph"]])
    missing = [path for path in required if not path.exists()]
    if missing:
        raise SystemExit("missing required marine GTDB-transfer inputs:\n" + "\n".join(map(str, missing)))

    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    taxid_to_gtdb, ambiguous_taxids = taxid_score.build_taxid_transfer(by_accession)
    taxmap = parse_species_taxmap(TAXMAP)

    quality_rows = []
    score_rows = []
    truth_rows = []
    for sample_id in all_gold_sample_ids(GOLD):
        ba_truth, scope = read_marine_ba_truth(sample_id)
        gtdb_truth, quality = taxid_score.transfer_truth_to_gtdb(
            sample_id,
            ba_truth,
            taxid_to_gtdb,
            ambiguous_taxids,
        )
        quality.update(scope)
        ba_mass = float(scope["gold_bacteria_archaea_species_mass_pct_all"])
        quality["truth_mass_mapped_pct_bacteria_archaea"] = (
            float(quality["truth_mass_mapped_pct_all"]) / ba_mass * 100.0 if ba_mass else 0.0
        )
        quality["release_grade_threshold_pct_bacteria_archaea"] = 95.0
        quality["release_grade_ready_bacteria_archaea"] = (
            quality["truth_mass_mapped_pct_bacteria_archaea"] >= 95.0
        )
        quality["scored_with_profiles"] = sample_id in SAMPLE_PATHS
        quality_rows.append(quality)

        if sample_id not in SAMPLE_PATHS:
            continue
        truth_rows.extend(gtdb_truth.to_dict(orient="records"))
        paths = SAMPLE_PATHS[sample_id]
        minco_pred, minco_extra = collapse_minco_prediction(
            paths["minco"],
            taxmap,
            taxid_to_gtdb,
            by_accession,
            by_core,
        )
        minco_extra.update({"profile": str(paths["minco"]), "reference_release_note": "cached_mixed_r226_paths_r232_taxmap"})
        sylph_pred, sylph_extra = taxid_score.load_sylph_predictions(paths["sylph"], by_accession, by_core)
        sylph_extra.update({"profile": str(paths["sylph"]), "reference_release_note": "cached_r226_paths_mapped_to_r232_metadata"})
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "minco_marine_s1000_unique_zipaaf_gtdb_transfer",
                minco_pred,
                gtdb_truth,
                minco_extra,
            )
        )
        score_rows.append(
            taxid_score.score_prediction(
                sample_id,
                "sylph_marine_r226_gtdb_transfer",
                sylph_pred,
                gtdb_truth,
                sylph_extra,
            )
        )

    quality_df = pd.DataFrame(quality_rows)
    score_df = pd.DataFrame(score_rows)
    summary_df = taxid_score.summarize(score_df)
    pd.DataFrame(truth_rows).to_csv(RESULTS / "marine_gtdb_taxid_transfer_truth.tsv", sep="\t", index=False)
    quality_df.to_csv(RESULTS / "marine_gtdb_taxid_transfer_quality.tsv", sep="\t", index=False)
    score_df.to_csv(RESULTS / "marine_gtdb_taxid_transfer_scores.tsv", sep="\t", index=False)
    summary_df.to_csv(RESULTS / "marine_gtdb_taxid_transfer_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
