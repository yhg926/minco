#!/usr/bin/env python3
"""Score CAMI3 samples3-5 using source-profile GTDB fallback truth.

The preferred CAMI3 extension route is still per-read source mapping. Those
files are not present locally for samples3-5. This diagnostic uses the CAMI
taxonomic profiles instead: strain/source rows provide source genome IDs,
species labels, and abundance percentages. Rows are transferred to GTDB species
only by deterministic metadata rules.

This does not replace source-readmap truth. It asks whether the locally
available source-profile truth is good enough to become an independent
diagnostic or a candidate release route after policy review.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


EXP = Path(__file__).resolve().parent
REPO_ROOT = EXP.parents[2]
RESULTS = EXP / "results"
sys.path.insert(0, str(REPO_ROOT))

from scripts import minco_profile_calibrated as wrapper  # noqa: E402

import audit_cami3_source_readmap_resolver_candidates as resolver  # noqa: E402
import score_cami3_gtdb_source_readmap as source_score  # noqa: E402
import score_cami3_gtdb_taxid_transfer as taxid_score  # noqa: E402


TAX_PROFILE_ROOT = Path("/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles")
RUN = Path("/tmp/cami3_toy_human_gut_20260626/run")
CANDIDATE_ROOT = Path("/tmp/minco_candidate_preset_replay_20260629")
TAXMAP = source_score.TAXMAP
SAMPLES = [3, 4, 5]
RELEASE_THRESHOLD = 95.0
METHOD_MINCO_BASE = "minco_candidate_source_profile"
METHOD_MINCO_REFINED = "minco_candidate_refined_allocator_source_profile"
METHOD_SYLPH = "sylph_source_profile"
ALLOCATOR_SWITCH = wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230

OUT_TRUTH = RESULTS / "cami3_source_profile_binomial_extension_truth.tsv"
OUT_SOURCE_ROWS = RESULTS / "cami3_source_profile_binomial_extension_source_rows.tsv"
OUT_QUALITY = RESULTS / "cami3_source_profile_binomial_extension_quality.tsv"
OUT_SCORES = RESULTS / "cami3_source_profile_binomial_extension_scores.tsv"
OUT_SUMMARY = RESULTS / "cami3_source_profile_binomial_extension_summary.tsv"
OUT_DELTA = RESULTS / "cami3_source_profile_binomial_extension_delta.tsv"
OUT_AUDIT = RESULTS / "cami3_source_profile_binomial_extension_audit.tsv"
OUT_MD = EXP / "CAMI3_SOURCE_PROFILE_BINOMIAL_EXTENSION.md"


def finite(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def source_profile_path(sample: int) -> Path:
    return TAX_PROFILE_ROOT / f"taxonomic_profile_{sample}.txt"


def candidate_profile_path(sample: int) -> Path:
    return CANDIDATE_ROOT / (
        f"cami3_toy_human_gut_gtdb_source_readmap.sample{sample}.candidate_preset.tsv"
    )


def raw_tables(sample: int) -> dict[str, Path]:
    return {
        "unique": RUN / f"sample{sample}_autoexact/minco.best_diff_unique.unfiltered.tsv",
        "split": RUN / f"sample{sample}_autoexact/minco.best_diff_split.unfiltered.tsv",
    }


def sylph_profile_path(sample: int) -> Path:
    return RUN / f"sylph_sample{sample}/profile.tsv"


def species_taxid_for_strain(row: Mapping[str, str]) -> str:
    path = str(row.get("TAXPATH", ""))
    parts = [part for part in path.split("|") if part]
    return parts[-2] if len(parts) >= 2 else str(row.get("TAXID", ""))


def species_label_for_strain(row: Mapping[str, str]) -> str:
    path = str(row.get("TAXPATHSN", ""))
    parts = [part for part in path.split("|") if part]
    return parts[-2] if len(parts) >= 2 else ""


def in_profile_scope(row: Mapping[str, str]) -> bool:
    path = str(row.get("TAXPATHSN", ""))
    return "Bacteria" in path or "Archaea" in path


def read_source_profile_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    header: list[str] | None = None
    with path.open() as handle:
        for line in handle:
            if line.startswith("@@"):
                header = line[2:].rstrip("\n").split("\t")
                continue
            if line.startswith("@") or not line.strip() or header is None:
                continue
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields + [""] * max(0, len(header) - len(fields))))
            if row.get("RANK") != "strain":
                continue
            source_id = str(row.get("_CAMI_genomeID", "")).strip()
            if not source_id:
                continue
            pct = finite(row.get("PERCENTAGE"))
            if pct <= 0.0:
                continue
            rows.append(
                {
                    "source_genome_id": source_id,
                    "source_profile_taxid": row.get("TAXID", ""),
                    "species_taxid": species_taxid_for_strain(row),
                    "species_label": species_label_for_strain(row),
                    "abundance_pct_all": pct,
                    "profile_scope": "gtdb_profile_scope" if in_profile_scope(row) else "out_of_scope",
                    "taxpathsn": row.get("TAXPATHSN", ""),
                }
            )
    return rows


def choose_resolution(
    source_row: Mapping[str, object],
    taxid_map: Mapping[str, set[str]],
    ncbi_name_map: Mapping[str, set[str]],
    gtdb_binomial_map: Mapping[str, set[str]],
) -> tuple[str, str, int, int, int]:
    taxid = str(source_row.get("species_taxid", "")).strip()
    label = source_score.normalize_name(source_row.get("species_label", ""))
    taxid_candidates = taxid_map.get(taxid, set())
    name_candidates = ncbi_name_map.get(label, set())
    binomial_candidates = gtdb_binomial_map.get(label, set())
    unique_taxid = resolver.one_or_empty(set(taxid_candidates))
    unique_name = resolver.one_or_empty(set(name_candidates))
    unique_binomial = resolver.one_or_empty(set(binomial_candidates))
    if unique_taxid:
        return unique_taxid, "unique_species_taxid", len(taxid_candidates), len(name_candidates), len(binomial_candidates)
    if unique_name:
        return unique_name, "unique_ncbi_species_name", len(taxid_candidates), len(name_candidates), len(binomial_candidates)
    if unique_binomial:
        return unique_binomial, "unique_gtdb_binomial", len(taxid_candidates), len(name_candidates), len(binomial_candidates)
    return "", "", len(taxid_candidates), len(name_candidates), len(binomial_candidates)


def build_truth() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    taxid_map, ncbi_name_map, gtdb_binomial_map = resolver.build_candidate_maps()
    source_rows: list[dict[str, object]] = []
    truth_acc: dict[int, Counter[str]] = defaultdict(Counter)
    truth_sources: dict[tuple[int, str], set[str]] = defaultdict(set)
    quality_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        rows = read_source_profile_rows(source_profile_path(sample))
        method_counts: Counter[str] = Counter()
        scope_total = 0.0
        mapped_total = 0.0
        ambiguous_or_unmapped = 0.0
        for row in rows:
            in_scope = row["profile_scope"] == "gtdb_profile_scope"
            species = ""
            method = "out_of_scope"
            taxid_n = name_n = binomial_n = 0
            if in_scope:
                scope_total += float(row["abundance_pct_all"])
                species, method, taxid_n, name_n, binomial_n = choose_resolution(
                    row,
                    taxid_map,
                    ncbi_name_map,
                    gtdb_binomial_map,
                )
                if species:
                    mapped_total += float(row["abundance_pct_all"])
                    truth_acc[sample][species] += float(row["abundance_pct_all"])
                    truth_sources[(sample, species)].add(str(row["source_genome_id"]))
                else:
                    method = "unresolved"
                    ambiguous_or_unmapped += float(row["abundance_pct_all"])
            method_counts[method] += 1
            source_rows.append(
                {
                    "sample": sample,
                    **row,
                    "recommended_gtdb_species": species,
                    "recommended_resolution_rule": method,
                    "taxid_gtdb_species_count": taxid_n,
                    "ncbi_name_gtdb_species_count": name_n,
                    "gtdb_binomial_species_count": binomial_n,
                }
            )
        quality_rows.append(
            {
                "sample": sample,
                "source_profile_rows": len(rows),
                "profile_scope_rows": sum(1 for row in rows if row["profile_scope"] == "gtdb_profile_scope"),
                "profile_scope_abundance_pct_all": scope_total,
                "mapped_profile_scope_abundance_pct_all": mapped_total,
                "unresolved_profile_scope_abundance_pct_all": ambiguous_or_unmapped,
                "mapped_profile_scope_pct": 100.0 * mapped_total / scope_total if scope_total else 0.0,
                "truth_gtdb_species": len(truth_acc[sample]),
                "unique_species_taxid_rows": method_counts.get("unique_species_taxid", 0),
                "unique_ncbi_species_name_rows": method_counts.get("unique_ncbi_species_name", 0),
                "unique_gtdb_binomial_rows": method_counts.get("unique_gtdb_binomial", 0),
                "unresolved_rows": method_counts.get("unresolved", 0),
                "out_of_scope_rows": method_counts.get("out_of_scope", 0),
                "release_threshold_pct": RELEASE_THRESHOLD,
                "profile_scope_release_ready": str(
                    scope_total > 0.0 and 100.0 * mapped_total / scope_total >= RELEASE_THRESHOLD
                ).lower(),
            }
        )

    truth_rows: list[dict[str, object]] = []
    for sample in sorted(truth_acc):
        total = sum(truth_acc[sample].values())
        for species, abundance_pct_all in sorted(truth_acc[sample].items()):
            truth_rows.append(
                {
                    "sample": sample,
                    "gtdb_species": species,
                    "abundance_pct_all": abundance_pct_all,
                    "truth_abundance": abundance_pct_all / total if total else 0.0,
                    "source_genomes": len(truth_sources[(sample, species)]),
                    "source_genome_ids": ",".join(sorted(truth_sources[(sample, species)])),
                }
            )
    return pd.DataFrame(truth_rows), pd.DataFrame(source_rows), pd.DataFrame(quality_rows)


def minco_prediction_from_profile(
    profile: pd.DataFrame,
    abundance_raw: np.ndarray | None,
    taxid_to_species: Mapping[str, str],
    name_to_species: Mapping[str, str],
    best_ref_species_by_taxid: Mapping[str, str],
    by_accession: Mapping[str, Mapping[str, str]],
    by_core: Mapping[str, list[Mapping[str, str]]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    call = profile["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"})
    selected = profile.loc[call].copy()
    if abundance_raw is None:
        abundance_values = pd.to_numeric(
            selected.get("calibrated_abundance", 0.0),
            errors="coerce",
        ).fillna(0.0).to_numpy(dtype=float)
    else:
        abundance_values = np.asarray(abundance_raw, dtype=float)[call.to_numpy()]
    rows = []
    unmapped = name_mapped = taxid_mapped = ref_mapped = 0
    for idx, row in enumerate(selected.itertuples(index=False)):
        taxid = str(getattr(row, "taxid", ""))
        species_name = source_score.normalize_name(getattr(row, "species_name", ""))
        gtdb_name = ""
        for col in ["s_best_accession", "u_best_accession", "candidate_surface_accession"]:
            if hasattr(row, col):
                gtdb_name = source_score.gtdb_from_accession(
                    getattr(row, col, ""),
                    by_accession,
                    by_core,
                )
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
        abundance = float(abundance_values[idx]) if idx < len(abundance_values) else 0.0
        ani = finite(getattr(row, "s_Ref_zip_aaf_ani_max", 0.0))
        rows.append((gtdb_name, abundance, ani))
    return taxid_score.collapse_prediction(rows), {
        "pred_rows_called": int(len(selected)),
        "pred_rows_unmapped": unmapped,
        "pred_rows_ref_mapped": ref_mapped,
        "pred_rows_taxid_mapped": taxid_mapped,
        "pred_rows_name_mapped": name_mapped,
    }


def score_profiles(truth_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    _wgs_to_species, taxid_to_species, name_to_species, _diag = source_score.build_transfer_maps()
    by_accession, by_core = source_score.truth.load_gtdb_metadata(source_score.truth.GTDB_METADATA)
    taxmap = source_score.parse_species_taxmap(TAXMAP)
    score_rows: list[dict[str, object]] = []
    delta_rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        best_ref_species, best_ref_diag = source_score.best_raw_ref_species_by_taxid(
            raw_tables(sample),
            taxmap,
            by_accession,
            by_core,
        )
        truth_one = truth_df.loc[truth_df["sample"].astype(str).eq(str(sample))].copy()
        profile = pd.read_csv(candidate_profile_path(sample), sep="\t")
        call = profile["calibrated_call"].astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()
        raw = pd.to_numeric(
            profile.get("calibrated_abundance_raw", profile.get("calibrated_abundance", 0.0)),
            errors="coerce",
        ).fillna(0.0).to_numpy(dtype=float)
        allocator_taxmap = wrapper.parse_species_taxmap(TAXMAP)
        adjusted, allocator_details = wrapper.guarded_feature_allocator_abundance_raw(
            profile,
            call,
            raw,
            ALLOCATOR_SWITCH,
            allocator_taxmap,
        )
        base_pred, base_meta = minco_prediction_from_profile(
            profile,
            None,
            taxid_to_species,
            name_to_species,
            best_ref_species,
            by_accession,
            by_core,
        )
        refined_pred, refined_meta = minco_prediction_from_profile(
            profile,
            adjusted,
            taxid_to_species,
            name_to_species,
            best_ref_species,
            by_accession,
            by_core,
        )
        sylph_pred, sylph_meta = taxid_score.load_sylph_predictions(
            sylph_profile_path(sample),
            by_accession,
            by_core,
        )
        base_score = taxid_score.score_prediction(
            sample,
            METHOD_MINCO_BASE,
            base_pred,
            truth_one,
            {**best_ref_diag, **base_meta},
        )
        refined_score = taxid_score.score_prediction(
            sample,
            METHOD_MINCO_REFINED,
            refined_pred,
            truth_one,
            {**best_ref_diag, **refined_meta, **allocator_details},
        )
        sylph_score = taxid_score.score_prediction(
            sample,
            METHOD_SYLPH,
            sylph_pred,
            truth_one,
            sylph_meta,
        )
        score_rows.extend([base_score, refined_score, sylph_score])
        delta_rows.append(
            {
                "sample": sample,
                "delta_refined_minus_base_F1": finite(refined_score["F1"]) - finite(base_score["F1"]),
                "delta_refined_minus_base_L1_union_pp": finite(refined_score["L1_union_pp"])
                - finite(base_score["L1_union_pp"]),
                "delta_refined_minus_base_Pearson_union": finite(refined_score["Pearson_union"])
                - finite(base_score["Pearson_union"]),
                "allocator_applied": allocator_details.get("abundance_feature_allocator_applied", False),
                "allocator_guard_passed": allocator_details.get("abundance_feature_allocator_guard_passed", False),
                "allocator_adjusted_rows_n": allocator_details.get("abundance_feature_allocator_adjusted_rows_n", 0),
                "allocator_base_mass_multi_genus_frac": allocator_details.get(
                    "abundance_feature_allocator_base_mass_multi_genus_frac",
                    0.0,
                ),
                "allocator_s_xny_median": allocator_details.get(
                    "abundance_feature_allocator_s_xny_median",
                    0.0,
                ),
            }
        )
    scores = pd.DataFrame(score_rows)
    summary = taxid_score.summarize(scores)
    return scores, summary, pd.DataFrame(delta_rows)


def build_audit(quality: pd.DataFrame, scores: pd.DataFrame, delta: pd.DataFrame) -> list[dict[str, object]]:
    ready = quality.loc[quality["profile_scope_release_ready"].astype(str).str.lower().eq("true")]
    min_mapped = float(quality["mapped_profile_scope_pct"].min()) if not quality.empty else 0.0
    base = scores.loc[scores["method"].eq(METHOD_MINCO_BASE)].set_index("sample")
    refined = scores.loc[scores["method"].eq(METHOD_MINCO_REFINED)].set_index("sample")
    sylph = scores.loc[scores["method"].eq(METHOD_SYLPH)].set_index("sample")
    refined_improved = int((delta["delta_refined_minus_base_L1_union_pp"] < -1e-12).sum())
    refined_worse = int((delta["delta_refined_minus_base_L1_union_pp"] > 1e-12).sum())
    mean_delta_l1 = float(delta["delta_refined_minus_base_L1_union_pp"].mean()) if not delta.empty else 0.0
    minco_f1_wins = sum(
        1 for sample in SAMPLES if finite(base.loc[sample, "F1"]) > finite(sylph.loc[sample, "F1"])
    )
    minco_l1_wins = sum(
        1
        for sample in SAMPLES
        if finite(base.loc[sample, "L1_union_pp"]) < finite(sylph.loc[sample, "L1_union_pp"])
    )
    return [
        {
            "metric": "source_profile_truth_scope",
            "value": (
                f"samples={','.join(map(str, SAMPLES))};"
                f"min_mapped_profile_scope_pct={min_mapped:.6f};"
                f"ready_samples={','.join(map(str, ready['sample'].astype(int).tolist()))}"
            ),
            "evidence": str(OUT_QUALITY.relative_to(EXP)),
            "decision": "source_profile_truth_mapping_passes_threshold"
            if len(ready) == len(SAMPLES)
            else "source_profile_truth_mapping_incomplete",
        },
        {
            "metric": "baseline_minco_vs_sylph",
            "value": f"minco_F1_wins={minco_f1_wins}/3;minco_L1_wins={minco_l1_wins}/3",
            "evidence": str(OUT_SCORES.relative_to(EXP)),
            "decision": "diagnostic_comparison_only",
        },
        {
            "metric": "refined_allocator_effect",
            "value": (
                f"improved={refined_improved};worsened={refined_worse};"
                f"mean_L1_delta_pp={mean_delta_l1:.9f}"
            ),
            "evidence": str(OUT_DELTA.relative_to(EXP)),
            "decision": "diagnostic_supports_refined_allocator"
            if refined_worse == 0 and refined_improved > 0
            else "do_not_promote_from_source_profile_diagnostic",
        },
        {
            "metric": "promotion_decision",
            "value": "source_profile_extension_diagnostic_not_source_readmap_release_truth",
            "evidence": str(OUT_AUDIT.relative_to(EXP)),
            "decision": "keep_refined_allocator_opt_in",
        },
    ]


def write_markdown(quality: pd.DataFrame, summary: pd.DataFrame, audit: list[dict[str, object]]) -> None:
    audit_by_metric = {str(row["metric"]): row for row in audit}
    lines = [
        "# CAMI3 Source-Profile Binomial Extension",
        "",
        "Date: 2026-06-29",
        "",
        "This diagnostic scores CAMI3 samples3-5 using local taxonomic-profile",
        "strain/source rows with deterministic GTDB transfer. It does not replace",
        "the missing per-read source-readmap truth files.",
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
            "## Truth Quality",
            "",
            "| Sample | In-scope abundance % | Mapped in-scope % | GTDB species | Ready |",
            "|---:|---:|---:|---:|---|",
        ]
    )
    for row in quality.to_dict("records"):
        lines.append(
            "| {sample} | {scope:.6f} | {mapped:.6f} | {species} | {ready} |".format(
                sample=int(row["sample"]),
                scope=float(row["profile_scope_abundance_pct_all"]),
                mapped=float(row["mapped_profile_scope_pct"]),
                species=int(row["truth_gtdb_species"]),
                ready=row["profile_scope_release_ready"],
            )
        )
    lines.extend(
        [
            "",
            "## Score Summary",
            "",
            "| Method | Samples | Mean F1 | Mean L1 union pp | Mean Pearson union |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in summary.to_dict("records"):
        lines.append(
            "| {method} | {samples} | {f1:.6f} | {l1:.6f} | {pearson:.6f} |".format(
                method=row["method"],
                samples=row["samples"],
                f1=float(row["mean_F1"]),
                l1=float(row["mean_L1_union_pp"]),
                pearson=float(row["mean_Pearson_union"]),
            )
        )
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Keep this as diagnostic evidence until the truth policy is reviewed.",
            "- The main release route remains recovering per-read source mapping for",
            "  CAMI3 samples3-5.",
            "- The refined allocator remains opt-in by default.",
            "",
            "## Outputs",
            "",
            f"- `{OUT_TRUTH.relative_to(EXP)}`",
            f"- `{OUT_SOURCE_ROWS.relative_to(EXP)}`",
            f"- `{OUT_QUALITY.relative_to(EXP)}`",
            f"- `{OUT_SCORES.relative_to(EXP)}`",
            f"- `{OUT_SUMMARY.relative_to(EXP)}`",
            f"- `{OUT_DELTA.relative_to(EXP)}`",
            f"- `{OUT_AUDIT.relative_to(EXP)}`",
            "",
            f"Release decision remains `{audit_by_metric['promotion_decision']['decision']}`.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    required = [TAXMAP]
    for sample in SAMPLES:
        required.extend(
            [
                source_profile_path(sample),
                candidate_profile_path(sample),
                raw_tables(sample)["unique"],
                raw_tables(sample)["split"],
                sylph_profile_path(sample),
            ]
        )
    missing = [str(path) for path in required if not Path(path).exists()]
    if missing:
        raise SystemExit("missing source-profile extension inputs:\n" + "\n".join(missing))

    truth_df, source_rows, quality = build_truth()
    scores, summary, delta = score_profiles(truth_df)
    audit = build_audit(quality, scores, delta)
    truth_df.to_csv(OUT_TRUTH, sep="\t", index=False)
    source_rows.to_csv(OUT_SOURCE_ROWS, sep="\t", index=False)
    quality.to_csv(OUT_QUALITY, sep="\t", index=False)
    scores.to_csv(OUT_SCORES, sep="\t", index=False)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    delta.to_csv(OUT_DELTA, sep="\t", index=False)
    write_tsv(OUT_AUDIT, audit, ["metric", "value", "evidence", "decision"])
    write_markdown(quality, summary, audit)
    print("\n".join(f"{row['metric']}\t{row['value']}\t{row['decision']}" for row in audit))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
