#!/usr/bin/env python3
"""Write a compact manifest for the current MinCO default strategy.

This is a release-facing recall artifact. It reads the already generated
decision TSVs and records the selected default, compatibility preset, evidence
size, score deltas, and claim boundary without rerunning any profiles.
"""

from __future__ import annotations

import csv
from pathlib import Path


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"
OUT = RESULTS / "current_default_strategy_manifest.tsv"
OUT_MD = EXP / "CURRENT_DEFAULT_STRATEGY.md"

RELEASE = RESULTS / "release_readiness.tsv"
CANDIDATE = RESULTS / "candidate_default_decision_audit.tsv"
PRESET_REPLAY = RESULTS / "candidate_preset_replay_audit.tsv"


def read_keyed_tsv(path: Path, key_col: str) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row[key_col]: row for row in reader}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def main() -> int:
    release = read_keyed_tsv(RELEASE, "gate")
    candidate = read_keyed_tsv(CANDIDATE, "metric")
    replay = read_keyed_tsv(PRESET_REPLAY, "metric")

    final = release.get("final_decision", {})
    preset = release.get("candidate_preset_ergonomics", {})
    abundance = release.get("abundance_default_not_replaced_by_local_candidate", {})
    ani = release.get("diagnostic_ani_reporting_available", {})

    require(
        preset.get("decision") == "candidate_preset_is_selected_default",
        "candidate preset is not the selected default in release_readiness.tsv",
    )
    require(
        final.get("status") == "pre_release_candidate",
        "release status changed; update this manifest script before regenerating",
    )
    require(
        "Use universal-auto-exact plus the candidate preset" in final.get("decision", ""),
        "final release decision no longer matches the documented default",
    )

    rows = [
        {
            "item": "selected_entrypoint",
            "value": "scripts/minco_profile",
            "evidence": "README.md; docs/USER_MANUAL.md; release_readiness.tsv",
            "status": release.get("single_default_entrypoint", {}).get("status", ""),
        },
        {
            "item": "selected_profile_preset",
            "value": "candidate",
            "evidence": preset.get("evidence", ""),
            "status": preset.get("status", ""),
        },
        {
            "item": "base_strategy",
            "value": "universal-auto-exact",
            "evidence": replay.get("preset_switch_values", {}).get("value", ""),
            "status": "selected_default",
        },
        {
            "item": "candidate_rescue_switch",
            "value": "emitted-ani90-xny100-br01-af70",
            "evidence": replay.get("preset_switch_values", {}).get("value", ""),
            "status": "selected_default",
        },
        {
            "item": "candidate_surface_switch",
            "value": "accession-ani90-xny100-br01-af70",
            "evidence": replay.get("preset_switch_values", {}).get("value", ""),
            "status": "selected_default",
        },
        {
            "item": "candidate_surface_taxmap_requirement",
            "value": "GTDB accession-level candidate_surface_taxmap.tsv/accession_species_taxmap.tsv sidecar or MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP",
            "evidence": "scripts/minco_profile_default.py; scripts/minco_profile_calibrated.py; hmp_gastrooral_sample1_sidecar_discovery_check",
            "status": "required_for_accession_surface_default",
        },
        {
            "item": "candidate_abundance_policy",
            "value": "normalized-depth-alpha2",
            "evidence": abundance.get("evidence", ""),
            "status": abundance.get("decision", ""),
        },
        {
            "item": "previous_default_reproduction",
            "value": "--profile-preset current or MINCO_PROFILE_PRESET=current",
            "evidence": "scripts/minco_profile_default.py; README.md; docs/USER_MANUAL.md",
            "status": "compatibility_path",
        },
        {
            "item": "default_user_command",
            "value": "scripts/minco_profile -r ref.minco --reads reads.fq.gz -p16 -o calibrated.profile.tsv",
            "evidence": "README.md; docs/USER_MANUAL.md",
            "status": "no_manual_strategy",
        },
        {
            "item": "default_preflight_command",
            "value": "scripts/minco_profile --check-ref -r ref.minco --reads reads.fq.gz",
            "evidence": "scripts/minco_profile_default.py; README.md; docs/USER_MANUAL.md",
            "status": "preflight_validates_packaged_sidecars",
        },
        {
            "item": "packaged_profile_defaults",
            "value": "optional minco_profile_defaults.tsv sidecar supports whitelisted domain defaults such as scope=bacteria",
            "evidence": "scripts/minco_profile_default.py; GTDBr232 S1000 packaged smoke",
            "status": "domain_default_without_user_options",
        },
        {
            "item": "scoring_priority",
            "value": "F1 first, then official abundance L1/Pearson, then ANI/reporting robustness",
            "evidence": "goal; release_readiness.tsv; candidate_default_decision_audit.tsv",
            "status": "documented_priority",
        },
        {
            "item": "candidate_vs_previous_minco_mean_delta",
            "value": candidate.get("candidate_vs_current_mean_delta", {}).get("value", ""),
            "evidence": candidate.get("candidate_vs_current_mean_delta", {}).get("evidence", ""),
            "status": candidate.get("candidate_vs_current_panel_safety", {}).get("decision", ""),
        },
        {
            "item": "candidate_preset_replay",
            "value": replay.get("max_delta_vs_prior_wrapper_candidate", {}).get("value", ""),
            "evidence": PRESET_REPLAY.relative_to(EXP).as_posix(),
            "status": replay.get("promotion_decision", {}).get("decision", ""),
        },
        {
            "item": "sylph_claim_boundary",
            "value": candidate.get("candidate_vs_sylph_panel_wins", {}).get("value", ""),
            "evidence": candidate.get("candidate_vs_sylph_panel_wins", {}).get("evidence", ""),
            "status": candidate.get("candidate_vs_sylph_panel_wins", {}).get("decision", ""),
        },
        {
            "item": "ani_reporting",
            "value": "reported_ani uses ZIP-AAF ANI; raw/emitted readwise ANI is diagnostic only",
            "evidence": ani.get("evidence", ""),
            "status": ani.get("decision", ""),
        },
        {
            "item": "release_status",
            "value": final.get("status", ""),
            "evidence": final.get("evidence", ""),
            "status": final.get("decision", ""),
        },
    ]

    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["item", "value", "evidence", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    with OUT_MD.open("w", newline="\n") as handle:
        handle.write("# Current Default Strategy\n\n")
        handle.write(
            "Generated by `write_current_default_strategy_manifest.py` from the "
            "release-readiness and candidate audit TSVs. It records the selected "
            "MinCO profiling default and its claim boundary without rerunning profiles.\n\n"
        )
        handle.write("| Item | Value | Status |\n")
        handle.write("| --- | --- | --- |\n")
        for row in rows:
            item = row["item"].replace("|", "\\|")
            value = row["value"].replace("|", "\\|")
            status = row["status"].replace("|", "\\|")
            handle.write(f"| `{item}` | {value} | {status} |\n")
        handle.write("\nThe full evidence paths are in `results/current_default_strategy_manifest.tsv`.\n")
    print(OUT)
    print(OUT_MD)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
