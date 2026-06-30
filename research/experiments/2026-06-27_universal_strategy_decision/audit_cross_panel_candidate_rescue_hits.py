#!/usr/bin/env python3
"""List species rescued by the best cross-panel zero-mass candidate rule."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

import audit_missed_truth_candidates as missed
import decompose_abundance_errors as decomp
import sweep_cross_panel_candidate_rescue as sweep


EXP = Path(__file__).resolve().parent
RESULTS = EXP / "results"

BEST_ANI_MIN = 0.90
BEST_XNY_MIN = 100.0
BEST_BREADTH_MIN = 0.01
BEST_REAL_AF_MIN = 0.70


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    log("initializing mapper")
    mapper = missed.CandidateMapper()
    log("initialized mapper")
    rows = []
    for panel in sweep.PANELS:
        cfg = decomp.PANELS[panel]
        for sample in cfg["samples"]:
            sample_id = int(sample)
            log(f"auditing {panel} sample{sample_id}")
            truth = decomp.normalize(
                decomp.load_truth(
                    Path(cfg["truth"](sample_id)),
                    sample_id,
                    str(cfg["truth_abundance_col"]),
                )
            )
            profile_path = Path(cfg["minco"](sample_id))
            if panel in {
                "hmp_airskin_gtdb_source_abundance",
                "hmp_gastrooral_gtdb_source_abundance",
            }:
                current_pred = sweep.load_hmp_current_pred(
                    panel, profile_path, str(cfg["minco_collapse"]), mapper
                )
                candidates = sweep.load_hmp_raw_candidates(panel, sample_id, set(current_pred))
            else:
                current_pred, candidates = sweep.load_profile_pred_and_candidates(
                    panel, sample_id, profile_path, str(cfg["minco_collapse"]), mapper
                )
            rescued = sweep.rescue_species(
                candidates,
                BEST_ANI_MIN,
                BEST_XNY_MIN,
                BEST_BREADTH_MIN,
                BEST_REAL_AF_MIN,
            )
            if not rescued:
                continue
            hit = candidates.loc[candidates["gtdb_species"].astype(str).isin(rescued)].copy()
            hit = hit.sort_values(["gtdb_species", "ANI", "XnY"], ascending=[True, False, False])
            hit = hit.drop_duplicates("gtdb_species", keep="first")
            for rec in hit.itertuples(index=False):
                species = str(rec.gtdb_species)
                truth_abundance = float(truth.get(species, 0.0))
                rows.append(
                    {
                        "panel": panel,
                        "sample": sample_id,
                        "gtdb_species": species,
                        "rescued_truth_status": "TP" if truth_abundance > 0.0 else "FP",
                        "truth_abundance_pct": truth_abundance * 100.0,
                        "candidate_source": str(rec.source),
                        "ANI": float(rec.ANI),
                        "XnY": float(rec.XnY),
                        "breadth": float(rec.breadth),
                        "real_af": float(rec.real_af),
                        "current_called_species": len(current_pred),
                        "profile_path": str(profile_path),
                    }
                )

    detail = pd.DataFrame(rows)
    if detail.empty:
        detail = pd.DataFrame(
            columns=[
                "panel",
                "sample",
                "gtdb_species",
                "rescued_truth_status",
                "truth_abundance_pct",
                "candidate_source",
                "ANI",
                "XnY",
                "breadth",
                "real_af",
                "current_called_species",
                "profile_path",
            ]
        )
    detail = detail.sort_values(
        ["rescued_truth_status", "truth_abundance_pct", "panel", "sample"],
        ascending=[False, False, True, True],
    )
    summary_rows = []
    if not detail.empty:
        for panel, sub in detail.groupby("panel", sort=True):
            summary_rows.append(
                {
                    "panel": panel,
                    "rescued_species": int(len(sub)),
                    "rescued_TP": int(sub["rescued_truth_status"].eq("TP").sum()),
                    "rescued_FP": int(sub["rescued_truth_status"].eq("FP").sum()),
                    "rescued_truth_mass_pct": sub.loc[
                        sub["rescued_truth_status"].eq("TP"), "truth_abundance_pct"
                    ].sum(),
                    "samples_with_rescue": int(sub["sample"].nunique()),
                    "candidate_sources": ",".join(sorted(set(sub["candidate_source"].astype(str)))),
                }
            )
    summary = pd.DataFrame(summary_rows)
    audit = pd.DataFrame(
        [
            {
                "metric": "selected_rule",
                "value": "cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass",
                "evidence": "ANI>=0.90; XnY>=100; breadth>=0.01; real_af>=0.70; zero rescued abundance mass",
                "decision": "best_cross_panel_candidate_rule",
            },
            {
                "metric": "rescued_species",
                "value": int(len(detail)),
                "evidence": f"TP={int(detail['rescued_truth_status'].eq('TP').sum()) if not detail.empty else 0};FP={int(detail['rescued_truth_status'].eq('FP').sum()) if not detail.empty else 0}",
                "decision": "detail_counts",
            },
            {
                "metric": "rescued_truth_mass_pct",
                "value": float(
                    detail.loc[detail["rescued_truth_status"].eq("TP"), "truth_abundance_pct"].sum()
                )
                if not detail.empty
                else 0.0,
                "evidence": "sum of truth abundance over rescued true-positive species across samples",
                "decision": "abundance_policy_target",
            },
            {
                "metric": "promotion_decision",
                "value": "diagnostic_only_not_default",
                "evidence": "detail explains rescued species; no wrapper implementation or abundance assignment policy",
                "decision": "use_for_wrapper_design",
            },
        ]
    )
    detail.to_csv(RESULTS / "cross_panel_candidate_rescue_best_hits.tsv", sep="\t", index=False)
    summary.to_csv(RESULTS / "cross_panel_candidate_rescue_best_hits_summary.tsv", sep="\t", index=False)
    audit.to_csv(RESULTS / "cross_panel_candidate_rescue_best_hits_audit.tsv", sep="\t", index=False)
    print(audit.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
