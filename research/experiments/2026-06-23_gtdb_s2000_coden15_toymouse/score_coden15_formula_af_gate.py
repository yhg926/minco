#!/usr/bin/env python3
"""Score coden15 gates whose AF floor is derived from ANI and context length."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
EXP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(EXP_DIR))

import optimize_coden15_params_fast as opt  # noqa: E402


def main() -> int:
    rows_by_sample = {sample: opt.load_cached_rows(sample, None, None) for sample in opt.SAMPLES}
    truth_by_sample = {sample: opt.load_truth(sample) for sample in opt.SAMPLES}
    base = {
        "xny_min": 15.0,
        "af_floor": 0.40,
        "vmr_min": 10.0,
        "delta_max": 0.03,
        "hit_mean_min": 3.0,
        "ani_min": 0.95,
        "abundance_exponent": 1.0,
        "median_cutoff": 20.0,
        "rescue_enabled": True,
        "rescue_ratio": 3.0,
        "rescue_zip_af_max": 0.25,
        "rescue_xny_min": 100.0,
        "rescue_effective_min": 5.0,
    }

    records = []
    for eff_len in [20, 22, 24, 26, 28, 30, 32]:
        converted = base["ani_min"] ** eff_len
        for scale in [1.0, 1.03, 1.05, 1.08, 1.10, 1.15]:
            af_floor = min(0.999, converted * scale)
            for xny_min in [10.0, 15.0, 20.0]:
                for vmr_min in [10.0, 15.0, 20.0, 50.0]:
                    for delta_max in [0.02, 0.03, 0.04]:
                        for median_cutoff in [10.0, 20.0]:
                            params = dict(base)
                            params.update(
                                {
                                    "xny_min": xny_min,
                                    "af_floor": af_floor,
                                    "vmr_min": vmr_min,
                                    "delta_max": delta_max,
                                    "median_cutoff": median_cutoff,
                                }
                            )
                            mean, _ = opt.score_config(rows_by_sample, truth_by_sample, params)
                            records.append(
                                {
                                    "eff_ctx_len": eff_len,
                                    "af_scale": scale,
                                    "af_floor": af_floor,
                                    "xny_min": xny_min,
                                    "vmr_min": vmr_min,
                                    "delta_max": delta_max,
                                    "median_cutoff": median_cutoff,
                                    **mean,
                                }
                            )

    grid = pd.DataFrame(records).sort_values("l1_pct_points")
    grid_path = EXP_DIR / "coden15_formula_af_grid.tsv"
    grid.to_csv(grid_path, sep="\t", index=False)

    best = grid.iloc[0].to_dict()
    best.update(
        {
            "hit_mean_min": 3.0,
            "ani_min": 0.95,
            "abundance_exponent": 1.0,
            "rescue_enabled": True,
            "rescue_ratio": 3.0,
            "rescue_zip_af_max": 0.25,
            "rescue_xny_min": 100.0,
            "rescue_effective_min": 5.0,
        }
    )
    _, sample_rows = opt.score_config(rows_by_sample, truth_by_sample, best)
    sample_df = pd.DataFrame(sample_rows)
    for key, value in best.items():
        if key not in sample_df.columns:
            sample_df[key] = value
    sample_path = EXP_DIR / "coden15_formula_af_best_sample_metrics.tsv"
    sample_df.to_csv(sample_path, sep="\t", index=False)

    print(grid.head(30).to_csv(sep="\t", index=False), end="")
    print(f"wrote {grid_path}")
    print(f"wrote {sample_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
