#!/usr/bin/env python3
"""Score defaked reliable-depth breadth gates on Toy Mouse sample0."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import build_and_score_gtdb_ground_truth as truth


EXP_DIR = Path(__file__).resolve().parent
MINCO_TSV = EXP_DIR / "toymouse_sample0_s2000_marker_split_naive_product_nb_reliable_abundance_p005.tsv"
DEFAULT_ANI_THRESHOLD = 0.95
DEFAULT_CTX_LENGTH = 22
DEFAULT_EFFECTIVE_CTX_LENGTH = DEFAULT_CTX_LENGTH + 2
DEFAULT_AF_FLOOR = math.exp((DEFAULT_ANI_THRESHOLD - 1.0) * DEFAULT_EFFECTIVE_CTX_LENGTH)


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


def ztp_conditional_variance(lam: float) -> float:
    if lam <= 0.0:
        return 0.0
    p_nonzero = 1.0 - math.exp(-lam)
    mean_pos = lam / p_nonzero
    return (lam + lam * lam) / p_nonzero - mean_pos * mean_pos


def ztp_adjusted_af(breadth: float, mean_pos: float) -> float:
    if not math.isfinite(breadth) or breadth <= 0.0:
        return 0.0
    lam = ztp_lambda_from_pos_mean(mean_pos)
    p_nonzero = 1.0 - math.exp(-lam) if lam > 0.0 else 0.0
    if p_nonzero <= 0.0:
        return 1.0
    return min(1.0, breadth / p_nonzero)


def nb_truncated_moments(mu: float, size: float):
    if mu <= 0.0 or size <= 0.0:
        return None
    log_p0 = size * (math.log(size) - math.log(size + mu))
    p0 = 0.0 if log_p0 < -745.0 else math.exp(log_p0)
    p_nonzero = 1.0 - p0
    if p_nonzero <= 0.0:
        return None
    var = mu + mu * mu / size
    mean_pos = mu / p_nonzero
    var_pos = (var + mu * mu) / p_nonzero - mean_pos * mean_pos
    if var_pos < 0.0 and var_pos > -1e-8:
        var_pos = 0.0
    return mean_pos, max(0.0, var_pos), p_nonzero


def ztnb_p_nonzero(mean_pos: float, var_pos: float, hit_ctx: int):
    if not math.isfinite(mean_pos) or mean_pos <= 0.0:
        return 0.0, "bad"

    lam = ztp_lambda_from_pos_mean(mean_pos)
    p_ztp = 1.0 - math.exp(-lam) if lam > 0.0 else 0.0
    var_ztp = ztp_conditional_variance(lam)
    if hit_ctx < 10 or (not math.isfinite(var_pos)) or var_pos <= var_ztp * 1.10 + 1e-9:
        return p_ztp, "ztp"

    eps = 1e-9

    def residual(log_params: np.ndarray) -> np.ndarray:
        mu = math.exp(float(log_params[0]))
        size = math.exp(float(log_params[1]))
        moments = nb_truncated_moments(mu, size)
        if moments is None:
            return np.array([1e3, 1e3])
        fit_mean, fit_var, _ = moments
        return np.array(
            [
                math.log((fit_mean + eps) / (mean_pos + eps)),
                math.log((fit_var + eps) / (var_pos + eps)),
            ]
        )

    starts = [
        (max(lam, 1e-6), 10.0),
        (max(lam, 1e-6), 1.0),
        (max(min(mean_pos, 100.0), 1e-6), 1.0),
        (max(mean_pos / 2.0, 1e-6), 0.2),
    ]
    best = None
    for mu0, size0 in starts:
        fit = least_squares(
            residual,
            np.log([mu0, size0]),
            bounds=([-30.0, -20.0], [20.0, 30.0]),
            max_nfev=200,
        )
        score = float(np.sum(fit.fun * fit.fun))
        if best is None or score < best[0]:
            best = (score, fit)

    if best is None or (not best[1].success and best[0] > 1e-4):
        return p_ztp, "fallback_ztp"
    mu = float(math.exp(best[1].x[0]))
    size = float(math.exp(best[1].x[1]))
    moments = nb_truncated_moments(mu, size)
    if moments is None:
        return p_ztp, "fallback_ztp"
    _, _, p_nonzero = moments
    if p_nonzero <= 0.0 or not math.isfinite(p_nonzero):
        return p_ztp, "fallback_ztp"
    return p_nonzero, "ztnb"


def score_pred(predicted, gold):
    pred = {str(x) for x in predicted if str(x)}
    tp = pred & gold
    fp = pred - gold
    fn = gold - pred
    precision = len(tp) / len(pred) if pred else 0.0
    recall = len(tp) / len(gold) if gold else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return pred, tp, fp, fn, precision, recall, f1


def main() -> int:
    by_accession, by_core = truth.load_gtdb_metadata(truth.GTDB_METADATA)
    profile = pd.read_csv(EXP_DIR / "mouse0_gtdb_species_profile.tsv", sep="\t")
    gold = set(profile["gtdb_species"].astype(str))
    abundance = dict(zip(profile["gtdb_species"].astype(str), profile["relative_abundance"]))

    rows = truth.load_minco_with_gtdb_species(MINCO_TSV, by_accession, by_core)
    rows = truth.add_naive_ani(rows)
    numeric_cols = [
        "XnY_ctx",
        "ANI_naive_calc",
        "Ref_zip_af",
        "Reliable_Ref_zip_af",
        "Reliable_Ref_breadth",
        "Reliable_Ref_hit_ctx",
        "Reliable_Ref_hit_mean_depth",
        "Reliable_Ref_hit_depth_variance",
    ]
    for col in numeric_cols:
        rows[col] = pd.to_numeric(rows[col], errors="coerce").fillna(0.0)

    ztp_af = []
    ztnb_af = []
    ztnb_model = []
    for breadth, mean_pos, var_pos, hit_ctx in zip(
        rows["Reliable_Ref_breadth"],
        rows["Reliable_Ref_hit_mean_depth"],
        rows["Reliable_Ref_hit_depth_variance"],
        rows["Reliable_Ref_hit_ctx"],
    ):
        ztp_af.append(ztp_adjusted_af(float(breadth), float(mean_pos)))
        p_nonzero, model = ztnb_p_nonzero(float(mean_pos), float(var_pos), int(hit_ctx))
        ztnb_model.append(model)
        if p_nonzero > 0.0:
            ztnb_af.append(min(1.0, float(breadth) / p_nonzero))
        else:
            ztnb_af.append(1.0 if float(breadth) > 0.0 else 0.0)
    rows["Reliable_ztp_af"] = ztp_af
    rows["Reliable_ztnb_af"] = ztnb_af
    rows["Reliable_depth_model"] = ztnb_model

    base = (
        (rows["XnY_ctx"] >= 10.0)
        & (rows["ANI_naive_calc"] > DEFAULT_ANI_THRESHOLD)
        & rows["gtdb_species"].astype(bool)
    )
    variants = {
        "raw_ref_zip": "Ref_zip_af",
        "reliable_zip": "Reliable_Ref_zip_af",
        "reliable_ztp": "Reliable_ztp_af",
        "reliable_ztnb_adaptive": "Reliable_ztnb_af",
    }

    score_rows = []
    detail_rows = []
    for method, col in variants.items():
        selected = rows.loc[base & (rows[col] >= DEFAULT_AF_FLOOR)].copy()
        pred, tp, fp, fn, precision, recall, f1 = score_pred(selected["gtdb_species"], gold)
        score_rows.append(
            {
                "method": method,
                "af_floor": DEFAULT_AF_FLOOR,
                "pred_taxa": len(pred),
                "TP": len(tp),
                "FP": len(fp),
                "FN": len(fn),
                "precision": precision,
                "recall": recall,
                "F1": f1,
            }
        )
        for kind, species_set in [("FP", fp), ("FN", fn)]:
            for species in sorted(species_set):
                rec = {
                    "method": method,
                    "kind": kind,
                    "gtdb_species": species,
                    "gold_abundance": abundance.get(species, ""),
                    "accession": "",
                    "ANI_naive_calc": "",
                    "XnY_ctx": "",
                    "Reliable_Ref_breadth": "",
                    "Reliable_ztp_af": "",
                    "Reliable_ztnb_af": "",
                    "Reliable_depth_model": "",
                }
                examples = selected.loc[selected["gtdb_species"] == species].copy()
                if not examples.empty:
                    best = examples.sort_values(["ANI_naive_calc", "XnY_ctx"], ascending=False).iloc[0]
                    rec.update(
                        {
                            "accession": best.get("accession", ""),
                            "ANI_naive_calc": best.get("ANI_naive_calc", ""),
                            "XnY_ctx": best.get("XnY_ctx", ""),
                            "Reliable_Ref_breadth": best.get("Reliable_Ref_breadth", ""),
                            "Reliable_ztp_af": best.get("Reliable_ztp_af", ""),
                            "Reliable_ztnb_af": best.get("Reliable_ztnb_af", ""),
                            "Reliable_depth_model": best.get("Reliable_depth_model", ""),
                        }
                    )
                detail_rows.append(rec)

    sweep_rows = []
    thresholds = sorted(set([DEFAULT_AF_FLOOR] + [x / 100.0 for x in range(20, 81, 5)]))
    for method, col in variants.items():
        for threshold in thresholds:
            selected = rows.loc[base & (rows[col] >= threshold)].copy()
            pred, tp, fp, fn, precision, recall, f1 = score_pred(selected["gtdb_species"], gold)
            sweep_rows.append(
                {
                    "method": method,
                    "af_floor": threshold,
                    "pred_taxa": len(pred),
                    "TP": len(tp),
                    "FP": len(fp),
                    "FN": len(fn),
                    "precision": precision,
                    "recall": recall,
                    "F1": f1,
                }
            )

    score_path = EXP_DIR / "gtdb_direct_scores_reliable_truncated_depth.tsv"
    detail_path = EXP_DIR / "gtdb_direct_score_reliable_truncated_depth_details.tsv"
    sweep_path = EXP_DIR / "gtdb_direct_score_reliable_truncated_depth_sweep.tsv"
    pd.DataFrame(score_rows).to_csv(score_path, sep="\t", index=False)
    pd.DataFrame(detail_rows).to_csv(detail_path, sep="\t", index=False)
    pd.DataFrame(sweep_rows).to_csv(sweep_path, sep="\t", index=False)

    print(f"af_floor={DEFAULT_AF_FLOOR:.12g}")
    print(pd.DataFrame(score_rows).to_csv(sep="\t", index=False), end="")
    print(f"wrote {score_path}")
    print(f"wrote {detail_path}")
    print(f"wrote {sweep_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
