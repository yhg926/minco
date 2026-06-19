#!/usr/bin/env python3
"""Build a stratified minco-vs-ANIm calibration pilot from the Nayfach table."""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ANI_BINS: Sequence[Tuple[str, float, float]] = (
    ("lt090", -math.inf, 0.90),
    ("090_094", 0.90, 0.94),
    ("094_095", 0.94, 0.95),
    ("095_097", 0.95, 0.97),
    ("097_099", 0.97, 0.99),
    ("gte099", 0.99, math.inf),
)

SAMPLE_HEADER = (
    "pair_id",
    "ani_bin",
    "ref_fna",
    "qry_fna",
    "anim_truth",
    "kssd3_t10_xny",
    "kssd3_t10_qaf",
    "kssd3_t10_raf",
    "kssd3_t10_n_diff_obj",
    "kssd3_t10_n_diff_obj_section",
    "kssd3_t10_n_mut2_ctx",
    "kssd3_t10_raw_ani",
)

MINCO_HEADER = (
    "Qry",
    "Ref",
    "ANI",
    "Distance",
    "Confidence",
    "Selected_metric",
    "XnY_ctx",
    "Qry_align_fraction",
    "blastn_Qry_align_fraction",
    "Ref_align_fraction",
    "blastn_Ref_align_fraction",
    "N_diff_obj",
    "N_diff_obj_section",
    "N_mut2_ctx",
    "Ref_annotation",
)

FEATURE_HEADER = SAMPLE_HEADER + tuple(f"minco_{c}" for c in MINCO_HEADER) + (
    "status",
    "seconds",
    "stderr_excerpt",
)

DERIVED_FEATURES = (
    "minco_ANI",
    "minco_Distance",
    "minco_XnY_ctx",
    "minco_Qry_align_fraction",
    "minco_blastn_Qry_align_fraction",
    "minco_Ref_align_fraction",
    "minco_blastn_Ref_align_fraction",
    "minco_N_diff_obj",
    "minco_N_diff_obj_section",
    "minco_N_mut2_ctx",
    "minco_min_af",
    "minco_max_af",
    "minco_af_delta",
    "minco_log_xny",
    "minco_mut2_per_xny",
    "minco_diffobj_per_xny",
    "minco_section_per_xny",
    "minco_raw_error",
)

DROPIN_HGB_FEATURES = (
    "raw_ani",
    "ref_af",
    "log1p_xny",
    "diff_obj_rate",
    "diff_section_rate",
    "mut2_rate",
    "ref_aaf_ani",
    "raw_minus_ref_aaf",
    "ref_af_lt_0_2",
    "ref_af_0_2_to_0_5",
    "ref_af_ge_0_5",
)


@dataclass(frozen=True)
class PairRow:
    pair_id: str
    ani_bin: str
    ref_fna: str
    qry_fna: str
    anim_truth: float
    kssd3_t10_xny: float
    kssd3_t10_qaf: float
    kssd3_t10_raf: float
    kssd3_t10_n_diff_obj: float
    kssd3_t10_n_diff_obj_section: float
    kssd3_t10_n_mut2_ctx: float
    kssd3_t10_raw_ani: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "pair_id": self.pair_id,
            "ani_bin": self.ani_bin,
            "ref_fna": self.ref_fna,
            "qry_fna": self.qry_fna,
            "anim_truth": self.anim_truth,
            "kssd3_t10_xny": self.kssd3_t10_xny,
            "kssd3_t10_qaf": self.kssd3_t10_qaf,
            "kssd3_t10_raf": self.kssd3_t10_raf,
            "kssd3_t10_n_diff_obj": self.kssd3_t10_n_diff_obj,
            "kssd3_t10_n_diff_obj_section": self.kssd3_t10_n_diff_obj_section,
            "kssd3_t10_n_mut2_ctx": self.kssd3_t10_n_mut2_ctx,
            "kssd3_t10_raw_ani": self.kssd3_t10_raw_ani,
        }


def ani_bin(value: float) -> str:
    for label, lo, hi in ANI_BINS:
        if lo <= value < hi:
            return label
    raise ValueError(f"ANI {value} did not fit configured bins")


def fna_path(raw: str, fna_root: Path) -> str:
    return str(fna_root / Path(raw).name)


def parse_anim_table_row(row: Sequence[str], fna_root: Path, pair_id: str) -> PairRow:
    if len(row) < 10:
        raise ValueError(f"expected at least 10 columns, saw {len(row)}")
    truth = float(row[9])
    return PairRow(
        pair_id=pair_id,
        ani_bin=ani_bin(truth),
        ref_fna=fna_path(row[0], fna_root),
        qry_fna=fna_path(row[1], fna_root),
        anim_truth=truth,
        kssd3_t10_xny=float(row[2]),
        kssd3_t10_qaf=float(row[3]),
        kssd3_t10_raf=float(row[4]),
        kssd3_t10_n_diff_obj=float(row[5]),
        kssd3_t10_n_diff_obj_section=float(row[6]),
        kssd3_t10_n_mut2_ctx=float(row[7]),
        kssd3_t10_raw_ani=float(row[8]),
    )


def sample_pairs(anim_table: Path, fna_root: Path, per_bin: int, seed: int) -> Tuple[List[PairRow], Dict[str, int]]:
    rng = random.Random(seed)
    reservoirs: Dict[str, List[PairRow]] = {label: [] for label, _, _ in ANI_BINS}
    counts: Dict[str, int] = {label: 0 for label, _, _ in ANI_BINS}
    total = 0
    with anim_table.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            total += 1
            pair = parse_anim_table_row(row, fna_root, f"P{total:09d}")
            label = pair.ani_bin
            counts[label] += 1
            bucket = reservoirs[label]
            if len(bucket) < per_bin:
                bucket.append(pair)
            else:
                j = rng.randrange(counts[label])
                if j < per_bin:
                    bucket[j] = pair
    sampled: List[PairRow] = []
    for label, _, _ in ANI_BINS:
        sampled.extend(reservoirs[label])
    sampled.sort(key=lambda x: (x.ani_bin, x.pair_id))
    counts["total_rows"] = total
    return sampled, counts


def write_sample_pairs(pairs: Sequence[PairRow], out_path: Path) -> None:
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_HEADER, delimiter="\t")
        writer.writeheader()
        for pair in pairs:
            writer.writerow(pair.as_dict())


def read_sample_pairs(path: Path) -> List[PairRow]:
    pairs: List[PairRow] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pairs.append(
                PairRow(
                    pair_id=row["pair_id"],
                    ani_bin=row["ani_bin"],
                    ref_fna=row["ref_fna"],
                    qry_fna=row["qry_fna"],
                    anim_truth=float(row["anim_truth"]),
                    kssd3_t10_xny=float(row["kssd3_t10_xny"]),
                    kssd3_t10_qaf=float(row["kssd3_t10_qaf"]),
                    kssd3_t10_raf=float(row["kssd3_t10_raf"]),
                    kssd3_t10_n_diff_obj=float(row["kssd3_t10_n_diff_obj"]),
                    kssd3_t10_n_diff_obj_section=float(row["kssd3_t10_n_diff_obj_section"]),
                    kssd3_t10_n_mut2_ctx=float(row["kssd3_t10_n_mut2_ctx"]),
                    kssd3_t10_raw_ani=float(row["kssd3_t10_raw_ani"]),
                )
            )
    return pairs


def parse_minco_output(path: Path) -> Dict[str, str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise RuntimeError("minco produced header but no comparison row")
    return rows[0]


def run_minco_pair(
    minco_bin: Path,
    pair: PairRow,
    sketch_size: int,
    minco_threads: int,
    timeout: int,
    tmp_dir: Path,
) -> Dict[str, object]:
    start = time.perf_counter()
    row = pair.as_dict()
    with tempfile.NamedTemporaryFile(
        prefix=f"{pair.pair_id}.",
        suffix=".tsv",
        dir=tmp_dir,
        delete=False,
    ) as tmp:
        out_path = Path(tmp.name)
    cmd = [
        str(minco_bin),
        "ani",
        "-S",
        str(sketch_size),
        "-p",
        str(minco_threads),
        "-f0",
        "-n0",
        "-t0",
        "--raw-output",
        "-s3",
        "-o",
        str(out_path),
        pair.ref_fna,
        pair.qry_fna,
    ]
    try:
        completed = subprocess.run(
            cmd,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
        row["seconds"] = f"{time.perf_counter() - start:.6f}"
        row["stderr_excerpt"] = " ".join(completed.stderr.split())[:240]
        if completed.returncode != 0:
            row["status"] = f"exit_{completed.returncode}"
            for key in MINCO_HEADER:
                row[f"minco_{key}"] = ""
            return row
        minco_row = parse_minco_output(out_path)
        for key in MINCO_HEADER:
            row[f"minco_{key}"] = minco_row.get(key, "")
        row["status"] = "ok"
        return row
    except Exception as exc:  # noqa: BLE001 - record per-pair failure and keep the run going.
        row["seconds"] = f"{time.perf_counter() - start:.6f}"
        row["stderr_excerpt"] = str(exc)[:240]
        row["status"] = "exception"
        for key in MINCO_HEADER:
            row[f"minco_{key}"] = ""
        return row
    finally:
        try:
            out_path.unlink()
        except FileNotFoundError:
            pass


def extract_features(
    pairs: Sequence[PairRow],
    minco_bin: Path,
    sketch_size: int,
    jobs: int,
    minco_threads: int,
    timeout: int,
    out_path: Path,
    tmp_dir: Path,
) -> None:
    tmp_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Dict[str, object]] = {}
    done = 0
    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(run_minco_pair, minco_bin, pair, sketch_size, minco_threads, timeout, tmp_dir): pair
            for pair in pairs
        }
        for future in as_completed(futures):
            pair = futures[future]
            results[pair.pair_id] = future.result()
            done += 1
            if done == 1 or done % 100 == 0 or done == len(pairs):
                elapsed = time.perf_counter() - started
                rate = done / elapsed if elapsed > 0 else 0.0
                print(
                    f"features {done}/{len(pairs)} pairs, {rate:.2f} pair/s",
                    file=sys.stderr,
                    flush=True,
                )
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FEATURE_HEADER, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for pair in pairs:
            writer.writerow(results[pair.pair_id])


def as_float(row: Dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "" or value is None:
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return float("nan")


def safe_div(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0:
        return float("nan")
    return num / den


def feature_vector(row: Dict[str, str]) -> Dict[str, float]:
    xny = as_float(row, "minco_XnY_ctx")
    qaf = as_float(row, "minco_Qry_align_fraction")
    raf = as_float(row, "minco_Ref_align_fraction")
    nmut = as_float(row, "minco_N_mut2_ctx")
    ndiff = as_float(row, "minco_N_diff_obj")
    nsection = as_float(row, "minco_N_diff_obj_section")
    raw_ani = as_float(row, "minco_ANI")
    values = {
        "minco_ANI": raw_ani,
        "minco_Distance": as_float(row, "minco_Distance"),
        "minco_XnY_ctx": xny,
        "minco_Qry_align_fraction": qaf,
        "minco_blastn_Qry_align_fraction": as_float(row, "minco_blastn_Qry_align_fraction"),
        "minco_Ref_align_fraction": raf,
        "minco_blastn_Ref_align_fraction": as_float(row, "minco_blastn_Ref_align_fraction"),
        "minco_N_diff_obj": ndiff,
        "minco_N_diff_obj_section": nsection,
        "minco_N_mut2_ctx": nmut,
        "minco_min_af": min(qaf, raf) if np.isfinite(qaf) and np.isfinite(raf) else float("nan"),
        "minco_max_af": max(qaf, raf) if np.isfinite(qaf) and np.isfinite(raf) else float("nan"),
        "minco_af_delta": abs(qaf - raf) if np.isfinite(qaf) and np.isfinite(raf) else float("nan"),
        "minco_log_xny": math.log1p(xny) if np.isfinite(xny) and xny >= 0 else float("nan"),
        "minco_mut2_per_xny": safe_div(nmut, xny),
        "minco_diffobj_per_xny": safe_div(ndiff, xny),
        "minco_section_per_xny": safe_div(nsection, xny),
        "minco_raw_error": 1.0 - raw_ani if np.isfinite(raw_ani) else float("nan"),
    }
    return values


def clamp01(value: float) -> float:
    if not np.isfinite(value):
        return value
    return min(1.0, max(0.0, value))


def dropin_hgb_feature_vector(row: Dict[str, str]) -> Dict[str, float]:
    raw_ani = clamp01(as_float(row, "minco_ANI"))
    ref_af = as_float(row, "minco_Ref_align_fraction")
    if not np.isfinite(ref_af) or ref_af < 1e-12:
        ref_af = 1e-12
    ref_af = min(1.0, ref_af)
    xny = as_float(row, "minco_XnY_ctx")
    nmut = as_float(row, "minco_N_mut2_ctx")
    ndiff = as_float(row, "minco_N_diff_obj")
    nsection = as_float(row, "minco_N_diff_obj_section")
    ref_aaf_ani = clamp01(1.0 + math.log(ref_af) / 11.0)
    return {
        "raw_ani": raw_ani,
        "ref_af": ref_af,
        "log1p_xny": math.log1p(xny) if np.isfinite(xny) and xny >= 0 else float("nan"),
        "diff_obj_rate": safe_div(ndiff, xny),
        "diff_section_rate": safe_div(nsection, xny),
        "mut2_rate": safe_div(nmut, xny),
        "ref_aaf_ani": ref_aaf_ani,
        "raw_minus_ref_aaf": raw_ani - ref_aaf_ani if np.isfinite(raw_ani) else float("nan"),
        "ref_af_lt_0_2": 1.0 if ref_af < 0.2 else 0.0,
        "ref_af_0_2_to_0_5": 1.0 if 0.2 <= ref_af < 0.5 else 0.0,
        "ref_af_ge_0_5": 1.0 if ref_af >= 0.5 else 0.0,
    }


def pearson_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2:
        return float("nan")
    y0 = y_true - y_true.mean()
    y1 = y_pred - y_pred.mean()
    den = math.sqrt(float((y0 * y0).sum() * (y1 * y1).sum()))
    if den == 0:
        return float("nan")
    return float((y0 * y1).sum() / den)


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    y_pred = np.clip(y_pred, 0.0, 1.0)
    err = y_pred - y_true
    truth_pos = y_true >= 0.95
    pred_pos = y_pred >= 0.95
    return {
        "n": float(len(y_true)),
        "bias": float(err.mean()),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "max_abs": float(np.max(np.abs(err))) if len(err) else float("nan"),
        "pearson_r": pearson_r(y_true, y_pred),
        "tp95": float(np.sum(truth_pos & pred_pos)),
        "fp95": float(np.sum(~truth_pos & pred_pos)),
        "fn95": float(np.sum(truth_pos & ~pred_pos)),
        "tn95": float(np.sum(~truth_pos & ~pred_pos)),
    }


def read_feature_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if row.get("status") == "ok"]


def c_float(value: float) -> str:
    return f"{float(value):.17g}"


def c_array(name: str, ctype: str, values: Sequence[object], per_line: int = 12) -> str:
    lines = [f"static const {ctype} {name}[{len(values)}] = {{"]
    for start in range(0, len(values), per_line):
        chunk = values[start : start + per_line]
        lines.append("    " + ", ".join(str(v) for v in chunk) + ",")
    lines.append("};")
    return "\n".join(lines)


def export_refaf_hgb_header(model: HistGradientBoostingRegressor, out_path: Path, note: str) -> None:
    baseline = float(np.asarray(model._baseline_prediction).ravel()[0])
    tree_offsets = [0]
    feature_idx: List[int] = []
    thresholds: List[str] = []
    values: List[str] = []
    left: List[int] = []
    right: List[int] = []
    is_leaf: List[int] = []
    missing_go_left: List[int] = []

    predictors = [stage[0] for stage in model._predictors]
    for predictor in predictors:
        nodes = predictor.nodes
        names = nodes.dtype.names or ()
        threshold_field = "num_threshold" if "num_threshold" in names else "threshold"
        missing_field = "missing_go_to_left" if "missing_go_to_left" in names else "missing_go_left"
        for node in nodes:
            feature_idx.append(int(node["feature_idx"]))
            thresholds.append(c_float(node[threshold_field]))
            values.append(c_float(node["value"]))
            left.append(int(node["left"]))
            right.append(int(node["right"]))
            is_leaf.append(int(node["is_leaf"]))
            missing_go_left.append(int(node[missing_field]))
        tree_offsets.append(tree_offsets[-1] + len(nodes))

    max_local_node = max(left + right + [0])
    if max_local_node > 65535:
        raise RuntimeError(f"local node index {max_local_node} exceeds unsigned short storage")
    n_nodes = tree_offsets[-1]
    lines = [
        "#ifndef MODEL_REFAF_HGB_H",
        "#define MODEL_REFAF_HGB_H",
        "",
        "#include <math.h>",
        "",
        "/* Generated by minco_nayfach_calibration_pilot.py.",
        f" * {note}",
        " * Feature order:",
        " *   0 raw_ani; 1 ref_af; 2 log1p_xny; 3 diff_obj_rate;",
        " *   4 diff_section_rate; 5 mut2_rate; 6 ref_aaf_ani;",
        " *   7 raw_minus_ref_aaf; 8 ref_af_lt_0_2;",
        " *   9 ref_af_0_2_to_0_5; 10 ref_af_ge_0_5.",
        " */",
        "#if NUM_CODENS == 11",
        "#define REFAF_HGB_N_FEATURES 11",
        f"#define REFAF_HGB_N_TREES {len(predictors)}",
        f"#define REFAF_HGB_N_NODES {n_nodes}",
        f"static const double refaf_hgb_baseline = {c_float(baseline)};",
        c_array("refaf_hgb_tree_offsets", "unsigned int", tree_offsets, per_line=12),
        c_array("refaf_hgb_feature_idx", "signed char", feature_idx, per_line=16),
        c_array("refaf_hgb_threshold", "double", thresholds, per_line=4),
        c_array("refaf_hgb_value", "double", values, per_line=4),
        c_array("refaf_hgb_left", "unsigned short", left, per_line=16),
        c_array("refaf_hgb_right", "unsigned short", right, per_line=16),
        c_array("refaf_hgb_is_leaf", "unsigned char", is_leaf, per_line=24),
        c_array("refaf_hgb_missing_go_left", "unsigned char", missing_go_left, per_line=24),
        "",
        "static inline double clamp01_refaf_hgb(double x)",
        "{",
        "    if (x < 0.0) return 0.0;",
        "    if (x > 1.0) return 1.0;",
        "    return x;",
        "}",
        "",
        "static inline double refaf_hgb_predict_from_values(const double x[REFAF_HGB_N_FEATURES])",
        "{",
        "    double pred = refaf_hgb_baseline;",
        "    for (unsigned int t = 0; t < REFAF_HGB_N_TREES; ++t) {",
        "        unsigned int node = refaf_hgb_tree_offsets[t];",
        "        const unsigned int begin = refaf_hgb_tree_offsets[t];",
        "        while (!refaf_hgb_is_leaf[node]) {",
        "            const int f = refaf_hgb_feature_idx[node];",
        "            const double v = x[f];",
        "            const unsigned int local_next = (isnan(v) ? refaf_hgb_missing_go_left[node] : (v <= refaf_hgb_threshold[node]))",
        "                ? refaf_hgb_left[node] : refaf_hgb_right[node];",
        "            node = begin + local_next;",
        "        }",
        "        pred += refaf_hgb_value[node];",
        "    }",
        "    return clamp01_refaf_hgb(pred);",
        "}",
        "",
        "static inline double refaf_hgb_predict_ani(double raw_ani, double ref_af, unsigned int xny,",
        "                                           unsigned int n_diff_obj, unsigned int n_diff_obj_section,",
        "                                           unsigned int n_mut2_ctx)",
        "{",
        "    if (xny == 0) return raw_ani;",
        "    raw_ani = clamp01_refaf_hgb(raw_ani);",
        "    if (ref_af < 1e-12) ref_af = 1e-12;",
        "    if (ref_af > 1.0) ref_af = 1.0;",
        "    const double xny_d = (double)xny;",
        "    const double xny_model_d = xny_d * ani_model_fold_scale();",
        "    const double ref_aaf_ani = clamp01_refaf_hgb(1.0 + log(ref_af) / 11.0);",
        "    double x[REFAF_HGB_N_FEATURES];",
        "    x[0] = raw_ani;",
        "    x[1] = ref_af;",
        "    x[2] = log1p(xny_model_d);",
        "    x[3] = (double)n_diff_obj / xny_d;",
        "    x[4] = (double)n_diff_obj_section / xny_d;",
        "    x[5] = (double)n_mut2_ctx / xny_d;",
        "    x[6] = ref_aaf_ani;",
        "    x[7] = raw_ani - ref_aaf_ani;",
        "    x[8] = ref_af < 0.2 ? 1.0 : 0.0;",
        "    x[9] = (ref_af >= 0.2 && ref_af < 0.5) ? 1.0 : 0.0;",
        "    x[10] = ref_af >= 0.5 ? 1.0 : 0.0;",
        "    return refaf_hgb_predict_from_values(x);",
        "}",
        "#else",
        "static inline double refaf_hgb_predict_ani(double raw_ani, double ref_af, unsigned int xny,",
        "                                           unsigned int n_diff_obj, unsigned int n_diff_obj_section,",
        "                                           unsigned int n_mut2_ctx)",
        "{",
        "    (void)ref_af; (void)xny; (void)n_diff_obj; (void)n_diff_obj_section; (void)n_mut2_ctx;",
        "    return raw_ani;",
        "}",
        "#endif",
        "",
        "#endif /* MODEL_REFAF_HGB_H */",
        "",
    ]
    out_path.write_text("\n".join(lines))


def train_and_evaluate(features_path: Path, outdir: Path, seed: int) -> None:
    rows = read_feature_rows(features_path)
    if len(rows) < 12:
        raise RuntimeError(f"need at least 12 successful feature rows, saw {len(rows)}")
    y = np.array([as_float(row, "anim_truth") for row in rows], dtype=float)
    bins = np.array([row["ani_bin"] for row in rows])
    raw = np.array([as_float(row, "minco_ANI") for row in rows], dtype=float)
    kssd3_t10 = np.array([as_float(row, "kssd3_t10_raw_ani") for row in rows], dtype=float)
    matrix_rows = [feature_vector(row) for row in rows]
    X = np.array([[fv[name] for name in DERIVED_FEATURES] for fv in matrix_rows], dtype=float)
    dropin_rows = [dropin_hgb_feature_vector(row) for row in rows]
    X_dropin = np.array([[fv[name] for name in DROPIN_HGB_FEATURES] for fv in dropin_rows], dtype=float)
    indices = np.arange(len(rows))
    stratify = bins if min(np.sum(bins == label) for label in set(bins)) >= 2 else None
    train_idx, test_idx = train_test_split(
        indices,
        test_size=0.30,
        random_state=seed,
        stratify=stratify,
    )

    numeric_preprocess = ColumnTransformer(
        transformers=[
            (
                "num",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                list(range(X.shape[1])),
            )
        ],
        remainder="drop",
    )
    ridge = Pipeline(
        steps=[
            ("preprocess", numeric_preprocess),
            ("model", RidgeCV(alphas=(0.001, 0.01, 0.1, 1.0, 10.0, 100.0))),
        ]
    )
    hgb = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            (
                "model",
                HistGradientBoostingRegressor(
                    max_iter=300,
                    learning_rate=0.05,
                    max_leaf_nodes=31,
                    l2_regularization=0.01,
                    random_state=seed,
                ),
            ),
        ]
    )
    hgb_dropin = HistGradientBoostingRegressor(
        max_iter=160,
        learning_rate=0.05,
        max_leaf_nodes=31,
        l2_regularization=0.01,
        random_state=seed,
    )

    models = {
        "minco_raw_ctxmoe": None,
        "kssd3_t10_raw_reference": None,
        "moe_style_ridge_full": ridge,
        "hgb_full": hgb,
        "hgb_refaf11_dropin": hgb_dropin,
    }
    predictions: Dict[str, np.ndarray] = {}
    for name, model in models.items():
        if name == "minco_raw_ctxmoe":
            predictions[name] = raw.copy()
        elif name == "kssd3_t10_raw_reference":
            predictions[name] = kssd3_t10.copy()
        elif name == "hgb_refaf11_dropin":
            model.fit(X_dropin[train_idx], y[train_idx])
            predictions[name] = np.clip(model.predict(X_dropin), 0.0, 1.0)
        else:
            model.fit(X[train_idx], y[train_idx])
            predictions[name] = np.clip(model.predict(X), 0.0, 1.0)

    metrics_path = outdir / "model_metrics.tsv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "model",
                "split",
                "n",
                "bias",
                "mae",
                "rmse",
                "max_abs",
                "pearson_r",
                "tp95",
                "fp95",
                "fn95",
                "tn95",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        split_masks = {
            "train": np.isin(indices, train_idx),
            "test": np.isin(indices, test_idx),
            "all": np.ones(len(indices), dtype=bool),
        }
        for model_name, pred in predictions.items():
            for split_name, mask in split_masks.items():
                row = {"model": model_name, "split": split_name}
                row.update(metrics(y[mask], pred[mask]))
                writer.writerow(row)

    by_bin_path = outdir / "model_metrics_by_bin.tsv"
    with by_bin_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("model", "split", "ani_bin", "n", "bias", "mae", "rmse", "max_abs", "pearson_r"),
            delimiter="\t",
        )
        writer.writeheader()
        for model_name, pred in predictions.items():
            for split_name, split_idx in (("train", train_idx), ("test", test_idx), ("all", indices)):
                for label in sorted(set(bins)):
                    mask = np.isin(indices, split_idx) & (bins == label)
                    if not np.any(mask):
                        continue
                    row = {"model": model_name, "split": split_name, "ani_bin": label}
                    row.update(metrics(y[mask], pred[mask]))
                    writer.writerow({key: row[key] for key in writer.fieldnames})

    pred_path = outdir / "model_predictions.tsv"
    split_lookup = {idx: "train" for idx in train_idx}
    split_lookup.update({idx: "test" for idx in test_idx})
    with pred_path.open("w", newline="") as handle:
        fieldnames = (
            "pair_id",
            "split",
            "ani_bin",
            "anim_truth",
            "minco_raw_ctxmoe",
            "kssd3_t10_raw_reference",
            "moe_style_ridge_full",
            "hgb_full",
            "hgb_refaf11_dropin",
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, row in enumerate(rows):
            writer.writerow(
                {
                    "pair_id": row["pair_id"],
                    "split": split_lookup[idx],
                    "ani_bin": row["ani_bin"],
                    "anim_truth": y[idx],
                    "minco_raw_ctxmoe": predictions["minco_raw_ctxmoe"][idx],
                    "kssd3_t10_raw_reference": predictions["kssd3_t10_raw_reference"][idx],
                    "moe_style_ridge_full": predictions["moe_style_ridge_full"][idx],
                    "hgb_full": predictions["hgb_full"][idx],
                    "hgb_refaf11_dropin": predictions["hgb_refaf11_dropin"][idx],
                }
            )

    ridge_model = ridge.named_steps["model"]
    coef_path = outdir / "ridge_coefficients.tsv"
    with coef_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("feature", "coefficient"))
        for name, coef in zip(DERIVED_FEATURES, ridge_model.coef_):
            writer.writerow((name, coef))
        writer.writerow(("intercept", ridge_model.intercept_))
        writer.writerow(("alpha", ridge_model.alpha_))

    export_refaf_hgb_header(
        hgb_dropin,
        outdir / "model_refaf_hgb.minco_generated.h",
        note=f"Trained on {len(train_idx)} stratified Nayfach pairs, tested on {len(test_idx)}; sketch-size 10000.",
    )


def write_counts(counts: Dict[str, int], out_path: Path) -> None:
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("ani_bin", "available_pairs"))
        for label, _, _ in ANI_BINS:
            writer.writerow((label, counts.get(label, 0)))
        writer.writerow(("total_rows", counts.get("total_rows", 0)))


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--anim-table",
        type=Path,
        default=Path("/mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv"),
    )
    parser.add_argument(
        "--fna-root",
        type=Path,
        default=Path("/mnt/new3T/skani_data/Nayfach_data/fna"),
    )
    parser.add_argument("--minco-bin", type=Path, default=Path("./bin/minco"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sketch-size", type=int, default=10000)
    parser.add_argument("--per-bin", type=int, default=200)
    parser.add_argument("--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--minco-threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument(
        "--steps",
        default="sample,features,train",
        help="Comma-separated subset of: sample, features, train.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    sample_path = args.outdir / "sample_pairs.tsv"
    counts_path = args.outdir / "sample_bin_counts.tsv"
    features_path = args.outdir / "minco_features.tsv"
    tmp_dir = args.outdir / "tmp"
    steps = {step.strip() for step in args.steps.split(",") if step.strip()}

    if "sample" in steps:
        pairs, counts = sample_pairs(args.anim_table, args.fna_root, args.per_bin, args.seed)
        write_sample_pairs(pairs, sample_path)
        write_counts(counts, counts_path)
        print(f"sampled {len(pairs)} pairs -> {sample_path}", file=sys.stderr)
    else:
        pairs = read_sample_pairs(sample_path)
        print(f"loaded {len(pairs)} pairs from {sample_path}", file=sys.stderr)

    if "features" in steps:
        if not args.minco_bin.exists():
            raise FileNotFoundError(args.minco_bin)
        extract_features(
            pairs=pairs,
            minco_bin=args.minco_bin,
            sketch_size=args.sketch_size,
            jobs=args.jobs,
            minco_threads=args.minco_threads,
            timeout=args.timeout,
            out_path=features_path,
            tmp_dir=tmp_dir,
        )
        print(f"features -> {features_path}", file=sys.stderr)

    if "train" in steps:
        train_and_evaluate(features_path, args.outdir, args.seed)
        print(f"model metrics -> {args.outdir / 'model_metrics.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
