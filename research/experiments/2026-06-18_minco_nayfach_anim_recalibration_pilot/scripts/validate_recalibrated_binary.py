#!/usr/bin/env python3
"""Validate compiled minco Recalibrated ANI against exported-model predictions."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error


def read_pairs(path: Path) -> Dict[str, Dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["pair_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def selected_prediction_rows(path: Path, per_bin: int) -> List[Dict[str, str]]:
    by_bin: Dict[str, List[Dict[str, str]]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("split") != "test":
                continue
            by_bin.setdefault(row["ani_bin"], []).append(row)
    rows: List[Dict[str, str]] = []
    for label in sorted(by_bin):
        rows.extend(by_bin[label][:per_bin])
    return rows


def run_pair(minco_bin: Path, pair: Dict[str, str], sketch_size: int, slmetrics: str, tmp_dir: Path, timeout: int) -> Dict[str, str]:
    start = time.perf_counter()
    with tempfile.NamedTemporaryFile(prefix=f"{pair['pair_id']}.", suffix=".tsv", dir=tmp_dir, delete=False) as tmp:
        out_path = Path(tmp.name)
    cmd = [
        str(minco_bin),
        "ani",
        "-S",
        str(sketch_size),
        "-p1",
        "-f0",
        "-n0",
        "-t0",
        f"-s{slmetrics}",
        "-o",
        str(out_path),
        pair["ref_fna"],
        pair["qry_fna"],
    ]
    out = dict(pair)
    try:
        completed = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
        out["seconds"] = f"{time.perf_counter() - start:.6f}"
        out["stderr_excerpt"] = " ".join(completed.stderr.split())[:200]
        if completed.returncode != 0:
            out["status"] = f"exit_{completed.returncode}"
            return out
        with out_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            out["status"] = "empty"
            return out
        row = rows[0]
        out["compiled_ani"] = row["ANI"]
        out["compiled_distance"] = row["Distance"]
        out["compiled_confidence"] = row["Confidence"]
        out["compiled_metric"] = row["Selected_metric"]
        out["status"] = "ok"
        return out
    except Exception as exc:  # noqa: BLE001 - validation should continue after a failed pair.
        out["status"] = "exception"
        out["stderr_excerpt"] = str(exc)[:200]
        out["seconds"] = f"{time.perf_counter() - start:.6f}"
        return out
    finally:
        try:
            out_path.unlink()
        except FileNotFoundError:
            pass


def metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    err = y_pred - y_true
    return {
        "n": float(len(y_true)),
        "bias": float(err.mean()),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "max_abs": float(np.max(np.abs(err))) if len(err) else float("nan"),
    }


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-pairs", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--minco-bin", type=Path, default=Path("./bin/minco"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--per-bin", type=int, default=10)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--sketch-size", type=int, default=10000)
    parser.add_argument("--slmetrics", default="2")
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args(argv)

    pairs = read_pairs(args.sample_pairs)
    pred_rows = selected_prediction_rows(args.predictions, args.per_bin)
    validation_pairs: List[Dict[str, str]] = []
    for pred in pred_rows:
        pair = dict(pairs[pred["pair_id"]])
        pair["split"] = pred["split"]
        pair["python_hgb_refaf11_dropin"] = pred["hgb_refaf11_dropin"]
        validation_pairs.append(pair)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = args.out.parent / "tmp_binary_validation"
    tmp_dir.mkdir(exist_ok=True)
    results: List[Dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(run_pair, args.minco_bin, pair, args.sketch_size, args.slmetrics, tmp_dir, args.timeout): pair
            for pair in validation_pairs
        }
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda row: (row["ani_bin"], row["pair_id"]))

    fieldnames = [
        "pair_id",
        "ani_bin",
        "split",
        "anim_truth",
        "python_hgb_refaf11_dropin",
        "compiled_ani",
        "compiled_distance",
        "compiled_confidence",
        "compiled_metric",
        "status",
        "seconds",
        "stderr_excerpt",
        "ref_fna",
        "qry_fna",
    ]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    ok = [row for row in results if row.get("status") == "ok"]
    y = np.array([float(row["anim_truth"]) for row in ok])
    c = np.array([float(row["compiled_ani"]) for row in ok])
    p = np.array([float(row["python_hgb_refaf11_dropin"]) for row in ok])
    metric_path = args.out.with_suffix(".metrics.tsv")
    with metric_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("comparison", "n", "bias", "mae", "rmse", "max_abs"), delimiter="\t")
        writer.writeheader()
        row = {"comparison": "compiled_vs_anim"}
        row.update(metrics(y, c))
        writer.writerow(row)
        row = {"comparison": "compiled_vs_python_dropin"}
        row.update(metrics(p, c))
        writer.writerow(row)
    print(f"validated {len(ok)}/{len(results)} pairs -> {args.out}")
    print(f"metrics -> {metric_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
