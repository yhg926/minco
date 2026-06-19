#!/usr/bin/env python3
"""Compare Mash -s 10000 ANI against ANIm and minco on held-out Nayfach pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def select_pairs(sample_pairs: Path, predictions: Path, per_bin: int) -> List[Dict[str, str]]:
    pair_by_id = {row["pair_id"]: row for row in read_tsv(sample_pairs)}
    by_bin: Dict[str, List[Dict[str, str]]] = {}
    for pred in read_tsv(predictions):
        if pred["split"] != "test":
            continue
        row = dict(pair_by_id[pred["pair_id"]])
        row["split"] = pred["split"]
        row["minco_raw_ctxmoe"] = pred["minco_raw_ctxmoe"]
        row["minco_best_hgb"] = pred["hgb_refaf11_dropin"]
        by_bin.setdefault(row["ani_bin"], []).append(row)
    selected: List[Dict[str, str]] = []
    for label in sorted(by_bin):
        selected.extend(by_bin[label][:per_bin])
    selected.sort(key=lambda row: (row["ani_bin"], row["pair_id"]))
    return selected


def sketch_prefix(path: str, sketch_dir: Path) -> Path:
    digest = hashlib.sha1(path.encode("utf-8")).hexdigest()[:16]
    stem = Path(path).name
    return sketch_dir / f"{digest}.{stem}"


def sketch_path(path: str, sketch_dir: Path) -> Path:
    return sketch_prefix(path, sketch_dir).with_suffix(sketch_prefix(path, sketch_dir).suffix + ".msh")


def run_mash_sketch(mash_bin: str, genome: str, sketch_dir: Path, sketch_size: int, timeout: int) -> Dict[str, str]:
    out_path = sketch_path(genome, sketch_dir)
    if out_path.exists() and out_path.stat().st_size > 0:
        return {"genome": genome, "sketch": str(out_path), "status": "cached", "seconds": "0"}
    prefix = sketch_prefix(genome, sketch_dir)
    start = time.perf_counter()
    cmd = [mash_bin, "sketch", "-s", str(sketch_size), "-o", str(prefix), genome]
    completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    seconds = time.perf_counter() - start
    status = "ok" if completed.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0 else f"exit_{completed.returncode}"
    return {
        "genome": genome,
        "sketch": str(out_path),
        "status": status,
        "seconds": f"{seconds:.6f}",
        "stderr_excerpt": " ".join(completed.stderr.split())[:240],
    }


def run_mash_dist(mash_bin: str, pair: Dict[str, str], sketch_dir: Path, timeout: int) -> Dict[str, object]:
    ref_msh = sketch_path(pair["ref_fna"], sketch_dir)
    qry_msh = sketch_path(pair["qry_fna"], sketch_dir)
    out = dict(pair)
    start = time.perf_counter()
    cmd = [mash_bin, "dist", str(ref_msh), str(qry_msh)]
    try:
        completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
        out["mash_seconds"] = f"{time.perf_counter() - start:.6f}"
        out["mash_stderr_excerpt"] = " ".join(completed.stderr.split())[:240]
        if completed.returncode != 0:
            out["status"] = f"dist_exit_{completed.returncode}"
            return out
        line = completed.stdout.strip().splitlines()[0]
        fields = line.split("\t")
        out["mash_ref"] = fields[0]
        out["mash_qry"] = fields[1]
        out["mash_distance"] = fields[2]
        out["mash_pvalue"] = fields[3] if len(fields) > 3 else ""
        out["mash_shared_hashes"] = fields[4] if len(fields) > 4 else ""
        out["mash_ani"] = f"{1.0 - float(fields[2]):.12g}"
        out["status"] = "ok"
    except Exception as exc:  # noqa: BLE001 - keep pair-level failure visible.
        out["mash_seconds"] = f"{time.perf_counter() - start:.6f}"
        out["mash_stderr_excerpt"] = str(exc)[:240]
        out["status"] = "exception"
    return out


def pearson_r(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 2:
        return float("nan")
    y0 = y_true - y_true.mean()
    y1 = y_pred - y_pred.mean()
    den = math.sqrt(float((y0 * y0).sum() * (y1 * y1).sum()))
    return float((y0 * y1).sum() / den) if den > 0 else float("nan")


def metric_row(rows: List[Dict[str, object]], model: str, pred_key: str, split: str, ani_bin: str = "all") -> Dict[str, object]:
    y = np.array([float(row["anim_truth"]) for row in rows], dtype=float)
    p = np.array([float(row[pred_key]) for row in rows], dtype=float)
    err = p - y
    truth_pos = y >= 0.95
    pred_pos = p >= 0.95
    return {
        "model": model,
        "split": split,
        "ani_bin": ani_bin,
        "n": len(rows),
        "bias": f"{float(err.mean()):.12g}",
        "mae": f"{float(np.mean(np.abs(err))):.12g}",
        "rmse": f"{float(math.sqrt(np.mean(err * err))):.12g}",
        "max_abs": f"{float(np.max(np.abs(err))):.12g}",
        "pearson_r": f"{pearson_r(y, p):.12g}",
        "tp95": int(np.sum(truth_pos & pred_pos)),
        "fp95": int(np.sum(~truth_pos & pred_pos)),
        "fn95": int(np.sum(truth_pos & ~pred_pos)),
        "tn95": int(np.sum(~truth_pos & ~pred_pos)),
    }


def write_metrics(pair_rows: List[Dict[str, object]], out_path: Path) -> None:
    ok = [row for row in pair_rows if row.get("status") == "ok"]
    models = (
        ("mash_ani_1_minus_dist", "mash_ani"),
        ("minco_raw_ctxmoe", "minco_raw_ctxmoe"),
        ("minco_best_hgb", "minco_best_hgb"),
    )
    rows: List[Dict[str, object]] = []
    for model, key in models:
        rows.append(metric_row(ok, model, key, "all", "all"))
        for label in sorted({row["ani_bin"] for row in ok}):
            sub = [row for row in ok if row["ani_bin"] == label]
            rows.append(metric_row(sub, model, key, "all", label))
    write_tsv(
        out_path,
        rows,
        ("model", "split", "ani_bin", "n", "bias", "mae", "rmse", "max_abs", "pearson_r", "tp95", "fp95", "fn95", "tn95"),
    )


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-pairs", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mash-bin", default="mash")
    parser.add_argument("--sketch-size", type=int, default=10000)
    parser.add_argument("--per-bin", type=int, default=100)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    sketch_dir = args.outdir / "mash_sketches_s10000"
    sketch_dir.mkdir(exist_ok=True)
    selected = select_pairs(args.sample_pairs, args.predictions, args.per_bin)
    write_tsv(
        args.outdir / "comparison_pairs.selected.tsv",
        selected,
        ("pair_id", "ani_bin", "split", "ref_fna", "qry_fna", "anim_truth", "minco_raw_ctxmoe", "minco_best_hgb"),
    )
    genomes = sorted({row["ref_fna"] for row in selected} | {row["qry_fna"] for row in selected})
    print(f"selected {len(selected)} pairs, {len(genomes)} unique genomes", file=sys.stderr, flush=True)

    sketch_rows: List[Dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_mash_sketch, args.mash_bin, genome, sketch_dir, args.sketch_size, args.timeout): genome for genome in genomes}
        for i, future in enumerate(as_completed(futures), 1):
            sketch_rows.append(future.result())
            if i == 1 or i % 100 == 0 or i == len(futures):
                print(f"mash sketch {i}/{len(futures)}", file=sys.stderr, flush=True)
    write_tsv(args.outdir / "mash_sketch_status.tsv", sketch_rows, ("genome", "sketch", "status", "seconds", "stderr_excerpt"))
    bad = [row for row in sketch_rows if row["status"] not in ("ok", "cached")]
    if bad:
        raise RuntimeError(f"{len(bad)} mash sketches failed")

    dist_rows: List[Dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_mash_dist, args.mash_bin, pair, sketch_dir, args.timeout): pair for pair in selected}
        for i, future in enumerate(as_completed(futures), 1):
            dist_rows.append(future.result())
            if i == 1 or i % 100 == 0 or i == len(futures):
                print(f"mash dist {i}/{len(futures)}", file=sys.stderr, flush=True)
    dist_rows.sort(key=lambda row: (row["ani_bin"], row["pair_id"]))
    write_tsv(
        args.outdir / "mash_minco_anim_pairs.tsv",
        dist_rows,
        (
            "pair_id", "ani_bin", "split", "anim_truth", "mash_ani", "mash_distance", "mash_pvalue",
            "mash_shared_hashes", "minco_raw_ctxmoe", "minco_best_hgb", "status", "mash_seconds",
            "ref_fna", "qry_fna", "mash_stderr_excerpt",
        ),
    )
    write_metrics(dist_rows, args.outdir / "mash_minco_anim_metrics.tsv")
    print(f"pair output -> {args.outdir / 'mash_minco_anim_pairs.tsv'}", file=sys.stderr)
    print(f"metrics -> {args.outdir / 'mash_minco_anim_metrics.tsv'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
