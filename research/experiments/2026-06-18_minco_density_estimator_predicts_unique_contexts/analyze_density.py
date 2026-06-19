#!/usr/bin/env python3
import argparse
import csv
import math
import statistics
import struct
from pathlib import Path


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pearson(xs, ys):
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    if den == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i + 1
        while j < len(vals) and vals[order[j]] == vals[order[i]]:
            j += 1
        avg = (i + j - 1) / 2 + 1
        for k in range(i, j):
            out[order[k]] = avg
        i = j
    return out


def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))


def metrics(est, truth):
    errs = [e - t for e, t in zip(est, truth)]
    rel = [(e - t) / t for e, t in zip(est, truth)]
    absrel = [abs(x) for x in rel]
    return {
        "pearson": pearson(est, truth),
        "spearman": spearman(est, truth),
        "mae": sum(abs(e) for e in errs) / len(errs),
        "rmse": math.sqrt(sum(e * e for e in errs) / len(errs)),
        "mape_pct": 100 * sum(absrel) / len(absrel),
        "mean_rel_bias_pct": 100 * sum(rel) / len(rel),
        "median_abs_pct": 100 * statistics.median(absrel),
        "max_abs_pct": 100 * max(absrel),
        "min_rel_pct": 100 * min(rel),
        "max_rel_pct": 100 * max(rel),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--std-ctxmeta", required=True)
    p.add_argument("--exact-ctxmeta", required=True)
    p.add_argument("--infilemeta", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--per-sample", required=True)
    args = p.parse_args()

    std = read_tsv(args.std_ctxmeta)
    exact = read_tsv(args.exact_ctxmeta)
    raw = Path(args.infilemeta).read_bytes()
    stride = len(raw) // len(std)
    lengths = [struct.unpack_from("<Q", raw, i * stride)[0] for i in range(len(std))]

    records = []
    for s, e, length_bp in zip(std, exact, lengths):
        exact_ctx = int(e["postconflict_observed_ctx"])
        est_pre = int(s["preconflict_estimated_unique_ctx"])
        est_post = int(s["postconflict_estimated_unique_ctx"])
        records.append(
            {
                "sample_id": s["sample_id"],
                "sample_path": s["sample_path"],
                "genome_bp": length_bp,
                "std_pre_est_total_ctx": est_pre,
                "std_post_est_conflictfree_ctx": est_post,
                "exact_total_unique_ctx": exact_ctx,
                "std_threshold": int(s["threshold"]),
                "std_pre_observed_ctx": int(s["preconflict_observed_ctx"]),
                "std_post_observed_ctx": int(s["postconflict_observed_ctx"]),
                "exact_sketch_entries_ctxobj": int(e["sketch_entries"]),
                "exact_threshold": int(e["threshold"]),
                "rel_err_pre_vs_exact": (est_pre - exact_ctx) / exact_ctx,
                "rel_err_post_vs_exact": (est_post - exact_ctx) / exact_ctx,
                "ctx_per_bp": exact_ctx / length_bp,
            }
        )

    per_sample_path = Path(args.per_sample)
    with per_sample_path.open("w", newline="") as f:
        fields = list(records[0].keys())
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(records)

    truth = [r["exact_total_unique_ctx"] for r in records]
    length = [r["genome_bp"] for r in records]
    est_pre = [r["std_pre_est_total_ctx"] for r in records]
    est_post = [r["std_post_est_conflictfree_ctx"] for r in records]

    rows = []
    for comparison, est, tru in [
        ("std_pre_est_total_ctx_vs_exact_total_unique_ctx", est_pre, truth),
        ("std_post_est_conflictfree_ctx_vs_exact_total_unique_ctx", est_post, truth),
        ("genome_bp_vs_exact_total_unique_ctx", length, truth),
        ("std_pre_est_total_ctx_vs_genome_bp", est_pre, length),
    ]:
        rows.append({"comparison": comparison, **metrics(est, tru)})

    summary_path = Path(args.summary)
    fields = [
        "comparison",
        "pearson",
        "spearman",
        "mae",
        "rmse",
        "mape_pct",
        "mean_rel_bias_pct",
        "median_abs_pct",
        "max_abs_pct",
        "min_rel_pct",
        "max_rel_pct",
    ]
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
