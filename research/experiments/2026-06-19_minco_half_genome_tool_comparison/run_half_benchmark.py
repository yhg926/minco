#!/usr/bin/env python3
"""Compare ANI estimators on a full same-species pair and query-genome halves."""

from __future__ import annotations

import csv
import gzip
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


PAIR_ID = "P000881414"
REF_FASTA = Path("/mnt/new3T/skani_data/Nayfach_data/fna/3300027414_1.fna")
QRY_FASTA = Path("/mnt/new3T/skani_data/Nayfach_data/fna/3300014912_1.fna")
NAYFACH_ANIM = 0.9890599554743553


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq_parts: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts).upper()))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts).upper()))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def write_single_fasta(path: Path, header: str, seq: str) -> None:
    with path.open("w") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            handle.write(seq[i : i + 80] + "\n")


def write_split_halves(source: Path, half1: Path, half2: Path) -> dict[str, int]:
    records = read_fasta(source)
    joined = "".join(seq for _, seq in records)
    midpoint = len(joined) // 2
    write_single_fasta(half1, f"{source.stem}_first_half", joined[:midpoint])
    write_single_fasta(half2, f"{source.stem}_second_half", joined[midpoint:])
    return {
        "source_records": len(records),
        "source_bases": len(joined),
        "half1_bases": midpoint,
        "half2_bases": len(joined) - midpoint,
    }


def run(cmd: list[str], cwd: Path, stdout_path: Path | None = None, stderr_path: Path | None = None) -> tuple[float, str, str]:
    start = time.perf_counter()
    completed = subprocess.run(
        cmd,
        cwd=str(cwd),
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    elapsed = time.perf_counter() - start
    if stdout_path is not None:
        stdout_path.write_text(completed.stdout)
    if stderr_path is not None:
        stderr_path.write_text(completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
        )
    return elapsed, completed.stdout, completed.stderr


def parse_minco(path: Path) -> dict[str, float | str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"expected one minco row in {path}, got {len(rows)}")
    row = rows[0]
    return {
        "ani": float(row["ANI"]),
        "distance": float(row["Distance"]),
        "xny": float(row["XnY_ctx"]),
        "qry_af": float(row["Qry_align_fraction"]),
        "ref_af": float(row["Ref_align_fraction"]),
        "confidence": row["Confidence"],
        "metric": row["Selected_metric"],
    }


def parse_skani(stdout: str) -> dict[str, float]:
    rows = [line for line in stdout.splitlines() if line and not line.startswith("[")]
    reader = csv.DictReader(rows, delimiter="\t")
    parsed = list(reader)
    if len(parsed) != 1:
        raise ValueError(f"expected one skani row, got {len(parsed)}: {stdout}")
    row = parsed[0]
    return {
        "ani": float(row["ANI"]) / 100.0,
        "ref_af": float(row["Align_fraction_ref"]) / 100.0,
        "qry_af": float(row["Align_fraction_query"]) / 100.0,
    }


def parse_mash(stdout: str) -> dict[str, float | str]:
    lines = [line for line in stdout.splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"expected one mash row, got {len(lines)}: {stdout}")
    fields = lines[0].split("\t")
    dist = float(fields[2])
    return {
        "ani": 1.0 - dist,
        "distance": dist,
        "pvalue": fields[3],
        "shared_hashes": fields[4],
    }


def parse_dnadiff_report(path: Path) -> dict[str, float]:
    text = path.read_text()
    aln = re.search(r"\[Alignments\](.*?)(?:\n\[|\Z)", text, flags=re.S)
    if not aln:
        raise ValueError(f"could not parse [Alignments] block from {path}")
    avg_identities = re.findall(r"AvgIdentity\s+([0-9.]+)\s+([0-9.]+)", aln.group(1))
    if not avg_identities:
        raise ValueError(f"could not parse AvgIdentity from {path}")
    # dnadiff reports 1-to-1 first, then M-to-M.
    ref_identity, qry_identity = avg_identities[0]
    aligned = re.search(
        r"AlignedBases\s+(\d+)\(([0-9.]+)%\)\s+(\d+)\(([0-9.]+)%\)",
        text,
    )
    out = {
        "ani": float(ref_identity) / 100.0,
        "ani_qry": float(qry_identity) / 100.0,
    }
    if aligned:
        out.update(
            {
                "ref_aligned_bases": float(aligned.group(1)),
                "ref_af": float(aligned.group(2)) / 100.0,
                "qry_aligned_bases": float(aligned.group(3)),
                "qry_af": float(aligned.group(4)) / 100.0,
            }
        )
    return out


def parse_fastani(path: Path) -> dict[str, float | str]:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    if len(lines) != 1:
        raise ValueError(f"expected one fastANI row in {path}, got {len(lines)}")
    fields = lines[0].split("\t")
    return {
        "ani": float(fields[2]) / 100.0,
        "mapped_fragments": float(fields[3]),
        "total_fragments": float(fields[4]),
        "fragment_fraction": float(fields[3]) / float(fields[4]) if float(fields[4]) else 0.0,
    }


def main() -> int:
    root = repo_root()
    out = Path(__file__).resolve().parent
    work = Path("/tmp/minco_half_genome_tool_comparison")
    tmp = work / "tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    raw = work / "raw_outputs"
    raw.mkdir(parents=True, exist_ok=True)

    half1 = tmp / "3300014912_1.first_half.fna"
    half2 = tmp / "3300014912_1.second_half.fna"
    split_meta = write_split_halves(QRY_FASTA, half1, half2)

    comparisons = [
        ("full_vs_full", QRY_FASTA, "complete query genome"),
        ("full_vs_half1", half1, "first half of query genome"),
        ("full_vs_half2", half2, "second half of query genome"),
    ]

    rows: list[dict[str, str | float]] = []
    tool_times: list[dict[str, str | float]] = []

    def display_path(path: Path) -> str:
        try:
            return str(path.resolve().relative_to(root))
        except ValueError:
            return str(path)

    for label, query, description in comparisons:
        dnadiff_prefix = raw / f"{label}.dnadiff"
        elapsed, _, _ = run(
            ["dnadiff", "-p", dnadiff_prefix.name, str(REF_FASTA), str(query)],
            cwd=raw,
            stdout_path=raw / f"{label}.dnadiff.stdout",
            stderr_path=raw / f"{label}.dnadiff.stderr",
        )
        tool_times.append({"comparison": label, "tool": "dnadiff", "seconds": elapsed})
        dnadiff = parse_dnadiff_report(raw / f"{label}.dnadiff.report")

        minco_best_out = raw / f"{label}.minco_best.tsv"
        elapsed, _, _ = run(
            [
                str(root / "bin/minco"),
                "ani",
                "-S",
                "10000",
                "-p2",
                "-f0",
                "-n0",
                "-t0",
                "-s1",
                "-o",
                str(minco_best_out),
                str(REF_FASTA),
                str(query),
            ],
            cwd=root,
            stderr_path=raw / f"{label}.minco_best.stderr",
        )
        tool_times.append({"comparison": label, "tool": "minco_best", "seconds": elapsed})
        minco_best = parse_minco(minco_best_out)

        minco_ctxmash_out = raw / f"{label}.minco_ctxmash.tsv"
        elapsed, _, _ = run(
            [
                str(root / "bin/minco"),
                "ani",
                "-S",
                "10000",
                "-p2",
                "-f0",
                "-n0",
                "-t0",
                "-s-5",
                "-o",
                str(minco_ctxmash_out),
                str(REF_FASTA),
                str(query),
            ],
            cwd=root,
            stderr_path=raw / f"{label}.minco_ctxmash.stderr",
        )
        tool_times.append({"comparison": label, "tool": "minco_ctxmash", "seconds": elapsed})
        minco_ctxmash = parse_minco(minco_ctxmash_out)

        elapsed, skani_stdout, skani_stderr = run(
            ["skani", "dist", "-t2", str(REF_FASTA), str(query)],
            cwd=raw,
            stdout_path=raw / f"{label}.skani.tsv",
            stderr_path=raw / f"{label}.skani.stderr",
        )
        tool_times.append({"comparison": label, "tool": "skani", "seconds": elapsed})
        skani = parse_skani(skani_stdout)

        mash_prefix = raw / f"{label}.ref_mash"
        for suffix in ["msh"]:
            path = Path(f"{mash_prefix}.{suffix}")
            if path.exists():
                path.unlink()
        elapsed, _, _ = run(
            ["mash", "sketch", "-s", "10000", "-k", "21", "-o", str(mash_prefix), str(REF_FASTA)],
            cwd=raw,
            stdout_path=raw / f"{label}.mash_sketch.stdout",
            stderr_path=raw / f"{label}.mash_sketch.stderr",
        )
        tool_times.append({"comparison": label, "tool": "mash_sketch_ref", "seconds": elapsed})
        elapsed, mash_stdout, _ = run(
            ["mash", "dist", f"{mash_prefix}.msh", str(query)],
            cwd=raw,
            stdout_path=raw / f"{label}.mash.tsv",
            stderr_path=raw / f"{label}.mash.stderr",
        )
        tool_times.append({"comparison": label, "tool": "mash_dist", "seconds": elapsed})
        mash = parse_mash(mash_stdout)

        fastani_out = raw / f"{label}.fastani.tsv"
        elapsed, _, _ = run(
            ["fastANI", "-q", str(REF_FASTA), "-r", str(query), "-o", str(fastani_out)],
            cwd=raw,
            stdout_path=raw / f"{label}.fastani.stdout",
            stderr_path=raw / f"{label}.fastani.stderr",
        )
        tool_times.append({"comparison": label, "tool": "fastANI", "seconds": elapsed})
        fastani = parse_fastani(fastani_out)

        truth = dnadiff["ani"]
        rows.append(
            {
                "pair_id": PAIR_ID,
                "comparison": label,
                "query": description,
                "ref_fasta": display_path(REF_FASTA),
                "query_fasta": display_path(query),
                "nayfach_full_anim": NAYFACH_ANIM if label == "full_vs_full" else "NA",
                "dnadiff_anim": truth,
                "dnadiff_ref_af": dnadiff.get("ref_af", "NA"),
                "dnadiff_qry_af": dnadiff.get("qry_af", "NA"),
                "fastani": fastani["ani"],
                "fastani_fragment_fraction": fastani["fragment_fraction"],
                "minco_best_ani": minco_best["ani"],
                "minco_best_error_vs_dnadiff": float(minco_best["ani"]) - truth,
                "minco_ctxmash_ani": minco_ctxmash["ani"],
                "minco_ctxmash_error_vs_dnadiff": float(minco_ctxmash["ani"]) - truth,
                "skani_ani": skani["ani"],
                "skani_error_vs_dnadiff": float(skani["ani"]) - truth,
                "mash_ani": mash["ani"],
                "mash_error_vs_dnadiff": float(mash["ani"]) - truth,
                "minco_best_xny": minco_best["xny"],
                "minco_best_qry_af": minco_best["qry_af"],
                "minco_best_ref_af": minco_best["ref_af"],
                "minco_ctxmash_xny": minco_ctxmash["xny"],
                "skani_ref_af": skani["ref_af"],
                "skani_qry_af": skani["qry_af"],
                "mash_distance": mash["distance"],
                "mash_shared_hashes": mash["shared_hashes"],
            }
        )

    with (out / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with (out / "tool_times.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=["comparison", "tool", "seconds"])
        writer.writeheader()
        writer.writerows(tool_times)

    with (out / "split_metadata.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in split_meta.items():
            writer.writerow([key, value])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
