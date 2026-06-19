#!/usr/bin/env python3
"""Build a minco VEuPathDB sketch from the existing VEuPathDB manifest."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


REQUIRED_PART_FILES = (
    "minco.ctxobj64",
    "minco.ctxobj64.offsets",
    "minco.stat",
    "minco.anno",
    "minco.infilemeta",
    "minco.ctxmeta.tsv",
    "minco.ctxsetmeta.tsv",
)

REQUIRED_FINAL_FILES = (
    "minco.ctxobj64",
    "minco.ctxobj64.offsets",
    "minco.stat",
    "minco.anno",
    "minco.infilemeta",
    "minco.ctxmeta.tsv",
    "minco.ctxsetmeta.tsv",
)

FILE_ALIASES = {
    "minco.ctxobj64": ("minco.ctxobj64", "comblco"),
    "minco.ctxobj64.offsets": ("minco.ctxobj64.offsets", "comblco.index"),
    "minco.stat": ("minco.stat", "lcofiles.stat"),
    "minco.anno": ("minco.anno", "lcofiles.anno"),
    "minco.infilemeta": ("minco.infilemeta", "lcofiles.infilemeta"),
    "minco.qc": ("minco.qc", "lcofiles.qc"),
    "minco.ctxobj64.abund": ("minco.ctxobj64.abund", "comblco.a"),
    "minco.ctxobj64.position": ("minco.ctxobj64.position", "comblco.position"),
    "minco.refindex.ctxgid64obj32": (
        "minco.refindex.ctxgid64obj32",
        "sortedcomb_ctxgid64obj32",
    ),
}


@dataclass(frozen=True)
class Genome:
    index: int
    site: str
    organism: str
    filename: str
    expected_bytes: int
    url: str

    @property
    def fasta_relpath(self) -> Path:
        return Path("reference_genomes") / self.site / self.organism / self.filename

    def part_dir(self, parts_dir: Path) -> Path:
        return parts_dir / self.site / self.organism


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def run(cmd: list[str], cwd: Path) -> None:
    log("+ " + " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def read_manifest(path: Path) -> list[Genome]:
    genomes: list[Genome] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(reader, 1):
            genomes.append(
                Genome(
                    index=index,
                    site=row["site"],
                    organism=row["organism"],
                    filename=row["filename"],
                    expected_bytes=int(row["bytes"]),
                    url=row["url"],
                )
            )
    return genomes


def apply_order_list(genomes: list[Genome], order_list: Path | None) -> list[Genome]:
    if order_list is None:
        return genomes
    by_path = {str(genome.fasta_relpath): genome for genome in genomes}
    ordered: list[Genome] = []
    missing: list[str] = []
    with order_list.open() as handle:
        for raw in handle:
            path = raw.strip()
            if not path:
                continue
            genome = by_path.get(path)
            if genome is None:
                missing.append(path)
            else:
                ordered.append(genome)
    if missing:
        preview = ", ".join(missing[:5])
        raise RuntimeError(f"{len(missing)} order-list entries are missing from manifest: {preview}")
    if not ordered:
        raise RuntimeError(f"order list is empty: {order_list}")
    return ordered


def has_required_file(path: Path, name: str) -> bool:
    for alias in FILE_ALIASES.get(name, (name,)):
        candidate = path / alias
        if candidate.is_file() and candidate.stat().st_size > 0:
            return True
    return False


def complete_dir(path: Path, required: tuple[str, ...]) -> bool:
    return all(has_required_file(path, name) for name in required)


def remove_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    else:
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def prune_empty_parents(path: Path, stop: Path) -> None:
    stop = stop.resolve()
    current = path.resolve()
    while current != stop and current.parent != current:
        try:
            current.rmdir()
        except OSError:
            return
        current = current.parent


def download_one(genome: Genome, root: Path) -> Path:
    fasta = root / genome.fasta_relpath
    tmp = fasta.with_suffix(fasta.suffix + ".part")
    fasta.parent.mkdir(parents=True, exist_ok=True)
    remove_path(tmp)
    remove_path(fasta)
    run(
        [
            "curl",
            "-fL",
            "--retry",
            "5",
            "--retry-delay",
            "10",
            "--retry-all-errors",
            "--connect-timeout",
            "30",
            "-o",
            str(tmp),
            genome.url,
        ],
        cwd=root,
    )
    size = tmp.stat().st_size
    if size != genome.expected_bytes:
        remove_path(tmp)
        raise RuntimeError(
            f"downloaded size mismatch for {genome.url}: got {size}, expected {genome.expected_bytes}"
        )
    tmp.replace(fasta)
    return fasta


def sketch_one(args: argparse.Namespace, genome: Genome) -> None:
    part_dir = genome.part_dir(args.parts_dir)
    if complete_dir(part_dir, REQUIRED_PART_FILES):
        log(f"SKIP {genome.index}: {part_dir}")
        return

    remove_path(part_dir)
    fasta: Path | None = None
    try:
        log(f"DOWNLOAD {genome.index}: {genome.site}/{genome.organism}")
        fasta = download_one(genome, args.root)
        log(f"SKETCH {genome.index}: {genome.fasta_relpath}")
        run(
            [
                str(args.minco),
                "sketch",
                "-T",
                "-S",
                str(args.sketch_size),
                "--anno",
                "--ctxmeta",
                "both",
                "-p",
                str(args.threads),
                "-o",
                str(part_dir),
                str(genome.fasta_relpath),
            ],
            cwd=args.root,
        )
        if not complete_dir(part_dir, REQUIRED_PART_FILES):
            raise RuntimeError(f"incomplete minco part sketch: {part_dir}")
        log(f"DONE {genome.index}: {part_dir}")
    finally:
        if fasta is not None:
            remove_path(fasta)
            prune_empty_parents(fasta.parent, args.root / "reference_genomes")


def merge_parts(args: argparse.Namespace, genomes: list[Genome]) -> None:
    if complete_dir(args.final_dir, REQUIRED_FINAL_FILES):
        log(f"SKIP merge: {args.final_dir}")
        return

    remove_path(args.final_dir)
    part_dirs = [str(genome.part_dir(args.parts_dir)) for genome in genomes]
    log(f"MERGE {len(part_dirs)} sketches into {args.final_dir}")
    run([str(args.minco), "sketch", "--append", "-o", str(args.final_dir), *part_dirs], cwd=args.root)
    if not complete_dir(args.final_dir, REQUIRED_FINAL_FILES):
        raise RuntimeError(f"incomplete merged minco sketch: {args.final_dir}")
    log(f"DONE merge: {args.final_dir}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/mnt/new3T/VEuPathDB"))
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("/mnt/new3T/VEuPathDB/metadata/veupathdb_reference_genomes_manifest.tsv"),
    )
    parser.add_argument(
        "--order-list",
        type=Path,
        default=Path("/mnt/new3T/VEuPathDB/genomepath.list"),
        help="Relative FASTA path list defining output sample order; use 'none' to keep manifest order.",
    )
    parser.add_argument(
        "--minco",
        type=Path,
        default=Path("/home/ubuntu/yihuiguang/tools/KSSD3mini/bin/minco"),
    )
    parser.add_argument(
        "--parts-dir",
        type=Path,
        default=Path("/mnt/new3T/VEuPathDB/VEuPathDB_minco_T_S10000_anno_parts"),
    )
    parser.add_argument(
        "--final-dir",
        type=Path,
        default=Path("/mnt/new3T/VEuPathDB/VEuPathDB_minco_T_S10000_anno"),
    )
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of genomes to download/sketch concurrently; each job uses --threads minco threads.",
    )
    parser.add_argument("--sketch-size", type=int, default=10000)
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--no-merge", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    if str(args.order_list).lower() == "none":
        args.order_list = None
    genomes = read_manifest(args.manifest)
    genomes = apply_order_list(genomes, args.order_list)
    if args.start < 1:
        raise RuntimeError("--start must be >= 1")
    if args.jobs < 1:
        raise RuntimeError("--jobs must be >= 1")
    selected = genomes[args.start - 1 :]
    if args.limit is not None:
        selected = selected[: args.limit]
    if not selected:
        raise RuntimeError("no genomes selected")
    log(
        f"selected={len(selected)} start={args.start} sketch_size={args.sketch_size} "
        f"threads={args.threads} jobs={args.jobs} parts={args.parts_dir} final={args.final_dir}"
    )
    if args.jobs == 1:
        for genome in selected:
            sketch_one(args, genome)
    else:
        failures: list[str] = []
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_to_genome = {executor.submit(sketch_one, args, genome): genome for genome in selected}
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001 - report all worker failures after shutdown.
                    failures.append(f"{genome.index}:{genome.site}/{genome.organism}: {exc}")
                    log(f"ERROR {genome.index}: {exc}")
        if failures:
            preview = "\n".join(failures[:10])
            raise RuntimeError(f"{len(failures)} genome sketches failed:\n{preview}")
    if not args.no_merge:
        merge_parts(args, selected)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
