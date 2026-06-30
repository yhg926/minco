#!/usr/bin/env python3
"""Preflight GTDB r232 Sylph database resources.

This script intentionally does not build the database. It checks local inputs,
disk headroom, completed chunked outputs, and writes the exact command needed
to reproduce same-release Sylph artifacts used by the release-readiness audit.
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
from pathlib import Path


NOTE_DIR = Path(__file__).resolve().parent
REPO_ROOT = NOTE_DIR.parents[2]
RESULTS = NOTE_DIR / "results"

SYLPH = Path("/home/ubuntu/yihuiguang/bin/sylph")
REP_LIST = Path("/mnt/new3T/gtdbr220/gtdb232/manifests/r232.reps_fullpath.list")
OUT_BASE = Path("/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1")
OUT_DB = OUT_BASE.with_suffix(".syldb")
R226_DB = Path("/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb")
THREADS = 16
C_VALUE = 200
K_VALUE = 31
CHUNK_SIZE = 10000
CHUNK_ROOT = Path("/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks")
CHUNK_LIST_DIR = CHUNK_ROOT / "lists"
CHUNK_DB_DIR = CHUNK_ROOT / "db"
CHUNK_PROFILE_INPUT_DIR = CHUNK_ROOT / "profile_inputs"
HMP_RUN = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run")
HMP_R232_OUT = Path("/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232")
HMP_SAMPLE_SKETCHES = {
    6: HMP_RUN / "sylph_sample6/airskinurogenital_sample6.nonzero.fastq.gz.sylsp",
    11: HMP_RUN / "sylph_sample11/airskinurogenital_sample11.nonzero.fastq.gz.sylsp",
}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def sylph_version() -> str:
    if not SYLPH.exists():
        return ""
    proc = subprocess.run([str(SYLPH), "--version"], check=False, text=True, capture_output=True)
    return (proc.stdout or proc.stderr).strip()


def read_paths(path: Path) -> list[Path]:
    if not path.exists():
        return []
    with path.open() as handle:
        return [Path(line.strip()) for line in handle if line.strip()]


def fmt_gib(value: float) -> str:
    return f"{value / (1024 ** 3):.3f}"


def write_script(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o700)


def main() -> int:
    paths = read_paths(REP_LIST)
    existing = []
    missing = []
    compressed_bytes = 0
    for path in paths:
        if path.exists():
            existing.append(path)
            compressed_bytes += path.stat().st_size
        else:
            missing.append(path)

    out_parent = OUT_DB.parent
    free_bytes = shutil.disk_usage(out_parent).free if out_parent.exists() else 0
    r226_bytes = R226_DB.stat().st_size if R226_DB.exists() else 0
    disk_heuristic_bytes = max(30 * 1024**3, 3 * r226_bytes)
    disk_heuristic_pass = free_bytes >= disk_heuristic_bytes
    inputs_complete = bool(paths) and not missing
    ready_to_build = SYLPH.exists() and REP_LIST.exists() and inputs_complete and disk_heuristic_pass and not OUT_DB.exists()
    chunk_count = math.ceil(len(paths) / CHUNK_SIZE) if paths else 0
    chunk_existing = 0
    for i in range(chunk_count):
        if (CHUNK_DB_DIR / f"chunk_{i:04d}.syldb").exists():
            chunk_existing += 1
    chunked_output_complete = bool(chunk_count) and chunk_existing == chunk_count
    chunk_list = CHUNK_ROOT / "gtdb-r232-c200-dbv1.chunk_syldb.list"
    chunk_list_exists = chunk_list.exists()
    chunked_ready_to_build = (
        SYLPH.exists()
        and REP_LIST.exists()
        and inputs_complete
        and disk_heuristic_pass
        and chunk_existing < chunk_count
    )

    command = (
        f"/usr/bin/time -v -o {OUT_BASE}.build.time.log "
        f"{SYLPH} sketch -t {THREADS} -c {C_VALUE} -k {K_VALUE} "
        f"--gl {REP_LIST} -o {OUT_BASE}"
    )
    chunked_build_script = f"""#!/usr/bin/env bash
set -euo pipefail

SYLPH={SYLPH}
REP_LIST={REP_LIST}
CHUNK_ROOT={CHUNK_ROOT}
CHUNK_LIST_DIR={CHUNK_LIST_DIR}
CHUNK_DB_DIR={CHUNK_DB_DIR}
CHUNK_SIZE={CHUNK_SIZE}
THREADS={THREADS}
C_VALUE={C_VALUE}
K_VALUE={K_VALUE}

mkdir -p "$CHUNK_LIST_DIR" "$CHUNK_DB_DIR"
if ! compgen -G "$CHUNK_LIST_DIR/chunk_*.list" > /dev/null; then
  split -l "$CHUNK_SIZE" -d -a 4 --additional-suffix=.list "$REP_LIST" "$CHUNK_LIST_DIR/chunk_"
fi

for list in "$CHUNK_LIST_DIR"/chunk_*.list; do
  name=$(basename "$list" .list)
  out="$CHUNK_DB_DIR/$name"
  if [[ -s "$out.syldb" ]]; then
    continue
  fi
  /usr/bin/time -v -o "$CHUNK_DB_DIR/$name.build.time.log" \\
    "$SYLPH" sketch -t "$THREADS" -c "$C_VALUE" -k "$K_VALUE" --gl "$list" -o "$out"
done

find "$CHUNK_DB_DIR" -maxdepth 1 -name 'chunk_*.syldb' | sort > "$CHUNK_ROOT/gtdb-r232-c200-dbv1.chunk_syldb.list"
wc -l "$CHUNK_ROOT/gtdb-r232-c200-dbv1.chunk_syldb.list"
"""
    profile_commands = []
    chunk_profile_commands = []
    for sample_id, sketch in HMP_SAMPLE_SKETCHES.items():
        out_dir = HMP_R232_OUT / f"sylph_sample{sample_id}"
        chunk_input = CHUNK_PROFILE_INPUT_DIR / f"sample{sample_id}.profile_inputs.list"
        chunk_out = out_dir / "profile.chunked.tsv"
        profile_commands.append(
            {
                "sample": sample_id,
                "sample_sketch": str(sketch),
                "sample_sketch_exists": str(sketch.exists()).lower(),
                "profile_path": str(out_dir / "profile.tsv"),
                "profile_command": (
                    f"mkdir -p {out_dir}\n"
                    f"/usr/bin/time -v -o {out_dir}/profile.time.log "
                    f"{SYLPH} profile -t {THREADS} {OUT_DB} {sketch} "
                    f"-o {out_dir}/profile.tsv"
                ),
            }
        )
        chunk_profile_commands.append(
            {
                "sample": sample_id,
                "sample_sketch": str(sketch),
                "sample_sketch_exists": str(sketch.exists()).lower(),
                "profile_input_list": str(chunk_input),
                "profile_path": str(chunk_out),
                "profile_command": (
                    f"mkdir -p {out_dir} {CHUNK_PROFILE_INPUT_DIR}\n"
                    f"cat {CHUNK_ROOT}/gtdb-r232-c200-dbv1.chunk_syldb.list > {chunk_input}\n"
                    f"echo {sketch} >> {chunk_input}\n"
                    f"/usr/bin/time -v -o {out_dir}/profile.chunked.time.log "
                    f"{SYLPH} profile -t {THREADS} -l {chunk_input} -o {chunk_out}"
                ),
            }
        )
    profile_commands_ready = all(
        Path(row["sample_sketch"]).exists() for row in profile_commands
    )

    summary = [
        {"metric": "sylph_binary_exists", "value": str(SYLPH.exists()).lower()},
        {"metric": "sylph_version", "value": sylph_version()},
        {"metric": "representative_list", "value": str(REP_LIST)},
        {"metric": "representative_list_exists", "value": str(REP_LIST.exists()).lower()},
        {"metric": "representative_paths_total", "value": len(paths)},
        {"metric": "representative_paths_existing", "value": len(existing)},
        {"metric": "representative_paths_missing", "value": len(missing)},
        {"metric": "representative_compressed_bytes", "value": compressed_bytes},
        {"metric": "representative_compressed_gib", "value": fmt_gib(compressed_bytes)},
        {"metric": "sylph_r226_db_exists", "value": str(R226_DB.exists()).lower()},
        {"metric": "sylph_r226_db_gib", "value": fmt_gib(r226_bytes)},
        {"metric": "sylph_r232_db_path", "value": str(OUT_DB)},
        {"metric": "sylph_r232_db_exists", "value": str(OUT_DB.exists()).lower()},
        {"metric": "output_parent_free_gib", "value": fmt_gib(free_bytes)},
        {"metric": "disk_heuristic_required_gib", "value": fmt_gib(disk_heuristic_bytes)},
        {"metric": "disk_heuristic_pass", "value": str(disk_heuristic_pass).lower()},
        {"metric": "ready_to_build", "value": str(ready_to_build).lower()},
        {"metric": "chunk_size", "value": CHUNK_SIZE},
        {"metric": "chunk_expected_count", "value": chunk_count},
        {"metric": "chunk_existing_syldb_count", "value": chunk_existing},
        {"metric": "chunk_root", "value": str(CHUNK_ROOT)},
        {"metric": "chunk_syldb_list", "value": str(chunk_list)},
        {"metric": "chunk_syldb_list_exists", "value": str(chunk_list_exists).lower()},
        {"metric": "chunked_output_complete", "value": str(chunked_output_complete).lower()},
        {"metric": "chunked_ready_to_build", "value": str(chunked_ready_to_build).lower()},
        {"metric": "hmp_sample_sketches_ready", "value": str(profile_commands_ready).lower()},
        {"metric": "hmp_r232_profile_output_dir", "value": str(HMP_R232_OUT)},
        {"metric": "build_command", "value": command},
        {
            "metric": "preflight_status",
            "value": "chunked_output_complete"
            if chunked_output_complete and chunk_list_exists
            else "ready_to_build_missing_output"
            if ready_to_build
            else "output_already_exists"
            if OUT_DB.exists()
            else "not_ready",
        },
    ]
    missing_rows = [{"path": str(path)} for path in missing[:100]]

    RESULTS.mkdir(parents=True, exist_ok=True)
    write_tsv(RESULTS / "gtdb232_sylph_db_preflight.tsv", summary, ["metric", "value"])
    write_tsv(RESULTS / "gtdb232_sylph_db_missing_paths.tsv", missing_rows, ["path"])
    write_tsv(
        RESULTS / "gtdb232_sylph_hmp_profile_commands.tsv",
        profile_commands,
        ["sample", "sample_sketch", "sample_sketch_exists", "profile_path", "profile_command"],
    )
    write_tsv(
        RESULTS / "gtdb232_sylph_chunked_hmp_profile_commands.tsv",
        chunk_profile_commands,
        ["sample", "sample_sketch", "sample_sketch_exists", "profile_input_list", "profile_path", "profile_command"],
    )
    write_script(RESULTS / "gtdb232_sylph_db_build_command.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"{command}\n",
    )
    write_script(RESULTS / "gtdb232_sylph_chunked_build_command.sh", chunked_build_script)
    profile_script_lines = ["#!/usr/bin/env bash", "set -euo pipefail", "", f"test -s {OUT_DB}", ""]
    for row in profile_commands:
        profile_script_lines.append(str(row["profile_command"]))
        profile_script_lines.append("")
    write_script(RESULTS / "gtdb232_sylph_hmp_profile_commands.sh", "\n".join(profile_script_lines))
    chunk_profile_script_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"test -s {CHUNK_ROOT}/gtdb-r232-c200-dbv1.chunk_syldb.list",
        "",
    ]
    for row in chunk_profile_commands:
        chunk_profile_script_lines.append(str(row["profile_command"]))
        chunk_profile_script_lines.append("")
    write_script(
        RESULTS / "gtdb232_sylph_chunked_hmp_profile_commands.sh",
        "\n".join(chunk_profile_script_lines),
    )
    scorer = NOTE_DIR / "score_hmp_gtdb_source_abundance.py"
    write_script(RESULTS / "gtdb232_sylph_hmp_rescore_command.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"cd {REPO_ROOT}\n"
        f"test -s {HMP_R232_OUT}/sylph_sample6/profile.tsv\n"
        f"test -s {HMP_R232_OUT}/sylph_sample11/profile.tsv\n"
        f"python3 -B {scorer} --sylph-source r232 --r232-run {HMP_R232_OUT}\n",
    )
    write_script(RESULTS / "gtdb232_sylph_chunked_hmp_rescore_command.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n\n"
        f"cd {REPO_ROOT}\n"
        f"test -s {HMP_R232_OUT}/sylph_sample6/profile.chunked.tsv\n"
        f"test -s {HMP_R232_OUT}/sylph_sample11/profile.chunked.tsv\n"
        f"cp {HMP_R232_OUT}/sylph_sample6/profile.chunked.tsv {HMP_R232_OUT}/sylph_sample6/profile.tsv\n"
        f"cp {HMP_R232_OUT}/sylph_sample11/profile.chunked.tsv {HMP_R232_OUT}/sylph_sample11/profile.tsv\n"
        f"python3 -B {scorer} --sylph-source r232 --r232-run {HMP_R232_OUT}\n",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
