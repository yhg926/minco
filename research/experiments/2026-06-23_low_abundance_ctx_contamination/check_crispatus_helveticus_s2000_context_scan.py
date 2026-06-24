#!/usr/bin/env python3
"""Scan whole genomes for S2000-compatible coden11 context IDs.

This avoids making a larger MinCO sketch, because MinCO's sketch ID includes
the target sketch size and S2000 context hashes are not comparable to S10000000.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import pandas as pd

import analyze_ctx_trace_sources as trace_sources
from check_origin_marker_membership import SketchCtxLookup


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-23_low_abundance_ctx_contamination"
SOURCE_REP = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_source_rep_ani.tsv"
S2000_FULL = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup")
S2000_MARKER = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker")

HEL_SOURCE_ACC = "GCF_001702095.1"
HEL_REP_ACC = "GCF_000160855.1"
CRISPATUS_REP_ACC = "GCF_018987235.1"

KLEN = 32
OBJ_BITS = 20
SEED = 0x9E3779B97F4A7C15
MASK64 = (1 << 64) - 1


def coden11_ctxmask() -> int:
    pattern = 0
    for i in range(11):
        pattern = (pattern << 6) | 0b111100
        if i == 10:
            pattern = (pattern << 4) | 0b1111
            break
    return pattern


CTXMASK = coden11_ctxmask()


BASEMAP = {
    ord("A"): 0,
    ord("a"): 0,
    ord("C"): 1,
    ord("c"): 1,
    ord("G"): 2,
    ord("g"): 2,
    ord("T"): 3,
    ord("t"): 3,
}


def mix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return (x ^ (x >> 31)) & MASK64


def hash_ctx(ctx: int) -> int:
    return mix64(ctx ^ SEED) >> OBJ_BITS


def coden11_ctx_obj(unituple: int) -> tuple[int, int]:
    high = 0
    low = 0
    u = unituple
    for i in range(10):
        high |= (u & 0xF) << (4 * i)
        u >>= 4
        low |= (u & 0x3) << (2 * i)
        u >>= 2
    high |= (u & 0xF) << 40
    return high, low


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def iter_fasta_sequences(path: Path):
    chunks: list[str] = []
    with open_text(path) as fh:
        for raw in fh:
            if raw.startswith(">"):
                if chunks:
                    yield "".join(chunks)
                    chunks = []
            else:
                chunks.append(raw.strip())
        if chunks:
            yield "".join(chunks)


def scan_watch_contexts(path: Path, watch: set[int]) -> dict[int, set[int]]:
    found: dict[int, set[int]] = {ctx: set() for ctx in watch}
    len_mv = 2 * KLEN - 2
    objmask = (1 << OBJ_BITS) - 1
    for seq in iter_fasta_sequences(path):
        tuple_fwd = 0
        tuple_rev = 0
        base_count = 0
        for ch in seq.encode("ascii", "ignore"):
            b = BASEMAP.get(ch)
            if b is None:
                tuple_fwd = 0
                tuple_rev = 0
                base_count = 0
                continue
            tuple_fwd = ((tuple_fwd << 2) | b) & MASK64
            tuple_rev = ((tuple_rev >> 2) | ((b ^ 3) << len_mv)) & MASK64
            base_count += 1
            if base_count < KLEN:
                continue
            use_fwd = (tuple_fwd & CTXMASK) < (tuple_rev & CTXMASK)
            unituple = tuple_fwd if use_fwd else tuple_rev
            ctx, obj = coden11_ctx_obj(unituple)
            hctx = hash_ctx(ctx)
            if hctx in watch:
                found[hctx].add(obj)
    return found


def load_paths() -> dict[str, Path]:
    df = pd.read_csv(SOURCE_REP, sep="\t", dtype=str)
    hel = df[df["source_accession"].eq(HEL_SOURCE_ACC)].iloc[0]
    crisp = df[df["gtdb_representative_accession"].eq(CRISPATUS_REP_ACC)].iloc[0]
    return {
        "helveticus_source": Path(hel["source_path"]),
        "helveticus_rep": Path(hel["representative_path"]),
        "crispatus_rep": Path(crisp["representative_path"]),
    }


def crispatus_helveticus_events() -> pd.DataFrame:
    cached = EXP / "crispatus_helveticus_wholegenome_membership_events.tsv"
    if cached.exists() and cached.stat().st_size > 0:
        df = pd.read_csv(cached, sep="\t", dtype={"read_id": str})
        df = df[["read_id", "qctx", "diff"]].copy()
        df["qctx"] = pd.to_numeric(df["qctx"], errors="coerce").astype("int64")
        df["diff"] = pd.to_numeric(df["diff"], errors="coerce")
        return df

    truth = pd.read_csv(trace_sources.SOURCE_TRUTH, sep="\t", dtype=str)
    target = next(t for t in trace_sources.TARGETS if t["target_label"] == "crispatus_low_bad")
    needed = trace_sources.collect_needed_read_ids()
    mapping = trace_sources.load_read_mapping(needed)
    trace = pd.read_csv(trace_sources.trace_path_for(target), sep="\t", dtype={"read_id": str})
    classified = trace_sources.classify_trace(
        trace,
        mapping,
        truth,
        target["gtdb_species"],
        set(truth[truth["gtdb_species"].eq(target["gtdb_species"])]["genome_id"].astype(str)),
    )
    out = classified[classified["source_accession"].eq(HEL_SOURCE_ACC)].copy()
    out["qctx"] = pd.to_numeric(out["qctx"], errors="coerce").astype("int64")
    out["diff"] = pd.to_numeric(out["diff"], errors="coerce")
    return out


def context_status(objset: set[int]) -> tuple[bool, bool, int]:
    raw_present = bool(objset)
    postconflict_present = len(objset) == 1
    return raw_present, postconflict_present, len(objset)


def make_rows(events: pd.DataFrame, scans: dict[str, dict[int, set[int]]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for _, event in events.iterrows():
        row = {
            "read_id": event["read_id"],
            "qctx": int(event["qctx"]),
            "diff": event["diff"],
        }
        for label, found in scans.items():
            raw, post, nobj = context_status(found[int(event["qctx"])])
            row[f"{label}_raw_present"] = raw
            row[f"{label}_postconflict_present"] = post
            row[f"{label}_object_count"] = nobj
        rows.append(row)
    detail = pd.DataFrame(rows)

    unique_rows = []
    for qctx, g in detail.groupby("qctx", sort=True):
        row = {"qctx": qctx, "event_n": len(g), "mean_diff": g["diff"].mean()}
        for label in scans:
            row[f"{label}_raw_present"] = bool(g[f"{label}_raw_present"].iloc[0])
            row[f"{label}_postconflict_present"] = bool(g[f"{label}_postconflict_present"].iloc[0])
            row[f"{label}_object_count"] = int(g[f"{label}_object_count"].iloc[0])
        unique_rows.append(row)
    return detail, pd.DataFrame(unique_rows)


def summarize(detail: pd.DataFrame, unique: pd.DataFrame) -> pd.DataFrame:
    rows = []
    labels = ["helveticus_source", "helveticus_rep", "crispatus_rep"]
    for label in labels:
        for mode in ["raw", "postconflict"]:
            col = f"{label}_{mode}_present"
            rows.append(
                {
                    "genome": label,
                    "mode": mode,
                    "event_present_n": int(detail[col].sum()),
                    "event_total": len(detail),
                    "event_present_frac": float(detail[col].mean()) if len(detail) else math.nan,
                    "unique_qctx_present_n": int(unique[col].sum()),
                    "unique_qctx_total": len(unique),
                    "unique_qctx_present_frac": float(unique[col].mean()) if len(unique) else math.nan,
                }
            )
    return pd.DataFrame(rows)


def validation_rows(scans: dict[str, dict[int, set[int]]]) -> pd.DataFrame:
    s2000_full = SketchCtxLookup(S2000_FULL)
    s2000_marker = SketchCtxLookup(S2000_MARKER)
    rows = []
    for label, acc in {
        "helveticus_rep": HEL_REP_ACC,
        "crispatus_rep": CRISPATUS_REP_ACC,
    }.items():
        for db_label, db in [("s2000_full", s2000_full), ("s2000_marker", s2000_marker)]:
            ctxs = db.ctx_set(acc)
            seen_raw = 0
            seen_post = 0
            # Only contexts that were also watched are included in scans; this
            # validation is for the suspicious set, not the full database.
            for qctx in scans[label]:
                if qctx in ctxs:
                    raw, post, _ = context_status(scans[label][qctx])
                    seen_raw += int(raw)
                    seen_post += int(post)
            rows.append(
                {
                    "genome": label,
                    "database": db_label,
                    "database_size": db.size(acc),
                    "watched_db_context_n": sum(1 for qctx in scans[label] if qctx in ctxs),
                    "watched_db_context_raw_present_n": seen_raw,
                    "watched_db_context_postconflict_present_n": seen_post,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    events = crispatus_helveticus_events()
    qctxs = set(int(x) for x in events["qctx"].unique())
    paths = load_paths()
    scans = {label: scan_watch_contexts(path, qctxs) for label, path in paths.items()}
    detail, unique = make_rows(events, scans)
    summary = summarize(detail, unique)
    validation = validation_rows(scans)

    detail.to_csv(EXP / "crispatus_helveticus_s2000_context_scan_events.tsv", sep="\t", index=False)
    unique.to_csv(EXP / "crispatus_helveticus_s2000_context_scan_unique_qctx.tsv", sep="\t", index=False)
    summary.to_csv(EXP / "crispatus_helveticus_s2000_context_scan_summary.tsv", sep="\t", index=False)
    validation.to_csv(EXP / "crispatus_helveticus_s2000_context_scan_validation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
