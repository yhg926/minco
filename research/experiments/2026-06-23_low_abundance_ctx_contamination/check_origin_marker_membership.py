#!/usr/bin/env python3
"""Check whether contaminant trace contexts exist in origin representative sketches."""

from __future__ import annotations

from pathlib import Path
import math
import struct

import numpy as np
import pandas as pd

import analyze_ctx_trace_sources as trace_sources


ROOT = Path("/home/ubuntu/yihuiguang/tools/KSSD3mini")
EXP = ROOT / "research/experiments/2026-06-23_low_abundance_ctx_contamination"
MARKER_REF = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker")
FULL_REF = Path("/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup")
S10000_REF = Path("/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno")
SOURCE_REP = ROOT / "research/experiments/2026-06-23_toymouse_source_rep_ani/toymouse_sample0_source_rep_ani.tsv"

OBJ_BITS = 20


def read_names(sketch_dir: Path) -> list[str]:
    stat_path = sketch_dir / "minco.stat"
    with stat_path.open("rb") as f:
        hdr = f.read(32)
        vals = struct.unpack("<IBBxxiiiiii", hdr)
        n = vals[-1]
        return [
            f.read(256).split(b"\0", 1)[0].decode("utf-8", "replace")
            for _ in range(n)
        ]


class SketchCtxLookup:
    def __init__(self, sketch_dir: Path):
        self.sketch_dir = sketch_dir
        self.names = read_names(sketch_dir)
        self.offsets = np.memmap(sketch_dir / "minco.ctxobj64.offsets", dtype="<u8", mode="r")
        self.ctxobj = np.memmap(sketch_dir / "minco.ctxobj64", dtype="<u8", mode="r")
        self.acc_to_gid = {}
        for i, name in enumerate(self.names):
            base = Path(name).name
            acc = base.split("_", 2)[0] + "_" + base.split("_", 2)[1] if base.startswith(("GCF_", "GCA_")) else ""
            if acc:
                self.acc_to_gid[acc] = i
        self._ctx_cache: dict[str, set[int]] = {}
        self._size_cache: dict[str, int] = {}

    def gid_for_accession(self, acc: str) -> int | None:
        return self.acc_to_gid.get(str(acc))

    def ctx_set(self, acc: str) -> set[int]:
        acc = str(acc)
        if acc in self._ctx_cache:
            return self._ctx_cache[acc]
        gid = self.gid_for_accession(acc)
        if gid is None:
            self._ctx_cache[acc] = set()
            self._size_cache[acc] = 0
            return self._ctx_cache[acc]
        begin = int(self.offsets[gid])
        end = int(self.offsets[gid + 1])
        vals = np.asarray(self.ctxobj[begin:end], dtype=np.uint64)
        ctxs = set((vals >> OBJ_BITS).astype(np.uint64).astype(object))
        self._ctx_cache[acc] = {int(x) for x in ctxs}
        self._size_cache[acc] = end - begin
        return self._ctx_cache[acc]

    def size(self, acc: str) -> int:
        if acc not in self._size_cache:
            self.ctx_set(acc)
        return self._size_cache.get(acc, 0)

    def contains(self, acc: str, qctx: int) -> bool:
        return int(qctx) in self.ctx_set(str(acc))


def frac(num: float, den: float) -> float:
    return float(num) / float(den) if den else math.nan


def load_classified_events() -> pd.DataFrame:
    truth = pd.read_csv(trace_sources.SOURCE_TRUTH, sep="\t", dtype=str)
    source_rep = pd.read_csv(SOURCE_REP, sep="\t", dtype=str)
    needed = trace_sources.collect_needed_read_ids()
    mapping = trace_sources.load_read_mapping(needed)
    rows = []
    for target in trace_sources.TARGETS:
        label = target["target_label"]
        species = target["gtdb_species"]
        target_rep = target["representative_accession"]
        species_truth = truth[truth["gtdb_species"].eq(species)].copy()
        target_source_ids = set(species_truth["genome_id"].astype(str))
        trace = pd.read_csv(trace_sources.trace_path_for(target), sep="\t", dtype={"read_id": str})
        trace["qctx"] = pd.to_numeric(trace["qctx"], errors="coerce").astype("Int64")
        classified = trace_sources.classify_trace(
            trace,
            mapping,
            truth,
            species,
            target_source_ids,
        )
        classified["target_label"] = label
        classified["target_species"] = species
        classified["target_rep"] = target_rep
        rows.append(classified)
    events = pd.concat(rows, ignore_index=True)
    rep_map = (
        source_rep.sort_values("source_abundance", ascending=False)
        .drop_duplicates("source_accession")
        .set_index("source_accession")["gtdb_representative_accession"]
        .to_dict()
    )
    events["origin_rep"] = events["source_accession"].map(rep_map)
    return events


def summarize(
    events: pd.DataFrame,
    marker: SketchCtxLookup,
    full: SketchCtxLookup,
    s10000: SketchCtxLookup,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    detail_rows = []
    events = events[events["origin_rep"].notna()].copy()
    events["qctx_int"] = events["qctx"].astype("int64")
    events["origin_marker_present"] = [
        marker.contains(acc, qctx) for acc, qctx in zip(events["origin_rep"], events["qctx_int"])
    ]
    events["origin_full_present"] = [
        full.contains(acc, qctx) for acc, qctx in zip(events["origin_rep"], events["qctx_int"])
    ]
    events["origin_s10000_present"] = [
        s10000.contains(acc, qctx) for acc, qctx in zip(events["origin_rep"], events["qctx_int"])
    ]
    events["target_marker_present"] = [
        marker.contains(acc, qctx) for acc, qctx in zip(events["target_rep"], events["qctx_int"])
    ]
    events["candidate_refs"] = pd.to_numeric(events["candidate_refs"], errors="coerce")

    other = events[events["category"].eq("other_or_unmapped")].copy()
    for label, g in other.groupby("target_label", sort=False):
        rows.append(
            {
                "target_label": label,
                "other_event_n": len(g),
                "other_unique_qctx_n": g["qctx_int"].nunique(),
                "candidate_refs_gt1_event_frac": frac((g["candidate_refs"] > 1).sum(), len(g)),
                "origin_marker_present_event_frac": frac(g["origin_marker_present"].sum(), len(g)),
                "origin_full_present_event_frac": frac(g["origin_full_present"].sum(), len(g)),
                "origin_s10000_present_event_frac": frac(g["origin_s10000_present"].sum(), len(g)),
                "full_present_but_marker_absent_event_frac": frac(
                    ((g["origin_full_present"]) & (~g["origin_marker_present"])).sum(), len(g)
                ),
                "s10000_present_but_s2000_absent_event_frac": frac(
                    ((g["origin_s10000_present"]) & (~g["origin_full_present"])).sum(), len(g)
                ),
                "origin_neither_event_frac": frac(
                    ((~g["origin_full_present"]) & (~g["origin_marker_present"])).sum(), len(g)
                ),
                "target_marker_present_event_frac": frac(g["target_marker_present"].sum(), len(g)),
                "origin_marker_present_unique_qctx_frac": frac(
                    g.drop_duplicates("qctx_int")["origin_marker_present"].sum(),
                    g["qctx_int"].nunique(),
                ),
                "origin_full_present_unique_qctx_frac": frac(
                    g.drop_duplicates("qctx_int")["origin_full_present"].sum(),
                    g["qctx_int"].nunique(),
                ),
                "origin_s10000_present_unique_qctx_frac": frac(
                    g.drop_duplicates("qctx_int")["origin_s10000_present"].sum(),
                    g["qctx_int"].nunique(),
                ),
                "full_present_but_marker_absent_unique_qctx_frac": frac(
                    (
                        g.drop_duplicates("qctx_int")["origin_full_present"]
                        & ~g.drop_duplicates("qctx_int")["origin_marker_present"]
                    ).sum(),
                    g["qctx_int"].nunique(),
                ),
                "s10000_present_but_s2000_absent_unique_qctx_frac": frac(
                    (
                        g.drop_duplicates("qctx_int")["origin_s10000_present"]
                        & ~g.drop_duplicates("qctx_int")["origin_full_present"]
                    ).sum(),
                    g["qctx_int"].nunique(),
                ),
                "origin_neither_unique_qctx_frac": frac(
                    (
                        ~g.drop_duplicates("qctx_int")["origin_full_present"]
                        & ~g.drop_duplicates("qctx_int")["origin_marker_present"]
                    ).sum(),
                    g["qctx_int"].nunique(),
                ),
            }
        )
    group_cols = [
        "target_label",
        "source_accession",
        "origin_rep",
        "gtdb_species",
        "category",
    ]
    for key, g in other.groupby(group_cols, dropna=False, sort=False):
        target_label, source_accession, origin_rep, gtdb_species, category = key
        detail_rows.append(
            {
                "target_label": target_label,
                "source_accession": source_accession,
                "origin_rep": origin_rep,
                "origin_marker_size": marker.size(origin_rep),
                "origin_full_size": full.size(origin_rep),
                "origin_s10000_size": s10000.size(origin_rep),
                "gtdb_species": gtdb_species,
                "event_n": len(g),
                "unique_qctx_n": g["qctx_int"].nunique(),
                "candidate_refs_gt1_event_frac": frac((g["candidate_refs"] > 1).sum(), len(g)),
                "origin_marker_present_event_frac": frac(g["origin_marker_present"].sum(), len(g)),
                "origin_full_present_event_frac": frac(g["origin_full_present"].sum(), len(g)),
                "origin_s10000_present_event_frac": frac(g["origin_s10000_present"].sum(), len(g)),
                "full_present_but_marker_absent_event_frac": frac(
                    ((g["origin_full_present"]) & (~g["origin_marker_present"])).sum(), len(g)
                ),
                "s10000_present_but_s2000_absent_event_frac": frac(
                    ((g["origin_s10000_present"]) & (~g["origin_full_present"])).sum(), len(g)
                ),
                "origin_neither_event_frac": frac(
                    ((~g["origin_full_present"]) & (~g["origin_marker_present"])).sum(), len(g)
                ),
                "target_marker_present_event_frac": frac(g["target_marker_present"].sum(), len(g)),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(detail_rows).sort_values(
        ["target_label", "event_n"], ascending=[True, False]
    )


def main() -> None:
    events = load_classified_events()
    marker = SketchCtxLookup(MARKER_REF)
    full = SketchCtxLookup(FULL_REF)
    s10000 = SketchCtxLookup(S10000_REF)
    summary, details = summarize(events, marker, full, s10000)
    summary.to_csv(EXP / "origin_membership_summary.tsv", sep="\t", index=False)
    details.to_csv(EXP / "origin_membership_by_source.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
