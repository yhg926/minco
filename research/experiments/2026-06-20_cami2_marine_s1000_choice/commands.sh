#!/usr/bin/env bash
set -euo pipefail

MINCO=./minco_core/bin/minco
REF_S1000=/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
GOLD=/tmp/gs_marine_short.profile
TAXMAP=/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv
OUT=/tmp/minco_s1000_zipaaf_20260620

mkdir -p "$OUT"

for SAMPLE in 0 1 2; do
  /usr/bin/time -v "$MINCO" ani -p16 -r "$REF_S1000" \
    --qraw "/tmp/cami_marine_sample${SAMPLE}_reads.fq.gz" \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-unique --readwise-ani zip-aaf \
    -m0 -f0.05 -n0.94 -t10 \
    -o "$OUT/s1000_gtdb_unique_zipaaf_sample${SAMPLE}_f0.05_n0.94_t10.tsv"

  python3 research/experiments/2026-06-20_cami2_marine_extra_samples/scripts/score_profiles.py \
    --sample-id "marmgCAMI2_short_read_sample_${SAMPLE}" \
    --minco "$OUT/s1000_gtdb_unique_zipaaf_sample${SAMPLE}_f0.05_n0.94_t10.tsv" \
    --gold "$GOLD" --taxmap "$TAXMAP" \
    --out "$OUT/sample${SAMPLE}_s1000_score.tsv"
done

python3 - <<'PY'
import pandas as pd

records = []
s1000_times = {
    "marmgCAMI2_short_read_sample_0": (130.96, 3.680080),
    "marmgCAMI2_short_read_sample_1": (129.45, 3.711520),
    "marmgCAMI2_short_read_sample_2": (113.04, 3.703028),
}
for sample_idx in range(3):
    sid = f"marmgCAMI2_short_read_sample_{sample_idx}"
    row = pd.read_csv(f"/tmp/minco_s1000_zipaaf_20260620/sample{sample_idx}_s1000_score.tsv", sep="\t").iloc[0].to_dict()
    row["label"] = "minco_s1000_gtdb_unique_zipaaf"
    row["wall_seconds"], row["peak_rss_gb"] = s1000_times[sid]
    records.append(row)
prev = pd.read_csv("research/experiments/2026-06-20_cami2_marine_extra_samples/summary.tsv", sep="\t")
prev = prev[prev["sample_id"].astype(str).str.startswith("marmgCAMI2_short_read_sample_")]
records.extend(prev.to_dict("records"))
out = pd.DataFrame(records)
keep = ["sample_id", "label", "gold_taxa", "pred_taxa", "TP", "FP", "FN", "precision", "recall", "F1", "tp_abundance_pearson", "tp_abundance_mae", "wall_seconds", "peak_rss_gb"]
out = out[keep].sort_values(["sample_id", "label"])
means = out.groupby("label", as_index=False)[["gold_taxa", "pred_taxa", "TP", "FP", "FN", "precision", "recall", "F1", "tp_abundance_pearson", "tp_abundance_mae", "wall_seconds", "peak_rss_gb"]].mean()
means.insert(0, "sample_id", "mean_samples_0_2")
pd.concat([out, means], ignore_index=True).to_csv("research/experiments/2026-06-20_cami2_marine_s1000_choice/summary.tsv", sep="\t", index=False)
PY
