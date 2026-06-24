#!/usr/bin/env bash
set -euo pipefail

EXP=research/experiments/2026-06-20_minco_readwise_correction_research
OUT=/tmp/minco_readwise_correction_20260620

python3 "$EXP/scripts/analyze_readwise_corrections.py" \
  --minco /tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_unfiltered.tsv \
  --sylph /tmp/sylph_marine_sample0/profile.tsv \
  --gold /tmp/gs_marine_short.profile \
  --taxmap /tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv \
  --sample-id marmgCAMI2_short_read_sample_0 \
  --ctx-k 11 \
  --outdir "$OUT"

cat "$OUT/manifest.tsv"
cat "$OUT/strategy_summary.tsv"
cat "$OUT/model_summary.tsv"
