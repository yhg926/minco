#!/usr/bin/env bash
set -euo pipefail

# Experiment: addback_abundance_normalization_search
# Date: 2026-06-24
# Project: KSSD3mini
#
# Run from the recorded working directory so relative paths resolve the same way.
cd /home/ubuntu/yihuiguang/tools/KSSD3mini

# Code repository: /home/ubuntu/yihuiguang/tools/KSSD3mini
# Git metadata directory: /home/ubuntu/yihuiguang/tools/KSSD3mini/.minco_git
# Commit: b55b4d3f4c986de29098d2f1a092ee28251008aa
# Ref/describe: main, b55b4d3-dirty
# Working tree: dirty
# Binary/script/tool: python3 research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py

python3 -m py_compile research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py

/usr/bin/time -v python3 research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py --max-rules 2000 > /tmp/minco_abundance_norm_2000.stdout 2> /tmp/minco_abundance_norm_2000.stderr

/usr/bin/time -v python3 research/experiments/2026-06-24_addback_abundance_normalization_search/search_abundance_normalization.py --max-rules 12000 > /tmp/minco_abundance_norm_blend12000.stdout 2> /tmp/minco_abundance_norm_blend12000.stderr

python3 - <<'PY'
from pathlib import Path

import pandas as pd

exp = Path("research/experiments/2026-06-24_addback_abundance_normalization_search")
summary = pd.read_csv(exp / "results/abundance_rule_summary.tsv", sep="\t")
sylph = pd.read_csv(exp / "results/sylph_baseline.tsv", sep="\t")

def minco_row(label, row):
    return {
        "method": label,
        "F1": row["all_F1"],
        "L1": row["all_L1"],
        "mouse_F1": row["mouse_F1"],
        "mouse_L1": row["mouse_L1"],
        "cami3_F1": row["cami3_F1"],
        "cami3_L1": row["cami3_L1"],
        "note": row["rule_id"],
    }

rows = []
s_all = sylph.loc[sylph["dataset"] == "all"].iloc[0]
s_mouse = sylph.loc[sylph["dataset"] == "mouse_gtdb"].iloc[0]
s_cami3 = sylph.loc[sylph["dataset"] == "cami3_ncbi"].iloc[0]
rows.append({
    "method": "Sylph",
    "F1": s_all["F1"],
    "L1": s_all["l1_pct_points"],
    "mouse_F1": s_mouse["F1"],
    "mouse_L1": s_mouse["l1_pct_points"],
    "cami3_F1": s_cami3["F1"],
    "cami3_L1": s_cami3["l1_pct_points"],
    "note": "baseline",
})
rows.append(minco_row("best_L1", summary.sort_values(["all_L1", "all_F1"], ascending=[True, False]).iloc[0]))
rows.append(minco_row("best_F1", summary.sort_values(["all_F1", "all_L1"], ascending=[False, True]).iloc[0]))
rows.append(minco_row("best_L1_F1ge078", summary.loc[summary["all_F1"] >= 0.78].sort_values(["all_L1", "all_F1"], ascending=[True, False]).iloc[0]))
pd.DataFrame(rows).to_csv(exp / "results/final_comparison.tsv", sep="\t", index=False)
PY
