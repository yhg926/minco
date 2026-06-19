# Artifacts

| Type | Path | Notes |
| --- | --- | --- |
| script | `research/experiments/2026-06-19_minco_half_genome_tool_comparison/run_half_benchmark.py` | Reproducible driver for splitting and tool runs. |
| commands | `research/experiments/2026-06-19_minco_half_genome_tool_comparison/commands.sh` | Single rerun entry point. |
| summary | `research/experiments/2026-06-19_minco_half_genome_tool_comparison/summary.tsv` | Final comparison table. |
| timings | `research/experiments/2026-06-19_minco_half_genome_tool_comparison/tool_times.tsv` | Per-tool wall-clock timings from Python. |
| split metadata | `research/experiments/2026-06-19_minco_half_genome_tool_comparison/split_metadata.tsv` | Query-genome split sizes. |
| raw outputs | `/tmp/minco_half_genome_tool_comparison/raw_outputs/` | Tool stdout/stderr and per-tool TSV outputs, including `*.minco_naive.tsv`. |
| half FASTA | `/tmp/minco_half_genome_tool_comparison/tmp/3300014912_1.first_half.fna` | First concatenated half of query genome. |
| half FASTA | `/tmp/minco_half_genome_tool_comparison/tmp/3300014912_1.second_half.fna` | Second concatenated half of query genome. |

External inputs:

| Type | Path | Notes |
| --- | --- | --- |
| reference genome | `/mnt/new3T/skani_data/Nayfach_data/fna/3300027414_1.fna` | Complete reference genome for pair `P000881414`. |
| query genome | `/mnt/new3T/skani_data/Nayfach_data/fna/3300014912_1.fna` | Complete query genome split into two halves. |
| ANIm label table | `research/experiments/2026-06-18_minco_nayfach_anim_large_joint/large_10k_perbin20000/sample_pairs.tsv` | Existing Nayfach ANIm label source for the full pair. |
