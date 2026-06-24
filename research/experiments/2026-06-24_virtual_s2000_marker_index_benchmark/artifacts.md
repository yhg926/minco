# Artifacts

| Kind | Path or Link | Source Path or URI | Description | Availability | Preserve? | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| script | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/benchmark_virtual_markers.py` | local workspace | Virtual marker-size benchmark from full sorted inverted index | available | yes |  |
| input-refdb | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup` | local temporary storage | Full S2000 dedup refdb | available | maybe | Temporary `/tmp` path |
| input-index | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup/minco.refindex.ctxgid64obj32` | local temporary storage | Full sorted inverted index scanned by benchmark | available | maybe | 4,812,648,000 bytes |
| input-markerdb | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker` | local temporary storage | Physical context markerdb equivalence target | available | maybe | Temporary `/tmp` path |
| input-psmp | `/tmp/gtdb232_s2000_dedup_marker.qKJofv/markerdb_ctx_s2000_dedup.psmp.tsv` | local temporary storage | Physical marker-size table equivalence target | available | maybe |  |
| command-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/commands.sh` | local workspace | Exact commands and working directory | available | yes | Rerunnable while `/tmp` data exist |
| parameter-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/parameters.tsv` | local workspace | Parameters and result-affecting constants | available | yes |  |
| code-provenance | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/provenance/code_status.txt` | local workspace | Commit, dirty status, and diffstat | available | yes |  |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results/virtual_marker_summary.tsv` | generated | Flat summary metrics | available | yes | Main result table |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results/virtual_marker_summary.json` | generated | Structured summary metrics | available | yes |  |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results/virtual_marker_sizes.tsv` | generated | Per-reference full/virtual/physical marker sizes | available | yes | 200,528 lines including header |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results/virtual_marker_mismatches.tsv` | generated | Mismatch table | available | yes | Header only after fix |
| output | `/home/ubuntu/yihuiguang/tools/KSSD3mini/research/experiments/2026-06-24_virtual_s2000_marker_index_benchmark/results/low_marker_refs_head200.tsv` | generated | First 200 low-marker refs sorted by marker size | available | yes | Diagnostic |
| temporary-log | `/tmp/virtual_s2000_marker_scan_fixed.stdout` | generated | Final scan stdout | temporary | no | Regenerable |
| temporary-log | `/tmp/virtual_s2000_marker_scan_fixed.stderr` | generated | Final scan timing log | temporary | no | Regenerable |
