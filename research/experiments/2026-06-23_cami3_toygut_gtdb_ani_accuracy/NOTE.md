# CAMI3 ToyGut ANI Accuracy With GTDB-Mapped Truth

Date: 2026-06-23

## Question

Using GTDB-mapped ground truth, how accurate are the ANI calls for the CAMI3 ToyGut comparison methods?

Methods scored:

- Sylph
- old ctx-marker robust rescue
- coden15 formula-AF gate robust rescue
- coden15 original gate robust rescue

## Truth Transfer

CAMI3 provides NCBI species taxid truth, not a direct GTDB species profile. I transferred each bacterial NCBI species taxid to GTDB r232 species through GTDB metadata.

Two truth views are used:

- `unique`: NCBI species taxid maps to exactly one GTDB species. This is strict and safe but excludes most CAMI3 truth.
- `possible`: NCBI species taxid maps to one or more GTDB species. A prediction is compatible if it lands in this possible GTDB set.

Truth transfer summary:

| sample | unique NCBI species | unique abundance | ambiguous NCBI species | ambiguous abundance | unmapped abundance |
|---:|---:|---:|---:|---:|---:|
| 0 | 34 | 0.5024 | 74 | 0.4776 | 0.0200 |
| 1 | 39 | 0.3181 | 76 | 0.6619 | 0.0200 |
| 2 | 32 | 0.2223 | 87 | 0.7577 | 0.0200 |

This means unique-only GTDB truth is too strict for final biological interpretation. The possible-GTDB view is better for asking whether a high-ANI call is at least compatible with the CAMI3 truth.

## Selected-Call ANI Result

Mean over samples 0-2:

| method | selected GTDB species | unique TP | unique recall | possible precision | outside-possible high-ANI calls | TP ANI min | TP ANI <0.95 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sylph | 72.00 | 22.33 | 0.6768 | 1.0000 | 0.00 | 0.9560 | 0.00 |
| coden15 formula-AF | 72.67 | 21.67 | 0.6553 | 0.9589 | 3.00 | 1.0000 | 0.00 |
| old ctx-marker robust rescue | 63.33 | 18.67 | 0.5655 | 0.9599 | 2.67 | 1.0000 | 0.00 |
| coden15 original | 56.00 | 17.67 | 0.5350 | 0.9879 | 0.67 | 1.0000 | 0.00 |

Under strict unique-GTDB truth, Sylph has the best selected-call F1 (`0.4245`), followed by coden15 formula-AF (`0.4082`), coden15 original (`0.3946`), and old ctx-marker (`0.3870`). This strict F1 is depressed because ambiguous GTDB truth is excluded.

Under possible-GTDB truth, Sylph is cleanest: every selected high-ANI GTDB species is compatible with CAMI3 truth. Coden15 formula-AF has higher recall than old/original but introduces about three high-ANI calls per sample outside possible GTDB truth.

Outside-possible selected high-ANI calls include:

- sample 0: `s__Anaerostipes sp019423795`, `s__Faecalibacterium sp902463235`
- sample 1: `s__Alistipes sp047233565`, `s__Anaerostipes sp019423795`, `s__Faecalibacterium sp902463235`, plus low-support `s__Lawsonibacter sp959022875` for coden15 formula-AF
- sample 2: `s__Blautia_A sp900751995`, `s__Enterococcus_D sp022711995`, `s__Veillonella sp902466275`

## Candidate ANI>=0.95 Result

Raw MinCO candidate ANI thresholding alone is not enough:

| method | candidate GTDB species | ANI>=0.95 species | strict precision | strict recall | possible precision |
|---|---:|---:|---:|---:|---:|
| coden15 candidates | 4038.67 | 1652.00 | 0.0193 | 0.9695 | 0.1574 |
| old ctx-marker candidates | 37908.00 | 1999.00 | 0.0162 | 0.9820 | 0.1283 |
| Sylph selected rows | 72.00 | 72.00 | 0.3093 | 0.6768 | 1.0000 |

This shows why the full gate is necessary: raw MinCO ANI `>=0.95` sees nearly all unique true species but also produces many high-ANI compatible-looking candidates outside truth.

## Conclusion

Using GTDB-compatible truth, Sylph has the cleanest selected ANI behavior on CAMI3 ToyGut. Coden15 formula-AF is the best MinCO variant for recall among the selected-call methods, but its relaxed gate creates a small number of high-ANI calls outside possible GTDB truth. Coden15 original is more conservative but misses more unique-GTDB true species. Old ctx-marker sits between them on outside-possible error, but has lower recall than formula-AF.

Main caveat: this is not true source-assembly ANI ground truth. CAMI3 source assembly IDs are not directly available in the current local files; this GTDB truth is transferred from CAMI3 NCBI species taxids through GTDB metadata, so ambiguous NCBI-to-GTDB species splits dominate the uncertainty.

## Artifacts

- Scorer: `score_gtdb_ani_accuracy.py`
- Truth transfer audit: `gtdb_truth_transfer_audit.tsv`
- Truth transfer summary: `gtdb_truth_transfer_summary.tsv`
- Selected-call mean metrics: `selected_call_gtdb_ani_mean_metrics.tsv`
- Selected-call details: `selected_call_gtdb_ani_details.tsv`
- Candidate ANI>=0.95 mean metrics: `candidate_gtdb_ani95_mean_metrics.tsv`
- Candidate details: `candidate_gtdb_ani95_details.tsv`
