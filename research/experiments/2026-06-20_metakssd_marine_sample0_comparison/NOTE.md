# MetaKSSD OPAL Marine Sample0 Comparison

Date: 2026-06-20

## Question

Compare current minco CAMI marine sample0 profiling performance with the MetaKSSD repository benchmark data for the same marine sample0 dataset.

## Sources

- MetaKSSD local repository: `/mnt/new3T/gtdbr220/smash_sylph_r226/MetaKSSD`
- MetaKSSD README benchmark section: `/mnt/new3T/gtdbr220/smash_sylph_r226/MetaKSSD/README.md`
- MetaKSSD marine OPAL report: `https://yhg926.github.io/KSSD2/OPAL/marine/`
- Fetched OPAL HTML: `/tmp/metakssd_opal_marine_index.html`
- MetaKSSD speed/memory screenshot: `/tmp/metakssd_speed_memory.png`
- Sylph DB: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`
- Sylph profile output: `/tmp/sylph_marine_sample0/profile.tsv`
- minco S1000 experiment: `research/experiments/2026-06-20_cami2_full_lineage_thresholds/NOTE.md`
- minco S10000 experiment: `research/experiments/2026-06-20_cami2_s10000_full_lineage_thresholds/NOTE.md`

## Important Caveats

- MetaKSSD/KSSD2 OPAL values are from the published OPAL report with GTDBr214 MarkerDB.
- minco values here use local GTDBr232-derived reference sketches. S1000 used GTDBr232 plus RefSeq virus; S10000 used GTDBr232 only.
- The minco TP/FP/FN values in `summary.tsv` use an OPAL-like species gold filter that drops three raw gold species entries: `unidentified plasmid`, `unidentified`, and `unidentified virus`. This matches the OPAL sample0 species gold count of 256.
- minco L1 and Bray-Curtis are local OPAL-like calculations from CAMI profiles. Weighted UniFrac was not recomputed locally.
- MetaKSSD speed/memory values are read from the README screenshot, not a machine-readable table.
- Sylph was run locally with `sylph 0.6.0`, default profile threshold, `-t 16`, and GTDB-r226 c200 DB.

## Results

Species-level sample0 detection:

- Best MetaKSSD/KSSD2 row in the OPAL table among KSSD2 variants was `KSSD2.L3K11S48`: F1 0.797, TP 208, FP 58, FN 48, L1 0.456, Bray-Curtis 0.246.
- Local Sylph with GTDB-r226 c200 DB gave F1 0.831, TP 222, FP 56, FN 34. This is better than the KSSD2/MetaKSSD rows by species F1 and close to the OPAL MetaPhlAn4 row.
- Best current minco row by F1 was `S10000_opt_combo_f0.05_n0.92_t10`: F1 0.731, TP 205, FP 100, FN 51, L1 0.752, Bray-Curtis 0.592.
- minco S10000 optimized recovered almost the same number of true species as KSSD2.L3K11S48 (205 vs 208), but with many more FP and substantially worse abundance reconstruction.
- minco S1000 optimized was lower recall: F1 0.686, TP 167, FP 64, FN 89.

Speed and memory:

- MetaKSSD README screenshot reports `marmgCAMI2 sample_0(10G)` runtime of 14 seconds and peak memory 0.5 GB.
- Sylph local run took 54.22 seconds to sketch and 159.40 seconds to profile, 213.62 seconds total. Peak RSS was 1.96 GB for sketching and 19.73 GB for profiling.
- minco S1000 plus-virus run on `/tmp/cami_marine_sample0_reads.fq.gz` took 223.60 seconds and 4.31 GB peak RSS.
- minco S10000 GTDB-only run took 187.76 seconds and 34.13 GB peak RSS.
- Relative to the MetaKSSD screenshot, minco was about 16.0x slower for S1000 and 13.4x slower for S10000 on sample0. S1000 used about 8.6x more peak memory, and S10000 used about 68x more peak memory.
- Relative to the same MetaKSSD screenshot, Sylph was about 15.3x slower end-to-end. If only counting the profile step after sketching, it was about 11.4x slower.

## Interpretation

MetaKSSD is still much faster and more memory efficient for CAMI marine sample0 profiling. Sylph gives the best local species-level detection among the tools compared here, but its runtime and memory are closer to minco than to the MetaKSSD screenshot. Current minco can approach MetaKSSD true-positive recovery with a permissive S10000 threshold, but it pays for that with higher false positives and worse abundance reconstruction.

This comparison is not fully apples-to-apples because the databases and scoring paths differ. The fair next step is to generate an official OPAL profile for minco using the same CAMI sample and upload/run it through OPAL, or run OPAL locally if available.

## Sylph TP Missing From minco

Question checked after the main comparison: since Sylph also uses GTDB, are Sylph true positives missing in minco because of GTDB-to-NCBI mapping, absent references, filtering, or ANI estimation?

Using the same OPAL-like 256-species gold set:

- Sylph TP: 222
- minco S10000 optimized-combo TP: 205
- Sylph TP absent from minco S10000 TP: 25
- minco S10000 TP absent from Sylph TP: 8

All 25 Sylph-only true-positive species had the exact Sylph reference accession present in both minco S10000 and S1000 sketches. They also had minco unfiltered candidate rows. Therefore this subset is not explained by missing references or taxmap absence.

The failure mode was ANI filtering:

- 24 of 25 failed only `ANI < 0.92` under the S10000 optimized-combo threshold.
- 1 of 25 failed `ANI < 0.92` plus `Ref_breadth < 0.05`.
- The best minco rows were usually the same accessions used by Sylph, but minco estimated ANI around 0.87-0.917 while Sylph estimated about 96.7-100.0 adjusted ANI.
- These were all very low-abundance gold species, mostly 0.0009-0.0183% in the CAMI profile.

Threshold tradeoff for S10000 with `Ref_breadth >= 0.05`, `XnY_ctx >= 10`, and varying ANI cutoff:

```text
ANI_cut  pred  TP   FP    FN   F1       recovered Sylph-only TP
0.92     305   205  100   51   0.7308   0/25
0.91     365   214  151   42   0.6892   8/25
0.90     440   225  215   31   0.6466   18/25
0.89     542   229  313   27   0.5739   22/25
0.88     823   233  590   23   0.4319   24/25
```

Interpretation: lowering minco's ANI cutoff recovers many Sylph-only TPs, but false positives rise faster than true positives. The main gap here is minco readwise ANI underestimation for low-abundance species, not GTDB reference absence.
