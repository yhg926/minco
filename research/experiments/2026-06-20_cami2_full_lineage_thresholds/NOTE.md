# CAMI II Marine Sample 0 Full-Lineage Taxmap and Threshold Tuning

Date: 2026-06-20

Code checkpoint before this experiment: `a9c98c1` (`Add CAMI lineage profile support`).

## Question

Can minco emit a full-lineage CAMI profile for the existing S1000 GTDBr232 + RefSeq-virus reference sketch, and what simple reporting thresholds improve CAMI sample-0 species/genus detection relative to the current default?

## Inputs

- Reference sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620`
- Cached CAMI II marine short-read sample 0 FASTQ: `/tmp/cami_marine_sample0_reads.fq.gz`
- Original species-only taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.tsv`
- NCBI taxdump used for lineage expansion: `/mnt/new3T/gtdbr220/gtdbr226/taxdump`
- Gold profile: `/tmp/gs_marine_short.profile`

## Taxmap Expansion

The original species-only taxmap had 214,952 rows and no duplicate reference keys. I expanded it to a multi-row CAMI taxmap:

- Output: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Rows: 1,147,634
- Unique reference keys: 214,952
- Maximum rows per key: 7
- Species rows: 214,952
- Superkingdom rows: 214,837
- Taxids with missing `taxidlineage.dmp` entries: 177
- Refs where source taxpath could recover Bacteria/Archaea/Viruses superkingdom: 89
- Refs where top-level fallback was root/unknown and therefore intentionally not emitted as superkingdom: 115

This keeps CAMI rank semantics: root is not reported as a `superkingdom`.

## Full-Lineage Profile Run

Command shape:

```bash
/usr/bin/time -v ./bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620 \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --cami-taxmap /tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv \
  --cami-profile /tmp/minco_cami_lineage_20260620/s1000_plusvirus_default.full_lineage.profile \
  --cami-sample-id marmgCAMI2_short_read_sample_0 \
  -m0 -o /tmp/minco_cami_lineage_20260620/s1000_plusvirus_default.full_lineage.tsv
```

Measured run:

- Reads processed: 33,294,790
- Wall time: 3:43.60
- Peak RSS: 4,306,792 KB
- Output rows: 104 detail rows in the post-hoc current-default profile

## Threshold Grid

I reran the current binary once with permissive filters:

```bash
/usr/bin/time -v ./bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620 \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 \
  -o /tmp/minco_cami_lineage_20260620/s1000_plusvirus_unfiltered.tsv
```

Measured run:

- Reads processed: 33,294,790
- Wall time: 3:41.75
- Peak RSS: 4,302,884 KB
- Candidate rows mapped to taxids: 208,749
- Unmapped candidate rows: 631

Grid dimensions:

- Support cut: `10,25,50,100,200,500,1000`
- Reference breadth: `0.05,0.1,0.2,0.3,0.5,0.7,0.9`
- ANI: `0.90,0.92,0.94,0.95,0.96,0.97,0.98,0.99`
- Relative depth: `0,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4`

Best simple species+genus tradeoff on this sample:

```text
-f0.1 -n0.94 -t10
```

This is not installed as a new default yet. It was selected from one CAMI sample and should be validated on more samples before changing CLI defaults.

## Results

Current default post-hoc filter:

```text
support >= 100, Ref_breadth >= 0.5, ANI >= 0.96, Relative_abundance_depth >= 1e-5
```

Species: TP 85, FP 14, FN 174, precision 0.8586, recall 0.3282, F1 0.4749.
Genus: TP 76, FP 7, FN 105, precision 0.9157, recall 0.4199, F1 0.5758.

Optimized sample-0 filter:

```text
support >= 10, Ref_breadth >= 0.1, ANI >= 0.94, Relative_abundance_depth >= 0
```

Species: TP 167, FP 64, FN 92, precision 0.7229, recall 0.6448, F1 0.6816.
Genus: TP 135, FP 20, FN 46, precision 0.8710, recall 0.7459, F1 0.8036.

The optimized filter reduces species FP+FN from 188 to 156 and genus FP+FN from 112 to 66 on CAMI marine sample 0. It trades precision for much higher recall.

## Caveats

- The CAMI example gold profile contains multiple sample sections; all scoring here uses only `marmgCAMI2_short_read_sample_0`.
- Some gold taxa are `unidentified`, plasmid, or otherwise not represented by the current refdb; threshold tuning cannot recover those.
- Taxdump is from the local GTDBr226 area, not necessarily identical to the CAMI/NCBI taxonomy used when the gold profile was generated.
- This is one sample and one sketch size. Do not change the minco default threshold from this alone.

## Next

1. Repeat threshold tuning for S10000 and more CAMI samples.
2. Separate "not in refdb" false negatives from "filtered out" false negatives by rank.
3. Upload the optimized full-lineage profile to CAMI/OPAL if website-side scoring is needed.
