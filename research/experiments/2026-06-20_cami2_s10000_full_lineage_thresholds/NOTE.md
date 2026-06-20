# CAMI II Marine Sample 0 S10000 Full-Lineage Threshold Tuning

Date: 2026-06-20

Code checkpoint: `4b3de14` (`Record CAMI full-lineage threshold experiment`).

## Question

Repeat the CAMI full-lineage profile and threshold tuning for the S10000 reference sketch, after the S1000 experiment showed a better sample-0 threshold than the current default.

## Important Difference From S1000

This experiment used the existing S10000 GTDB-only refdb:

```text
/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
```

I did not find an S10000 GTDB+RefSeq-virus appended refdb. Therefore S10000 cannot recover viral taxa in this run, while the previous S1000 experiment used:

```text
sketch_T_S1000_anno_plus_refseq_virus_20260620
```

Interpret superkingdom and viral false negatives with that database difference in mind.

## Inputs

- Reference sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`
- Cached CAMI II marine short-read sample 0 FASTQ: `/tmp/cami_marine_sample0_reads.fq.gz`
- Full-lineage taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Gold profile: `/tmp/gs_marine_short.profile`

The same full-lineage taxmap can be used because matching is accession-based; unused virus keys are harmless.

## Runs

Default filtered profile:

```bash
/usr/bin/time -v ./bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --cami-taxmap /tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv \
  --cami-profile /tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_default.full_lineage.profile \
  --cami-sample-id marmgCAMI2_short_read_sample_0 \
  -m0 -o /tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_default.full_lineage.tsv
```

Unfiltered candidate table:

```bash
/usr/bin/time -v ./bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 \
  -o /tmp/minco_cami_lineage_s10000_20260620/s10000_gtdb_unfiltered.tsv
```

## Performance

Default run:

- Reads processed: 33,294,790
- Wall time: 3:07.76
- Peak RSS: 34,133,416 KB

Unfiltered run:

- Reads processed: 33,294,790
- Wall time: 3:06.04
- Peak RSS: 34,143,756 KB
- Candidate rows mapped to taxids: 199,923
- Unmapped candidate rows: 786

## Threshold Grid

Grid dimensions:

- Support cut: `10,25,50,100,200,500,1000,2000,5000,10000`
- Reference breadth: `0.05,0.1,0.2,0.3,0.5,0.7,0.9`
- ANI: `0.90,0.92,0.94,0.95,0.96,0.97,0.98,0.99`
- Relative depth: `0,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4`

Best species-only FP+FN:

```text
-f0.05 -n0.94 -t10
```

Best combined species+genus FP+FN:

```text
-f0.05 -n0.92 -t10
```

## Results

Current default post-hoc filter:

```text
support >= 100, Ref_breadth >= 0.5, ANI >= 0.96, Relative_abundance_depth >= 1e-5
```

- Species: TP 89, FP 13, FN 170, precision 0.8725, recall 0.3436, F1 0.4931
- Genus: TP 79, FP 6, FN 102, precision 0.9294, recall 0.4365, F1 0.5940

Optimized combined species+genus filter:

```text
support >= 10, Ref_breadth >= 0.05, ANI >= 0.92, Relative_abundance_depth >= 0
```

- Species: TP 205, FP 100, FN 54, precision 0.6721, recall 0.7915, F1 0.7270
- Genus: TP 156, FP 23, FN 25, precision 0.8715, recall 0.8619, F1 0.8667

Optimized species-only filter:

```text
support >= 10, Ref_breadth >= 0.05, ANI >= 0.94, Relative_abundance_depth >= 0
```

- Species: TP 171, FP 62, FN 88, precision 0.7339, recall 0.6602, F1 0.6951
- Genus: TP 137, FP 18, FN 44, precision 0.8839, recall 0.7569, F1 0.8155

## Interpretation

S10000 improves the best sample-0 tradeoff over S1000:

- S1000 optimized species F1 was 0.6816; S10000 combined optimized species F1 is 0.7270.
- S1000 optimized genus F1 was 0.8036; S10000 combined optimized genus F1 is 0.8667.

The best S10000 filter is more permissive than the current default. It sharply improves recall but increases species FP from 13 to 100. This is appropriate for an exploratory detection profile but still too aggressive to install as the universal default without validating more samples and an S10000 plus-virus refdb.

## Caveats

- S10000 here is GTDB-only; virus recall cannot be fairly compared to S1000 plus-virus.
- The CAMI example gold profile contains multiple sample sections; scoring used only `marmgCAMI2_short_read_sample_0`.
- Some CAMI gold taxa are unidentified/plasmid taxa not expected to be in the GTDB refdb.
- This is still one sample. Do not change the default threshold from this alone.

## Next

1. Build S10000 GTDB+RefSeq-virus to make the comparison with S1000 fair.
2. Run the same threshold grid on more CAMI samples.
3. Consider exposing a `--profile-sensitive` preset rather than changing the default.
