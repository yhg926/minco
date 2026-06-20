# CAMI II Marine Sample 0 Minco Profile

Date: 2026-06-20
Author/agent: Codex
Project: minco
Code commit: b92eed2 plus working-tree CAMI/profile-only changes
Tool version: minco 0.1

## Question

Can minco generate a CAMI-format taxonomic profile for CAMI II marine short-read sample 0 and complete OPAL evaluation through the CAMI website?

## Dataset

- CAMI marine short-read sample: `marmgCAMI2_sample_0_reads.tar.gz`
- Streamed member: `simulation_short_read/2018.08.15_09.49.32_sample_0/reads/anonymous_reads.fq.gz`
- Sample ID: `marmgCAMI2_short_read_sample_0`
- Gold profile used for local sanity checks: `/tmp/gs_marine_short.profile`
- Reference sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620`
- CAMI taxmap generated from GTDB metadata plus RefSeq viral taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.tsv`

## Methods

The original exact readwise mode was attempted first and was killed by memory pressure:

```text
reads processed before kill: 21,000,000
exit status: 137
max RSS: 58,808,256 KB
wall time: 8:55.48
```

I added `--readwise-profile-only`, which preserves reference breadth/depth abundance and CAMI profile output but skips exact global query-context and query-ref-context sets. The final full-sample run used:

```text
minco ani -p16 -r REF --qraw - --query-density ref --abundance-est depth \
  --readwise-profile-only --cami-taxmap TAXMAP --cami-profile PROFILE \
  --cami-sample-id marmgCAMI2_short_read_sample_0 -m0 -o TSV
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
Full sample reads processed: 33,294,790
Runtime: 6:20.24 wall
Max RSS: 4,174,092 KB
CAMI profile rows: 17 species rows plus header
CAMI upload: /submission/9089740f62614e468422/
OPAL run: 471
OPAL software: OPAL 1.0.14
```

OPAL species-rank result:

```text
Predicted taxa: 17
True positives: 16
False positives: 1
False negatives: 243
Completeness: 0.062
Purity: 0.941
F1: 0.116
L1 norm error: 1.746
Bray-Curtis distance: 0.887
```

Local species-level sanity check against sample 0 of the gold profile:

```text
Predicted species: 17
Gold positive species: 259
Overlap: 16
Pearson on rooted species after excluding unrooted placeholder taxa: 0.873
```

### S=1000 False-Negative Diagnosis And Threshold Optimization

To understand whether false negatives were missing from the reference DB or
filtered after candidate generation, I generated a full permissive candidate
table with:

```text
minco ani -p16 -r REF --qraw - --query-density ref --abundance-est depth \
  --readwise-profile-only -m0 -f0 -n0 -t0 -o UNFILTERED_TSV
```

Permissive run performance:

```text
Candidate rows: 209,380
Runtime: 5:59.01 wall
Max RSS: 4,173,212 KB
```

False-negative split against gold-positive species:

```text
Gold-positive species: 259
False-negative species: 243
FN not in current ref DB/taxmap: 17 taxa, 70.2333% abundance
FN present in current ref DB/taxmap: 226 species, 15.4642% abundance
```

Among the 226 represented false negatives, the primary blocker under the
current default low-abundance rule was:

```text
Primary blocker       species  gold abundance
ANI < 0.95            73       12.0255%
breadth < 0.5         149       3.2937%
XnY < 1000             4        0.1450%
```

Counting overlapping failed metrics, 225/226 represented false negatives had
best-candidate ANI below 0.95. The high-abundance missed species usually had
enough breadth and XnY support but were rejected by the 0.95 ANI cutoff.

The best measured threshold tradeoff for minimizing species-count `FP + FN` on
the S=1000 candidate table was:

```text
ANI >= 0.84
Ref_breadth >= 0.55
Relative_abundance_depth >= 0
XnY_ctx >= 100
Unique_ref_ctx_hit >= 10
```

Result:

```text
TP: 82
FP: 18
FN: 177
FP + FN: 195
Predicted species: 100
Recovered gold abundance: 22.8142%
```

A stricter near-optimal variant with stronger XnY support was:

```text
ANI >= 0.84
Ref_breadth >= 0.55
Relative_abundance_depth >= 0
XnY_ctx >= 1000
Unique_ref_ctx_hit >= 10
```

Result:

```text
TP: 81
FP: 18
FN: 178
FP + FN: 196
Recovered gold abundance: 22.7710%
```

For lower false positives, the best tested `FP <= 10` point was:

```text
ANI >= 0.84
Ref_breadth >= 0.70
Relative_abundance_depth >= 0
XnY_ctx >= 100
Unique_ref_ctx_hit >= 10

TP: 58
FP: 6
FN: 201
FP + FN: 207
Recovered gold abundance: 21.4994%
```

The structured threshold table is saved in
`threshold_optimization_s1000.tsv`.

### S=10000 GTDB-Only Follow-Up

I then repeated the same profiling and threshold-optimization workflow using
the existing indexed S=10000 GTDB reference:

```text
/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
```

This is not the same reference composition as the S=1000 appended DB: it has
200,709 GTDB samples and does not include the appended RefSeq viral references.
The S=10000 sketch-set metadata reports:

```text
universal_density: 0.07022037931096747
sample_count: 200,709
```

The original network streaming pipeline failed for S=10000 because minco had to
load the 23 GB reference index before reading stdin, causing the remote stream
to close early. I therefore cached the inner CAMI FASTQ gzip first:

```text
/tmp/cami_marine_sample0_reads.fq.gz
```

Default S=10000 GTDB-only run:

```text
Reads processed: 33,294,790
Runtime: 2:55.79 wall
Max RSS: 34,142,552 KB
Species rows: 8
TP: 8
FP: 0
FN: 251
FP + FN: 251
Recovered gold abundance: 10.2535%
```

The permissive S=10000 candidate run used:

```text
minco ani -p16 -r REF_S10000 --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  -m0 -f0 -n0 -t0 -o UNFILTERED_TSV_S10000
```

Permissive S=10000 run performance:

```text
Candidate rows: 200,709
Runtime: 2:55.99 wall
Max RSS: 34,136,604 KB
```

False-negative split against the actual S=10000 GTDB candidate universe:

```text
FN not in S=10000 GTDB candidates: 17 taxa, 70.2333% abundance
FN present in S=10000 GTDB candidates: 234 species, 16.2840% abundance
```

Among the 234 represented false negatives, the primary blocker under the
current default low-abundance rule was:

```text
Primary blocker       species  gold abundance
ANI < 0.95            86       13.1877%
breadth < 0.5         148       3.0963%
```

The exact selected-threshold verification showed that the S=1000 optimized
rule transfers to S=10000 with the same species-count objective:

```text
ANI >= 0.84
Ref_breadth >= 0.55
Relative_abundance_depth >= 0
XnY_ctx >= 100
Unique_ref_ctx_hit >= 100

TP: 82
FP: 18
FN: 177
FP + FN: 195
Predicted species: 100
Recovered gold abundance: 22.8142%
```

Using `XnY_ctx >= 1000` with the same other cutoffs gave the same result in
the exact row-level check. A more sensitive breadth-only point increased
abundance recall but added false positives:

```text
ANI >= 0.50
Ref_breadth >= 0.50
Relative_abundance_depth >= 0
XnY_ctx >= 0
Unique_ref_ctx_hit >= 0

TP: 95
FP: 32
FN: 164
FP + FN: 196
Recovered gold abundance: 23.4757%
```

The structured S=10000 table is saved in
`threshold_optimization_s10000_gtdb.tsv`.

### Why S=10000 Did Not Improve Default Recall

The first apparent S=1000 versus S=10000 comparison was not clean because the
S=1000 result used a GTDB+RefSeq-virus appended DB with universal density 1.0,
while the S=10000 result used a GTDB-only DB with universal density 0.0702.
The clean GTDB-only metadata are:

```text
S=1000 GTDB-only universal_density:  0.00689210289544917
S=10000 GTDB-only universal_density: 0.07022037931096747
```

Thus S=10000 really does sample about 10x more query context than S=1000 when
the reference composition is held constant.

I ran the S=1000 GTDB-only default profile on the same cached CAMI FASTQ:

```text
Runtime: 2:10.32 wall
Max RSS: 3,776,160 KB
Species rows: 17
TP: 17
FP: 0
FN: 242
FP + FN: 242
Recovered gold abundance: 11.1251%
```

Compared with S=10000 GTDB-only default:

```text
Species rows: 8
TP: 8
FP: 0
FN: 251
FP + FN: 251
Recovered gold abundance: 10.2535%
```

The nine species reported at S=1000 but lost at S=10000 did not lose breadth or
support. Their support increased about 10x, but their ANI estimates moved below
0.95. Every loss failed only the ANI cutoff:

```text
taxid    species                                  S1000 ANI  S10000 ANI
379547   Candidatus Aciduliprofundum boonei       0.950309   0.926618
430914   Halorhabdus tiamatea                     0.963928   0.925783
1124597  Magnetococcus marinus                    0.961168   0.906928
1006     Marivirga tractuosa                      0.967657   0.923252
160232   Nanoarchaeum equitans                    0.972191   0.945976
1737403  Nanohaloarchaea archaeon SG9             0.952765   0.944033
1564114  Rhodococcus sp. B7740                    0.956459   0.911901
1292     Staphylococcus warneri                   0.967957   0.938138
53633    Sulfobacillus acidophilus                0.963095   0.912721
```

The detailed row is saved in
`s1000_vs_s10000_gtdb_default_lost_taxa.tsv`.

Interpretation: increasing sketch size improved support and likely reduced
sampling noise, but the default `ANI >= 0.95` species-detection cutoff became
more conservative. In the S=1000 run, sampling noise appears to push several
borderline or non-identical references above 0.95. In S=10000, the ANI estimate
is more stable and lower, so those taxa are filtered. This supports using
breadth/depth as the primary metagenomic detection signal and a lower ANI gate
for profiling, rather than the assembled-genome species ANI cutoff.

### Density-Block Hypothesis Test

One possible mechanism was block pooling: with sparse S=1000 density and
`--density-block-ctx 100`, each pseudo-read block contains many reads, and a
block could potentially choose the closest object among multiple query objects.
To test this, I reran S=1000 GTDB-only permissive profiling with smaller blocks:

```text
--density-block-ctx 10
--density-block-ctx 0
```

The block counts were:

```text
setting          total blocks   reads/block
S1000 block=100      272,335    122.26
S1000 block=10     2,627,436     12.67
S1000 block=0     18,695,514      1.78
S10000 block=100   2,671,932     12.46
```

For the nine borderline taxa, ANI was essentially unchanged when moving S=1000
from `block=100` to `block=10` or exact per-read `block=0`:

```text
species                         S1000 b100  S1000 b10  S1000 b0   S10000 b100
Magnetococcus marinus           0.961168    0.961188   0.960838   0.906928
Marivirga tractuosa             0.967657    0.967575   0.967108   0.923252
Staphylococcus warneri          0.967957    0.967893   0.967684   0.938138
Nanoarchaeum equitans           0.972191    0.972161   0.971673   0.945976
Sulfobacillus acidophilus       0.963095    0.963116   0.963018   0.912721
```

This rejects block pooling as the main explanation. The detailed comparison is
saved in `s1000_density_block_borderline_ani.tsv`.

Across 242 gold species with candidates in both S=1000 block=0 and S=10000,
the S=10000 ANI was not globally lower:

```text
mean(S10000 - S1000 ANI):   -0.003198
median(S10000 - S1000 ANI):  0.0007645
S10000 lower:  109 species
S10000 higher: 133 species
```

Interpretation update: the apparent S=10000 ANI underestimation is mostly a
threshold-selection/regression-to-the-mean effect. We focused on taxa that
passed `ANI >= 0.95` at S=1000 but failed at S=10000, so they are enriched for
cases where the smaller S=1000 sample was optimistic. The denser S=10000 sample
does not generally lower ANI across all represented gold species.

### Naive ANI Feature Decomposition For The Nine Lost Taxa

Readwise profiling currently forces `ani_opt->v = true`, so these ANI values
come from the Naive distance:

```text
D = N_diff_obj
S = N_diff_obj_section
X = XnY_ctx
ratio = S / D
dist0 = D / (X + D)
final_dist = 1 - (1 - dist0)^ratio
naive_dist = 0.1544286 * final_dist + 0.0007133
ANI = 1 - naive_dist
```

Thus `N_mut2_ctx` is printed but does not directly enter the Naive ANI. The two
effective inputs are:

```text
D / XnY_ctx                    fraction of matched contexts with object mismatch
N_diff_obj_section / N_diff_obj average differing object sections per mismatching context
```

Using S=1000 `--density-block-ctx 0` as the baseline, I decomposed the S=10000
distance increase with a two-feature Shapley-style swap of `D/XnY_ctx` and
`N_diff_obj_section/D`.

Aggregate over the nine taxa:

```text
total distance increase:                  0.317043
D/XnY_ctx contribution:                   0.101820  (32.1%)
N_diff_obj_section/N_diff_obj contribution: 0.215223  (67.9%)
```

For every one of the nine taxa, the dominant driver was
`N_diff_obj_section / N_diff_obj`, not `D / XnY_ctx`. Example:

```text
Magnetococcus marinus:
  S1000 D/XnY:      0.190415
  S10000 D/XnY:     0.317505
  S1000 section/D:  1.642617
  S10000 section/D: 3.305567
  ANI:              0.960838 -> 0.906928
  contribution:     39.4% D/XnY, 60.6% section/D
```

The full decomposition is saved in
`s1000_s10000_naive_ani_feature_decomposition.tsv`.

## Validation

- `make test` passed after adding profile-only mode.
- Bounded 100k-read smoke test with `--readwise-profile-only` completed and produced a valid CAMI profile.
- CAMI website accepted the profile as taxonomic profiling for the Marine Illumina HiSeq dataset and sample 0.
- OPAL evaluation completed successfully.

## Caveats

- The output profile is species-rank only; higher taxonomic ranks are absent and scored as zero by OPAL.
- The GTDB+RefSeq-virus reference DB does not represent the dominant CAMI placeholder taxa `45202` (`unidentified plasmid`) and `32644` (`unidentified`), which together dominate sample 0. This strongly limits completeness and abundance accuracy.
- `--readwise-profile-only` reports `Unique_query_ctx` as 0 and approximates query AF from reference breadth; it is an abundance profiling mode, not a replacement for exact readwise ANI/AF.
- The S=1000/S=10000 threshold optimum is exploratory and currently based on
  one CAMI marine sample. The S=10000 run used a GTDB-only reference, while the
  S=1000 run used GTDB plus appended RefSeq viral references.

## Next Experiment

Build or append an S=10000 RefSeq-virus sketch if we want exact reference
composition parity with the S=1000 appended DB. Then repeat the same threshold
table and compare whether the densest viral references force universal density
back to 1.0.
