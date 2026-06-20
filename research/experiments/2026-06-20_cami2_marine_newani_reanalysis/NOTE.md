# CAMI Marine New-ANI Reanalysis

Date: 2026-06-20
Code base: minco, `b92eed2` plus working-tree unique-best readwise ANI changes
Binary: `bin/minco_stage3_native`

## Question

After changing readwise abundance ANI to use the best observed object difference
per unique reference context entry, how do the previous CAMI marine sample 0
results change?

## Dataset And References

- Query: `/tmp/cami_marine_sample0_reads.fq.gz`
- Gold profile: `/tmp/gs_marine_short.profile`
- Gold species baseline: first occurrence per species taxid in that file, which
  matches the previous 259-species sample-0 denominator.
- S1000 GTDB reference:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
- S10000 GTDB reference:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`
- S1000 GTDB plus RefSeq virus reference:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620`

## Methods

I reran default-filtered and permissive `-f0 -n0 -t0` readwise abundance
profiling for the GTDB-only S1000 and S10000 references. I also reran the
default-filtered S1000 GTDB plus RefSeq virus profile used by the earlier CAMI
upload. Commands are in `commands.sh`.

The compiled binary was checked for hardware popcount instructions:

```text
objdump -d bin/minco_stage3_native | rg -n "popcnt|vpopcnt"
```

The binary contains multiple `popcnt` instructions, so the object-difference
count path is using the intrinsic route inherited from KSSD3.

## Results

Default GTDB-only comparison:

```text
Reference   Rows  Pred species  TP  FP  FN   Recovered gold abundance
S1000       11    10            10  0   249  13.2237%
S10000      129   121           93  28  166  23.3463%
```

Compared with the previous occurrence-weighted ANI:

```text
Reference   Old TP/FP/FN       New TP/FP/FN
S1000       17 / 0 / 242       10 / 0 / 249
S10000      8 / 0 / 251        93 / 28 / 166
```

The old problem where S10000 lost S1000 default calls due to ANI
underestimation disappeared:

```text
S1000 default species lost at S10000 default: 0
S10000 default species gained vs S1000:       111
Common default species:                       10
Mean S10000-S1000 ANI on common species:      0.0000246
```

The S1000 plus RefSeq virus default run produced the same species-level calls as
S1000 GTDB-only:

```text
Rows 11, predicted species 10, TP 10, FP 0, FN 249
Runtime 3:40.53 wall, max RSS 4,292,184 KB
```

Best tested FP+FN grid points on permissive candidate tables:

```text
Reference   TP   FP  FN  FP+FN  ANI  breadth  rel_depth  XnY  ctxhit
S1000       167  63  92  155    .94  0        0          0    10
S10000      171  62  88  150    .94  0        0          0    0
```

Low-FP grid examples:

```text
Reference   FP cap  TP  FP  FN   FP+FN  rule
S1000       10      80  9   179  188    ANI>=0.97, ctxhit>=10
S10000      10      67  10  192  202    breadth>=0.65
S10000      20      109 14  150  164    ANI>=0.96
```

## Interpretation

The new unique-best ANI fixes the specific S10000 underestimation symptom: no
S1000 default species are lost at S10000. However, it also makes readwise ANI
much more permissive. S10000 default recall increased substantially, but 28
false positives appeared under the current default rule.

This means the unique-best rule is not a safe standalone default for metagenomic
profiling. It is better than occurrence-weighted mutation counts for avoiding
coverage-driven ANI underestimation, but it can overestimate ANI because a
single exact object observation can dominate a reference context entry.

## Caveats

- The sample-0 gold set is reconstructed from the first species occurrence in
  the CAMI example gold profile because that file repeats species across
  samples under one header.
- The threshold grid is exploratory and species-count optimized, not abundance
  optimized.
- These results are for one CAMI marine short-read sample.

## Next Step

Keep the packed coverage plus best-diff data structure, but expose or test a
less permissive mutation summary for readwise metagenomes, such as coverage-
weighted best-diff with a minimum evidence count, median/min-after-QC diff, or
separate exact-hit fraction. The default filters should be recalibrated after
choosing that mutation summary.
