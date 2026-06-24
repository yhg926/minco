# CAMI II Marine S1000 Choice

Date: 2026-06-20
Author/agent: Codex
Project: minco
Code checkpoint: `b55b4d3` plus working-tree readwise assignment and ZIP AAF changes

## Question

Can the practical S1000 GTDB-only reference sketch match the S10000 recipe on
CAMI II marine samples 0-2 after the `best-diff-unique` and ZIP AAF fixes?

## Method

Used the same detection recipe as the S10000 runs:

```bash
./minco_core/bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno \
  --qraw FASTQ.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o OUT.tsv
```

Scoring reused the local species-level taxid scorer from
`2026-06-20_cami2_marine_extra_samples`.

## Results

Per-sample detection:

```text
sample  method   TP   FP  FN  precision  recall  F1
0       S1000    223  36  33  0.861      0.871   0.866
0       S10000   221  37  35  0.857      0.863   0.860
0       Sylph    222  56  34  0.799      0.867   0.831

1       S1000    261  39  39  0.870      0.870   0.870
1       S10000   263  40  37  0.868      0.877   0.872
1       Sylph    264  59  36  0.817      0.880   0.848

2       S1000    236  46  38  0.837      0.861   0.849
2       S10000   236  47  38  0.834      0.861   0.847
2       Sylph    236  67  38  0.779      0.861   0.818
```

Mean over samples 0-2:

```text
method   TP     FP     FN     precision  recall  F1      time_s  RSS_GB
S1000    240.0  40.3   36.7   0.856      0.867   0.862   124.5   3.70
S10000   240.0  41.3   36.7   0.853      0.867   0.860   169.8   33.65
Sylph    240.7  60.7   36.0   0.798      0.870   0.832   212.6   20.00
```

S1000 had nearly identical sensitivity to S10000 and slightly fewer false
positives on average. It also reduced mean runtime by about `26.7%` and peak
memory by about `89%` relative to S10000.

## Interpretation

For CAMI marine samples 0-2, there is no evidence that S10000 is needed for the
current readwise profiling recipe. S1000 is the better practical default for
this benchmark: it matches or slightly exceeds S10000 F1, keeps the Sylph F1
advantage, and uses much less memory.

This also means the useful signal in the new method is not simply denser
sampling. The `best-diff-unique` assignment and ZIP AAF model are carrying most
of the gain; S1000 has enough sampled contexts for these samples.

## Caveats

- Marine samples 0-2 are related benchmark replicates; a different environment
  or host-associated dataset is needed before making S1000 the universal
  default.
- This used the GTDB-only S1000 sketch for fair comparison to the GTDB-only
  S10000 run. The plus-virus S1000 database may have different FP behavior.
- This is local OPAL-like taxid scoring, not a CAMI web OPAL evaluation.

## Next Step

Use S1000 as the default benchmark target on a different dataset. If S1000 holds
there, prioritize memory/index cleanup and official OPAL profile generation
over further S10000 tuning.
