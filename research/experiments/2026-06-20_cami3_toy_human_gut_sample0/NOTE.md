# CAMI III Toy Human Gut Sample0

Date: 2026-06-20
Project: minco
Status: preliminary non-marine benchmark

## Question

Does the current `S=1000` readwise profiling recipe generalize from CAMI marine
to a different dataset: CAMI III Toy Human Gut sample0?

## Dataset

CAMI III Toy Human Gut provides simulated short- and long-read human gut
metagenomes. The page reports 20 Illumina HiSeq short-read samples totaling
100 Gbp, with sample IDs `cami3_toy_human_gut_short_read_sample_[0..19]`.

This note tested short-read sample0:

- Read archive URL:
  `https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_0_reads.tar.gz`
- Read archive size by HTTP HEAD: `4,866,106,604` bytes
- Archive members:
  - `sample_0_reads/reads_mapping.tsv.gz`
  - `sample_0_reads/anonymous_reads.fq.gz`
- Gold profile:
  `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt`

Disk check before download:

```text
/tmp      9.8G free, 99% used
/mnt/new3T 441G free, 77% used
```

Therefore all large files were cached under
`/mnt/new3T/minco_cami3_toygut_20260620`, not `/tmp`.

Gold species composition:

```text
bacteria 118 species, 49.6574%
virus    135 species, 48.2064%
```

## Methods

Minco used the existing S1000 plus-virus reference:

```text
/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620
```

Current marine-tuned command:

```bash
tar -xOzf sample_0_reads.tar.gz sample_0_reads/anonymous_reads.fq.gz |
  minco ani -p16 -r REF --qraw - \
    --query-density ref --abundance-est depth --readwise-profile-only \
    --readwise-assign best-diff-unique \
    --readwise-ani zip-aaf \
    -m0 -f0.05 -n0.94 -t10 -o current.tsv
```

An unfiltered minco candidate table was also generated with `-f0 -n0 -t0`.

Sylph used the local GTDB-r226 c200 database:

```text
/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb
```

This Sylph baseline is not all-species comparable because the local Sylph DB is
GTDB-only and cannot represent the viral half of the sample.

A strict GTDB-only minco run was also added after the initial plus-virus run:

```text
/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno
```

An S10000 GTDB-only run was later added for the same sample:

```text
/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno
```

## Results

Runtime and memory:

```text
method                  wall time  peak RSS
minco current            4:12.37   3.80 GB
minco unfiltered         3:49.80   3.80 GB
minco S1000 GTDB-only    2:09.80   3.31 GB
minco S10000 current     5:40.83   29.65 GB
minco S10000 unfiltered  3:00.15   29.66 GB
Sylph sketch+profile     1:36.46   19.32 GB
```

Species-level scores:

```text
label            scope     TP   FP  FN   precision  recall  F1
minco_current    all      108  15  145   0.878      0.427   0.574
minco_current    bacteria  32   9   86   0.780      0.271   0.403
minco_current    virus     76   6   59   0.927      0.563   0.700

minco_best_grid  all      177  58   76   0.753      0.700   0.725
minco_best_grid  bacteria  63  21   55   0.750      0.534   0.624
minco_best_grid  virus    114  37   21   0.755      0.844   0.797

Sylph_default    all       55  13  198   0.809      0.217   0.343
Sylph_default    bacteria  55  13   63   0.809      0.466   0.591
Sylph_default    virus      0   0  135   0.000      0.000   0.000
```

Strict GTDB-only bacterial comparison:

```text
label                   TP  FP  FN  precision  recall  F1
S1000 current           32  10  86  0.762      0.271   0.400
S1000 best              64  22  54  0.744      0.542   0.627
S10000 current          33   9  85  0.786      0.280   0.413
S10000 best             60  17  58  0.779      0.508   0.615
Sylph_default           55  13  63  0.809      0.466   0.591
```

So minco did not beat Sylph with the current marine threshold on the GTDB-only
bacterial comparison. It did beat Sylph after applying the same post-hoc
toy-gut relaxed threshold. S10000 did not improve the result: its best
post-hoc bacterial F1 was `0.615`, below S1000's `0.627`, while using about
9x more memory.

The best simple post-hoc minco threshold in the unfiltered table was:

```text
XnY_ctx >= 1
Ref_breadth >= 0.005
ANI >= 0.90
Relative_abundance_depth >= 0
```

This improved all-species F1 from `0.574` to `0.725`, showing that the marine
thresholds are too strict for this dataset.

## Interpretation

This benchmark did not reproduce the CAMI marine result directly, but it
identified useful next improvements.

Key findings:

- The current marine-tuned threshold has high precision but poor recall on toy
  human gut.
- Viral detection is already reasonably strong after threshold relaxation:
  viral F1 `0.797`.
- Bacterial detection improves over Sylph under the best simple minco threshold
  (`0.624` vs `0.591` bacterial F1), but current thresholds miss many bacterial
  species.
- In strict GTDB-only bacterial scoring, minco also beats Sylph only after
  threshold relaxation (`0.627` vs `0.591` F1); the marine threshold is worse
  (`0.400` F1).
- Testing S10000 did not solve the toy-gut threshold/domain issue. S10000's
  best bacterial F1 was `0.615`, slightly worse than S1000's `0.627`, and its
  memory was much higher (`29.66 GB` vs `3.31 GB`).
- Minco uses much less memory than Sylph (`3.8 GB` vs `19.3 GB`) but is slower
  on this plus-virus mixed-reference run.

## Improvement Directions

1. Add dataset/domain-aware default thresholds for readwise profiling.
   The marine setting `-f0.05 -n0.94 -t10` is too strict here. A lower breadth
   and ANI threshold performs much better on toy-gut.
2. Add optional rank/domain filters or domain-specific reports. Viral and
   bacterial taxa behave differently and should not necessarily share one
   reporting threshold.
3. Improve plus-virus reference lookup speed. Toy-gut S1000 plus-virus was
   about 2x slower than marine S1000 GTDB-only despite similar read count.
4. Add a broader reference or separate eukaryotic/fungal mode if scoring all
   CAMI species, because toy-gut includes non-bacterial, non-viral entries.
5. Use unfiltered candidate diagnostics to distinguish true reference absence
   from threshold-filtered false negatives.

## Caveats

- This is only sample0.
- The best minco threshold is post-hoc and must not be treated as a general
  default until tested on more samples.
- Sylph comparison is database-limited because the local Sylph database is
  GTDB-only.
- Official CAMI/OPAL evaluation was not run.
