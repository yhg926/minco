# CAMI II Marine Extra Samples

Date: 2026-06-20
Author/agent: Codex
Project: minco
Code checkpoint: `b55b4d3` plus working-tree readwise assignment and ZIP AAF changes

## Question

Does the sample0 minco S10000 recipe generalize to additional CAMI II marine
short-read samples, and does it still beat local Sylph species-level F1?

## Recipe

Minco command shape:

```bash
./minco_core/bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno \
  --qraw FASTQ.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o OUT.tsv
```

Sylph baseline:

```bash
/home/ubuntu/yihuiguang/bin/sylph sketch -t 16 -r FASTQ.fq.gz -d SYLPH_OUT
/home/ubuntu/yihuiguang/bin/sylph profile -t 16 \
  /mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb \
  SYLPH_OUT/FASTQ.fq.gz.sylsp -o profile.tsv
```

## Dataset

- Gold profile: `/tmp/gs_marine_short.profile`
- Tested sample IDs:
  - `marmgCAMI2_short_read_sample_0`
  - `marmgCAMI2_short_read_sample_1`
  - `marmgCAMI2_short_read_sample_2`
- Sample1 and sample2 were cached from the CAMI II marine short-read archives:
  - `marmgCAMI2_sample_1_reads.tar.gz`
  - `marmgCAMI2_sample_2_reads.tar.gz`
- Scoring used local species-level taxid matching with
  `research/experiments/2026-06-20_cami2_marine_extra_samples/scripts/score_profiles.py`.

## Results

Species-level detection:

```text
sample  method   gold  pred  TP   FP  FN  precision  recall  F1
0       minco    256   258   221  37  35  0.857      0.863   0.860
0       sylph    256   278   222  56  34  0.799      0.867   0.831
1       minco    300   303   263  40  37  0.868      0.877   0.872
1       sylph    300   323   264  59  36  0.817      0.880   0.848
2       minco    274   283   236  47  38  0.834      0.861   0.847
2       sylph    274   303   236  67  38  0.779      0.861   0.818
```

Mean over samples 0-2:

```text
method  precision  recall  F1     TP     FP     FN
minco   0.8528     0.8671  0.8599 240.0  41.3   36.7
sylph   0.7983     0.8695  0.8323 240.7  60.7   36.0
```

Timing and memory:

```text
sample  method  wall_seconds  peak_RSS_GB
0       minco   170.86        33.41
0       sylph   213.62        19.73
1       minco   170.43        33.80
1       sylph   216.07        20.13
2       minco   168.09        33.73
2       sylph   208.25        20.13
```

The minco advantage is mainly lower false positives at nearly identical recall.
Across samples 0-2, minco averaged `41.3` FP per sample versus Sylph `60.7`,
while Sylph averaged only `0.7` more TP per sample.

TP-only abundance correlation remained high for both tools. Minco Pearson was
0.981, 0.980, and 0.990 for samples 0-2; Sylph Pearson was 0.985, 0.991, and
0.996. Sylph had slightly lower TP-only abundance MAE on these local scores.

## Validation

- `score_profiles.py` passed `python3 -m py_compile`.
- Sample1 minco completed normally on 33,301,262 reads.
- Sample2 minco completed normally on 33,298,832 reads.
- Sylph sample1 and sample2 sketch/profile commands completed normally.

## Interpretation

The sample0 improvement replicated on two additional CAMI II marine samples.
The current minco S10000 recipe beat local Sylph F1 on all three tested samples:

- Sample0: `0.8599` vs `0.8315`
- Sample1: `0.8723` vs `0.8475`
- Sample2: `0.8474` vs `0.8180`

The result supports the conclusion that `best-diff-unique` shared-context
assignment is controlling false positives better than Sylph on this local
species-taxid scoring path, while ZIP AAF keeps recall close to Sylph.

## Caveats

- This is still local OPAL-like species-taxid scoring, not a fresh CAMI web OPAL
  evaluation.
- Minco and Sylph use different GTDB releases/databases:
  minco GTDBr232-derived S10000, Sylph GTDB-r226 c200.
- Minco is faster than Sylph on these runs, but uses substantially more memory:
  about 33.7 GB versus about 20.0 GB peak RSS.
- Only samples 0-2 were tested. CAMI marine samples 3-9 remain useful for a
  broader holdout check.

## Paper-Relevant Claim

On CAMI II marine short-read samples 0-2, minco S10000 with unique best-diff
readwise assignment and ZIP AAF ANI achieved higher local species-level F1 than
Sylph while running faster, mainly by reducing false positives. This claim still
needs official OPAL validation and more samples.
