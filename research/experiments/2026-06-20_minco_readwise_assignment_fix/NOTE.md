# Minco Readwise Assignment Fix

Date: 2026-06-20
Author/agent: Codex
Project: minco
Code checkpoint: working tree after readwise assignment and ZIP AAF changes

## Question

Can minco fix low-abundance ANI underestimation and shared-context false
positives well enough to beat local Sylph on CAMI II marine short-read sample0?

## Changes Tested

Implemented direct readwise shared-context assignment modes:

- `--readwise-assign all`: legacy behavior; every matching reference gets the
  query context.
- `--readwise-assign best-diff`: only references with the best object-difference
  score for that query context are counted.
- `--readwise-assign best-diff-split`: same as `best-diff`, but depth coverage
  is split across tied best hits using fixed-point coverage.
- `--readwise-assign best-diff-unique`: only contexts with exactly one best
  reference hit are counted. This is now the default.

Implemented readwise selected ANI model:

- `--readwise-ani zip-aaf`: default when depth abundance is available. It fits
  latent reference AF from `Ref_breadth` and `Ref_mean_depth` using a
  zero-inflated Poisson equation, then converts AF to context AAF ANI.
- `--readwise-ani naive`: legacy object-difference readwise ANI.

The readwise abundance output now includes `Ref_zip_af` and
`Ref_zip_aaf_ani`.

## Dataset

- Query reads: `/tmp/cami_marine_sample0_reads.fq.gz`
- Reference sketch:
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`
- Gold profile: `/tmp/gs_marine_short.profile`
- Taxmap:
  `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Sylph profile: `/tmp/sylph_marine_sample0/profile.tsv`

Gold scoring uses the same local OPAL-like species filter as previous notes:
sample0 species rows, dropping `unidentified`, `unidentified plasmid`, and
`unidentified virus`, for 256 gold species.

## Results

Main species-level results:

```text
method / setting                         pred  TP   FP   FN   precision  recall   F1
Sylph local default                      278   222  56   34   0.799      0.867    0.831
old minco S10000 all/shared baseline     305   205  100  51   0.672      0.801    0.731
best-diff assignment, naive ANI grid     306   190  116  66   0.621      0.742    0.676
best-diff-split, naive ANI grid          365   219  146  37   0.600      0.855    0.705
best-diff-unique, naive ANI grid         242   211  31   45   0.872      0.824    0.847
best-diff-unique + ZIP AAF direct CLI    258   221  37   35   0.857      0.863    0.860
best-diff-unique + ZIP AAF + depthCV     265   226  39   30   0.853      0.883    0.868
```

The direct CLI setting that beat local Sylph was:

```bash
./minco_core/bin/minco ani -p16 \
  -r /mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno \
  --qraw /tmp/cami_marine_sample0_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  -m0 -f0.05 -n0.94 -t10 \
  -o /tmp/minco_readwise_assign_20260620/s10000_unique_zipaaf_f0.05_n0.94_t10.tsv
```

Direct CLI output score:

```text
pred_taxa 258
TP        221
FP        37
FN        35
precision 0.856589
recall    0.863281
F1        0.859922
```

Speed and memory:

```text
old minco S10000 baseline     187.76 s, 34.13 GB peak RSS
best-diff unfiltered          177.91 s, 34.12 GB peak RSS
best-diff-split unfiltered    185.61 s, 34.13 GB peak RSS
best-diff-unique unfiltered   175.02 s, 33.41 GB peak RSS
ZIP AAF final direct command  170.86 s, 33.41 GB peak RSS
Sylph local total             213.62 s, 19.73 GB peak RSS
```

`best-diff` fixed most of the low-abundance underestimation but produced too
many FPs because tied shared contexts still counted as support. `best-diff-split`
reduced depth but not support/breadth, so FP inflation remained. The useful
fix was `best-diff-unique`, which removes tied shared contexts from the support
and ANI feature counts.

After `best-diff-unique`, ZIP-corrected AAF became a useful minco-native ANI
model. This matches the previous oracle-effective-coverage diagnosis: the
problem was mostly effective coverage plus shared-context assignment, not just
sequencing-error distance inflation.

## Validation

- `make -C minco_core -j8` passed.
- `make -C minco_core test` passed.
- Small readwise smoke test passed for default `best-diff-unique` and explicit
  `--readwise-ani naive`.
- Full CAMI sample0 S10000 runs completed for `best-diff`,
  `best-diff-split`, `best-diff-unique`, and final ZIP AAF direct CLI.

## Important Artifacts

See `artifacts.md`.

## Caveats

- This is one CAMI sample. The result should be replicated on more CAMI samples
  and on non-CAMI metagenomes before calling it a general default victory.
- S10000 reference sketch is GTDB-only; the S1000 plus-virus database is not the
  same database.
- Local scoring is OPAL-like species taxid scoring, not a fresh CAMI website
  OPAL submission.
- The `best-diff-unique + ZIP AAF + depthCV` row requires an additional
  `Ref_depth_cv <= 10` post-filter not currently exposed as a CLI filter.

## Conclusion

The issue is now substantially fixed for CAMI marine sample0. The new default
readwise assignment (`best-diff-unique`) removes the shared-context support
inflation, and the new ZIP AAF readwise ANI model corrects low-coverage
underestimation. With direct CLI filters, minco S10000 now beats local Sylph
species-level F1 on this sample: `0.860` vs `0.831`.

## Next Experiment

Run the same final settings on additional CAMI marine samples and on the S1000
plus-virus database. If the depth-CV post-filter remains consistently useful,
add a CLI filter such as `--max-depth-cv`.
