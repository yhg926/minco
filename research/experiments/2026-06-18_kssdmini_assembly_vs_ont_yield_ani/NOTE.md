# KSSDmini assembly-vs-ONT-yield ANI for 2009K-1545

Date: 2026-06-18
Project: KSSD3mini

## Question

Estimate ANI between the matched `Salmonella enterica` assembly for sample
`2009K-1545` and ONT read subsets from `SRR19787631` at 5m to 250m yield.

## Dataset

- Reference assembly: `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz`
- Assembly source: ATB/Shovill-SPAdes assembly for Illumina BioSample `SAMN29253066`
- ONT run: `SRR19787631`, ONT BioSample `SAMN29253116`
- Species in metadata: `Salmonella enterica`
- Read subsets: 5m, 10m, 20m, 50m, 75m, 100m, and 250m

## Methods

Reference assembly was sketched once with:

```text
bin/kssd3mini_stage3_native sketch -p 4 --ctxmeta both
```

Read subsets were sketched two ways:

```text
No QC:   sketch -p 4 --conflict --ctxmeta both
readsQC: sketch -p 4 --conflict --readsQC --ctxmeta both
```

A no-conflict follow-up was also run by omitting `--conflict` for both modes.

ANI was run as raw-read query ANI:

```text
ani -r assembly_sketch --qraw read_sketch -n0 -f0 -t0 -s4 -m0
```

`-n0 -f0 -t0` disables filtering so low-overlap rows remain visible. `-s4`
selects the naive metric, which is the active raw-read metric for this mode.

## Results

Summary metrics are in `summary.tsv`.

Key readsQC ANI estimates:

```text
5m:   1.000000, shared contexts 499,  ref AF 0.0499
10m:  0.999131, shared contexts 1982, ref AF 0.1982
20m:  0.999054, shared contexts 4647, ref AF 0.4647
50m:  0.998820, shared contexts 989,  ref AF 0.0989
75m:  0.998749, shared contexts 4010, ref AF 0.4010
100m: 0.999238, shared contexts 6361, ref AF 0.6361
250m: 0.999270, shared contexts 9045, ref AF 0.9045
```

No-QC high-yield rows saturated at ANI 1.0 while shared context fraction dropped
to 5-15%, consistent with error-context inflation in raw ONT sketches. The
readsQC 100m and 250m rows have much stronger overlap and are the most credible
KSSDmini estimates from this run.

No-conflict follow-up:

```text
yield  conflict_readsQC_ANI  noconflict_readsQC_ANI  noconflict_readsQC_refAF
5m     1.000000              1.000000                 0.0499
10m    0.999131              0.999131                 0.1981
20m    0.999054              0.999054                 0.4642
50m    0.998820              0.998820                 0.0989
75m    0.998749              0.998748                 0.4001
100m   0.999238              0.999238                 0.6359
250m   0.999270              0.999270                 0.9040
```

For this dataset, removing conflicting context-objects from the readsQC sketches
changed ANI by at most 0.000001 and did not materially change the usable
high-yield estimates. No-conflict no-QC still saturated at ANI 1.0 for 75m+,
with even lower 250m reference AF than conflict-kept no-QC.

KSSD3A-sampled readsQC patch:

After changing KSSDmini to infer `--readsQC` ranges from a KSSD3A-style
fixed-reduction abundance sample, then applying that range before final
context-minhash selection, the inferred ranges matched KSSD3A on these FASTQs:

```text
5m:   mode=2  lower=2  upper=4
10m:  mode=2  lower=2  upper=3
20m:  mode=3  lower=2  upper=10
50m:  mode=6  lower=4  upper=16
75m:  mode=10 lower=5  upper=22
100m: mode=13 lower=6  upper=29
250m: mode=33 lower=15 upper=58
```

The patched readsQC sketch filled 10,000 entries for every yield. ANI against
the matched assembly became:

```text
yield  ANI       queryAF   refAF    shared_ctx
5m     0.999168  0.129804  0.1297   1297
10m    0.999094  0.320252  0.3199   3199
20m    0.999134  0.710101  0.7079   7079
50m    0.999269  0.853678  0.8518   8518
75m    0.999237  0.929480  0.9279   9279
100m   0.999271  0.960020  0.9581   9581
250m   1.000000  0.970246  0.9685   9685
```

This fixes the earlier 50m readsQC failure, where refAF had dropped to 0.0989.

## Caveats

- The assembly is matched by sample name but comes from the Illumina BioSample,
  not directly from the ONT BioSample.
- Low-yield rows have low shared context counts, so ANI can be unstable.
- KSSDmini raw-read ANI currently uses the naive raw-read metric.
- No external ground-truth ANI was recomputed in this run.
