# minco Cross-Read ANI Against Assembly Truth

Date: 2026-06-18
Project: minco

## Question

Use assembly-to-assembly ANI as ground truth for two same-species
`Salmonella enterica` samples, then test whether minco raw-read ANI is close
to that truth for:

- reads from one sample against the other sample's assembly
- reads from one sample against the other sample's reads

## Dataset

- Sample A: `2009K-1545`, ONT run `SRR19787631`
- Sample B: `2010K-1947`, ONT run `SRR19787632`
- Assembly A: `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz`
- Assembly B: `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2010K-1947.fa.gz`
- ONT yield subsets: 5m, 10m, 20m, 50m, 75m, 100m, and 250m

Ground truth from
`/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/metadata/Salmonella_ATB_asm_ani.Tf4.clean.tsv`:

```text
2010K-1947 assembly vs 2009K-1545 assembly ANI = 0.999280
```

## Methods

minco binary:

```text
/home/ubuntu/yihuiguang/tools/minco/bin/minco_stage3_native
```

Read sketches used patched `--readsQC`, where QC count ranges are inferred from
a KSSD3A-style fixed-reduction abundance sample and then applied before final
fixed-size context-minhash selection.

ANI was run with:

```text
ani -r REF_SKETCH --qraw READ_SKETCH -n0 -f0 -t0 -s4 -m0
```

`-n0 -f0 -t0` keeps all rows visible. `-s4` selects the naive raw-read metric.

## Results

Main summary is in `summary.tsv`; aggregate errors are in `metrics.tsv`.

High-yield rows, 50m and above:

```text
comparison                         MAE       max_abs_error
2009 reads -> 2010 assembly         0.000081  0.000110
2010 reads -> 2009 assembly         0.000030  0.000045
2009 reads -> 2010 reads, matched   0.000039  0.000047
2010 reads -> 2009 reads, matched   0.000039  0.000047
```

Matched-yield read-to-read ANI was symmetric in ANI; the reverse direction only
swapped query/reference alignment fractions.

## Interpretation

For this close same-species pair, patched minco raw-read ANI tracks the
assembly-to-assembly truth well once yield reaches 50m. Low yields are noisier,
especially read-to-read at 5m where the error was 0.001066 ANI.

## Caveats

- Ground truth here is the available ATB assembly-to-assembly ANI table, not a
  newly recomputed ANIm run.
- This is one close Salmonella pair. More pairs and wider ANI ranges are needed
  before treating the error rates as general.
- Read-to-read results reported here use matched yield rows from the 7x7
  read-to-read output matrix.
