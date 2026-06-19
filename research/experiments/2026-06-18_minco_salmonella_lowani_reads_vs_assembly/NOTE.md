# minco Salmonella Low-ANI Read-to-Assembly Check

Date: 2026-06-18

## Question

For a same-species Salmonella enterica pair below 99% ANI, does a minco sketch from real ONT reads estimate ANI to the other assembly close to the assembly-vs-assembly truth?

## Inputs

- Assembly with ONT reads: `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies/2009K-1545.fa.gz`
- ONT reads: `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/{5m,10m,20m,50m,75m,100m,250m}/2009K-1545__SRR19787631__*.fastq.gz`
- Other same-species assembly: `/mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_001272115.1_Salmonella_enterica_CVM_N46840_v1.0_genomic.fna.gz`

The FASTQ sketches were the existing 2009K-1545 minco sketches made with `--conflict --readsQC --ctxmeta both`.

## Ground Truth

`fastANI` assembly-vs-assembly ANI:

- `2009K-1545` assembly vs `GCF_001272115`: `98.7777%` (`0.987777`)

## Results

minco assembly-vs-assembly ANI was `0.987683`, only `0.000094` below the `fastANI` truth.

Read-vs-assembly estimates approached the assembly truth as yield increased:

- `5m`: `0.985297`, error `0.002480`
- `10m`: `0.986880`, error `0.000897`
- `20m`: `0.987113`, error `0.000664`
- `50m`: `0.987192`, error `0.000585`
- `75m`: `0.987452`, error `0.000325`
- `100m`: `0.987447`, error `0.000330`
- `250m`: `0.987502`, error `0.000275`

Direction check (`assembly -> reads` vs `reads -> assembly`) produced the same ANI values and swapped query/ref AF as expected.

## Interpretation

For this below-99% same-species Salmonella pair, minco read-to-assembly estimates are close to assembly truth after moderate yield. The remaining high-yield bias is about `0.00028` ANI against `fastANI`, or `0.00018` ANI against minco's own assembly-vs-assembly estimate.

## Caveats

- Only one side has real ONT reads in this test.
- The assembly truth here is `fastANI`, not full ANIm.
- This uses minco `--readsQC`; no no-QC FASTQ control was rerun for this pair.
