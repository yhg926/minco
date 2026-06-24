# CAMI III Toy Gut Source-Genome ANI Diagnosis

Date: 2026-06-20
Project: minco
Code commit: b55b4d3
Status: diagnostic

## Question

For CAMI III Toy Human Gut sample0, several bacterial taxa are true positives
for Sylph but fail minco's current default because minco reports ANI below the
default cutoff. Are those references truly distant from the simulated source
genomes, or is minco's readwise ANI model underestimating ANI?

## Hypothesis

Minco's current readwise `zip-aaf` ANI is underestimating real source-vs-ref ANI
because it converts sparse context breadth directly to ANI. In high-depth gut
taxa, missing context breadth is not mainly a zero-inflated coverage problem.

## Dataset

CAMI III Toy Human Gut sample0:

```text
/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz
/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt
```

New provenance files used in this diagnosis:

```text
/mnt/new3T/minco_cami3_toygut_20260620/sample_0_reads_mapping.tsv.gz
/mnt/new3T/minco_cami3_toygut_20260620/gsa_pooled_mapping.tsv.gz
/mnt/new3T/minco_cami3_toygut_20260620/source_genomes.tar.gz
```

The sample read mapping directly provides:

```text
anonymous_read_id    genome_id    tax_id    read_id
S0R1000/1            ASV507.1     853       NZ_QVFB01000020.1-965/1
```

The pooled GSA mapping uses the same source genome IDs as `BINID`, allowing
ASV IDs to be matched to source FASTA files by contig accession.

## Methods

I focused on the largest Sylph-TP/minco-current-missing bacterial taxa where
the missing reason was ANI-related:

```text
853     Faecalibacterium prausnitzii
626940  Phascolarctobacterium succinatutens
40520   Blautia obeum
818     Bacteroides thetaiotaomicron
```

For these taxa, minco and Sylph selected the same GTDB reference accession:

```text
F. prausnitzii        GCF_003324185.1
P. succinatutens      GCF_000188175.1
B. obeum              GCF_025147765.1
B. thetaiotaomicron   GCF_000011065.1
```

That makes the comparison direct: source genome vs selected reference.

Source-genome counts from the read mapping, one row per read mate:

```text
818     ASV102.10  1432884
818     ASV153.8    437894
818     ASV99.3      13590
853     ASV507.1    960068
853     ASV507.0    593842
853     ASV507.3    407460
40520   ASV406.4   1053872
626940  ASV579.0    763090
```

The exact source FASTA mappings were identified by contig accession:

```text
ASV507.1   Faecalibacterium_prausnitzii_AM37_13AC_strain_genomic.fna
ASV507.0   Faecalibacterium_prausnitzii_2789STDY5834970_genomic.fna
ASV507.3   Faecalibacterium_prausnitzii_AM100_B16A_strain_genomic.fna
ASV579.0   Phascolarctobacterium_succinatutens_NB2A_12_FMU_strain_genomic.fna
ASV406.4   Blautia_obeum_AM37_4AC_strain_genomic.fna
ASV153.8   Bacteroides_thetaiotaomicron_E8_m1001262Bd2_200113_strain_genomic.fna
ASV102.10  Bacteroides_thetaiotaomicron_1001095st1_D10_1001095H_141210_strain_genomic.fna
ASV99.3    Bacteroides_thetaiotaomicron_D52t1_170925_B7_strain_genomic.fna
```

Truth ANI was measured with `fastANI`; `skani dist` was run as a second check.
Commands are recorded in `commands.sh`.

## Results

Detailed rows are in `summary.tsv`.

Read-count weighted fastANI truth, compared with minco and Sylph:

```text
species                       truth ANI  minco ANI  Sylph ANI
Faecalibacterium prausnitzii   95.97%     89.97%     98.22%
Phascolarctobacterium succ.    98.70%     93.40%     99.32%
Blautia obeum                  98.36%     92.13%     98.58%
Bacteroides thetaiotaomicron   97.94%     91.35%     97.75%
```

Skani confirmed the same direction. For example:

```text
P. succinatutens source vs GCF_000188175.1: skani ANI 99.13%
B. obeum source vs GCF_025147765.1:       skani ANI 98.60%
B. thetaiotaomicron source ASV102.10 vs
  GCF_000011065.1:                         skani ANI 98.36%
F. prausnitzii source ASV507.0 vs
  GCF_003324185.1:                         skani ANI 94.77%
```

For three of the four taxa, Sylph is close to source-genome truth and minco is
clearly too low. For `F. prausnitzii`, the source strains are mixed and one
major source strain is below 95% ANI to the selected reference, so Sylph's
98.22% is probably high; minco is still too low.

Minco feature table for the same selected references:

```text
species/ref                   minco ANI  XnY  breadth  mean depth  hit depth  diff_obj
P. succinatutens               0.933980  234  0.234    8.256       35.28      5
B. obeum                       0.921291  177  0.177    5.891       33.28      9
B. thetaiotaomicron            0.913463  149  0.149    3.270       21.95      5
F. prausnitzii                 0.899669  110  0.110    3.864       35.13      4
```

The current ZIP model did not correct AF here:

```text
Ref_zip_af == observed Ref_breadth
```

## Interpretation

Minco failed the ANI cutoff mainly because `zip-aaf` is using observed context
breadth as if it were nucleotide containment/aligned fraction. That is not true
for these readwise context sketches.

The zero-inflated coverage model cannot fix these rows because average depth is
already high. With observed breadth `b` and mean depth `d`, the current model
solves:

```text
observed_breadth = latent_af * (1 - exp(-d / latent_af))
```

When `d` is much larger than `b`, the model's maximum-likelihood answer is
effectively `latent_af = observed_breadth`; increasing latent AF would predict
many more hit contexts than observed. Therefore the missing breadth is not
being interpreted as incomplete sampling; it is being interpreted as real
absence.

The small object-difference counts show a separate signal:

```text
P. succinatutens:      5 diff objects / 234 hits
B. obeum:             9 / 177
B. thetaiotaomicron:  5 / 149
F. prausnitzii:       4 / 110
```

So the per-hit object identity signal is high. The failure is mostly in the
conversion from sparse context breadth to ANI, not in the object-difference
counting itself.

## Readwise Controls

Additional controls were run after the initial source-genome truth check.

First, direct assembled-genome minco ANI was run for the same source/ref pairs
with `-S1000 -m0 -f0 -n0 -t0`. This did not reproduce the low readwise ANI:

```text
species                       fastANI truth  minco assembled ANI
F. prausnitzii ASV507.1        96.44%         98.65%
P. succinatutens ASV579.0      98.70%         99.41%
B. obeum ASV406.4             98.36%         99.10%
B. thetaiotaomicron ASV102.10 97.99%         99.04%
```

Thus the observed low minco ANI is not an inherent low bias in the
assembled source-vs-reference minco comparison. If anything, the assembled
default model overestimates several of these pairs.

Second, an error-free perfect FASTQ was generated from
`P. succinatutens ASV579.0` by tiling 150 bp reads every 50 bp and repeating the
tiling 10 times. This removes sequencing errors while keeping fragmented,
high-depth readwise input.

The same readwise density path reproduced the low ANI:

```text
input/model                         ANI       Ref_breadth  Ref_zip_af
real Toy Gut reads, zip-aaf          0.933980  0.234000     0.234000
perfect 150 bp tiled reads, zip-aaf  0.933198  0.230000     0.230006
perfect 150 bp tiled reads, naive    0.996001  0.230000     0.230006
```

So sequencing error is not the main cause of the low ANI in this case. The
readwise `zip-aaf` selected ANI is biased low because high-depth short reads
observe only sparse reference context breadth, and the ZIP model does not
reinterpret that missing breadth as high-ANI sequence divergence. The naive
object-difference signal remains high on the same perfect reads.

## Conclusion

The truth-source diagnosis confirms a real minco ANI-model problem in CAMI Toy
Gut readwise profiling. For several high-abundance bacterial taxa, true
source-vs-reference ANI is about 98-99%, but minco reports about 91-93% because
the current `zip-aaf` model treats context breadth as aligned fraction.

The deeper answer is: the low bias is readwise-model specific, not a
source-vs-reference assembled minco low bias, and not primarily caused by
sequencing errors.

Sylph is closer for most of these taxa, but not perfect: for mixed-strain
`F. prausnitzii`, source-genome truth is about 96% weighted ANI, while Sylph
reports 98.22%.

## Caveats

- This diagnosis covers four high-impact ANI-missing bacterial taxa, not every
  missed taxon in Toy Gut sample0.
- fastANI and skani are assembly-to-assembly references; they do not model read
  sampling directly.
- Species-level truth may be a mixture of multiple source strains. Weighted
  truth should be preferred over a single-source comparison when a species has
  multiple source ASVs.

## Next Experiment

Train or fit a readwise ANI calibration model using source-genome truth for all
Toy Gut source/ref hits. Predict real ANI from minco features such as observed
breadth, mean depth, hit mean depth, depth CV, XnY, and object-difference
features. The model should not rely on raw context AAF alone.
