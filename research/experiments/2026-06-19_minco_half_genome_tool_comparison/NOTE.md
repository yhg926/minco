# minco Half-Genome Tool Comparison

Date: 2026-06-19
Author/agent: Codex
Project: minco
Code commit: 763a226
minco binary/tool version: minco 0.1

## Question

For a same-species pair, how do minco default ANI, minco context-Mash ANI,
skani, and Mash behave when one genome is replaced by its first or second
half?

## Hypothesis

Methods that explicitly account for aligned fraction should keep ANI close to
regional ANIm after truncation, while pure Mash-style fixed-denominator
similarity should drop because only about half the genome is present.

## Dataset

- Input data: Nayfach genome pair `P000881414`.
- Reference genome: `/mnt/new3T/skani_data/Nayfach_data/fna/3300027414_1.fna`.
- Query genome: `/mnt/new3T/skani_data/Nayfach_data/fna/3300014912_1.fna`.
- Full-pair ANIm label: `0.9890599554743553` from the Nayfach ANIm table.
- Query split: concatenate all 453 query records, then split into two equal
  sequence halves: 3,242,967 bp and 3,242,968 bp.
- Storage location: `research/experiments/2026-06-19_minco_half_genome_tool_comparison`.

## Methods

Key parameters:

```text
ANIm truth: dnadiff 1-to-1 AvgIdentity
minco default: minco ani -S 10000 -p2 -f0 -n0 -t0 -s1
minco context-Mash: minco ani -S 10000 -p2 -f0 -n0 -t0 -s-5
skani: skani dist -t2
Mash: mash sketch -s 10000 -k 21, then mash dist; Mash ANI = 1 - distance
fastANI: run as a secondary check, not the main truth
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
comparison      ANIm/dnadiff  minco_best  minco_ctxmash  skani    Mash
full_vs_full    0.9890        0.988175    0.986221       0.9894   0.986859
full_vs_half1   0.9890        0.986711    0.957802       0.9892   0.968856
full_vs_half2   0.9887        0.985193    0.954844       0.9886   0.966481
```

Errors versus dnadiff ANIm:

```text
comparison      minco_best  minco_ctxmash  skani     Mash
full_vs_full    -0.000825   -0.002779      +0.000400 -0.002141
full_vs_half1   -0.002289   -0.031198      +0.000200 -0.020145
full_vs_half2   -0.003507   -0.033856      -0.000100 -0.022219
```

Aligned-fraction behavior:

```text
comparison      dnadiff_ref_af  dnadiff_qry_af  minco_best_af  skani_ref_af  skani_qry_af
full_vs_full    0.9040          0.8808          0.7385         0.8715        0.8360
full_vs_half1   0.4823          0.9039          0.3952         0.8981        0.4308
full_vs_half2   0.4286          0.8577          0.3703         0.8447        0.4052
```

## Validation

- Checks run: generated equal-size halves; ran dnadiff, minco, skani, Mash, and
  fastANI for all three comparisons.
- Expected behavior: full-vs-full dnadiff ANIm should agree with the precomputed
  Nayfach full-pair ANIm label.
- Observed behavior: dnadiff full-pair ANIm was `0.9890`, matching the Nayfach
  label `0.9890599554743553` after rounding.
- Failure modes checked: command exits are checked by the driver; parsers require
  exactly one row per tool output.

## Important Artifacts

See `artifacts.md` for paths to generated files and external inputs.

## Conclusion

On this one high-ANI same-species pair, truncating one genome to half length
barely changes dnadiff regional ANIm and skani ANI, while minco default drops
modestly by about 0.2-0.35 ANI percentage points. minco context-Mash and Mash
drop much more strongly under truncation because the shared-sketch fraction is
affected by missing genome content.

## Paper-Relevant Claim

Tentative only: minco default ANI is more robust than context-Mash-style ANI to
large one-sided genome truncation, but the current default still has some
coverage-related downward bias compared with skani on this pair.

## Caveats

- This is one same-species pair only.
- The half-genome truth is dnadiff regional identity, not a precomputed Nayfach
  label.
- The query halves are artificial concatenated sequence halves, not real MAG
  incompleteness with random contig loss.
- Mash ANI is reported as `1 - Mash distance`, matching prior quick-comparison
  convention.

## Next Experiment

Repeat on 50-100 same-species pairs with random 25%, 50%, and 75% genome
subsampling to separate truncation position effects from general completeness
effects.
