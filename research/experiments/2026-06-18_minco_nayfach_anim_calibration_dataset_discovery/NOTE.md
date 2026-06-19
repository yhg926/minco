# minco Nayfach ANIm calibration dataset discovery

Date: 2026-06-18
Author/agent: Codex
Project: minco
Code commit: NA
Tool version: `bin/minco 0.1`

## Question

Find a local Nayfach ANIm ground-truth dataset that can be reused to recalibrate
minco ANI models, specifically the raw MoE-style distance and the HGB/Best
calibration layer.

## Hypothesis

Existing Nayfach ANIm artifacts from prior KSSD3 model work are present on
`/mnt/new3T` and contain enough pair-level labels and feature columns to seed a
minco-specific calibration experiment.

## Dataset

- Primary ANIm label table:
  `/mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv`
- Alternate T11 label/feature table:
  `/mnt/new3T/skani_data/Nayfach_data/Nayfach52k.kssd3_codenpatternT11_vs_ANIm`
- FASTA directory:
  `/mnt/new3T/skani_data/Nayfach_data/fna`
- Metadata:
  `/mnt/new3T/skani_data/Nayfach_data/genome_metadata.tsv`

The primary table has no header. Observed columns are:

```text
1  query/ref path A, relative to Nayfach_data
2  query/ref path B, relative to Nayfach_data
3  XnY_ctx
4  query alignment fraction
5  reference alignment fraction
6  N_diff_obj
7  N_diff_obj_section
8  N_mut2_ctx
9  KSSD/KSSD3 raw ANI for that feature set
10 ANIm ground-truth ANI
```

## Methods

Commands are recorded in `commands.sh`.

Checks performed:

```text
rg/find for Nayfach, ANIm, and ground-truth artifacts
head/wc/awk schema and range checks on candidate tables
FASTA and metadata existence checks
```

## Results

Key metrics are recorded in `summary.tsv`.

```text
Primary table rows: 1,075,493
Primary table columns: 10
ANIm truth range: 0.840905409902203 to 1.0
Primary T10 raw ANI range: 0.857929 to 1.0
XnY_ctx range: 148 to 41,283
Nayfach FASTA count: 52,515 .fna files
```

Prior KSSD3 model artifacts in the same folder include fitted ridge/correction
models and summaries. These are useful as references, but they are not a
substitute for minco recalibration because minco uses fixed-size bottom-k
context minhash rather than the older fixed reduction-fold sketch distribution.

## Validation

- Checked the primary label table has 10 columns and 1,075,493 rows.
- Confirmed the alternate T11 table also has 1,075,493 rows.
- Confirmed `/mnt/new3T/skani_data/Nayfach_data/fna` exists and contains 52,515
  FASTA files.
- Confirmed `genome_metadata.tsv` exists and includes genome IDs, length,
  completeness, contamination, quality, OTU, and taxonomy/ecosystem fields.

## Important Artifacts

See `artifacts.md`.

## Conclusion

The Nayfach ANIm ground truth needed for recalibration is available locally.
For true minco model training, use this ANIm table as the label source but
recompute pair-level minco features on a sampled or full pair set with the
target minco parameters, especially `--sketch-size`.

## Paper-Relevant Claim

None yet. This note only records dataset discovery and suitability.

## Caveats

- The primary table's feature columns are KSSD3 T10 features, not minco
  features.
- The table can train historical KSSD3 models directly, but minco needs its own
  feature extraction before fitting MoE/HGB calibration.
- Full 1.075M-pair minco feature extraction may be expensive; start with a
  stratified sample and then scale to the full table if the model improves.

## Next Experiment

Create a minco Nayfach feature table:

1. Sample or stratify pair rows from the primary ANIm table.
2. Recompute minco `ani -m0 -f0 -n0 -t0 --raw-output` features for those pairs,
   using fixed `--sketch-size` values such as 10k and 20k.
3. Fit minco-specific MoE/ridge and HGB calibration models against ANIm truth.
4. Validate by held-out genome/OTU split, not only random row split.
