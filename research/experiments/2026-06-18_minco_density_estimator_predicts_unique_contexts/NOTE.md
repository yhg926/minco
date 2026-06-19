# minco density estimator predicts unique contexts

Date: 2026-06-18
Author/agent: ubuntu
Project: minco
Code commit: NA
minco binary/tool version: minco_stage3_native + minco_exact20m

## Question

Can the minco density estimator from the final bottom-k threshold predict
total unique genome contexts, and how well does it track assembly genome size?

## Hypothesis

Because context hashes are approximately uniform, `observed_ctx_under_threshold
* 2^44 / (threshold + 1)` should estimate total unique context hashes with
sampling error near the expected bottom-k scale. Genome size should also be
strongly correlated, but less exact because invalid bases, contig boundaries,
and repeated contexts reduce unique contexts per base.

## Dataset

- Input data: first 50 genomes from `/tmp/kssd3a_based_minco_stages/genomes_50.list`.
- Input sketch or matrix: none; sketches generated for this experiment.
- Sample count: 50.
- Selection criteria: existing GTDB stage-test genome list, first 50 entries.
- Storage location: source genomes under `/mnt/new3T/gtdbr220/GTDBr226_genomes`.

## Methods

Key parameters:

```text
Normal estimator:
  bin/minco_stage3_native sketch -p 8 --ctxmeta both
  minco fixed sketch size = 10000 context hashes

Exact total unique context proxy:
  bin/minco_exact20m sketch -p 2 --conflict --ctxmeta both
  MINCO_SKETCH_SIZE=20000000
  --conflict keeps all context-object conflicts.
  Since max exact sketch_entries was 8,363,072, no sample hit the 20M cap.

Genome size:
  Parsed total_length_bp from lcofiles.infilemeta.

Exact total unique contexts:
  postconflict_observed_ctx from the exact --conflict run. With --conflict and
  no cap hit, this is the exact number of distinct retained 44-bit context
  hashes before conflict removal.
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

```text
n = 50
Exact total unique context range: 242,746 to 8,344,574; mean 3,419,856.92.
Genome size range: 243,933 bp to 8,491,258 bp; mean 3,509,350.34 bp.
Mean exact unique contexts per bp: 0.9708; range 0.8114 to 0.9988.

Normal 10k preconflict estimator vs exact total unique contexts:
  Pearson = 0.999846
  Spearman = 0.999328
  MAPE = 0.813%
  mean relative bias = +0.185%
  max absolute relative error = 2.642%

Genome bp vs exact total unique contexts:
  Pearson = 0.999066
  MAPE = 3.172%
```

## Validation

- Checks run: exact-count debug sketch used a 20M cap and `--conflict`.
- Expected behavior: exact sketch entry counts should remain below the 20M cap.
- Observed behavior: max exact sketch_entries was 8,363,072, so no sample hit the cap.
- Failure modes checked: estimator sidecar was generated for both standard and exact runs; `lcofiles.infilemeta` size was consistent with 50 records at 32 bytes each.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, matrices, trees, figures, and sketches.

## Conclusion

The density estimator is a good predictor of total unique context count on this
50-genome set: about 0.81% MAPE and 0.99985 Pearson correlation against exact
unique context counts. It also tracks genome size very strongly, but genome size
is a less exact target because unique contexts per bp varied from 0.811 to
0.999 in this sample.

## Paper-Relevant Claim

Exploratory: minco can estimate per-genome total unique context cardinality
from the 10,000-context bottom-k threshold with sub-1% mean absolute percentage
error on this 50-genome GTDB test set.

## Caveats

- Tested on only 50 GTDB genomes from an existing stage benchmark list.
- Exact total context count is for distinct 44-bit context hashes, not raw un-hashed context strings.
- Needs replication across more taxa, small/large genomes, fragmented assemblies, and high-N assemblies.
- The exact-count debug sketch is large and intended only for validation, not routine use.

## Next Experiment

Repeat on 1,000 genomes with exact counting implemented as a counter-only mode
instead of writing all exact sketch entries, then stratify error by genome size,
GC content, contig count, N fraction, and taxonomic group.
