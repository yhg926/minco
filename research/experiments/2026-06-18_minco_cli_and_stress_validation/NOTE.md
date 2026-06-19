# minco CLI and stress validation

Date: 2026-06-18
Author/agent: ubuntu
Project: minco
Code commit: NA
Tool version: minco 0.1

## Question

After converting the repo to pure C and keeping only the public `minco`
subcommands, do the shown commands and currently involved sketch/ANI/matrix
functions work on both synthetic coverage tests and larger real data?

## Hypothesis

The pure C minco binary should pass the public CLI surface, sketch modification
functions, read-QC paths, and large genome/FASTQ stress tests without reverting
to the removed C++ prototype code.

## Dataset

- Synthetic FASTA/FASTQ generated inside `tests/smoke.sh` and `tests/full_cli.sh`.
- First 1000 paths from
  `/mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list`.
- Four 250m ONT FASTQ files from
  `/home/ubuntu/kssd3_prjna851853_salmonella/ONT_yield_curve_alltools/fastq_by_yield/250m`.
- Four matching Salmonella assemblies from
  `/home/ubuntu/kssd3_prjna851853_salmonella/atb_assemblies`.

## Methods

Key parameters:

```text
binary: bin/minco
version: minco 0.1
functional tests: make test
genome stress: 1000 genomes, sketch -p8 --ctxmeta both, index, ANI triangle, matrix triangle
FASTQ stress: 4 x 250m ONT FASTQ.gz, sketch --conflict --readsQC --ctxmeta both -p8
qraw check: reads_sketch vs matching assembly sketch, -m0 -f0 -n0 -t0 -p8
```

Commands are recorded in `commands.sh`.

## Results

Key metrics are recorded in `summary.tsv`.

Important measured results:

```text
make test: pass, 0.33 s, 5.4 MB RSS
1000-genome sketch: pass, 4.24 s, 372 MB RSS, sketch dir 192 MB
1000-genome index: pass, 0.58 s, 198 MB RSS
1000-genome ANI triangle: pass, 34.24 s, 82 MB RSS
1000-genome matrix triangle: pass, 0.67 s, 201 MB RSS
4 x 250m FASTQ readsQC sketch: pass, 26.94 s, 690 MB RSS
FASTQ qraw ANI vs 4 assemblies: pass, 16 comparisons
```

## Validation

- `make test` now runs both smoke and expanded CLI tests.
- Expanded tests cover shown public subcommands: `sketch`, `ani`, `matrix`,
  `examples`, and `doctor`.
- Expanded `sketch` coverage includes list/positional/stdin/pipecmd inputs,
  `--ctxmeta` modes, `--anno`, `--position`, `--splitmfa`, `--asone`,
  `--psmp`, `--psketch`, `--pindex`, `--ppos`, `-i`, `--keep`, `--remove`,
  `--append`, `--dedup`, `--dedup-index`, `--readsQC`, `-A`, and
  `--sketchQC`.
- Expanded `ani` coverage includes sketch-vs-sketch detail, self matrix,
  self triangle, direct FASTA input, sequence query auto-sketch, `--qraw`, and
  `--reflist/--qrylist`.
- Expanded `matrix` coverage includes full, triangle, rectangular, edges,
  clusters, dedup-plan, and PHYLIP full matrix output.

## Debug/Fixes

The expanded test found two real bottom-k gaps:

- `--position` sketches aborted because the minhash bottom-k stage only handled
  key-only records. Fixed by preserving key/position alignment through sorting,
  deduplication, conflict removal, and truncation.
- Abundance sketches used by `--sketchQC` aborted because the same bottom-k
  stage rejected value arrays. Fixed by preserving key/count and
  key/count/position alignment.

## Conclusion

The shown public CLI surface and the currently involved sketch, ANI, matrix,
read-QC, metadata, and sketch-modification paths passed functional tests and
larger real-data stress tests. The implementation remains pure C.

The remaining efficiency observation is that `ani -m2` self-triangle took
34.24 s for 1000 genomes while `matrix --format triangle` took 0.67 s on the
same sketch. This is not a correctness failure, but it is a useful target for a
future ANI-specific optimization pass.

## Caveats

- Stress output was written under `/tmp/minco_validation_20260618`; preserve or
  move it if these exact artifacts are needed later.
- The 1000-genome ANI triangle used default ANI filters, so many cross-species
  off-diagonal values are reported as zero/no-call.
- The FASTQ stress test used four 250m files, not the full 12 GB 250m directory.

## Next Experiment

Optimize the self `ani -m2` object-difference path or add an intentional fast
dispatch for matrix-compatible metrics, then rerun the 1000-genome timing.
