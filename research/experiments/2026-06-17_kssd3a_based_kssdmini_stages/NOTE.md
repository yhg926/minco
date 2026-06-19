# KSSD3A-Based KSSDmini Staged Port Plan

Date: 2026-06-17

## Goal

Build a KSSDmini sketcher from the KSSD3A sketch path, then replace only
KSSD3A's fixed-rate context selection with fixed-size context MinHash selection.
Measure speed, memory, sketch size, and output changes at each stage.

## Rationale

The current C++ KSSDmini compact BMI2 vector/radix backend improved 1000-genome
sketch time from 9.76 s to 4.93 s, but KSSD3A `sketch -f8` remains faster at
2.86 s on the same input. The next implementation should start from KSSD3A's
known-fast sketch architecture and make the context-selection change in small
controlled steps.

## Fixed Baseline

Dataset:

- `/tmp/kssdmini_1000_speed/genomes.list`
- 1000 GTDBr226 genome paths from
  `/mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list`

Baseline commands:

- KSSD3A: `/home/ubuntu/yihuiguang/tools/KSSD3codex/bin/kssd3a sketch -f8 -p8`
- Current KSSDmini: `bin/kssdmini_compact_bmi2 sketch -p 8 -s 10000`

Known results:

| stage | sketch wall | max RSS | sketch size | note |
|---|---:|---:|---:|---|
| KSSD3A `-f8` | 2.86 s | 237.9 MB | 108.50 MB raw sketch; 270.80 MB after m2 index | fixed-rate context selection |
| KSSDmini compact old | 9.76 s | 195.4 MB | 80.16 MB | C++ map/heap |
| KSSDmini compact BMI2 vector/radix | 4.93 s | 349.7 MB | 80.16 MB | fixed-size MinHash, byte-identical to old compact |
| copied KSSD3A stage0 `-f8` | 2.87 s | 243.7 MB | 108.50 MB raw sketch | copied KSSD3A code path |
| stage2 exact post-collection MinHash | 59.72 s | 27.84 GB | 80.30 MB | exact but too slow/high-memory |
| stage3 native streaming MinHash `m2` | 4.10 s | 369.1 MB | 80.30 MB | exact and selected speed candidate |

## Stages

### Stage 0: Exact KSSD3A Sketch Skeleton

Copy the KSSD3A sketch hot path into this repo as a new isolated target, without
changing context selection.

Implementation target:

- copied KSSD3A source subtree under `kssd3a_stage0/`
- binary: `bin/kssd3mini_stage0`

Selection rule:

- keep KSSD3A fixed-rate rule:
  `SKETCH_HASH(unictx) <= FILTER`, with `FILTER = UINT32_MAX >> drfold`

Expected invariant:

- Same number of retained ctxobj entries as KSSD3A `-f8` for each sample, or a
  documented difference caused only by output-wrapper metadata.
- Speed should be close to KSSD3A. If it is not, the copied skeleton is not close
  enough and should be fixed before adding MinHash.

Measured result:

- Command:
  `bin/kssd3mini_stage0 sketch -f8 -p8 -l /tmp/kssdmini_1000_speed/genomes.list -o /tmp/kssd3a_based_kssdmini_stages/stage0_kssd3a_f8`
- Wall time: 2.87 s
- User time: 21.68 s
- Max RSS: 243.7 MB
- Raw output size: 108,499,888 bytes
- File check versus native KSSD3A baseline:
  `comblco`, `comblco.index`, and `lcofiles.infilemeta` are byte-identical;
  `lcofiles.stat` differs.
- Gate status: passed. The copied skeleton is close enough to native KSSD3A to
  proceed with isolated stage changes.

### Stage 1: KSSD3A Filter Plus Hash Bottom-k

Keep Stage 0 `-f8` filtering unchanged, but transform retained exact
`ctxobj` keys to KSSDmini-style `(ctxhash44 << 20) | obj` keys and cap each
sample at 10,000 entries. This keeps KSSD3A directory output so output-format
cost is not mixed into the selection experiment.

Selection rule:

- unchanged KSSD3A fixed-rate prefilter, then bottom-k cap

Measured result:

- Wall time: 2.90 s
- Max RSS: 243.5 MB
- Mean entries/genome: 8,801.02
- Interpretation: hashing/capping the already-filtered set is cheap, but this
  is not a valid fixed-size sketch because KSSD3A `-f8` often leaves fewer than
  10,000 contexts.

### Stage 2: Exact Fixed-size Bottom-k After Full Collection

Disable the KSSD3A fixed-rate prefilter, collect all per-genome contexts, then
transform to `(ctxhash44 << 20) | obj`, deduplicate, remove conflicting
contexts, and keep the lowest 10,000 context hashes.

Selection rule:

- exact fixed-size MinHash after full collection

Measured result:

- Wall time: 59.72 s
- Max RSS: 27.84 GB
- Mean entries/genome: 10,000
- Byte-identical to current KSSDmini compact output.
- Interpretation: the major slowdown is collecting/sorting all contexts before
  selection.

### Stage 3: Streaming/Threshold Bottom-k Optimization

Move bottom-k rejection into the KSSD3A hot loop. The selected implementation
hashes the compressed context first, rejects contexts above the current
threshold before extracting the object, keeps an oversampled buffer
(`keep=2*s`, `trigger=4*s`), and defers conflict removal until the final pass so
conflict evidence is not forgotten.

Selection rule:

- exact fixed-size MinHash with streaming early rejection

Measured result:

- `bin/kssd3mini_stage3_native`
- Wall time: 4.10 s
- User time: 31.30 s
- Max RSS: 369.1 MB
- Mean entries/genome: 10,000
- Byte-identical to Stage 2 exact post-collection output and to current
  KSSDmini compact output.
- Tuning: `keep=1*s, trigger=2*s` was faster (3.64 s) but not exact; `keep=8*s,
  trigger=16*s` was exact but slower/larger.

### Stage 5: Integrated ANI Check

Use the stage selected as production candidate to sketch same-species benchmark
sets and compare ANI accuracy with the current KSSDmini compact output, KSSD3A,
skani, and Mash.

Expected invariant:

- If the MinHash selection is semantically the same as current KSSDmini compact,
  same-species ANI should match current KSSDmini compact.

## Results Summary

The speed regression comes from removing KSSD3A's fixed-rate prefilter and
collecting all contexts before selection. That exact post-collection control
took 59.72 s and 27.84 GB RSS.

Streaming bottom-k fixes the regression. The best valid staged build is
`bin/kssd3mini_stage3_native`, which uses the copied KSSD3A parser/data path,
context-first hashing, early threshold rejection, and deferred final conflict
removal. On the 1000-genome benchmark it took 4.10 s wall and 369.1 MB RSS,
with exactly 10,000 entries per genome. Its `comblco` payload is byte-identical
to the exact Stage 2 output and to the current KSSDmini compact sketch.

Compared with the previous best C++ KSSDmini compact BMI2/vector/radix sketcher
at 4.93 s, the selected staged C/native build is about 1.20x faster for sketching
on this dataset. It is still slower than KSSD3A `-f8` at 2.86 s, which is
expected because fixed-size MinHash must inspect more contexts than fixed-rate
KSSD3A sampling.

## Metrics

For each stage, record:

- exact command
- git/source snapshot or copied-source manifest
- wall time, user time, and max RSS from `/usr/bin/time -v`
- output sketch size
- per-sample retained entry count summary: min, mean, max
- byte identity when expected
- semantic identity when byte identity is not expected
- indexed ANI time if `.kmini` output is produced
- accuracy only after Stage 3 or later

## Decision Rules

- Do not proceed past Stage 0 until the copied KSSD3A skeleton is within about
  10-20% of KSSD3A sketch speed on the 1000-genome benchmark.
- Do not mix output-format changes with selection changes.
- Do not tune calibration during sketch-speed stages.
- Treat memory increases as real cost; report speed/RSS together.

## Caveats

- KSSD3A `-f8` and KSSDmini fixed-size `-s10000` do not retain the same number
  of contexts, so speed is not a perfect apples-to-apples comparison until Stage
  3 isolates selection cost.
- KSSDmini compact stores context hashes rather than exact contexts, so rare hash
  collisions remain possible.
