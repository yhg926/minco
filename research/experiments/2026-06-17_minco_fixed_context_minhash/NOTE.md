# Fixed-Size Context MinHash minco V1

Date: 2026-06-17

## Question

Can a fixed-size bottom-k context sketch combine KSSD3 context/object scoring
with Mash-style genome-length-independent sketch size, and how does the first
implementation behave on same-species ANI benchmarks?

## Implementation

Implemented `bin/minco` in this repo.

Core design:

- default coden layout compatible with KSSD3's 32-mer, 44-bit context and
  20-bit object split;
- bottom-k MinHash over context hash values, default `-s 10000`;
- object differences counted only on shared contexts;
- reports `ANI_hybrid`, `ANI_obj`, and `ANI_ctx_mash`;
- low-overlap guard marks object-based ANI as `low_overlap` unless
  `XnY_ctx >= 1000` and max alignment fraction is at least 0.10;
- bounded candidate heap is default; `--full-map` keeps exact all-context
  accounting for debugging only;
- `sketch` and `ani` support `-p/--threads`;
- speed pass added an early bottom-k threshold reject before hash-table lookup,
  precomputed hash-sorted sketch views for bottom-k Jaccard, buffered ANI row
  formatting, and a direct `fread` parser for non-gzip FASTA/plain inputs;
- compact comparison binary `bin/minco_compact` stores each retained item as
  one `uint64_t`: high 44 bits are the context hash key and low 20 bits are the
  object. This trades exact stored context for a small hash-collision risk;
- compact builds use a 256-entry lookup table for per-base decoding;
- compact builds use the same MIT `kseq.h` FASTA/FASTQ reader used by KSSD3A
  for gzip inputs;
- compact builds use a KSSD3A-style vector/radix candidate buffer for bottom-k
  context MinHash instead of the previous `unordered_map` + heap store;
- BMI2 compact comparison binary `bin/minco_compact_bmi2` writes the same
  compact sketch format but also uses KSSD3A-style `_pext_u64` context/object
  extraction on BMI2-capable x86 CPUs;
- rolling-hash comparison binary `bin/minco_compact_roll` uses the same
  compact `uint64_t` layout, but replaces `mix64(ctx)` with an ntHash-style
  rotated-base hash over the 22 selected coden-context bases.

Build/test command:

```bash
make test
```

## Data

Same-species ANIm-grounded benchmark tables from:

- Neisseria meningitidis n=20:
  `/mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_neisseria_n20`
- Five-species n=20 each:
  `/mnt/new3T/gtdbr220/ani_same_species_word/same_species_ani_benchmark_all_n20`

The five species are Brucella melitensis, Escherichia coli, Kingella kingae,
Neisseria meningitidis, and Pseudomonas_E paracarnis.

## Results

Neisseria n=20:

- minco hybrid: MAE 0.496870 ANI percentage points, RMSE 0.508351,
  Pearson r 0.989031.
- Existing baseline summary: KSSD3A best MAE 0.498071; skani MAE 0.374781.
- Optimized minco `-p8` timing: sketch 0.16 s / 41.4 MB RSS, compare
  0.05 s / 22.1 MB RSS.
- Local timing: KSSD3A sketch 0.03 s / 40.9 MB RSS, ANI 0.05 s / 6.6 MB RSS;
  skani triangle 0.12 s / 46.2 MB RSS; Mash triangle 0.19 s / 17.1 MB RSS.

Five-species n=20 each:

- minco hybrid overall: MAE 0.249972, RMSE 0.333985, Pearson r 0.969672.
- Existing baseline summary: KSSD3A best MAE 0.162096; skani MAE 0.177801.
- minco hybrid by species:
  Brucella 0.056981 MAE, E. coli 0.324455, Kingella 0.053644,
  Neisseria 0.496870, Pseudomonas_E 0.317910.
- Optimized minco `-p8` timing for 100 genomes: sketch 1.05 s / 87.5 MB
  RSS, all-vs-all compare 0.40 s / 98.1 MB RSS.
- Compact minco `-p8` timing for 100 genomes: sketch 0.98 s / 58.2 MB RSS,
  all-vs-all compare 0.21 s / 36.4 MB RSS. Same-species accuracy summary was
  identical to the full representation. Sketch file size was 8.02 MB versus
  20.02 MB for the full ctx+hash+obj representation.
- Rolling compact minco `-p8` timing for 100 genomes: sketch 1.22 s /
  59.6 MB RSS, all-vs-all compare 0.21 s / 36.5 MB RSS. Same-species hybrid
  MAE improved to 0.226773, RMSE 0.302028, Pearson r 0.974402.
- Baseline timing on the same 100 genomes: KSSD3A sketch 0.28 s / 109.6 MB
  RSS plus ANI 2.05 s / 29.5 MB RSS; skani triangle 1.48 s / 240.3 MB RSS;
  Mash triangle 1.95 s / 106.9 MB RSS.

Cross-species guard check from the same 100-genome all-vs-all output:

- 950/950 same-species pairs pass the default overlap guard.
- 0/4000 cross-species pairs pass the default overlap guard.

GTDBr226 first 1000 genomes speed check:

- Input list: first 1000 paths from
  `/mnt/new3T/gtdbr220/GTDBr226_kssd3a_Tf8_anno_20260604/GTDBr226_genomes.fna_gz.list`.
- minco full: sketch 10.59 s / 498.7 MB RSS, all-vs-all ANI 17.19 s /
  1095.5 MB RSS, total 27.78 s.
- minco compact mix64: sketch 9.76 s / 195.4 MB RSS, all-vs-all ANI
  12.29 s / 479.2 MB RSS, total 22.05 s.
- minco compact mix64 after lookup-table base decoding: sketch 9.27 s /
  193.9 MB RSS. The `.kmini` file was byte-identical to the regular compact
  mix64 baseline (`cmp` exit 0).
- minco compact BMI2 after lookup-table base decoding: sketch 7.37 s /
  194.0 MB RSS. The `.kmini` file was byte-identical to the regular compact
  mix64 baseline (`cmp` exit 0). A repeated regular compact run after restoring
  heap pruning and switching final selection to `nth_element` was still about
  10.05 s before the base-decoding lookup, so the largest measured sketch speedup
  comes from `_pext_u64` extraction rather than final sorting.
- minco compact vector/radix after porting KSSD3A-style gzip reading and
  candidate buffering: sketch 6.32 s / 352.3 MB RSS. The `.kmini` file was
  byte-identical to the regular compact mix64 baseline (`cmp` exit 0).
- minco compact BMI2 vector/radix: sketch 4.93 s / 349.7 MB RSS. The
  `.kmini` file was byte-identical to the regular compact mix64 baseline
  (`cmp` exit 0). This is still slower than KSSD3A `-f8` sketching on the same
  list (2.86 s / 237.9 MB RSS), but the gap narrowed from 3.4x to 1.7x for
  sketching.
- minco compact mix64 indexed detail: sketch 9.76 s / 195.4 MB RSS,
  all-vs-all detail ANI 1.81 s / 265.4 MB RSS, total 11.57 s. Object-based
  columns and `ANI_hybrid` matched the pairwise detail output exactly on
  1,000,000 rows; `ANI_ctx_mash` uses set-Jaccard in indexed detail.
- minco compact mix64 indexed triangle: sketch 9.76 s / 195.4 MB RSS,
  triangle ANI 1.54 s / 247.8 MB RSS, total 11.30 s.
- minco compact rolling hash: sketch 11.68 s / 194.8 MB RSS, all-vs-all
  ANI 12.37 s / 478.5 MB RSS, total 24.05 s.
- Baseline timings on the same list: KSSD3A `-f8` with corrected self-triangle
  command total 4.08 s when forcing the indexed matrix route
  (`KSSD3A_ANI_MATRIX_DIRECT_THRESHOLD=0 kssd3a ani -q sketch -m2 -s -1 -d`);
  skani triangle 16.49 s; Mash `-s10000` triangle 17.19 s.
- KSSD3A caveat: default `kssd3a ani -q sketch -m2 -s -1 -d` took 58.41 s
  total on exactly 1000 samples because the documented default direct/index
  threshold is 1000. Forcing the indexed route produced byte-identical
  triangle output in 4.08 s total.
- Output formats differ: minco writes a 1,000,000-pair full matrix plus
  header; Mash writes 499,500 triangle pairs; skani writes a triangle matrix;
  corrected KSSD3A `-m2` writes 1000 lower-triangle matrix rows.

## Interpretation

The algorithm works as intended as a fixed-size context sketch and gives usable
same-species ranking signal. The first raw hybrid metric is competitive with
KSSD3A on the Neisseria n=20 benchmark, but across five species it is less
accurate than KSSD3A best and skani.

The dominant issue is calibration, not sketch mechanics. Species-specific bias
is visible: Neisseria is biased high, while E. coli and Pseudomonas_E are biased
low. The bias-removed MAE in the full benchmark is 0.238552, so a learned or
species-agnostic calibration model should be tested next.

The default bounded bottom-k implementation is much more memory efficient than
the initial exact full-map prototype. On the 100-genome benchmark, exact full
map took 3:41.77 and 8.9 GB RSS; bounded bottom-k took 14.33 s and 54.1 MB RSS
for sketching with identical same-species summary metrics in this run.

The speed pass preserved byte-identical ANI output relative to the previous
threaded implementation on the 100-genome benchmark. Wall time improved from
2.01 s sketch + 1.43 s compare to 1.05 s sketch + 0.40 s compare. Compared with
the original single-thread bounded implementation, total time dropped from
24.35 s to 1.45 s on the 100-genome benchmark.

The compact representation changed 24/10000 all-vs-all rows relative to the
full representation, all cross-species low-overlap rows. It changed 0
same-species benchmark rows and did not change any `Pass_default` calls in this
dataset.

The ntHash-style rolling compact variant did not improve sketch speed on this
implementation: it was 0.24 s slower than compact `mix64` on the 100-genome
benchmark. It did improve same-species hybrid MAE on this dataset, likely
because the different hash permutation selected a different bottom-k context
sample. Treat this as an accuracy/sampling observation, not proof that rolling
hash is generally better.

The 1000-genome sketching bottleneck was not final sorting. Replacing full final
hash sort with `nth_element` preserved bytes but did not explain the speed gap.
Lookup-table base decoding reduced compact sketch wall time from 9.76 s to
9.27 s. KSSD3A-style `kseq` gzip reading plus vector/radix bottom-k buffering
reduced portable compact sketching to 6.32 s, and adding BMI2 `_pext_u64`
context/object extraction reduced the same compact sketch to 4.93 s on the same
1000-genome list. The speed gain trades memory: compact BMI2 vector/radix used
349.7 MB RSS versus 195.4 MB for the old compact sketcher.

## Caveats

- ANIm truth was reused from existing benchmark artifacts; it was not
  regenerated in this run.
- Timing was measured on one local machine and should be repeated for larger
  genome counts before making broad throughput claims.
- The 1000-genome speed check used a broad GTDB path prefix, not an ANIm truth
  dataset; it is speed-only.
- The optimized ANI path stores an additional hash-sorted copy of each sketch,
  increasing compare RSS versus the earlier threaded version.
- The rolling-hash result is one benchmark run on 100 genomes. It needs
  replicate seeds and more species before interpreting the accuracy difference.
- Low-overlap cross-species pairs are no-calls for object-based ANI; raw
  object ANI should not be interpreted without the guard columns.

## Next Experiment

Train a fold/sketch-size-aware calibration model on multiple species using
features already emitted by minco: `ANI_hybrid`, `ANI_obj`, `ANI_ctx_mash`,
`Jaccard_ctx`, `XnY_ctx`, alignment fractions, and object-difference counts.
Evaluate on held-out species and include KSSD3A/skani/Mash timings on the same
full benchmark. For speed, profile sketching on larger datasets and test a
packed C-style candidate table/radix selection path before considering a full C
rewrite.
