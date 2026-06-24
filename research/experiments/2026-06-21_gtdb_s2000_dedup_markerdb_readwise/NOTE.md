# GTDB S2000 Dedup Context-Markerdb Readwise Test

Date: 2026-06-21
Author/agent: Codex
Project: MinCO / KSSD3mini
Code commit: unavailable; `/home/ubuntu/yihuiguang/tools/KSSD3mini` is not a valid git repository in this session
KSSD3A binary/tool version: `minco 0.1`

## Question

Does building an S2000 GTDB MinCO sketch from the existing S10000 sketch, removing near-duplicate references with an AAF `0.03` dedup plan, and then building a context markerdb improve readwise ANI behavior on CAMI II Toy Mouse Gut sample0?

## Hypothesis

Keeping only reference-specific contexts should remove many fake contexts from readwise assignment. A dedup pass before markerdb construction should reduce redundant references and make the context markerdb less sparse.

## Dataset

- Source sketch: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S10000_anno_20260619/sketch_T_S10000_anno`
- Source sample count: 200,709
- Source sketch size: S10000, 2,007,090,000 total entries
- Downsampled sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_anno`
- Deduped sketch: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup`
- Context markerdb: `/tmp/gtdb232_s2000_dedup_marker.qKJofv/sketch_T_S2000_aaf003_dedup_ctxmarker`
- Readwise sample: `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz`
- Toy Mouse positive source genomes and abundances: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/positive_source_genomes.tsv`
- GTDB r232 metadata: `/mnt/new3T/gtdbr220/gtdb232/metadata/bac120_metadata_r232.tsv.gz` and `/mnt/new3T/gtdbr220/gtdb232/metadata/ar53_metadata_r232.tsv.gz`
- GTDB species ground truth generated here: `mouse0_gtdb_species_profile.tsv`
- Source-aware taxmap: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/mouse0_sourceaware_taxmap_ani95_minaf50.tsv`

## Methods

Key parameters:

```text
downsample S=2000
dedup-plan metric=aaf cut=0.03 threads=8
markerdb mode=--uniq_union --markerdb-ctx
markerdb warning threshold=500 default
readwise assignment=best-diff-unique and best-diff-split
readwise ANI=zip-aaf
direct scoring rule=XnY_ctx >= 10, ANI >= 0.94, Real_min_align_fraction >= 0.05
benchmark truth=GTDB r232 species profile built from CAMISIM source genomes
always-naive variant=same direct rule but ANI replaced by ANI_naive_calc from XnY/N_diff_obj/N_diff_obj_section; tested ANI >= 0.94 and strict ANI > 0.95
```

Commands are recorded in `commands.sh`. GTDB ground truth and authoritative scoring were regenerated with `build_and_score_gtdb_ground_truth.py`; compact outputs are `mouse0_gtdb_species_profile.tsv`, `gtdb_direct_scores.tsv`, `gtdb_direct_scores_naive_ani.tsv`, `gtdb_direct_scores_naive_ani_gt095.tsv`, `gtdb_direct_scores_naive_ani_gt095_breadth_depth.tsv`, and `mouse0_gtdb_truth_summary.tsv`. Source-aware scoring was also regenerated with `score_sourceaware_direct.py`, but it is diagnostic only.

## GTDB Ground Truth

Benchmarking should use GTDB species truth before comparing methods. I transferred each CAMISIM positive source genome from its NCBI-labeled assembly to GTDB r232 species using the GTDB metadata accession fields:

```text
source genomes: 75
mapped source genomes: 75
exact accession matches: 63
unique GCF/GCA assembly-core matches: 12
unmapped source genomes: 0
aggregated GTDB species: 65
raw CAMISIM abundance sum: 5335
normalized GTDB abundance sum: 1.0
```

The per-source mapping is `mouse0_gtdb_source_species_ground_truth.tsv`; the benchmark truth profile is `mouse0_gtdb_species_profile.tsv`. The GCF/GCA assembly-core fallback only accepts unique assembly number/version matches and records `mapping_method=unique_assembly_core` per source genome.

For all method outputs scored in `gtdb_direct_scores.tsv`, selected direct-call rows had `selected_unmapped_rows=0`; any unmapped rows were below the direct-call thresholds and did not enter the benchmark call sets.

## Results

Key metrics are recorded in `summary.tsv`.

```text
S10000 -> S2000 downsample completed: 200,709 refs, 401,418,000 entries.
AAF 0.03 dedup removed only 182 / 200,709 refs, leaving 200,527 refs.
Dedup-stage markerdb prediction exactly matched the actual markerdb:
  total_entries=239,132,175 min=11 max=1934 mean=1192.52 zero=0 below500=9936
Context markerdb construction completed with peak RSS 9,755,392 KB.
Toy Mouse GTDB-species direct F1, authoritative benchmark:
  S1000 unique: 0.862069
  S1000 split: 0.883333
  S2000 dedup ctx-marker unique: 0.878049
  S2000 dedup ctx-marker split: 0.878049
Toy Mouse source-aware direct F1, diagnostic only:
  S1000 unique: 0.869565
  S1000 split: 0.890756
  S2000 dedup ctx-marker unique: 0.885246
  S2000 dedup ctx-marker split: 0.885246
```

I also tested an always-naive direct-scoring variant: use the same GTDB species truth and the same direct thresholds, but replace each row's output `ANI` with `ANI_naive_calc` from `XnY_ctx`, `N_diff_obj`, and `N_diff_obj_section`. This is a negative result:

```text
S1000 unique always-naive:              F1=0.813008 TP=50 FP=8  FN=15
S1000 split always-naive:               F1=0.577114 TP=58 FP=78 FN=7
S2000 marker unique always-naive:       F1=0.846715 TP=58 FP=14 FN=7
S2000 marker split always-naive:        F1=0.846715 TP=58 FP=14 FN=7
S2000 marker split poisson-depth naive: F1=0.834532 TP=58 FP=16 FN=7
S2000 marker split poisson-product:     F1=0.822695 TP=58 FP=18 FN=7
```

The naive formula improves recall for some low-abundance taxa, but direct use at `ANI >= 0.94` lets many related-reference false positives through. The S1000 split case is the clearest failure: predicted GTDB species rose from 55 to 136 and FP rose from 2 to 78.

I then tightened the always-naive filter to strict `ANI_naive_calc > 0.95`. This improved precision slightly but did not restore the original ZIP-AAF direct F1:

```text
S1000 unique naive >0.95:              F1=0.819672 TP=50 FP=7  FN=15
S1000 split naive >0.95:               F1=0.594872 TP=58 FP=72 FN=7
S2000 marker unique naive >0.95:       F1=0.850746 TP=57 FP=12 FN=8
S2000 marker split naive >0.95:        F1=0.850746 TP=57 FP=12 FN=8
S2000 marker split poisson-depth >0.95 F1=0.846715 TP=58 FP=14 FN=7
S2000 marker split product >0.95:      F1=0.828571 TP=58 FP=17 FN=7
```

Compared with the authoritative ZIP-AAF S2000 marker split result, strict always-naive still has lower F1: `0.850746` versus `0.878049`.

I then added the breadth-depth consistency rule suggested by the FP inspection:

```text
base call: XnY_ctx >= 10, ANI_naive_calc > 0.95, Real_min_align_fraction >= 0.05
reject if: Ref_hit_mean_depth >= 4 and Ref_breadth < 0.5
```

This improves the S2000 markerdb naive path:

```text
S2000 marker split naive >0.95:                  F1=0.850746 TP=57 FP=12 FN=8
S2000 marker split naive >0.95 + breadth-depth:  F1=0.883721 TP=57 FP=7  FN=8
S2000 marker split poisson-depth + breadth-depth F1=0.892308 TP=58 FP=7  FN=7
```

For `s2000_marker_split_direct`, the breadth-depth rule rejected 5 candidate rows, all false positives in this GTDB truth. A total-depth Poisson/binomial concentration filter was also tested as a low-depth rule, but it was too aggressive in this sample and cut recall sharply; it is not included in the retained scoring table.

I also tested the stricter adjusted-breadth rule:

```text
base call: XnY_ctx >= 10, ANI_naive_calc > 0.95, Real_min_align_fraction >= 0.05
require:   Ref_zip_af >= 0.5
```

This removes every false positive, but it is too strict for Toy Mouse low-abundance truth:

```text
S2000 marker split poisson-depth + breadth-depth: F1=0.892308 TP=58 FP=7 FN=7
S2000 marker split poisson-depth + Ref_zip_af>=0.5: F1=0.869565 TP=50 FP=0 FN=15
```

The adjusted-breadth cutoff removes 22 of 72 candidate rows in the Poisson-depth output. It eliminates all remaining FP, but also loses eight additional TP, including `Lactobacillus crispatus`.

### Candidate ANI-Coupled Adjusted-Breadth Gate

The fixed adjusted-breadth floor `Ref_zip_af >= 0.5` is too conservative, so the next gate template should tie the breadth floor to the same ANI threshold used by the direct-call ANI filter.

Parameters:

```text
ANI_threshold = 0.95 by default
effective_ctx_length = ctx_length + 2
adjusted_breadth_floor = exp((ANI_threshold - 1) * effective_ctx_length)
```

Template:

```text
1. readwise assign contexts to references
2. estimate per-reference hit-context depth
3. defake contexts before computing naive ANI
4. compute ANI_naive_calc from non-fake contexts only
5. compute adjusted_breadth from the low-tail depth model
6. call if:
     XnY_ctx_nonfake >= 10
     ANI_naive_calc_nonfake > ANI_threshold
     adjusted_breadth >= adjusted_breadth_floor
```

This deliberately uses the same `ANI_threshold` in both the naive-ANI gate and the breadth-derived AF-to-ANI gate, instead of repeating a hard-coded numeric threshold in two places.

The defake model should target high-depth nonzero-diff contexts before `ANI_naive_calc` is computed:

```text
for each reference:
  product_i = best_diff_i * depth_i
  use only product_i > 0 for the product distribution

  preferred fit:
    robust negative-binomial high-tail model
    trim the largest 1-5% product_i values, or at least the single largest value
    fit mean and overdispersion from the trimmed product_i values

  fallback:
    if there are too few nonzero products, do not defake
    if variance <= mean, use a Poisson high-tail model with trimmed mean

  fake if:
    best_diff_i > 0
    high_tail_p(product_i) * ref_marker_ctx_count <= alpha
```

Default:

```text
alpha = 0.05
fake_threshold = -log10(alpha) = 1.30103
minimum nonzero product contexts for NB fit = 10
```

Zero-diff contexts are never removed by the defake model. The model is meant to remove suspicious high-depth, nonzero-diff contexts before the naive ANI formula receives `XnY_ctx`, `N_diff_obj`, and `N_diff_obj_section`.

Naive ANI should be zero-distance based after defaking:

```text
raw_dist = N_diff_obj_nonfake / (XnY_ctx_nonfake + N_diff_obj_nonfake)
section_ratio = (N_diff_obj_section_nonfake + epsilon) / (N_diff_obj_nonfake + epsilon)
final_dist = 1 - (1 - raw_dist) ^ section_ratio
naive_dist = final_dist * 0.1544286
ANI_naive_calc_nonfake = 1 - naive_dist
```

The breadth adjustment model should estimate the breadth expected if read depth were sufficient:

```text
adjusted_breadth = Ref_breadth / (1 - exp(-lambda))
```

where `lambda` is the fitted normal per-context read depth. For the current low-tail version:

```text
1. collect hit reference-context depths
2. use the low-tail depths from 1..75th percentile to avoid high-depth contaminated contexts
3. extrapolate a full Poisson lambda from that low-tail distribution
4. compute adjusted_breadth = Ref_breadth / (1 - exp(-lambda))
5. cap adjusted_breadth at 1.0
```

On the current Toy Mouse S2000 markerdb Poisson-depth table, with `ctx_length=22`, this gives:

```text
effective_ctx_length = 24
adjusted_breadth_floor = exp((0.95 - 1) * 24) = 0.301194
```

Scoring with no `Real_min_align_fraction >= 0.05` gate:

```text
p1..p75 extrapolated adjusted breadth, floor 0.301194: F1=0.920635 TP=58 FP=3 FN=7
builtin Ref_zip_af, same floor:                         F1=0.912000 TP=57 FP=3 FN=8
```

The nearby sweep was:

```text
effective k=22 floor 0.332871: F1=0.912000 TP=57 FP=3 FN=8
effective k=24 floor 0.301194: F1=0.920635 TP=58 FP=3 FN=7
effective k=25 floor 0.286505: F1=0.920635 TP=58 FP=3 FN=7
effective k=27 floor 0.259240: F1=0.906250 TP=58 FP=5 FN=7
effective k=32 floor 0.201897: F1=0.907692 TP=59 FP=6 FN=6
```

Interpretation: tying the adjusted-breadth floor to `ctx_length + 2` and the same `ANI_threshold` is the cleanest current template. The p1..p75 adjusted breadth is modestly better than builtin `Ref_zip_af` at the same floor, but the calibration is still exploratory and single-sample.

For the known low-abundance `Lactobacillus crispatus` target `GCF_018987235.1`, markerdb improved zip-aaf ANI but did not lift the row above the direct call threshold. The old S1000 split readwise naive ANI was already higher than zip-aaf:

```text
S1000 unique zip-aaf=0.878808 naive=0.931216 XnY=68 weak
S1000 split  zip-aaf=0.897978 naive=0.934744 XnY=170 weak
S2000 marker unique/split zip-aaf=0.901884 naive=0.954948 XnY=100 weak
```

## L. crispatus Read-Origin Trace

I reran the S2000 context-markerdb split readwise pass with `--density-block-ctx 0` and the hidden trace hook:

```text
MINCO_READWISE_TRACE_REF=GCF_018987235
MINCO_READWISE_TRACE_OUT=l_crispatus_s2000_marker_split_readwise_assignment_trace.tsv
```

The trace produced 303 selected assignment events from 301 unique reads for the target reference. Joining those read IDs to CAMISIM `reads_mapping.tsv.gz`, then mapping `genome_id` through `positive_source_genomes.tsv`, shows the remaining signal is still dominated by related-species contamination:

```text
Lactobacillus helveticus  203 events, 203 unique reads, 67.0% events
Lactobacillus crispatus    99 events,  97 unique reads, 32.7% events
Lactobacillus johnsonii     1 event,    1 unique read,  0.3% events
```

The top contaminating source genome is `GCF_001702095.1` / genome ID `463794.1`, labeled `Lactobacillus helveticus`. The true `L. crispatus` source genome is `GCF_002218965.1` / genome ID `denovo11386.0`.

This is much cleaner than the old S1000 trace, where `L. helveticus` contributed about 95% of the selected events, but the new markerdb trace still has a majority non-`L. crispatus` contribution. That explains why the S2000 markerdb readwise naive ANI improves to `0.954948` but remains below the assembly-to-ref naive ANI near `0.986`.

## Poisson Depth Chomp Test

I added an explicit experimental readwise context filter:

```text
--readwise-ctx-filter poisson-depth
--readwise-fake-threshold 1.30103
```

For this mode, the threshold is `-log10(Poisson high-tail p * markerdb_size)`, so `1.30103` corresponds to Bonferroni-adjusted p-value `<= 0.05`. The filter only removes nonzero-diff contexts whose depth is a high-tail Poisson outlier.

On the `GCF_018987235.1` trace, the offline version removed exactly one context: the `L. helveticus` context with 187 reads, `best_diff=4`, and no true `L. crispatus` reads. Recomputed naive features changed:

```text
before: XnY=100 Ndiff=21 Nsection=38 Nmut2=7 ANI=0.954948
after:  XnY=99  Ndiff=20 Nsection=34 Nmut2=6 ANI=0.958519
```

Whole-sample Toy Mouse runs using the new filter completed successfully. The filter rejected 8,280 nonzero-diff contexts across 6,976 reported rows. It reduced output rows from 65,414 to 54,626 because some references lost all effective naive contexts.

GTDB-species direct scoring did not improve:

```text
s2000 marker split zip-aaf unfiltered:        F1=0.878049 TP=54 FP=4  FN=11
s2000 marker split zip-aaf poisson-depth:     F1=0.878049 TP=54 FP=4  FN=11
s2000 marker split naive unfiltered:          F1=0.846715 TP=58 FP=14 FN=7
s2000 marker split naive poisson-depth:       F1=0.834532 TP=58 FP=16 FN=7
```

Interpretation: the Poisson-depth idea is real for the diagnostic target and cleanly catches the high-depth contaminated `L. helveticus` context. Alone, at this threshold, it is not enough to improve species-level Toy Mouse GTDB F1. It should remain a diagnostic feature or be combined with a broader call model.

I also tested a stronger product outlier variant:

```text
--readwise-ctx-filter poisson-product
--readwise-fake-threshold 1.30103
```

This fits a robust per-reference Poisson lambda to `best_diff * depth` after dropping the largest product outlier, then applies the same Bonferroni-adjusted high-tail test. On the `L. crispatus` trace, this is much closer to the oracle contamination removal:

```text
observed:        XnY=100 Ndiff=21 Nsection=38 Nmut2=7 ANI=0.954948
poisson-depth:   XnY=99  Ndiff=20 Nsection=34 Nmut2=6 ANI=0.958519
poisson-product: XnY=96  Ndiff=17 Nsection=22 Nmut2=3 ANI=0.970625
oracle clean:    XnY=91  Ndiff=15 Nsection=17 Nmut2=2 ANI=0.975477
```

However, whole-sample GTDB direct scoring got worse with selected naive ANI:

```text
s2000 marker split naive poisson-product: F1=0.822695 TP=58 FP=18 FN=7
```

The product rule is therefore a strong target-level diagnostic, but too aggressive as a standalone direct-call filter in this sample.

I then added an experimental negative-binomial product filter:

```text
--readwise-ctx-filter product-nb
--readwise-fake-threshold 1.30103
```

This uses positive `best_diff * depth` values, drops the largest product outlier per reference before fitting, requires at least 10 positive product contexts, fits a negative-binomial high-tail model when the trimmed variance exceeds the trimmed mean, and falls back to Poisson otherwise. Zero-diff contexts are never rejected.

The full Toy Mouse run completed successfully:

```text
output rows excluding header: 65414
wall time: 1:58.45
peak RSS: 4,289,920 KB
rejected contexts: 33
rows with rejected contexts: 28
```

Scoring the current candidate gate template on this output:

```text
Gate:
  XnY_ctx_nonfake >= 10
  ANI_naive_calc_nonfake > 0.95
  adjusted_breadth >= exp((0.95 - 1) * (22 + 2)) = 0.301194
  no Real_min_align_fraction gate

product-nb + p1..p75 adjusted breadth: F1=0.920635 TP=58 FP=3 FN=7
product-nb + builtin Ref_zip_af:       F1=0.912000 TP=57 FP=3 FN=8
```

Interpretation: product-NB is much more conservative than Poisson-product and does not damage the full-sample call set, but at the current threshold it is nearly neutral. The useful part of the current template is still the ANI-coupled adjusted-breadth floor rather than NB defaking.

## Defaked Truncated Depth Trial

I added experimental readwise abundance columns that are recomputed after the same product-NB fake-context decisions used for the naive ANI features:

```text
Reliable_Ref_breadth
Reliable_Ref_mean_depth
Reliable_Ref_hit_ctx
Reliable_Ref_hit_mean_depth
Reliable_Ref_hit_depth_variance
Reliable_Ref_zip_af
```

This lets the adjusted-breadth model use only non-fake hit contexts. The rerun completed successfully:

```text
output rows excluding header: 65414
wall time: 1:55.75
peak RSS: 4,295,408 KB
```

I then scored these direct-call gates:

```text
shared base:
  XnY_ctx_nonfake >= 10
  ANI_naive_calc_nonfake > 0.95
  no Real_min_align_fraction gate
  adjusted-breadth floor = exp((0.95 - 1) * (22 + 2)) = 0.301194

raw Ref_zip_af floor:          F1=0.912000 TP=57 FP=3 FN=8
reliable ZIP floor:            F1=0.920635 TP=58 FP=3 FN=7
reliable zero-truncated Poisson: F1=0.920635 TP=58 FP=3 FN=7
reliable adaptive zero-truncated NB: F1=0.906250 TP=58 FP=5 FN=7
```

The zero-truncated Poisson result is effectively the same correction as `Reliable_Ref_zip_af` when using the full non-fake positive-depth mean. It ties the current best default-floor result but does not improve it.

The adaptive zero-truncated NB branch is too permissive here. High overdispersion in a few non-gold Lactobacillus rows drives `Reliable_ztnb_af` toward 1.0 and adds two FPs:

```text
s__Lactobacillus sp945980025
s__Lactobacillus sp947575965
```

A threshold sweep showed that the default ANI-coupled floor is not the best floor for the defaked ZTP/ZIP metric on this sample:

```text
Reliable ZTP/ZIP floor 0.301194: F1=0.920635 TP=58 FP=3 FN=7
Reliable ZTP/ZIP floor 0.400000: F1=0.926829 TP=57 FP=1 FN=8
Reliable ZTP/ZIP floor 0.450000: F1=0.907563 TP=54 FP=0 FN=11
```

At floor `0.40`, the remaining FP is `s__Lactobacillus mulieris_A`; the extra FN relative to the default floor is `s__Allobosea robiniae`.

Interpretation: yes, defake before depth modeling is the right data flow. For this sample, full-depth zero-truncated Poisson/reliable ZIP is stable, but zero-truncated NB should not be a hard gate without additional regularization or a higher floor because it inflates adjusted breadth for overdispersed FP rows.

## Product Top-Fraction Trial

I then tested the alternative suggested after inspecting `best_diff * depth` distributions: instead of a hard product cutoff such as `product <= 1`, drop the highest per-reference fraction of hit contexts ranked by `best_diff * depth`.

Command parameters:

```text
--readwise-ctx-filter product-topfrac
--readwise-fake-threshold 0.25
```

Implementation detail: the threshold is computed from all hit contexts for each reference, but only positive-product contexts can be rejected. Zero-diff contexts are preserved.

For the difficult `s__Limosilactobacillus kinnaridis` row, the local behavior is exactly what we wanted. The row had 42 hit contexts with product distribution:

```text
0:22, 1:9, 2:3, 4:2, 6:1, 8:2, 12:1, 14:1, 21:1
```

Dropping the top 25% removes 11 contexts, equivalent in this row to removing contexts with `product >= 2`. This rescues the local direct-call gate:

```text
unfiltered: ANI=0.907408 XnY=42 Ndiff=20 adjusted_AF=0.102492
top25:      ANI=0.965254 XnY=31 Ndiff=9  Reliable_ztp_af=0.423233
```

The whole-sample benchmark is a negative result, however. The run completed successfully:

```text
output rows excluding header: 16411
wall time: 1:56.44
peak RSS: 4,294,328 KB
rows with rejected contexts: 15923
rejected contexts: 20139
```

GTDB species scoring:

```text
default floor 0.301194, reliable ZTP/ZIP: F1=0.892562 TP=54 FP=2 FN=11
floor 0.40, reliable ZTP/ZIP:             F1=0.888889 TP=52 FP=0 FN=13
floor 0.40, raw Ref_zip_af:               F1=0.918033 TP=56 FP=1 FN=9
```

Compared with the current best product-NB reliable-depth result (`F1=0.926829`, TP=57, FP=1, FN=8 at floor `0.40`), top-25% product filtering is too aggressive as a global hard filter. It is useful as a targeted diagnostic for rows like `kinnaridis`, but it drops enough contexts across the sample to lose recall.

## Product Top-Fraction Median Replacement Trial

I then tested the less destructive version of the same idea: identify the top per-reference `best_diff * depth` fraction, but do not remove those contexts. Instead, keep each context and replace its `best_diff` and depth with the per-reference lower median among hit contexts before computing downstream naive ANI and reliable abundance statistics.

Command parameters:

```text
--readwise-ctx-filter product-topfrac-median
--readwise-fake-threshold 0.25
```

The full run completed successfully:

```text
output rows excluding header: 65414
wall time: 1:56.09
peak RSS: 4,294,900 KB
rows with adjusted contexts: 64926
adjusted contexts: 70733
```

This fixed the main weakness of hard top-25% removal: support and breadth were preserved. Against product-NB at the same gates, median replacement gained `s__Limosilactobacillus kinnaridis` and did not lose any previous TP, but added three very-low-support FPs:

```text
s__Bacteroides caecimuris_A       XnY=12
s__NK4A144 sp000701625            XnY=13
s__Parabacteroides distasonis_A   XnY=12
```

Main GTDB species scores:

```text
default floor 0.301194, reliable ZTP/ZIP: F1=0.907692 TP=59 FP=6 FN=6
floor 0.40, reliable ZTP/ZIP:             F1=0.913386 TP=58 FP=4 FN=7
```

The added FPs are small enough that a modest support guard helps. A context-support sweep showed:

```text
XnY_ctx >= 15, floor 0.40, reliable ZTP/ZIP: F1=0.935484 TP=58 FP=1 FN=7
XnY_ctx >= 20, floor 0.40, reliable ZTP/ZIP: F1=0.935484 TP=58 FP=1 FN=7
XnY_ctx >= 30, floor 0.40, reliable ZTP/ZIP: F1=0.926829 TP=57 FP=1 FN=8
```

At `XnY_ctx >= 15` and floor `0.40`, the remaining FP is the same high-support Lactobacillus row:

```text
s__Lactobacillus mulieris_A  ANI=1.0 XnY=461 Reliable_ztp_af=0.422952
```

Interpretation: median replacement is a better form of the top-product idea than hard removal because it preserves breadth and rescues `kinnaridis`. The price is that replacing high-product contexts with median diff/depth can make very low-support rows look perfectly clean. A support guard around `XnY_ctx >= 15` is therefore required for this sample.

## Product Plus-One Top-Fraction Median Trial

The remaining FP `s__Lactobacillus mulieris_A` showed a high-depth exact-match tail. Because `best_diff * depth` is zero for exact-match contexts, I tested the modified ranking score:

```text
product1 = (best_diff + 1) * depth
--readwise-ctx-filter product1-topfrac-median
--readwise-fake-threshold 0.25
```

This mode keeps the same top-fraction median replacement rule, but high-depth exact-match contexts can now be selected and median-replaced too.

The full run completed successfully:

```text
output rows excluding header: 65414
wall time: 1:56.12
peak RSS: 4,294,580 KB
rows with adjusted contexts: 65414
adjusted contexts: 79320
```

For `s__Lactobacillus mulieris_A`, the modified score did what it was meant to do mechanically:

```text
previous median top25:
  ANI=1.000000 XnY=461 adjusted_ctx=44  Ndiff=0  Reliable_hit_var=664.284 Reliable_ztp_af=0.422952

product1 median top25:
  ANI=0.981402 XnY=461 adjusted_ctx=122 Ndiff=29 Reliable_hit_var=4.019   Reliable_ztp_af=0.425142
```

So `(best_diff + 1) * depth` removes the suspicious high-depth variance signal, but the row still passes the current `Reliable_ztp_af >= 0.40` gate. Worse, globally this variant is too permissive because many low-support related rows become perfect or near-perfect after median replacement:

```text
ctxcut 10, floor 0.40: F1=0.832215 TP=62 FP=22 FN=3
ctxcut 15, floor 0.40: F1=0.884058 TP=61 FP=12 FN=4
ctxcut 30, floor 0.40: F1=0.906250 TP=58 FP=5  FN=7
```

Interpretation: adding `+1` to the product score is useful as a diagnostic for high-depth exact-match tails, but it is not a better hard gate in this form. It over-corrects exact-match high-depth contexts across the sample and increases FP substantially.

## ANI-AF Delta Gate Correction

The `s__Lactobacillus mulieris_A` value `ANI=0.981402` is the product-adjusted naive ANI row, not the standard product-NB row. In the files on disk this value is in `toymouse_sample0_s2000_marker_split_naive_product1_topfrac_median025.tsv`; the non-plus-one `product_topfrac_median025` row for the same accession still has `ANI=1.0`.

Using the AF-to-ANI conversion coupled to the same template threshold,

```text
ANI_from_AF = 1 + ln(Reliable_ztp_af) / (ctx_length + 2)
```

with `ctx_length + 2 = 24`, the product-adjusted `L. mulieris_A` row is:

```text
ANI_naive_calc=0.981401665
Reliable_ztp_af=0.425142412
ANI_from_AF=0.964361205
delta=ANI_naive_calc - ANI_from_AF = 0.017040460
```

Therefore it fails a strict delta gate of `ANI_naive_calc - ANI_from_AF < 0.016`. The earlier interpretation that this FP passed came from using the standard product-NB `ANI_naive_calc=0.967678`, not the product-adjusted value. With MinCO context length `22` instead of effective length `24`, the delta is even larger (`0.020280350`), so it also fails.

Whole-sample scoring showed that the delta rule is useful for identifying this FP but is not safe as a global hard gate in this form:

```text
product1 top25 median, ctx>=15, Reliable_ztp_af>=0.40, delta<0.016, k_eff=24:
  F1=0.882353 TP=60 FP=11 FN=5
  Previous without delta: F1=0.884058 TP=61 FP=12 FN=4
```

The delta gate removes `L. mulieris_A`, but also removes the true-positive `s__Enterocloster clostridioformis`, so F1 decreases slightly.

For the non-plus-one product median table, applying the same delta rule at the previous best gate is much too aggressive:

```text
product top25 median, ctx>=15, Reliable_ztp_af>=0.40, delta<0.016, k_eff=24:
  F1=0.818182 TP=45 FP=0 FN=20
  Previous without delta: F1=0.935484 TP=58 FP=1 FN=7
```

This removes the remaining `L. mulieris_A` FP, but it also removes many real low-breadth species with adjusted `ANI_naive_calc=1.0`. Interpretation: the ANI-AF delta is a good diagnostic feature for suspicious clean-but-narrow rows, but it should not be used as a standalone hard gate without an abundance/support rescue model.

A looser product0 delta gate at `0.03` was also tested:

```text
product top25 median, ctx>=15, Reliable_ztp_af>=0.40, delta<0.03, k_eff=24:
  F1=0.888889 TP=52 FP=0 FN=13
  Previous without delta: F1=0.935484 TP=58 FP=1 FN=7
```

This gives perfect precision and removes `s__Lactobacillus mulieris_A`, but loses six true positives that have clean product0-adjusted `ANI_naive_calc=1.0` and low-to-moderate reliable adjusted breadth:

```text
s__Bacteroides fragilis              delta=0.033745
s__Blautia hansenii                  delta=0.030545
s__Enterocloster clostridioformis    delta=0.030881
s__Lactobacillus crispatus           delta=0.032148
s__Pseudobutyrivibrio ruminis        delta=0.032666
s__Salmonella enterica               delta=0.037497
```

Therefore `delta < 0.03` should be treated as a high-precision diagnostic variant, not the best-F1 active product0 rule.

The better form is to apply the `0.03` delta only to high-depth, high-variance rows. Here the high-variance statistic is the variance-to-mean ratio, not conventional CV:

```text
if product0 adjusted median coverage > 3 and adjusted depth variance / mean > 50:
  require ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03
```

With the product0 `ctx>=15`, `ANI_naive_calc>0.95`, `Reliable_ztp_af>=0.40` gate, this conditional rule gives:

```text
F1=0.943089 TP=58 FP=0 FN=7
```

It removes only `s__Lactobacillus mulieris_A`: product0 adjusted median coverage is `6.0`, adjusted hit-depth variance/mean is `65.41`, and delta is `0.035854`. The six true positives lost by the unconditional `delta<0.03` rule all have adjusted hit mean depth at most `1.86`, so their hit median cannot exceed `3`; they do not trigger the conditional delta check.

Using the current TSV fields, the inclusive proxy `adjusted hit mean depth > 3` plus `variance/mean > 50` triggers three selected rows. Two are true positives and already pass the delta check:

```text
s__Lacticaseibacillus rhamnosus    variance/mean=286.93 delta=0.021387
s__Lacticaseibacillus paracasei    variance/mean=51.31  delta=0.008111
s__Lactobacillus mulieris_A        variance/mean=65.41  delta=0.035854
```

Therefore the exact median coverage output would not change this selected set unless a new row both fails delta and has mean depth above 3 but median coverage at or below 3.

## Current Working Gate Strategy

This is the active product0-adjusted direct-call strategy for the Toy Mouse sample0 GTDB-species benchmark.

Step 1: product0 context adjustment.

```text
product0_score = best_diff * depth
mode = top-25% product median replacement
```

For each reference, rank hit contexts by `best_diff * depth`. The top `25%` product tail is not removed. Instead, each selected context is replaced with that reference's lower-median hit-context `diff` and lower-median hit-context coverage. In product0 mode, exact-match contexts with `best_diff=0` are not adjusted; product1 `(best_diff+1)*depth` remains diagnostic only.

Step 2: compute adjusted naive ANI.

Use the product0-adjusted `XnY_ctx`, `N_diff_obj`, and `N_diff_obj_section`:

```text
dist0 = N_diff_obj / (XnY_ctx + N_diff_obj)
ratio = (N_diff_obj_section + eps) / (N_diff_obj + eps)
final_dist = 1 - (1 - dist0)^ratio
ANI_naive_calc = 1 - final_dist * 0.1544286
```

When `final_dist <= 0`, `ANI_naive_calc` is `1.0`. This keeps the zero-intercept naive distance form: `naive_dist = final_dist * 0.1544286`.

Step 3: compute post-adjustment reliable breadth.

After product0 median replacement, compute reliable hit-context depth statistics:

```text
Reliable_Ref_breadth
Reliable_Ref_hit_mean_depth
Reliable_Ref_hit_depth_variance
product0 adjusted median coverage
```

Fit a zero-truncated Poisson model to positive adjusted hit-context depths by solving:

```text
Reliable_Ref_hit_mean_depth = lambda / (1 - exp(-lambda))
P_observed = 1 - exp(-lambda)
```

Then adjust breadth:

```text
Reliable_ztp_af = min(1, Reliable_Ref_breadth / P_observed)
```

This is the breadth-adjusted AF used by the call gate.

Step 4: base direct-call gate.

```text
XnY_ctx >= 15
ANI_naive_calc > 0.95
Reliable_ztp_af >= 0.40
```

The active `0.40` AF floor is empirical for this sample and is stricter than the ANI-derived floor for `ANI_threshold=0.95` and `ctx_length+2=24`:

```text
exp((0.95 - 1) * 24) = 0.301194
```

Step 5: conditional high-depth variance consistency gate.

Convert breadth-adjusted AF back to ANI with the effective context length:

```text
effective_ctx_length = ctx_length + 2 = 24
ANI_from_Reliable_ztp_af = 1 + ln(Reliable_ztp_af) / effective_ctx_length
delta = ANI_naive_calc - ANI_from_Reliable_ztp_af
```

Only apply the delta requirement to high-depth, high-variance rows:

```text
if product0 adjusted median coverage > 3
and Reliable_Ref_hit_depth_variance / Reliable_Ref_hit_mean_depth > 50:
    require delta < 0.03
```

The statistic above is variance-to-mean ratio, not conventional coefficient of variation. The exact gate should use the product0 median coverage already computed internally by the median-replacement filter; it should be emitted as an output column before this rule is made a CLI/default option.

Step 6: current benchmark result.

```text
F1=0.943089
TP=58
FP=0
FN=7
```

This is the current best MinCO strategy in this Toy Mouse sample0 GTDB-species benchmark, now slightly above the Sylph GTDB-species result while preserving the `s__Limosilactobacillus kinnaridis` rescue. The product1 branch should be kept as a diagnostic for high-depth exact-match tails, not as the current direct-call strategy.

As of the 2026-06-22 pairwise markerdb retest, this old global context markerdb result remains the current best MinCO S2000 markerdb result:

```text
current best MinCO S2000 markerdb: old active global ctx markerdb
F1=0.943089
TP=58
FP=0
FN=7
```

## Abundance Correlation Against GTDB Truth

I compared the current product0 gate against the GTDB species truth abundance profile and Sylph sample0 output:

```text
truth: mouse0_gtdb_species_profile.tsv
MinCO: toymouse_sample0_s2000_marker_split_naive_product_topfrac_median025.tsv
Sylph: /mnt/new3T/minco_cami2_toymouse_20260621/sylph_sample0/profile.tsv
```

Both MinCO and Sylph refs were mapped to GTDB species. For species-level abundance, rows mapping to the same GTDB species were collapsed by maximum abundance, matching the earlier profile-scorer convention and avoiding duplicate representative inflation. Missing truth taxa were assigned predicted abundance `0`.

Call-set comparison:

```text
MinCO current product0 gate: F1=0.943089 TP=58 FP=0 FN=7
Sylph default:               F1=0.936508 TP=59 FP=2 FN=6
```

Abundance comparison on all 65 GTDB truth species, with missed taxa set to zero:

```text
method                  predicted abundance sum  Pearson   Spearman  MAE percentage points  L1 percentage points
MinCO product0 current  0.705454                 0.998262  0.888085  0.461136               29.973822
Sylph default           0.999522                 0.999997  0.979936  0.009365                0.608751
```

If the MinCO selected abundance vector is renormalized to sum to `1`, Pearson and Spearman are unchanged, but absolute error improves:

```text
MinCO renormalized: Pearson=0.998262 Spearman=0.888085 MAE=0.155189 percentage points L1=10.087294 percentage points
Sylph renormalized: Pearson=0.999997 Spearman=0.979936 MAE=0.009364 percentage points L1=0.608649 percentage points
```

Largest MinCO abundance errors:

```text
s__Lacticaseibacillus paracasei       truth=0.573008 MinCO=0.375515 Sylph=0.574453
s__Lactobacillus helveticus           truth=0.237316 MinCO=0.169297 Sylph=0.236623
s__Secundilactobacillus pentosiphilus truth=0.080037 MinCO=0.073119 Sylph=0.079825
s__Lactobacillus sp910589675          truth=0.005248 MinCO=0.000000 Sylph=0.005306
```

Interpretation: the current MinCO gate slightly beats Sylph on species presence/absence for this GTDB truth, but not on abundance calibration. Sylph's abundance estimates are essentially on target for the dominant species. MinCO's abundance ranking is still highly correlated because the largest species dominate Pearson, but the selected MinCO abundance mass is only `0.705454` and the two dominant Lactobacillaceae are under-estimated.

## Abundance Quantification Debug

The MinCO readwise abundance code currently computes:

```text
Ref_mean_depth = sum(per-context coverage) / ref_ctx_total
Relative_abundance_depth = Ref_mean_depth / sum(Ref_mean_depth over refs)
Normalized_abundance_depth = Relative_abundance_depth renormalized over emitted survivor rows
```

The key implementation detail is that `Normalized_abundance_depth` is filled from the raw `abundance_stats` path, while product0-adjusted depth is reported only through `Reliable_Ref_*` fields. Product0 median replacement changes `ANI_naive_calc` and reliable depth/breadth columns, but it does not currently create a `Reliable_Normalized_abundance_depth` column.

Because the current product0 gate is applied offline from the unfiltered TSV, the selected rows keep abundance values normalized over the broad unfiltered output. That is why the current selected species sum to only `0.705454`.

Comparison against non-markerdb S1000 MinCO, using the same GTDB species truth and max abundance per GTDB species:

```text
method                         call F1   pred_sum  Pearson   Spearman  L1 pct-pt
S1000 split direct             0.883333  0.615619  0.998566  0.875246  39.867100
S1000 unique direct            0.862069  0.713732  0.998025  0.875784  29.660957
S2000 marker product0 current  0.943089  0.705454  0.998262  0.888085  29.973822
```

So the current S2000 product0 gate is better than S1000 split for abundance error, but not better than S1000 unique if using the raw `Normalized_abundance_depth` field from the unfiltered output.

If the called species are renormalized to sum to 1:

```text
method                                Pearson   Spearman  L1 pct-pt
S1000 split direct, raw depth          0.998566  0.875246  10.558131
S1000 unique direct, raw depth         0.998025  0.875784  10.018768
S2000 product0 current, raw depth      0.998262  0.888085  10.087294
S2000 product0 current, reliable mean  0.998609  0.892020   9.109039
```

The product0 reliable mean depth gives the best MinCO abundance among these tested variants, but it still trails Sylph strongly (`L1=0.608649` percentage points after equivalent normalization).

Dominant-species check after selected-species renormalization:

```text
species                              truth     S2000 product0 reliable  S1000 split raw  S1000 unique raw
s__Lacticaseibacillus paracasei      0.573008  0.537653                 0.536850         0.559756
s__Lactobacillus helveticus          0.237316  0.240749                 0.224249         0.211508
s__Secundilactobacillus pentosiphilus 0.080037 0.101192                 0.104961         0.108401
s__Lactobacillus sp910589675         0.005248  0.000000                 0.000000         0.000000
```

Interpretation: markerdb/product0 improves species calling and improves abundance over S1000 split, but the main raw abundance errors are not fully fixed. The abundance path needs a product0-aware output, likely `Reliable_Relative_abundance_depth` and `Reliable_Normalized_abundance_depth`, computed after the final gate rather than from raw all-row abundance.

## ZIP Lambda Abundance Test

Because Sylph abundance is coverage/lambda based, I tested a MinCO abundance estimate that uses the same zero-truncated Poisson lambda already needed for reliable breadth adjustment.

For each selected product0 row:

```text
positive_mean = Reliable_Ref_hit_mean_depth
solve lambda from:
  positive_mean = lambda / (1 - exp(-lambda))
abundance_weight = lambda
```

Then collapse rows to GTDB species by max abundance weight and renormalize over the called species. I also tested controls using raw mean depth, reliable mean depth, hit mean depth, breadth-derived lambda, and lambda multiplied by AF/breadth.

Results against GTDB species truth:

```text
method                                      Pearson   Spearman  MAE pct-pt  L1 pct-pt
ZTP lambda from reliable hit mean           0.999649  0.904805  0.079286    5.153596
ZTP lambda * adjusted AF                    0.998608  0.895820  0.139468    9.065434
Reliable_Ref_mean_depth                     0.998609  0.892020  0.140139    9.109039
ZTP lambda * raw breadth                    0.998587  0.902015  0.143388    9.320231
Raw Ref_mean_depth                          0.998262  0.888088  0.155199   10.087936
Reliable_Ref_hit_mean_depth                 0.999667  0.904805  0.170101   11.056533
Raw Normalized_abundance_depth as-is        0.998262  0.888085  0.461136   29.973822
Lambda from breadth zero fraction only      0.071367  0.870449  2.370309  154.070102
```

Dominant species under the best MinCO ZIP-lambda abundance:

```text
species                                truth     ZIP-lambda  reliable mean  raw normalized
s__Lacticaseibacillus paracasei        0.573008  0.554401    0.537653       0.375515
s__Lactobacillus helveticus            0.237316  0.242318    0.240749       0.169297
s__Secundilactobacillus pentosiphilus  0.080037  0.086070    0.101192       0.073119
s__Lactobacillus sp910589675           0.005248  0.000000    0.000000       0.000000
s__Butyrivibrio sp000423945            0.011809  0.012807    0.010796       0.007756
```

This is a large improvement over reliable mean depth (`L1=9.109039`) and raw normalized abundance (`L1=29.973822`), but it still trails Sylph (`L1=0.608649` in the same all-truth-species view). The remaining largest failure is `s__Lactobacillus sp910589675`, which is still a false negative under the current call gate and therefore has abundance zero.

Conclusion for implementation: add a product0/reliable abundance output that reports ZTP lambda abundance, normalized after final called species selection:

```text
Reliable_ZTP_lambda
Reliable_ZTP_lambda_relative_abundance
Reliable_ZTP_lambda_normalized_abundance
```

## Validation

- `bin/minco set --downsample -S 2000` completed with exit status 0.
- `bin/minco sketch -i` succeeded for both S2000 and S2000 context markerdb.
- `bin/minco matrix --format dedup-plan -m aaf --cut 0.03` completed with exit status 0 and emitted the new markerdb prediction warning.
- The markerdb prediction from `matrix` matched the actual post-markerdb psmp summary.
- `bin/minco ani` completed for both `best-diff-unique` and `best-diff-split` readwise modes.
- `build_and_score_gtdb_ground_truth.py` mapped all 75 Toy Mouse source genomes to GTDB r232 species and regenerated `gtdb_direct_scores.tsv`, `gtdb_direct_scores_naive_ani.tsv`, `gtdb_direct_scores_naive_ani_gt095.tsv`, and `gtdb_direct_scores_naive_ani_gt095_breadth_depth.tsv`.
- `score_sourceaware_direct.py` regenerated the direct source-aware comparison from the old S1000 outputs and the new S2000 markerdb outputs.
- The S2000 markerdb `GCF_018987235.1` read-origin trace was joined against CAMISIM `reads_mapping.tsv.gz`; all 301 unique traced reads were assigned to source genomes.
- `--readwise-ctx-filter poisson-depth` was added and tested on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- `--readwise-ctx-filter product-nb` was added, built, and tested on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- `--readwise-ctx-filter product-topfrac` was added, built, and tested with a top-25% product cutoff on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- `--readwise-ctx-filter product-topfrac-median` was added, built, and tested with top-25% product median replacement on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- `--readwise-ctx-filter product1-topfrac-median` was added, built, and tested with `(best_diff+1)*depth` top-25% product median replacement on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- Experimental post-defake reliable abundance columns were added, built, and tested on the full Toy Mouse sample0 S2000 markerdb split readwise run.
- `score_reliable_truncated_depth_gate.py` scored reliable ZIP, reliable zero-truncated Poisson, and adaptive reliable zero-truncated NB gates against GTDB r232 species truth.

## Important Artifacts

See `artifacts.md` for paths to generated files, logs, sketches, and scoring outputs.

## Conclusion

The plan is technically sound and ran without memory risk on this machine: peak preprocessing RSS stayed below 10 GB, and readwise RSS stayed below 4.3 GB.

The benchmark result should be read against the GTDB species truth. Context markerdb improved the S1000 unique baseline from F1 `0.862069` to `0.878049`, mainly by increasing recall. It did not improve over the existing S1000 split baseline, which was F1 `0.883333`; the S2000 markerdb split result was slightly lower at `0.878049` because false positives increased from 2 to 4 and false negatives stayed near the same level.

AAF `0.03` was conservative on this full GTDB S2000 sketch: only 182 refs were removed. That means most of the markerdb change came from context uniqueness filtering, not from deduplication.

For the `L. crispatus` diagnostic, markerdb helped substantially but did not eliminate contamination. The target-ref trace is still mostly `L. helveticus` reads, although the contaminating fraction dropped from about 95% in the old S1000 trace to about 67% in this S2000 markerdb trace.

The Poisson-depth chomp rule is mechanistically promising for detecting narrow, over-deep contaminated contexts, but this first whole-sample GTDB test is neutral for zip-aaf calls and worse for selected naive calls. The `best_diff * depth` product rule moves the target row much closer to the oracle clean ANI, but worsens whole-sample GTDB precision.

Always using naive ANI alone is not a good direct-call default under the tested thresholds. At `ANI >= 0.94`, it recovers a few additional true low-abundance species but the precision loss is larger than the recall gain, especially for S1000 split where FP increases from 2 to 78. Tightening to strict `ANI > 0.95` reduces FP slightly but still trails ZIP-AAF direct scoring.

The breadth-depth consistency rule is the first variant here that improves the S2000 markerdb path over the S1000 split baseline on GTDB species truth: `0.892308` F1 for S2000 marker split with Poisson-depth context filtering plus breadth-depth direct scoring, versus `0.883333` for S1000 split ZIP-AAF. This is still a single-sample diagnostic result and should be validated before becoming a default.

Requiring adjusted breadth `Ref_zip_af >= 0.5` is too conservative here. It gives perfect precision but lowers recall enough to drop F1 to `0.869565`, below the retained breadth-depth result.

The post-defake reliable-depth trial supports the ordering of operations: remove fake contexts first, then model depth. However, full-depth zero-truncated NB is not yet safe as a hard gate. It lowered F1 from `0.920635` to `0.906250` at the default floor by adding two overdispersed Lactobacillus FPs. Defaked reliable ZIP/ZTP at a stricter floor `0.40` gave the best local score from this trial, `0.926829`, but still trails Sylph GTDB-species F1 `0.936508` and needs multi-sample calibration.

The product top-25% filter rescued `s__Limosilactobacillus kinnaridis` locally, but it is too aggressive as a global hard filter in this sample. It lowered reliable ZTP/ZIP F1 to `0.888889` at floor `0.40` by increasing FN from 8 to 13.

Replacing top-product contexts with per-reference median diff/depth is much better than removing them. With an additional support guard `XnY_ctx >= 15`, reliable ZTP/ZIP at floor `0.40` reached `F1=0.935484` (TP=58, FP=1, FN=7), essentially matching the Sylph GTDB-species score from this sample while preserving the `kinnaridis` rescue.

Changing the top-product score to `(best_diff + 1) * depth` catches high-depth exact-match tails, but as a hard median-replacement gate it is worse here. It lowered the best support-swept F1 to `0.906250` and added many low-support FPs.

The later pairwise context markerdb retest improved marker retention but did not improve Toy Mouse GTDB-species calling. The old active product0 gate collapsed on the pairwise markerdb because low-depth, tiny-breadth rows saturated the ZTP breadth adjustment, giving `F1=0.256881` with 315 FPs. The best small retuning tested here recovered precision but only reached `F1=0.925620`, still below the old global context markerdb active gate (`0.943089`) and Sylph (`0.936508`).

Therefore the active MinCO S2000 markerdb benchmark should continue to use the old global context markerdb product0 active result: `F1=0.943089`, `TP=58`, `FP=0`, `FN=7`.

## Paper-Relevant Claim

Diagnostic only: context markerdb can reduce one known low-abundance Toy Mouse ANI underestimation. The plain S2000 markerdb ZIP-AAF direct call improves S1000 unique but not S1000 split; adding naive ANI with a high-depth breadth floor is promising and beats S1000 split on this sample, but needs replication before changing the default.

## Caveats

- This is one CAMI II Toy Mouse sample0 run.
- The authoritative benchmark truth is now `mouse0_gtdb_species_profile.tsv`; source-aware NCBI-species scoring is retained only for diagnostics.
- The GTDB truth uses unique GCF/GCA assembly-core fallback for 12/75 source genomes. These are recorded explicitly in `mouse0_gtdb_source_species_ground_truth.tsv`.
- The always-naive test is a scoring rerun from existing unfiltered MinCO output columns. Because the original commands used `-m0 -f0 -n0 -t0`, this should be equivalent to rerunning with `--readwise-ani naive` for the direct-call decision, without repeating the large readwise jobs.
- The breadth-depth rule uses empirical cutoffs from this sample: `Ref_hit_mean_depth >= 4` and `Ref_breadth < 0.5`. It needs cross-sample calibration before it should become a CLI default.
- The adjusted-breadth floor `Ref_zip_af >= 0.5` removed all false positives but also removed low-abundance true positives; do not use it as the primary direct-call rule without a low-abundance rescue path.
- The ANI-coupled adjusted-breadth template is exploratory. The current best local form uses `adjusted_breadth_floor = exp((ANI_threshold - 1) * (ctx_length + 2))`, but the `+2`, the p1..p75 lambda fit, and the absence of `Real_min_align_fraction` need multi-sample validation.
- The product-NB defake model is implemented and much less aggressive than Poisson-product here, but with `alpha=0.05` it rejected only 33 contexts and did not improve beyond the adjusted-breadth gate. Tune or model it as a feature before treating it as a default hard filter.
- The post-defake reliable abundance columns are experimental output fields. The adaptive zero-truncated NB fit used moment matching from aggregate positive-depth mean/variance, not the full depth histogram, and over-corrected several overdispersed FP rows.
- The product top-fraction filter is rank-based and does not estimate a statistical tail probability. At `0.25` it removed 20,139 contexts across 15,923 rows and should be treated as a diagnostic, not a default hard filter.
- The product top-fraction median replacement mode is also rank-based. It preserves support but can make very low-support rows look too clean after median replacement; in this sample it needs a support guard around `XnY_ctx >= 15`.
- The `(best_diff + 1) * depth` top-fraction median mode should be treated as a diagnostic only. It reduces `L. mulieris_A` depth variance but does not remove that FP and increases FP count globally.
- The AAF `0.03` dedup plan removed too few GTDB representatives to test a strong representative-collapse hypothesis.
- Generated sketches and large TSVs are under `/tmp/gtdb232_s2000_dedup_marker.qKJofv`; preserve them elsewhere if this result needs long-term storage.
- `poisson-depth` currently filters the naive/object-difference context features and effective support. It does not recompute ZIP-AAF breadth/depth from the filtered context set, so selected zip-aaf ANI does not change for rows that remain reportable.

## Next Experiment

Try a stronger dedup sweep, for example AAF cuts `0.04`, `0.05`, and possibly a species-aware representative keep list, then rerun markerdb prediction before building the markerdb. The decision point should be whether the warning set below 500 markers and the Toy Mouse GTDB-species F1 improve without a large FP increase.

For the Poisson-depth filter, the next specific test should recompute ZIP-AAF abundance from the filtered context set or add the rejection fraction/maximum-depth-outlier score as a feature in a call model rather than using it as a hard direct-call filter alone.

For the product rule, test it as a model feature first: maximum product outlier, rejected-product fraction, post-product-filter naive ANI, and a rank-tail top-fraction diagnostic. Do not use it directly as a default hard filter before multi-sample validation.

For the reliable-depth model, test a regularized NB or an empirical floor around `0.40` across more CAMI samples. Do not use unrestricted zero-truncated NB as a default gate from this single-sample result.

## Pairwise Context MarkerDB Retest (2026-06-22)

The pairwise context markerdb keeps a context unless it is shared by a sufficiently close reference pair. This was intended to avoid losing markers because of remote, biologically unrelated references. The tested markerdb was:

```text
/tmp/gtdb232_s2000_pairctx_af005.YHvD9c/sketch_T_S2000_aaf003_dedup_pairctx_af005
```

Build parameters were the S2000 AAF-0.03 deduped GTDB sketch, `--markerdb-ctx-pairwise`, default `--markerdb-ctx-min-af=0.05`, and `--markerdb-ctx-min-xny=10`. Marker retention improved relative to the global context markerdb:

```text
global context markerdb:   9,936 / 200,527 refs below 500 markers; mean size 1,192.52
pairwise context markerdb: 3,283 / 200,527 refs below 500 markers; mean size 1,656.16
```

The Toy Mouse sample0 product0 active readwise run completed cleanly:

```text
bin/minco ani -p16 \
  -r /tmp/gtdb232_s2000_pairctx_af005.YHvD9c/sketch_T_S2000_aaf003_dedup_pairctx_af005 \
  --qraw /mnt/new3T/minco_cami2_toymouse_20260621/sample_0/2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0 \
  -o /tmp/gtdb232_s2000_pairctx_af005.YHvD9c/toymouse_sample0_pairctx_af005_split_naive_product_topfrac_median025.tsv
```

Runtime was `2:20.48`, peak RSS was `5,813,920 KB`, and exit status was `0`.

The same active product0 gate was first rechecked on the old global context markerdb output and reproduced the previous result exactly:

```text
old global ctx markerdb active gate: F1=0.943089 TP=58 FP=0 FN=7
```

On the new pairwise markerdb, the same gate failed badly:

```text
pairwise ctx markerdb active gate: F1=0.256881 TP=56 FP=315 FN=9
```

The failure mode is that many false positives have extremely small raw reliable breadth but positive hit depth near or below one, so the zero-truncated Poisson breadth adjustment saturates to `Reliable_ztp_af=1.0`. Across the 315 FPs, median `Reliable_Ref_breadth` was `0.002104`, median `XnY_ctx` was `21`, and median normalized abundance was `0.000004`.

A small retuning sweep added raw support floors. The best tested point was:

```text
XnY_ctx >= 30
ANI_naive_calc > 0.95
Reliable_ztp_af >= exp((0.95 - 1) * 24) = 0.301194
Ref_zip_af >= 0.20
active conditional high-depth variance delta rule unchanged
```

This removed all FPs but still did not beat the old active markerdb:

```text
pairwise best retuned gate: F1=0.925620 TP=56 FP=0 FN=9
```

Compared with the old active gate, this retuned pairwise gate rescued `s__Allobosea robiniae`, but newly missed `s__Lactobacillus crispatus`, `s__Ligilactobacillus agilis`, and `s__Limosilactobacillus kinnaridis`. Abundance was also worse than the old active MinCO result: selected-species max abundance had predicted sum `0.562562`, Pearson `0.998042`, Spearman `0.890205`, MAE `0.677435` percentage points, and L1 `44.033258` percentage points; after renormalizing selected MinCO abundance to sum to one, MAE was `0.168013` and L1 was `10.920820`.

Interpretation: pairwise markerdb does solve a marker-retention problem, but for this Toy Mouse readwise benchmark it introduces many low-depth remote-support rows that break the current ZTP-adjusted breadth gate. Under the current gate and the best small retuning tried here, Toy Mouse GTDB species performance did not improve.
