# CAMI II Toy Mouse Gut Sample0 External Test

Date: 2026-06-21

## Question

Test the current minco S1000 readwise profiler on a new CAMI dataset, Toy Mouse Gut sample0, and compare species-level calls with Sylph under the same local GTDB-to-NCBI taxid mapping.

This run also checked whether the dataset has usable ground truth under `distributions/` or `source_genomes/`.

## Code And Tools

- minco repo: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- minco git short hash: `b55b4d3`
- minco version: `minco 0.1`
- Sylph version: `sylph 0.6.0`
- minco reference: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
- Sylph reference: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`
- Species taxmap: `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv`
- Source-aware taxmap: `/mnt/new3T/minco_cami2_toymouse_20260621/source_aware_crosswalk/mouse0_sourceaware_taxmap_ani95_minaf50.tsv`

## Ground Truth

The sample tar contains the clean official gold profile:

- `/mnt/new3T/minco_cami2_toymouse_20260621/sample_0/taxonomic_profile_0.txt`

The setup archive also contains source-genome-level inputs:

- `distributions/distribution_0.txt`: genome ID abundance weights
- `internal/meta_data.tsv`: genome ID to NCBI taxid
- `internal/genome_locations.tsv`: genome ID to `source_genomes/*.fa`
- `source_genomes/*.fa`: source assemblies inside the setup archive

For scoring, I used the official species profile because it already contains CAMI taxonomic paths and percentages. It has 64 positive bacterial species, summing to 99.9991%.

Important caveat: only 52 of the 64 official gold species are represented in the current GTDB/taxmap mapping. Twelve gold species are absent from the reference/taxmap bridge, so exact taxid recall is capped before considering algorithm behavior.

Because CAMISIM provides source genomes for each sample, I also built a
sample-specific source-aware taxmap. It compares selected GTDB representatives
against this sample's 75 positive source genomes using skani, then remaps a
called representative to the source genome's NCBI species when the nearest
source match passes `ANI >= 95` and `min(refAF, queryAF) >= 50`.

Source-aware crosswalk summary:

| Metric | Value |
|---|---:|
| Selected called representatives | 124 |
| Representatives with skani source hit | 102 |
| Source-aware overrides | 68 |
| Best hits with ANI >= 95 | 77 |
| Best hits with ANI >= 94 | 88 |

## Runtime

| Tool | Mode | Wall time | Peak RSS | CPU |
|---|---:|---:|---:|---:|
| minco | S1000 unique readwise | 1:51.00 | 3.33 GB | 438% |
| minco | S1000 split readwise | 1:51.19 | 3.62 GB | 437% |
| Sylph | sketch | 0:53.69 | 0.99 GB | 99% |
| Sylph | profile | 0:52.12 | 19.44 GB | 210% |

Sylph total wall time was about 1:46 when sketch and profile are counted. minco needs one pass per assignment mode in this experiment; production default would not normally run both.

## Species Call Results

Exact species-taxid scoring against the official CAMI profile:

| Method | Pred | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sylph default | 60 | 47 | 13 | 17 | 0.783 | 0.734 | 0.758 |
| minco train9 RF/HGB avg | 58 | 45 | 13 | 19 | 0.776 | 0.703 | 0.738 |
| minco train9 RF | 54 | 43 | 11 | 21 | 0.796 | 0.672 | 0.729 |
| minco unique direct | 50 | 39 | 11 | 25 | 0.780 | 0.609 | 0.684 |
| minco split direct | 54 | 40 | 14 | 24 | 0.741 | 0.625 | 0.678 |

On this independent Toy Mouse sample0, the current train9 minco model does not beat Sylph. This is consistent with the earlier external strain result: the current trained rule is not ready as a general default.

## Source-Aware Species Scoring

The sample-specific source-aware taxmap corrects GTDB representative labels to
the NCBI species of the nearest positive CAMISIM source genome. This is a more
appropriate benchmark view when the simulator source genomes are known.

| Method | Original F1 | Source-aware F1 | Source-aware TP | Source-aware FP | Source-aware FN |
|---|---:|---:|---:|---:|---:|
| Sylph default | 0.758 | 0.968 | 60 | 0 | 4 |
| minco train9 RF/HGB avg | 0.738 | 0.921 | 58 | 4 | 6 |
| minco train9 RF | 0.729 | 0.897 | 52 | 0 | 12 |
| minco split direct | 0.678 | 0.891 | 53 | 2 | 11 |
| minco unique direct | 0.684 | 0.870 | 50 | 1 | 14 |
| minco unique relaxed | 0.600 | 0.763 | 58 | 30 | 6 |

This confirms the user's hypothesis: the earlier exact-taxid score was strongly
confounded by GTDB representative taxonomy versus CAMISIM source-genome NCBI
taxonomy. However, after this correction, Sylph is still ahead on mouse sample0.

## Abundance

For species that both a method calls and the gold contains, abundance correlation is high:

| Method | TP Pearson | TP Spearman | TP MAE percentage points |
|---|---:|---:|---:|
| Sylph default | 0.99999 | 0.9745 | 0.0090 |
| minco split direct | 0.9930 | 0.9494 | 0.1475 |
| minco train9 RF/HGB avg | 0.9928 | 0.9419 | 0.1418 |
| minco unique direct | 0.9880 | 0.9367 | 0.1935 |

The all-gold zero-for-missing correlation is much lower because some high-abundance source genomes are assigned to nearby GTDB representative species with different exact NCBI species taxids.

## Taxonomy Caveat

The dominant gold genome is:

- genome ID `135956.0`
- source accession `GCF_001066565.1`
- NCBI species taxid `1582`
- gold abundance 57.3008%

The GTDB representative databases instead report the closest dominant call as:

- accession `GCF_000829035.1`
- species taxid `1597`
- species name `Lacticaseibacillus paracasei`

This `1582` versus `1597` mismatch is counted as one FN and one FP for both Sylph and minco under exact taxid scoring. It is not a minco-specific failure.

The source-aware crosswalk fixes this dominant mismatch:

| Called representative | Source accession | Source species taxid | ANI to source | refAF | queryAF |
|---|---|---:|---:|---:|---:|
| `GCF_000829035.1` | `GCF_001066565.1` | 1582 | 98.97 | 86.60 | 79.78 |
| `GCF_000160855.1` | `GCF_001702095.1` | 1587 | 99.05 | 76.12 | 77.54 |
| `GCF_002217985.1` | `GCF_002217985.1` | 1714682 | 100.00 | 100.00 | 100.00 |

## Conclusion

This dataset is useful as an external benchmark, but exact species-taxid scoring is partly limited by the GTDB representative/reference mapping. Source-aware scoring is the better CAMISIM benchmark view for this sample and changes the interpretation substantially: Sylph improves from F1 `0.758` to `0.968`, and minco train9 RF/HGB improves from `0.738` to `0.921`.

The correction validates that taxonomy mismatch was a major cause of the apparent mouse failure. It does not make minco beat Sylph on mouse sample0. The next improvement should focus on minco recall under source-aware scoring: direct minco remains high precision but misses more source species than Sylph.

## Why Sylph Still Leads

Against the best minco source-aware result (`train9_rf_hgb_avg`), Sylph has
three true-positive source species that minco misses:

| taxid | species | gold % | source accession(s) | Sylph adjusted ANI | Sylph naive ANI | Sylph eff cov | minco probability | best minco split ANI | best minco split AF | interpretation |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| 1598 | Lactobacillus reuteri | 0.0326 | GCF_000159615.1,GCF_000410995.1 | 98.60 | 92.38 | 0.139 | 0.269 | 0.939 | 0.103 | minco has support, but split ANI is just below 0.94 and model probability is just below 0.30 |
| 47770 | Lactobacillus crispatus | 0.0187 | GCF_002218965.1 | 96.69 | 92.14 | 0.243 | 0.025 | 0.898 | 0.170 | minco support/AF are adequate, but minco ANI is strongly underestimated |
| 28116 | Bacteroides ovatus | 0.0187 | GCF_001578575.1 | 97.52 | 91.13 | 0.128 | 0.056 | 0.873 | 0.048 | minco signal is sparse and just below the AF threshold |

This points to Sylph's main advantage here: its effective-coverage/adjusted-ANI
model rescues low-abundance species whose naive ANI is around 91-92. minco still
behaves closer to the uncorrected object-difference ANI in these sparse cases.

minco RF/HGB also has four source-aware false positives. Three are very low
support model-threshold artifacts; one (`Lactobacillus mulieris`) has substantial
support and likely reflects a related-Lactobacillus collision:

| taxid | species | probability | best split ANI | split XnY | split AF |
|---:|---|---:|---:|---:|---:|
| 2508708 | Lactobacillus mulieris | 0.434 | 0.957 | 518 | 0.518 |
| 462198 | uncultured Oribacterium sp. | 0.370 | 1.000 | 24 | 0.024 |
| 2600149 | Streptococcus sp. sy004 | 0.316 | 0.735 | 6 | 0.006 |
| 46127 | Staphylococcus felis | 0.315 | 0.736 | 6 | 0.006 |

A simple mouse-only sweep suggests a hard support/AF guard plus a slightly lower
RF/HGB threshold can improve minco from F1 `0.921` to `0.937`, but this remains
below Sylph's `0.968` and is not a deployable conclusion without multi-sample
validation.

## ANIm Truth For Sylph-Only Cases

I computed pyani-style ANIm for the source-aware source/reference pairs:
`nucmer --mum`, `delta-filter -1`, then identity from the filtered delta
similarity-error counts.

| taxid | species | source | called ref | role | ANIm ANI | minco ZIP-AAF | minco naive | Sylph adjusted | Sylph naive |
|---:|---|---|---|---|---:|---:|---:|---:|---:|
| 1598 | Limosilactobacillus reuteri | GCF_000410995.1 | GCF_000016825.1 | minco best split | 0.9667 | 0.9392 | 0.9574 | NA | NA |
| 1598 | Limosilactobacillus reuteri | GCF_000159615.1 | GCF_036621895.1 | Sylph best | 0.9694 | 0.9118 | 0.9669 | 0.9860 | 0.9238 |
| 28116 | Bacteroides ovatus | GCF_001578575.1 | GCF_001314995.1 | both | 0.9692 | 0.8725 | 0.9679 | 0.9752 | 0.9113 |
| 47770 | Lactobacillus crispatus | GCF_002218965.1 | GCF_018987235.1 | both | 0.9810 | 0.8980 | 0.9340 | 0.9669 | 0.9214 |

This changes the interpretation of the missed cases. minco ZIP-AAF is clearly
too low here. minco naive ANI is near ANIm for `L. reuteri` and `B. ovatus`,
while Sylph adjusted ANI is closest for `L. crispatus`.

## L. crispatus Read-Origin Trace

I traced the final readwise assignments for the `Lactobacillus crispatus`
target reference `GCF_018987235.1` in the S=1000 GTDB readwise run, using
per-read density mode (`--density-block-ctx 0`) so every selected assignment
could be joined exactly to CAMISIM `reads_mapping.tsv.gz`.

The per-read run reproduced the same final naive features as the block run:
`XnY_ctx=170`, `N_diff_obj=52`, `N_diff_obj_section=107`, naive ANI
`0.934031`. This validates the trace as a diagnostic for the biased readwise
ANI.

Read-origin summary for selected `GCF_018987235.1` assignment events:

| origin genome | source accession | source species | events | unique reads | fraction |
|---|---|---|---:|---:|---:|
| `463794.1` | `GCF_001702095.1` | Lactobacillus helveticus | 5676 | 5624 | 0.950 |
| `denovo11386.0` | `GCF_002218965.1` | Lactobacillus crispatus true source | 119 | 115 | 0.020 |
| `133892.0` | `GCF_000498675.1` | Lactobacillus johnsonii | 55 | 55 | 0.009 |
| `851733.0` | `GCF_002217985.1` | Lactobacillus pentosiphilus | 49 | 49 | 0.008 |

Final-context summary:

- 170 final reference contexts were retained for naive ANI.
- 106 / 170 contexts had at least one best-diff read from the true
  `L. crispatus` source.
- 64 / 170 contexts had no true-source read among the best-diff contributors.
- The true-source-supported contexts alone imply naive ANI `0.972246`.
- Contexts without true-source best reads imply naive ANI `0.898391`.

This shows the `L. crispatus` readwise ANI failure is mainly an assignment
contamination problem from abundant related Lactobacillus reads, especially
`L. helveticus` (`GCF_001702095.1`), not an intrinsic source/reference
distance problem. Complete source-vs-reference minco naive ANI remains near
ANIm (`0.984944` vs ANIm `0.980997`).

## Experimental Fake-Context Filter

I added an experimental readwise context filter:

```text
--readwise-ctx-filter poisson-diff
--readwise-fake-threshold FLOAT
```

The model is applied only to the context set used for readwise
naive/object-difference ANI. It keeps `best_diff=0` contexts in this first
version, then rejects nonzero-diff contexts whose score exceeds the threshold:

```text
score = best_diff
      + 0.30 * -log10(PoissonTail(depth >= observed_depth | lambda=-ln(zero_fraction)))
      + 0.25 * log1p(mean_depth / -ln(zero_fraction))
```

With `--readwise-fake-threshold 3.5`, the exact per-read S=1000 mouse run gives
for `L. crispatus` target ref `GCF_018987235.1`:

| metric | before | after |
|---|---:|---:|
| readwise naive ANI | 0.934031 | 0.962932 |
| XnY_ctx | 170 | 154 |
| N_diff_obj | 52 | 36 |
| N_diff_obj_section | 107 | 46 |
| rejected contexts | 0 | 16 |
| fake context fraction | 0 | 0.0941 |

This rescues the `L. crispatus` ANI above 0.95 and is closer to ANIm
(`0.980997`), but it is not ready as a global selected call metric. When the
filtered naive ANI was used as the selected ANI for the whole mouse sample,
source-aware direct calls over-predicted many taxa:

| method | precision | recall | F1 | FP | FN |
|---|---:|---:|---:|---:|---:|
| Sylph default | 1.000 | 0.938 | 0.968 | 0 | 4 |
| minco split direct, filtered naive selected | 0.487 | 0.906 | 0.634 | 61 | 6 |

The filter is therefore best treated as a new diagnostic/model feature. The
next model should use filtered naive ANI together with raw ANI, ZIP-AAF,
fake-context fraction, breadth/depth inflation, and support to avoid converting
contaminated weak refs into false positives.

## Experimental Fake-Context Probability Weighting

I added a second experimental mode:

```text
--readwise-ctx-filter fake-prob
--readwise-fake-threshold FLOAT
```

This uses the same breadth/depth/object-difference score, converts it to an
estimated fake-context probability with a logistic transform, and soft-weights
nonzero-diff contexts by `1 - P(fake)`. Zero-diff contexts keep full weight for
the naive ANI denominator in this first implementation. Output columns now also
include `Fake_ctx_prob_mean` and `Fake_ctx_prob_weighted`.

Focused S=1000 mouse sample0 result for `L. crispatus` target ref
`GCF_018987235.1`:

| mode | threshold | ANI | XnY_ctx | N_diff_obj | N_diff_obj_section | Raw_XnY_ctx | fake/reject signal |
|---|---:|---:|---:|---:|---:|---:|---:|
| raw naive | NA | 0.934031 | 170 | 52 | 107 | 170 | NA |
| hard `poisson-diff` | 3.5 | 0.962932 | 154 | 36 | 46 | 170 | rejected 0.094118 |
| soft `fake-prob` | 3.0 | 0.966331 | 145 | 27 | 38 | 170 | mean P(fake) 0.222310 |
| soft `fake-prob` | 3.5 | 0.961411 | 150 | 32 | 46 | 170 | mean P(fake) 0.167735 |

Runtime for the two soft-probability runs was essentially unchanged from the
hard-filter run: `1:51.57` and `1:51.77` wall time with about `3.63 GB` peak
RSS on 16 threads. Threshold `3.0` is the better focused setting here, but it
still underestimates the ANIm source/reference truth (`0.980997`). This should
therefore remain a candidate feature/correction path, not a default global call
rule, until validated against full-sample FP/FN behavior.

Source-aware sample-level scoring confirms that warning. With the current
direct rule (`ANI >= 0.94`, `XnY_ctx >= 10`, `Real_min_align_fraction >= 0.05`),
`fake-prob` threshold `3.0` has the same outcome as the hard-filtered naive
path: `TP=58`, `FP=61`, `FN=6`, `F1=0.633880`. A small threshold sweep could
only recover `F1=0.755245` for `fake-prob`, still below the original split
ZIP/default direct call (`F1=0.890756`) and its tuned operating point
(`F1=0.928000`). Thus the probability mode improves the targeted ANI failure,
but it does not improve whole-sample F1 as a direct selected metric.

## Validation

- `python3 -m py_compile` passed for the updated gold parser, Toy Mouse profile builder, source-aware taxmap builder, and ANIm helper.
- `make -C minco_core test` passed.
- `make -C minco_core -j4` passed after adding the hidden readwise trace hook.
- `make -C minco_core test` passed after adding `--readwise-ctx-filter`.
