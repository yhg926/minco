# Low-Abundance Context Source Contamination Diagnostic

## Question

For low-abundance true-positive references in CAMI II Toy Mouse Gut sample0, is MinCO readwise naive ANI bias mainly caused by too few hit contexts, or by a higher fraction of contexts assigned from contaminating source genomes?

## Method

I used the existing MinCO readwise trace hook:

- `MINCO_READWISE_TRACE_REF=<GTDB representative accession>`
- `MINCO_READWISE_TRACE_OUT=<trace.tsv>`

The trace records assigned context events for one reference, including `read_id`, `qctx`, `diff`, `best_diff`, number of candidate references, number of selected references, split coverage increment, and reference name.

The trace rows were joined to CAMISIM sample0 `reads_mapping.tsv.gz`, then to the GTDB-transferred source truth table. A trace event is classified as:

- `target_source`: read source genome is one of the source genomes mapped to the target GTDB species;
- `same_gtdb_species`: read source genome maps to the same GTDB species but is not necessarily the same source accession;
- `other_or_unmapped_source`: read source genome does not map to the target GTDB species.

The main contamination fraction is `other_or_unmapped_source_event_frac`.

Important implementation detail: the first trace attempt used the default `--density-block-ctx 100`, but density-block traces leave `read_id` empty. The source-contamination summaries therefore use per-read traces rerun with `--density-block-ctx 0`.

## Results

Per-read trace summary:

| target | raw source abundance | ANIm | MinCO readwise naive | error | XnY ctx | trace events | unique qctx | same GTDB event frac | other-source event frac | same GTDB diff mean | other-source diff mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `kinnaridis_low_bad` | 1.000000 | 0.980758 | 0.907408 | -0.073350 | 42 | 70 | 42 | 0.428571 | 0.571429 | 0.666667 | 3.775000 |
| `johnsonii_low_bad` | 1.259272 | 0.979937 | 0.946397 | -0.033540 | 29 | 59 | 29 | 0.322034 | 0.677966 | 0.421053 | 0.625000 |
| `crispatus_low_bad` | 1.000000 | 0.981000 | 0.954948 | -0.026052 | 100 | 303 | 100 | 0.326733 | 0.673267 | 0.242424 | 4.799020 |
| `rhamnosus_low_multisource` | 11.375254 | 0.975527 | 0.971288 | -0.004239 | 923 | 5885 | 923 | 0.299745 | 0.700255 | 0.323696 | 1.428537 |
| `paracasei_high_control` | 3057.000000 | 0.989315 | 0.991954 | +0.002639 | 1247 | 518863 | 1247 | 0.999958 | 0.000042 | 0.276366 | 1.909091 |

Top contributing sources:

- `kinnaridis_low_bad`: true `s__Limosilactobacillus kinnaridis` gives 30/70 events. Other events come mostly from `s__Limosilactobacillus reuteri_L` 12, `s__Companilactobacillus crustorum` 10, `s__Lactobacillus mulieris` 7, and `s__Lacticaseibacillus paracasei` 6.
- `johnsonii_low_bad`: true `s__Lactobacillus johnsonii` gives 19/59 events. `s__Lactobacillus sp910589675` gives 37/59 events and `s__Lactobacillus helveticus` gives 3/59.
- `crispatus_low_bad`: true `s__Lactobacillus crispatus` gives 99/303 events. `s__Lactobacillus helveticus` gives 203/303 events, from only 10 unique `qctx`, with high mean diff 4.803.
- `rhamnosus_low_multisource`: same-GTDB `s__Lacticaseibacillus rhamnosus` sources together give 1764/5885 events. `s__Lacticaseibacillus paracasei` gives 4121/5885 events, from 63 unique `qctx`.
- `paracasei_high_control`: true `s__Lacticaseibacillus paracasei` gives 518841/518863 events; only 22 events are from other sources.

## Interpretation

The low-abundance readwise ANI bias is not purely caused by too few hit contexts. Sparse true support is part of the failure mode, but the low-abundance biased rows also have large other-source event fractions:

- 57.1% for `kinnaridis`;
- 67.8% for `johnsonii`;
- 67.3% for `crispatus`;
- 70.0% for `rhamnosus`.

The high-abundance control has only 0.0042% other-source events. This contrast is strong evidence that low abundance makes the trace vulnerable to contaminating/shared contexts from related sources.

The contaminating contexts are not equally damaging in all rows. For `crispatus`, `helveticus` contamination is a clear depth-weighted artifact: only 10 unique contaminant `qctx` produce 203 read-level events, and the contaminant mean diff is 4.80 versus 0.24 for true `crispatus` events. For `johnsonii`, the largest contaminant source has lower diff, so sparse true support plus related-species cross-assignment is more important than high-diff contamination alone.

Algorithm implication: a good defake rule should use source-like evidence available internally: repeated high-depth nonzero-diff contexts, best_diff*depth outliers, and maybe genus-local winner assignment. A simple low-XnY penalty alone would miss the fact that some false context mass is high-depth and source-skewed.

## Origin Membership Check

I then checked whether each other-source trace context could have been absorbed by its origin GTDB representative in the S2000 and S10000 sketches.

For the S2000 ctx-marker run:

- `candidate_refs > 1` fraction was `0.0`;
- `origin_marker_present` fraction was `0.0`;
- `origin_full_S2000_present` fraction was `0.0`;
- `target_marker_present` fraction was `1.0`.

For S10000 origin representatives, the answer was target-dependent:

| target | other events | other unique qctx | origin S10000 event frac | origin S10000 unique qctx frac |
|---|---:|---:|---:|---:|
| `kinnaridis_low_bad` | 40 | 18 | 0.800000 | 0.833333 |
| `johnsonii_low_bad` | 40 | 10 | 0.100000 | 0.100000 |
| `crispatus_low_bad` | 204 | 11 | 0.000000 | 0.000000 |
| `rhamnosus_low_multisource` | 4121 | 63 | 0.000000 | 0.000000 |
| `paracasei_high_control` | 22 | 7 | 0.000000 | 0.000000 |

So the earlier S2000 conclusion needs qualification. These contaminant contexts were absent from the origin representative's S2000 sketch, but some are rescued by S10000. The clearest case is `kinnaridis_low_bad`: 80% of other-source events and 83.3% of other-source unique qctx are present in the origin S10000 sketch but absent from origin S2000. In contrast, the `crispatus` helveticus-contaminant qctx and the `rhamnosus` paracasei-contaminant qctx are absent even from the origin S10000 representative.

In the active S2000 ctx-marker run, MinCO assigned them to the traced target because the target markerdb contained that `qctx`, while the origin S2000 representative was not a candidate at all.

## Crispatus <- Helveticus Whole-Genome Check

For the `crispatus_low_bad` case, the top other-source contributor is CAMISIM source genome `GCF_001702095.1` (`s__Lactobacillus helveticus`) against target representative `GCF_018987235.1` (`s__Lactobacillus crispatus`). The GTDB representative for that helveticus source is `GCF_000160855.1`.

Important sketch-identity detail: a larger MinCO sketch such as `-S 10000000` is not directly comparable to S2000, because MinCO includes target sketch size in `sketch_id`. The valid check therefore scanned the whole genomes and computed S2000-compatible coden11 context hashes directly.

Results for the 203 helveticus-labeled trace events / 10 unique `qctx`:

| genome checked | mode | event present | unique qctx present |
|---|---|---:|---:|
| helveticus source `GCF_001702095.1` | raw/postconflict | 187/203 | 1/10 |
| helveticus representative `GCF_000160855.1` | raw/postconflict | 0/203 | 0/10 |
| crispatus representative `GCF_018987235.1` | raw/postconflict | 203/203 | 10/10 |

The dominant contaminating context, `qctx=1674286752`, accounts for 187/203 events. It is present in the helveticus source genome with one object, absent from the helveticus GTDB representative, and present in the crispatus representative. The other 9 helveticus-labeled `qctx` are absent from both the helveticus source genome and helveticus representative but present in the crispatus representative.

Whole-genome ANIm for the same example:

| pair | ANIm | aligned bp | ref aligned fraction | query aligned fraction |
|---|---:|---:|---:|---:|
| helveticus source vs helveticus representative | 0.992292 | 1,569,340 | 0.762438 | 0.776677 |
| helveticus source vs crispatus representative | 0.851865 | 654,435 | 0.317946 | 0.292277 |
| helveticus source vs crispatus source | 0.852187 | 663,624 | 0.322411 | 0.323275 |
| helveticus representative vs crispatus representative | 0.850095 | 641,385 | 0.317426 | 0.286449 |

Interpretation for this example: `representative/source mismatch` means the exact CAMISIM source genome and the GTDB representative for that source species are not identical, so a real source context can be absent from the representative used by the database. This explains the dominant `qctx=1674286752`. The remaining 9 low-count helveticus-labeled contexts are not source-representative mismatch; they look more like read-level artifacts, source-label/crosswalk edge cases, or rare hash/context coincidences because they are not present in the helveticus source whole-genome scan.

## Caveats

- This is a five-reference diagnostic on one sample, not a full benchmark.
- The trace follows MinCO's selected readwise assignments after the active best-diff-split rule. It does not enumerate every possible read/reference match.
- Same GTDB species but different source strain is not counted as contamination here, because the profiled row is GTDB-representative/species-centric.
- `trace events` are read-level selected context assignments. `unique qctx` is the number of distinct context hashes observed in the trace. These should not be confused with every possible reference marker.
