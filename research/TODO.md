# TODO

## Patent

- Draft a Chinese patent disclosure before broad public disclosure of the
  method. Keep the repo and detailed benchmark notes private until filing
  strategy is decided.
- Invention focus:
  fixed-size context-minhash readwise profiling, reference-density extraction,
  unique best object-difference assignment for shared contexts, ZIP/effective
  coverage context AAF ANI, and context-coverage abundance estimation.
- Prepare materials for a patent professional:
  algorithm description, flow diagram, prior-art comparison, implementation
  evidence, and CAMI marine benchmark results.

## Immediate Benchmarking

- Continue testing the current `S=1000` readwise recipe on independent datasets
  beyond the completed CAMI marine + toy gut + plant + strain-madness panel:

```bash
minco ani -p16 -r REF_S1000 --qraw reads.fq.gz \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-unique \
  --readwise-ani zip-aaf \
  -m0 -f0.05 -n0.94 -t10 \
  -o out.tsv
```

- Preferred next datasets:
  Toy Mouse Gut, Toy Human Microbiome Project, pathogen-rich samples,
  viral-spiked samples, or low-complexity mock
  communities with known truth.
- Compare against Sylph, Kraken2/Bracken when applicable, and existing
  MetaKSSD/KSSD baselines when the database mapping is fair.
- Record speed, peak RSS, species TP/FP/FN/F1, abundance correlation, and
  failure cases.

## Method Improvements

- Validate the adaptive edge-EM trigger on a larger multi-domain panel before
  promoting it. Current pilot rule:
  `marker_raw_targets>=150 -> beta=0.02`,
  `marker_raw_targets<=25 -> beta=0.0015`, otherwise `beta=0`. It improves
  marine/strain L1 versus S2000 edge-marker beta 0 and avoids CAMI3 regression,
  but it has only been tested as a three-sample spot check and does not improve
  F1.
- Add optional `Ref_depth_cv` reporting filter if it improves across datasets,
  not only sample0.
- Validate whether `S=1000` remains the best practical default on non-marine
  datasets; otherwise derive S-dependent defaults.
- Avoid dataset-label features in call models. The third-domain plant extension
  showed that dataset one-hot features inflated the six-sample result; use only
  deployable minco evidence such as breadth, depth CV, XnY support, ANI,
  abundance, and strict/relaxed rule evidence.
- Do not install the current unique/split hybrid threshold rule as a default.
  The 2026-06-21 two-sample experiment showed useful split-mode recall, but the
  joint hybrid rule did not beat Sylph or the best unique baseline. Next step is
  multi-sample calibration, not more same-sample threshold tuning.
- Do not export the current RF/HGB ensemble to pure C yet. External
  strain-madness holdout testing showed it does not generalize beyond Sylph:
  mean F1 0.466 versus 0.527 for local Sylph and 0.463 for minco split-direct.
- Prefer a simpler deployable filter first: direct split/unique evidence plus a
  small single model or distilled rule. Candidate features: unique and split
  XnY, breadth, depth, depth CV, ZIP AAF ANI, naive ANI, unique-vs-split deltas,
  estimated abundance, and strict/relaxed rule evidence.
- Keep CAMI strain-madness samples 0-2 as a locked external holdout while
  improving filters on other training/development data. Do not tune thresholds
  directly on strain-madness unless the result is clearly marked diagnostic.
- Evaluate Toy Mouse Gut only after building a reliable species-level gold
  profile from the CAMISIM setup. Its download directory exposes distributions
  and source genomes, not ready `taxonomic_profile_*.txt` files, so a fair
  minco/Sylph comparison needs source-genome-to-taxid/profile reconstruction
  before downloading many 5.1 GB read tarballs.
- Keep `u_or_s_direct_validtax` as a candidate optional reporting filter. It
  improved external strain F1 from 0.477 to 0.497 and FP+FN from 17.3 to 16.0,
  but lowered the nine-sample training F1 from 0.620 to 0.603.
- Support domain-specific reporting or scoring views for mixed references
  because viral and bacterial taxa showed different recall/precision behavior.
- Improve reference index memory layout so large references are closer to
  Sylph memory while preserving readwise speed.
- Improve plus-virus reference lookup speed; CAMI III Toy Human Gut sample0 was
  slower than CAMI marine S1000 despite similar read counts.
- Test plus-virus S1000 database separately from GTDB-only S1000 to measure
  viral sensitivity and FP behavior.
- Improve CAMI/NCBI/GTDB lineage mapping and generate official OPAL profiles
  for public benchmark comparison after patent filing strategy is settled.

## Release Gate

- Do not present the algorithm as a general default until it passes at least one
  substantially different dataset.
- After non-marine validation, update README defaults and command-line help to
  recommend `S=1000` for readwise profiling.
