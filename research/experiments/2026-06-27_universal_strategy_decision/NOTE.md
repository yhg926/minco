# Universal MinCO Profile Strategy Decision

Date: 2026-06-27
Author/agent: Codex
Project: MinCO / KSSD3mini

## Code Provenance

- Code repository: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Commit: `b10e1104a9bfd02eb08f58d0c0c656de5e0a9cfc`
- Ref/describe: `b10e110-dirty`
- Working tree: dirty; this decision note summarizes uncommitted wrapper,
  documentation, and benchmark-note work. Use `commands.sh` to reproduce the
  status checks.
- Binary/script/tool: `scripts/minco_profile_calibrated.py` and
  `scripts/minco_profile_default.py`
- Provenance files: `build_strategy_decision_summary.py`, `summary.tsv`,
  `results/decision_checks.tsv`, `results/objective_audit.tsv`,
  `parameters.tsv`, `artifacts.md`, `commands.sh`

## Command and Parameter Provenance

- Working directory: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Exact commands: `commands.sh`
- Parameters, defaults, thresholds, and source result files: `parameters.tsv`

## Question

Given the current multi-panel benchmark evidence, what is the most stable
single MinCO profiling strategy to expose as the default, without requiring
users to choose a dataset-specific mode?

## Hypothesis

A guarded F1-priority `universal-auto-exact` wrapper is a better default than
legacy probability-only gates because it improves detection on the mixed
readiness panel while avoiding the unsafe broad split/near-split rescues that
increase false positives. It should be described as MinCO's best current
default, not as a universal Sylph-beating strategy.

## Dataset

- Input data: existing benchmark summaries from CAMI II/III marine, human gut,
  mouse gut, plant-associated, strainmadness, and diagnostic ANI experiments.
- Source path or URI: see `artifacts.md`.
- Local path or soft link: source result files under `research/experiments/`.
- Availability status: available for small repo-tracked notes and summaries;
  large read/sketch artifacts are external or temporary and referenced by the
  source notes.
- Last verified: 2026-06-28.
- Retention/deletion risk: low for repo note files; external `/tmp` and
  `/mnt/new3T` data are not guaranteed by this note.
- Sample count: source panels include 26-sample mixed readiness, CAMI III Toy
  Human Gut samples0-5, CAMI II Toy Mouse Gut samples5-7, and earlier
  dataset-specific panels.
- Selection criteria: prioritize F1, then abundance L1/Pearson, then
  ANI/reporting robustness, with one automatic strategy for users.

## Methods

I reviewed the current wrapper implementation and the existing benchmark
summaries, then added release-grade raw-read HMP airskin
sample0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28
profiling and corrected aggregate HMP rescoring after fixing stale
sample6/sample11 scorer defaults. The decision is a stage summary over
existing evidence plus lightweight rescoring/decomposition from cached
truth/profile tables and the new same-release HMP raw-read runs.

`score_cami3_gtdb_taxid_transfer.py` adds a conservative CAMI III Toy Human
Gut GTDB-species transfer view: a CAMI/NCBI species taxid is converted to GTDB
only when GTDB metadata maps that taxid to exactly one GTDB species. Ambiguous
or unmapped truth mass is reported and excluded from the transferred truth set.
This gives a comparable MinCO/Sylph GTDB-species diagnostic, but it is not as
strong as source-genome accession truth.

`score_cami3_gtdb_source_readmap.py` adds a stronger CAMI III samples0-2
source-readmap transfer view. It maps source read rows to GTDB species through
unique GTDB WGS-prefix matches from source contigs, then falls back to unique
taxid transfer. The scorer also recovers MinCO best-reference accessions from
raw unique/split tables so MinCO is not penalized for taxid-only calibrated
profile output.

`score_hmp_gtdb_source_abundance.py` adds a CAMI II HMP airskin source-genome
abundance transfer view, initially for samples6 and 11 and now expanded through
samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28. It maps CAMI source genome
accessions to GTDB r232 species by exact accession or unique assembly core, then
scores MinCO and Sylph in the same GTDB species namespace. The scorer now uses
`--sylph-source auto` by default: it consumes r232 Sylph profiles from
`/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232` only when both sample6
and sample11 profiles exist, otherwise it falls back to cached r226 diagnostic
profiles. It also accepts dynamic airskin sample IDs using standard current-run
paths under `/tmp/minco_current_code_hmp_airskin${sample}_20260627` and
`/tmp/cami2_hmp_unseen_transfer_20260626/run_sylph_r232`, while still honoring
explicit per-sample path overrides for alternate-root runs. The mapped truth
abundance is high across the scored HMP airskin samples; sample11 is the main
imperfect case at 96.9321%, while sample18 and sample21 map 100.0%. The r232
chunked Sylph DB and scored r232 HMP profiles now exist, so the current HMP
rows are same-release GTDB r232 evidence.

`score_hmp_gastrooral_gtdb_source_abundance.py` adds the matching CAMI II HMP
gastrooral source-genome abundance transfer view for samples0 and 6. It uses
the current universal table-mode MinCO outputs, existing unique/split raw
tables for best-reference accession recovery, and new chunked GTDB r232 Sylph
profiles generated from the cached gastrooral `.sylsp` sketches. The mapped
truth abundance is high, 99.1736% and 99.2685%, so this is now a release-grade
GTDB source-abundance counterexample panel.

`cross_validate_cami3_loose_split_rescue.py` tests the most promising
CAMI3 source-readmap high-depth split rescue against the broader 26-sample
readiness panel. It applies the rule without truth labels and is used only to
decide whether the rescue is promotable.

The fixed-call abundance cross-check retests the CAMI3 source-readmap-favored
`s_mean / zip_af^0.25` formula on the broader 32-sample fixed-call abundance
panel spanning Toy Mouse, HMP airskin/gastrooral, and CAMI3 source-readmap
truth views. Calls are unchanged, so F1
is constant and only L1/Pearson are evaluated.

`decompose_abundance_errors.py` breaks MinCO/Sylph abundance L1 into
matched-true-positive absolute error, missing truth mass, and extra predicted
mass for Toy Mouse samples5-7, same-release HMP airskin samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28, HMP
gastrooral samples0/6, and CAMI3 source-readmap samples0-2. It reuses the
official GTDB truth/profile scoring helpers and validates against the official
score TSVs before drawing conclusions.

`audit_abundance_oracle_bounds.py` adds a fixed-call abundance upper-bound
diagnostic on the same 32 GTDB panel samples. It keeps the current MinCO call
set for F1 accounting, then assigns true abundance to the detected true
species under two bounds: unnormalized detected truth mass and truth mass
renormalized over detected species. This asks whether perfect abundance
allocation alone could close the Sylph gap without changing calls.

`audit_missed_truth_candidates.py` checks the next question raised by the
oracle: for truth species that current MinCO misses, is there any mapped
uncalled row in the emitted profile that could be rescued by a safer gate? It
uses the same GTDB mapping helpers as the official scorers. The audit is
profile-surface limited: many HMP outputs emit only final rows, so absence from
this audit does not prove absence from raw unique/split evidence.

`audit_hmp_missed_truth_raw_tables.py` follows up the HMP profile-surface
limitation. For high-abundance HMP truth species absent from emitted profile
rows, it loads the corresponding raw unique/split/exact tables, maps reference
accessions to GTDB species with the same scorer helpers, and checks whether the
missed species exists before final calibrated emission.

`sweep_cross_panel_abundance_variants.py` holds calls fixed and sweeps
candidate abundance formulas on those four GTDB panels. It preserves the
panel-specific official collapse rules: Toy Mouse uses max abundance per GTDB
species and truth-only L1, while HMP/CAMI3 source-transfer panels use summed
species abundance and union L1. The current abundance row is validated against
`abundance_error_decomposition_summary.tsv` before interpreting variants.
Genus-local reallocation candidates use the profile's own species names, not
external scorer-only GTDB remappings, so the sweep reflects what the wrapper can
actually apply at output time.

`summarize_abundance_variant_safety.py` reduces that sweep into a promotion
audit. It excludes trivial current-abundance equivalents, counts tested
nontrivial variants, and records whether any variant improves all evaluated
samples or all panels. This is stricter than comparing only mean L1.

`audit_adaptive_abundance_switch.py` tests whether a simple automatic switch
can choose between current abundance and top fixed-call variants using only
MinCO output-derived sample features. Truth is used only for scoring candidate
rules. The audit reports in-sample safety counts and a leave-one-panel-out
selection test to separate real generalization from panel overfit.

`sweep_cross_panel_call_filters.py` tests whether simple output-only call
filters can prune false positives without dataset-specific choices. It keeps
the same four GTDB panels and official collapse/scoring rules, starts from
current calibrated calls, and uses only profile columns such as probability,
XnY, breadth, real align fraction, ZIP AF, depth, and ANI. This is post-hoc
diagnostic evidence; a filter is promotable only if F1 improves without
sample/panel tradeoffs.
The same script now has an opt-in `--include-hmp-omitted` mode that adds the
restored HMP omitted samples2/8/26 and writes separate
`cross_panel_call_filter_with_hmp_omitted_*` result files, leaving the original
32-sample audit outputs unchanged.

`audit_adaptive_call_filter_switch.py` tests whether a simple automatic switch
can choose among top output-only call filters using only MinCO output-derived
sample features. Candidate filters are selected from the fixed call-filter
sweep, truth is used only for scoring rules, and leave-one-panel-out validation
is used to separate an in-panel signal from a default-ready rule.
The same script now has an opt-in `--include-hmp-omitted` mode that repeats
the adaptive switch audit after adding restored HMP samples2/8/26, using the
separate `cross_panel_call_filter_with_hmp_omitted_*` sweep outputs.

`--adaptive-call-filter-switch lopo-min-xny25` was added only as an explicit
experimental wrapper option. Its default is `off`, so the selected default call
gate is unchanged. `run_adaptive_call_filter_wrapper_validation.sh` and
`score_adaptive_call_filter_wrapper_validation.py` replay the wrapper on 10
cached-table profiles, then check that the wrapper-emitted call mask matches
the offline adaptive switch audit. Abundance parity versus the offline post-hoc
tables is recorded as informational because this option's implementation
contract is the call mask, not an abundance formula replacement.

`audit_adaptive_call_filter_external_exactsplit.py` then stress-tests the same
`lopo-min-xny25` rule on cached plant-associated and strainmadness exact-split
profiles outside the four GTDB panels that selected the switch. The test uses
the same local bacteria-scope scorers as the existing exact-split summaries and
validates the current-call baseline against the cached score TSVs before
interpreting the filtered result.

`audit_supervised_abundance_calibrator.py` tests row-level supervised abundance
models trained only from numeric MinCO output columns, with no species, ref,
sample, or panel identifiers as features. It evaluates five model specs across
four current-output blend values by leave-one-panel-out. Truth is used for
training targets and held-out scoring, so this is a candidate-discovery audit,
not a default promotion by itself.

`audit_supervised_abundance_sample28_holdout.py` tests the predefined RF
abundance candidate more strictly by excluding HMP airskin sample28 from
training, then scoring sample28 as a post-hoc held-out sample. This evaluates
whether the candidate survives the newly restored same-release sample, but it is
still a same-panel check rather than a new independent dataset.

`audit_supervised_abundance_external_exactsplit.py` tests the same predefined
RF abundance model with blend strengths 0.25, 0.5, and 0.75 on cached
plant-associated and strainmadness exact-split profiles. These panels are
outside the four-panel supervised abundance training set, but they are still
diagnostic/non-release because they use cached `/tmp` profiles, local scorers,
and mixed reference-release baselines.

`--abundance-genus-xny-blend-alpha` was added only as an explicit experimental
wrapper option. Its default is 0, so the selected default abundance is
unchanged. `run_abundance_blend_wrapper_validation.sh` and
`score_abundance_blend_wrapper_validation.py` validate the alpha 0.25
genus-XnY blend as a wrapper-realistic experiment; it is not promoted because
it still has panel tradeoffs.

For ANI reporting, the calibrated wrapper now exposes `reported_ani` and
`reported_ani_source` near the front of the output. `reported_ani` uses
ZIP-AAF ANI (`Ref_zip_aaf_ani`) rather than the raw/emitted ANI field, because
cached source/ref ANI diagnostics showed that raw/emitted readwise ANI can
saturate in metagenome profiles.

`build_strategy_decision_summary.py` rebuilds `summary.tsv`,
`results/decision_checks.tsv`, and `results/objective_audit.tsv` from source
benchmark TSVs. This keeps the stage decision tied to measurable source rows
rather than chat memory.

`build_release_readiness.py` applies stricter release gates over the generated
decision checks and holdout-bundle audit. It separates "best current MinCO
default candidate" from "release-grade broad universal claim".

`audit_holdout_resources.py` records local resources needed to promote
diagnostic panels to release grade. For the HMP airskin source-abundance
panel, it confirms cached MinCO inputs, the chunked GTDB r232 Sylph DB list,
and sample6/sample11 r232 Sylph profiles are present before marking HMP
same-release scoring ready.

`plan_next_release_grade_actions.py` combines the HMP airskin candidate audit,
release-readiness gates, abundance-error decomposition, and runtime-readiness
gates to pick the next concrete evidence action. After completing sample8, it
now selects HMP airskin sample0 as the next release-grade candidate, but marks
the run as blocked by scratch-space headroom until more space is freed.
The same planner writes the tracked `NEXT_RELEASE_GRADE_ACTIONS.md` recall
file so the next actions remain visible without committing ignored result TSVs.

`audit_tmp_headroom_for_next_action.py` checks whether the next release-grade
run has enough `/tmp` space before starting any network-backed restore/profile
command. It estimates the required footprint from observed HMP FASTQ, Sylph
sketch, and MinCO workdir sizes, adds a 2 GiB safety buffer, and writes a
manual-review cleanup candidate list. It is non-destructive.

`audit_hmp_airskin_sample8_extra_calls.py` is a focused post-run diagnostic
for the sample8 HMP holdout. It uses the same source-accession GTDB mapping
policy as `score_hmp_gtdb_source_abundance.py`, labels each called MinCO row
as TP/FP in the mapped truth namespace, records the single FN, and evaluates a
small set of already-studied output-only filters on sample8. It is diagnostic
only and does not change the default strategy.

`preflight_gtdb232_sylph_db.py` checks GTDB r232 Sylph DB resources. It
verifies the representative path list, missing FASTA count, compressed input
size, output free space, chunked DB completion, and exact monolithic/chunked
commands. Because this machine has 61 GiB RAM and no swap, the completed r232
database is the chunked 20-part build rather than a single 199,923-genome
monolithic sketch.

The user-facing command boundary was tightened after the HMP source-abundance
diagnostic: `scripts/minco_profile_default.py --help`, README, and the user
manual now put the calibrated species default first and explicitly reserve the
C `minco profile` subcommand for conservative direct species/AMR/virus/gene or
mixed-domain profiling.

The default launcher now also supports packaged species databases: when
`--taxmap`, `--train-features`, and their environment-variable equivalents are
absent, it discovers `species_taxmap.tsv` and `joined_feature_training/`
sidecars beside `--ref`. This lets a user run the selected default strategy
without choosing a strategy or wiring calibration paths by hand.

The launcher now also supports fitted RF/HGB model caches. `--model-cache`,
`MINCO_PROFILE_MODEL_CACHE`, or a compatible
`minco_profile_rf_hgb*.joblib` sidecar beside the ref/training directory lets
the calibrated wrapper skip loading training tables and refitting RF/HGB for
each sample. On CAMI II Toy Mouse sample6, precomputed-table scoring dropped
from 39.97 s uncached to 15.85 s cached.

For speed rather than F1-priority default behavior, I tested cached
`--strategy universal` with normal calls-only output. This skips the
`universal-auto-exact` split rerun. On CAMI II Toy Mouse sample6 it was faster
than Sylph sketch+profile, but Sylph retained better F1 and abundance accuracy
on the same read file. Therefore this is documented as a speed-priority mode,
not as a replacement for the current F1-priority default.

I then retested the exact-skip choice on CAMI II Toy Mouse samples5-7 using
current-code table mode: same unique table, either normal block split or exact
split, cached RF/HGB model, and `--strategy universal`. Exact split rescued one
true low-abundance sample5 species and changed no calls in samples6-7. The
aggregate difference is small but in the requested priority order it favors
keeping exact available for the F1-priority default: exact mean F1 0.893839
versus block mean F1 0.892792. Block-only has slightly better L1/Pearson, so it
remains the speed-priority mode.

Parameter summary:

```text
strategy = universal-auto-exact
threshold = 0.35
exact split trigger = probability_extra_mass_ratio <= 0.10
exact split guard = joined_base_median_uaf < 0.35 OR raw_unique_to_base_ratio >= 0.80
exact split low-extra mode = skip exact rerun when block low-extra split rescue added candidates
high-extra tail rescue = joined_base_median_uaf < 0.35 and calibrated_probability >= 0.25
low-extra split rescue = base mode, probability_extra_mass_ratio < 0.005,
                         joined_base_median_uaf < 0.45,
                         top 1/genus with P >= 0.20 and strict split evidence
strict split evidence = s_XnY_ctx >= 50, s_ANI >= 0.95,
                        s_Real_min_align_fraction >= 0.05,
                        s_Ref_breadth >= 0.05
rejected loose CAMI3 rescue = P >= 0.02, s_XnY >= 300, s_ANI >= 0.93,
                              s_Real_min_align_fraction >= 0.30,
                              s_Ref_breadth >= 0.20, s_mean_depth >= 1,
                              top 1/genus
abundance default = max(split mean depth / split zip AF,
                        unique mean depth / unique zip AF) for normal rows;
                    split mean depth for tail-rescue-added rows
rejected abundance alternative = split mean depth / zip AF^0.25
reported ANI = split Ref_zip_aaf_ani by default;
               unique Ref_zip_aaf_ani when raw-unique fallback is active
```

Commands are recorded in `commands.sh`; detailed parameters are recorded in
`parameters.tsv`.

## Results

Key metrics are recorded in `summary.tsv`; machine-readable decision checks are
recorded in `results/decision_checks.tsv`; requirement-level status is recorded
in `results/objective_audit.tsv` and the tracked
`OBJECTIVE_COMPLETION_AUDIT.md`; release-readiness gates are recorded in
`results/release_readiness.tsv`.

The current best stable MinCO default candidate is:

```text
universal-auto-exact
+ guarded exact split rerun
+ guarded high-extra low-uAF probability-tail rescue
+ guarded low-extra split rescue
```

Important measured evidence:

| panel | method | samples | F1 | L1 | Pearson | interpretation |
|---|---|---:|---:|---:|---:|---|
| mixed readiness | MinCO current default | 26 | 0.709348 | 0.623300 | 0.939155 | best documented MinCO F1-priority default over legacy gates |
| mixed readiness | legacy probability gate | 26 | 0.696214 | 0.615277 | 0.938123 | slightly lower F1 but slightly better L1 |
| CAMI III Toy Human Gut | MinCO current default | 6 | 0.764745 | 0.654447 | 0.960564 | beats Sylph on F1/FP+FN, not abundance |
| CAMI III Toy Human Gut | Sylph | 6 | 0.684913 | 0.471147 | 0.996573 | stronger abundance baseline |
| CAMI III GTDB taxid transfer | MinCO current default | 6 | 0.701329 | 74.674523 pp | 0.685932 | beats Sylph in the partial unique-taxid transfer view |
| CAMI III GTDB taxid transfer | Sylph | 6 | 0.404228 | 144.731757 pp | 0.427660 | many FPs after GTDB transfer in this diagnostic view |
| CAMI III GTDB source-readmap | MinCO current default | 3 | 0.815957 | 62.732029 pp | 0.823915 | beats Sylph on F1 after accession-aware MinCO transfer |
| CAMI III GTDB source-readmap | Sylph | 3 | 0.741521 | 32.632841 pp | 0.945708 | stronger abundance L1/Pearson |
| CAMI II Toy Mouse Gut | MinCO current-code default refresh | 3 | 0.893839 | 13.767837 pp | 0.987977 | refreshed with current wrapper schema; improves MinCO versus previous mouse default |
| CAMI II Toy Mouse Gut | Sylph | 3 | 0.952806 | 6.731523 pp | 0.994040 | clearly stronger on this panel |
| CAMI II Toy Mouse exact-skip audit | MinCO cached exact split `universal` | 3 | 0.893839 | 13.767837 pp | 0.987977 | one extra TP versus block-only; supports exact for F1-priority default |
| CAMI II Toy Mouse exact-skip audit | MinCO cached block split `universal` | 3 | 0.892792 | 13.619065 pp | 0.988226 | faster mode; slightly lower F1, slightly better abundance |
| CAMI II Toy Mouse sample6 speed | MinCO current `universal-auto-exact` low-extra skip | 1 | 0.854749 | 19.831577 pp | 0.974331 | 1:47.90, 3.43 GiB RSS; faster than Sylph here but lower accuracy |
| CAMI II Toy Mouse sample6 speed | MinCO exact sidecar allow | 1 | 0.854749 | 20.062956 pp | 0.973755 | 2:28.51, 4.50 GiB RSS; faster than old exact rerun but slower than Sylph |
| CAMI II Toy Mouse sample6 speed | Sylph sketch+profile | 1 | 0.939086 | 11.759075 pp | 0.984876 | 1:50.98, 18.81 GiB RSS; higher F1/abundance accuracy |
| CAMI II marine species-taxid | MinCO S1000 unique ZIP-AAF domain recipe | 6 | 0.847097 | MAE 0.002482 | 0.983377 | useful MinCO domain evidence, but not the current universal default |
| CAMI II marine species-taxid | Sylph | 6 | 0.829047 | MAE 0.001977 | 0.992436 | better abundance correlation/MAE |
| CAMI II marine GTDB taxid transfer | MinCO S1000 unique ZIP-AAF | 0,3,4,5 | 0.730691 | 50.894507 pp | 0.941545 | beats Sylph in a partial GTDB transfer namespace, but not current universal-default evidence |
| CAMI II marine GTDB taxid transfer | Sylph r226 | 0,3,4,5 | 0.703504 | 73.864971 pp | 0.894244 | lower F1 and union-L1 in this partial transfer view |
| CAMI II marine GTDB transfer quality | truth audit | 0-9 | NA | scored B/A mapped 71.3731-77.0927% | NA | not release-grade; only about 74.4% of scored Bacteria/Archaea species mass maps uniquely to GTDB |
| CAMI II plant bacteria-scope | MinCO RF/HGB P>=0.35 train12 | 3 | 0.615015 | 1.050342 | 0.733767 | slightly beats Sylph, but local bacteria-scope only |
| CAMI II plant bacteria-scope | Sylph | 3 | 0.599993 | 1.070001 | 0.611248 | close external baseline |
| CAMI II plant bacteria-scope | C `minco profile` direct | 3 | 0.175519 | 1.236095 | 0.997935 | too conservative for calibrated species default; documented as direct/domain path |
| CAMI II plant current exact-split | MinCO current universal exact-split | 3 | 0.612016 | 0.822986 | 0.922354 | competitive current-gate evidence, but non-release because it uses cached `/tmp` exact-split outputs and the plant bacteria-scope scorer |
| CAMI II plant current exact-split | Sylph | 3 | 0.599993 | 1.070001 | 0.611248 | close external baseline in the same non-release scoring namespace |
| CAMI II plant GTDB transfer feasibility | source-abundance truth audit | 3 | NA | mapped 19.0835-22.2414% | NA | not release-grade; most abundance is RNODE/no-taxid or genus/family-only CAMI source truth |
| CAMI II strainmadness current exact-split | MinCO current universal exact-split | 3 | 0.628737 | 0.759779 | 0.823822 | higher F1 than Sylph, but cached/train-leaky non-release evidence and weaker abundance |
| CAMI II strainmadness current exact-split | Sylph | 3 | 0.527374 | 0.203000 | 0.994006 | lower F1 but much stronger abundance |
| CAMI II strainmadness GTDB transfer feasibility | source-abundance truth audit | 3 | NA | mapped 0.4775-2.9454% | NA | not release-grade; most abundance is new_strain/new_genus/new_order truth that does not map uniquely to GTDB species |
| CAMI II HMP gastrooral pilot | MinCO calibrated train12 | 2 | 0.779956 | 0.821302 | 0.820224 | counterexample for always-calibrated gate |
| CAMI II HMP gastrooral pilot | MinCO unique direct ANI>=0.95 | 2 | 0.832996 | 0.986546 | 0.841353 | best MinCO signal here, still below Sylph |
| CAMI II HMP gastrooral pilot | Sylph | 2 | 0.901587 | 0.311628 | 0.996839 | clear winner here |
| CAMI II HMP gastrooral current check | MinCO current universal table-mode | 2 | 0.861646 | 0.766191 | 0.867024 | current universal fallback improves over old calibrated and unique-direct MinCO |
| CAMI II HMP gastrooral current check | Sylph | 2 | 0.901587 | 0.311628 | 0.996839 | still stronger on this panel |
| CAMI II HMP gastrooral source-abundance | MinCO current universal table-mode r232 | 0,6 | 0.919568 | 68.950069 pp | 0.858401 | release-grade GTDB source-abundance counterexample; superseded for current-default runtime by raw rerun below |
| CAMI II HMP gastrooral raw default source-abundance | MinCO current raw default r232 | 0,6 | 0.919568 | 48.986991 pp | 0.882277 | raw-wrapper rerun improves MinCO L1 by 19.963078 pp, but Sylph remains much better |
| CAMI II HMP gastrooral source-abundance | Sylph r232 chunked | 0,6 | 0.955300 | 5.063129 pp | 0.999125 | stronger F1 and abundance |
| CAMI II HMP gastrooral sample0 speed | MinCO exact sidecar | 0 | panel F1 0.919568 | panel L1 48.986991 pp | panel Pearson 0.882277 | 2:25.93, 4.63 GiB RSS; same calls/abundance as exact rerun and faster than Sylph sample0 |
| CAMI II HMP gastrooral sample0 speed | MinCO exact rerun | 0 | panel F1 0.919568 | panel L1 48.986991 pp | panel Pearson 0.882277 | 4:00.22, 3.80 GiB RSS; current default exact-needed path |
| CAMI II HMP gastrooral sample0 speed | Sylph r232 chunked | 0 | panel F1 0.955300 | panel L1 5.063129 pp | panel Pearson 0.999125 | 3:27.18, 25.92 GiB RSS; more accurate but slower than sidecar here |
| exact-sidecar policy audit | six current instances | 6 | exact candidates 3 | avoid-exact rows 3 | NA | reject unconditional sidecar default; needs preflight or cheaper conditional sidecar |
| prefix exact-preflight audit | HMP gastrooral first50k/first200k | 4 preflights | raw mismatches 2 | support-aware decisive 0 | NA | reject prefix preflight default; small prefixes are too sparse |
| candidate-restricted exact audit | HMP sample0 and Toy Mouse sample6 | 2 | final-taxid exact keep 0.140-0.154% | post-hoc exact-row removal 99.846-99.860% | NA | row filtering after full exact scan is large, but candidate-only exact rerun is rejected because aggregate tables lack per-context competitor closure |
| CAMI II HMP airskin source-abundance | MinCO current-code refresh exact6 | 6,11 | 0.961058 | 22.780570 pp | 0.925614 | same-release r232 counterexample for abundance after restoring sample6 FASTQ and rerunning exact split |
| CAMI II HMP airskin source-abundance | Sylph r232 chunked | 6,11 | 0.970588 | 5.133659 pp | 0.996572 | clear abundance winner and slightly higher F1 |
| CAMI II HMP airskin sample28 source-abundance | MinCO current default | 28 | 0.888889 | 18.352634 pp | 0.989034 | new same-release r232 holdout; all 12 true species detected but 3 FPs |
| CAMI II HMP airskin sample28 source-abundance | Sylph r232 chunked | 28 | 0.923077 | 1.803063 pp | 0.999933 | stronger F1 and abundance |
| CAMI II HMP airskin sample22 source-abundance | MinCO current default | 22 | 0.644068 | 20.648591 pp | 0.988592 | new same-release r232 counterexample; 19 TPs, 19 FPs, 2 FNs |
| CAMI II HMP airskin sample22 source-abundance | Sylph r232 chunked | 22 | 0.909091 | 6.043692 pp | 0.997079 | stronger F1 and abundance |
| CAMI II HMP airskin sample5 source-abundance | MinCO current default | 5 | 0.961538 | 10.337736 pp | 0.998382 | same-release r232 holdout; MinCO has higher F1 but worse abundance |
| CAMI II HMP airskin sample5 source-abundance | Sylph r232 chunked | 5 | 0.943396 | 5.965327 pp | 0.998278 | lower F1 than MinCO here, stronger abundance |
| CAMI II HMP airskin sample0 source-abundance | MinCO current default | 0 | 0.840580 | 18.780365 pp | 0.989504 | same-release r232 counterexample; 29 TPs, 7 FPs, 4 FNs |
| CAMI II HMP airskin sample0 source-abundance | Sylph r232 chunked | 0 | 0.903226 | 2.510444 pp | 0.999779 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code exact6 plus sample28 | 6,11,28 | 0.937001 | 21.304591 pp | 0.946754 | expanded same-release panel confirms MinCO abundance gap |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 6,11,28 | 0.954751 | 4.023461 pp | 0.997693 | stronger F1 and much stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples6,11,22,28 | 6,11,22,28 | 0.863768 | 21.140591 pp | 0.957213 | sample22 expands the release-grade counterexample; F1 and L1 both worse than Sylph |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 6,11,22,28 | 0.943336 | 4.528518 pp | 0.997539 | stronger F1 and much stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples5,6,11,22,28 | 5,6,11,22,28 | 0.883322 | 18.980020 pp | 0.965447 | sample5 improves MinCO single-sample F1 but the five-sample panel still favors Sylph |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 5,6,11,22,28 | 0.943348 | 4.815880 pp | 0.997687 | stronger five-sample F1 and much stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,5,6,11,22,28 | 0,5,6,11,22,28 | 0.876198 | 18.946744 pp | 0.969457 | sample0 strengthens the same-release HMP counterexample |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,5,6,11,22,28 | 0.936661 | 4.431641 pp | 0.998036 | stronger six-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample21 source-abundance | MinCO current default | 21 | 0.880000 | 34.744967 pp | 0.962690 | same-release r232 counterexample restored under `/mnt/new3T`; 33 TPs, 8 FPs, 1 FN |
| CAMI II HMP airskin sample21 source-abundance | Sylph r232 chunked | 21 | 0.970588 | 2.197983 pp | 0.999889 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,5,6,11,21,22,28 | 0,5,6,11,21,22,28 | 0.876741 | 21.203633 pp | 0.968490 | sample21 further strengthens the same-release HMP counterexample |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,5,6,11,21,22,28 | 0.941508 | 4.112547 pp | 0.998300 | stronger seven-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample18 source-abundance | MinCO current default | 18 | 0.943396 | 12.023035 pp | 0.994733 | same-release r232 counterexample restored under `/mnt/new3T`; 25 TPs, 1 FP, 2 FNs |
| CAMI II HMP airskin sample18 source-abundance | Sylph r232 chunked | 18 | 1.000000 | 1.882765 pp | 0.999930 | perfect F1 here and stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,5,6,11,18,21,22,28 | 0,5,6,11,18,21,22,28 | 0.885073 | 20.056058 pp | 0.971770 | corrected eight-sample panel after sample18; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,5,6,11,18,21,22,28 | 0.948819 | 3.833824 pp | 0.998504 | stronger eight-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample13 source-abundance | MinCO current default | 13 | 0.925000 | 54.747984 pp | 0.910524 | same-release r232 counterexample restored under `/mnt/new3T`; 37 TPs, 3 FPs, 3 FNs |
| CAMI II HMP airskin sample13 source-abundance | Sylph r232 chunked | 13 | 0.987654 | 8.344359 pp | 0.997403 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,5,6,11,13,18,21,22,28 | 0,5,6,11,13,18,21,22,28 | 0.889510 | 23.910717 pp | 0.964965 | corrected nine-sample panel after sample13; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,5,6,11,13,18,21,22,28 | 0.953134 | 4.334995 pp | 0.998382 | stronger nine-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample25 source-abundance | MinCO current default | 25 | 0.720000 | 23.876830 pp | 0.984476 | same-release r232 counterexample restored under `/mnt/new3T`; 45 TPs, 31 FPs, 4 FNs |
| CAMI II HMP airskin sample25 source-abundance | Sylph r232 chunked | 25 | 0.969072 | 4.836648 pp | 0.998867 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,5,6,11,13,18,21,22,25,28 | 0,5,6,11,13,18,21,22,25,28 | 0.872559 | 23.907328 pp | 0.966916 | corrected ten-sample panel after sample25; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,5,6,11,13,18,21,22,25,28 | 0.954728 | 4.385160 pp | 0.998430 | stronger ten-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample3 source-abundance | MinCO current default | 3 | 0.957447 | 12.670667 pp | 0.997025 | same-release r232 counterexample restored under `/mnt/new3T`; 45 TPs, 2 FPs, 2 FNs |
| CAMI II HMP airskin sample3 source-abundance | Sylph r232 chunked | 3 | 0.989247 | 2.961929 pp | 0.999926 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,3,5,6,11,13,18,21,22,25,28 | 0,3,5,6,11,13,18,21,22,25,28 | 0.880276 | 22.885814 pp | 0.969653 | corrected eleven-sample panel after sample3; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,3,5,6,11,13,18,21,22,25,28 | 0.957866 | 4.255775 pp | 0.998566 | stronger eleven-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample1 source-abundance | MinCO current default | 1 | 0.680000 | 18.105884 pp | 0.997394 | same-release r232 counterexample restored under `/mnt/new3T`; 51 TPs, 48 FPs, 0 FNs |
| CAMI II HMP airskin sample1 source-abundance | Sylph r232 chunked | 1 | 0.990291 | 2.178406 pp | 0.999907 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,5,6,11,13,18,21,22,25,28 | 0,1,3,5,6,11,13,18,21,22,25,28 | 0.863586 | 22.487486 pp | 0.971965 | corrected twelve-sample panel after sample1; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,5,6,11,13,18,21,22,25,28 | 0.960568 | 4.082661 pp | 0.998678 | stronger twelve-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample17 source-abundance | MinCO current default | 17 | 0.773109 | 26.525181 pp | 0.996010 | same-release r232 counterexample restored under `/mnt/new3T`; 46 TPs, 27 FPs, 0 FNs |
| CAMI II HMP airskin sample17 source-abundance | Sylph r232 chunked | 17 | 0.978723 | 1.872577 pp | 0.999967 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,5,6,11,13,17,18,21,22,25,28 | 0,1,3,5,6,11,13,17,18,21,22,25,28 | 0.856626 | 22.798078 pp | 0.973815 | corrected thirteen-sample panel after sample17; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,5,6,11,13,17,18,21,22,25,28 | 0.961965 | 3.912655 pp | 0.998777 | stronger thirteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample24 source-abundance | MinCO current default | 24 | 0.909091 | 13.776347 pp | 0.994696 | same-release r232 counterexample restored under `/mnt/new3T`; 25 TPs, 1 FP, 4 FNs |
| CAMI II HMP airskin sample24 source-abundance | Sylph r232 chunked | 24 | 0.915254 | 4.668799 pp | 0.999353 | slightly stronger F1 and much stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,5,6,11,13,17,18,21,22,24,25,28 | 0,1,3,5,6,11,13,17,18,21,22,24,25,28 | 0.860374 | 22.153669 pp | 0.975306 | corrected fourteen-sample panel after sample24; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,5,6,11,13,17,18,21,22,24,25,28 | 0.958628 | 3.966665 pp | 0.998818 | stronger fourteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample16 source-abundance | MinCO current default | 16 | 1.000000 | 6.914590 pp | 0.999835 | same-release r232 holdout restored under `/mnt/new3T`; 10 TPs, 0 FPs, 0 FNs |
| CAMI II HMP airskin sample16 source-abundance | Sylph r232 chunked | 16 | 1.000000 | 1.288589 pp | 0.999960 | ties F1 and remains much stronger on abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,5,6,11,13,16,17,18,21,22,24,25,28 | 0,1,3,5,6,11,13,16,17,18,21,22,24,25,28 | 0.869682 | 21.137730 pp | 0.976941 | corrected fifteen-sample panel after sample16; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,5,6,11,13,16,17,18,21,22,24,25,28 | 0.961386 | 3.788127 pp | 0.998894 | stronger fifteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample15 source-abundance | MinCO current default | 15 | 0.913043 | 8.896280 pp | 0.999525 | same-release r232 counterexample restored under `/mnt/new3T`; 21 TPs, 1 FP, 3 FNs |
| CAMI II HMP airskin sample15 source-abundance | Sylph r232 chunked | 15 | 0.938776 | 0.978464 pp | 0.999984 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0,1,3,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0.872392 | 20.372639 pp | 0.978353 | corrected sixteen-sample panel after sample15; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0.959973 | 3.612523 pp | 0.998962 | stronger sixteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample4 source-abundance | MinCO current default | 4 | 0.839506 | 28.726741 pp | 0.988378 | same-release r232 counterexample restored under `/mnt/new3T`; 102 TPs, 35 FPs, 4 FNs |
| CAMI II HMP airskin sample4 source-abundance | Sylph r232 chunked | 4 | 0.971963 | 4.592948 pp | 0.999771 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0,1,3,4,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0.870458 | 20.864057 pp | 0.978943 | corrected seventeen-sample panel after sample4; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,11,13,15,16,17,18,21,22,24,25,28 | 0.960679 | 3.670195 pp | 0.999010 | stronger seventeen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample23 source-abundance | MinCO current default | 23 | 0.898551 | 28.419232 pp | 0.984447 | same-release r232 counterexample restored under `/mnt/new3T`; 31 TPs, 5 FPs, 2 FNs |
| CAMI II HMP airskin sample23 source-abundance | Sylph r232 chunked | 23 | 0.969697 | 20.895611 pp | 0.983133 | stronger F1 and L1, while MinCO has slightly higher Pearson |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,11,13,15,16,17,18,21,22,23,24,25,28 | 0,1,3,4,5,6,11,13,15,16,17,18,21,22,23,24,25,28 | 0.872019 | 21.283789 pp | 0.979248 | corrected eighteen-sample panel after sample23; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,11,13,15,16,17,18,21,22,23,24,25,28 | 0.961180 | 4.627162 pp | 0.998128 | stronger eighteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample20 source-abundance | MinCO current default | 20 | 0.821622 | 25.952986 pp | 0.990482 | same-release r232 counterexample restored under `/mnt/new3T`; 76 TPs, 30 FPs, 3 FNs |
| CAMI II HMP airskin sample20 source-abundance | Sylph r232 chunked | 20 | 0.968153 | 4.085099 pp | 0.999816 | same TP/FN counts with far fewer FPs and much stronger abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0,1,3,4,5,6,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0.869366 | 21.529536 pp | 0.979840 | corrected nineteen-sample panel after sample20; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0.961547 | 4.598633 pp | 0.998217 | stronger nineteen-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample7 source-abundance | MinCO current default | 7 | 0.841463 | 27.609400 pp | 0.987554 | same-release r232 counterexample restored under `/mnt/new3T`; 69 TPs, 20 FPs, 6 FNs |
| CAMI II HMP airskin sample7 source-abundance | Sylph r232 chunked | 7 | 0.960526 | 5.205718 pp | 0.999444 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0.867971 | 21.833529 pp | 0.980225 | corrected twenty-sample panel after sample7; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,11,13,15,16,17,18,20,21,22,23,24,25,28 | 0.961496 | 4.628987 pp | 0.998278 | stronger twenty-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample14 source-abundance | MinCO current default | 14 | 0.862069 | 14.824168 pp | 0.998636 | same-release r232 counterexample restored under `/mnt/new3T`; 25 TPs, 5 FPs, 3 FNs |
| CAMI II HMP airskin sample14 source-abundance | Sylph r232 chunked | 14 | 0.949153 | 4.072864 pp | 0.999842 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.867690 | 21.499750 pp | 0.981102 | corrected 21-sample panel after sample14; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.960908 | 4.602505 pp | 0.998353 | stronger 21-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample9 source-abundance | MinCO current default | 9 | 0.929293 | 31.526930 pp | 0.986667 | same-release r232 counterexample restored under `/mnt/new3T`; 92 TPs, 10 FPs, 4 FNs |
| CAMI II HMP airskin sample9 source-abundance | Sylph r232 chunked | 9 | 0.973822 | 7.687184 pp | 0.998098 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,9,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,9,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.870490 | 21.955531 pp | 0.981355 | corrected twenty-two-sample panel after sample9; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,9,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.961495 | 4.742718 pp | 0.998341 | stronger twenty-two-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample10 source-abundance | MinCO current default | 10 | 0.913043 | 20.729208 pp | 0.988950 | same-release r232 counterexample restored under `/mnt/new3T`; 63 TPs, 6 FPs, 6 FNs |
| CAMI II HMP airskin sample10 source-abundance | Sylph r232 chunked | 10 | 0.941176 | 11.172738 pp | 0.993135 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.872340 | 21.902213 pp | 0.981685 | corrected 23-sample panel after sample10; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,20,21,22,23,24,25,28 | 0.960611 | 5.022284 pp | 0.998115 | stronger 23-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample19 source-abundance | MinCO current default | 19 | 0.806202 | 37.173740 pp | 0.963409 | same-release r232 counterexample restored under `/mnt/new3T`; 104 TPs, 43 FPs, 7 FNs |
| CAMI II HMP airskin sample19 source-abundance | Sylph r232 chunked | 19 | 0.968326 | 13.049354 pp | 0.994156 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.869584 | 22.538526 pp | 0.980924 | corrected 24-sample panel after sample19; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.960933 | 5.356745 pp | 0.997950 | stronger 24-sample F1 and much stronger abundance |
| CAMI II HMP airskin sample26 source-abundance | MinCO current default | 26 | 0.882353 | 31.580382 pp | 0.973380 | omitted-sample holdout; another negative case for broad external-baseline superiority |
| CAMI II HMP airskin sample26 source-abundance | Sylph r232 chunked | 26 | 0.926471 | 20.141523 pp | 0.983203 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28 | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28 | 0.870095 | 22.900201 pp | 0.980622 | 25-sample aggregate after sample26; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28 | 0.959554 | 5.948136 pp | 0.997360 | stronger 25-sample F1 and abundance |
| CAMI II HMP airskin sample2 source-abundance | MinCO current default | 2 | 0.872727 | 35.247753 pp | 0.940035 | omitted-sample holdout; negative abundance and F1 case |
| CAMI II HMP airskin sample2 source-abundance | Sylph r232 chunked | 2 | 0.892857 | 6.585788 pp | 0.997264 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.869710 | 23.046896 pp | 0.979288 | 25-sample aggregate after sample2; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,2,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.958210 | 5.405907 pp | 0.997922 | stronger 25-sample F1 and abundance |
| CAMI II HMP airskin sample8 source-abundance | MinCO current default | 8 | 0.808824 | 19.710972 pp | 0.975710 | omitted-sample holdout; high recall but 25 extra calls |
| CAMI II HMP airskin sample8 source-abundance | Sylph r232 chunked | 8 | 0.982456 | 7.967710 pp | 0.992947 | stronger F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | MinCO current-code samples0,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.867154 | 22.425424 pp | 0.980715 | 25-sample aggregate after sample8; still worse than Sylph on F1 and abundance |
| CAMI II HMP airskin source-abundance expanded | Sylph r232 chunked | 0,1,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,28 | 0.961794 | 5.461184 pp | 0.997750 | stronger 25-sample F1 and abundance |
| CAMI II HMP airskin sample8 FP diagnostic | MinCO current calls | 8 | 0.808824 | 19.710972 pp | 0.975710 | 55 TP, 25 FP, 1 FN; FP median split XnY=4 and split breadth=0.004 |
| CAMI II HMP airskin sample8 FP diagnostic | `min_breadth_ge_0.05` post-hoc | 8 | 0.964912 | 19.557355 pp | 0.974871 | removes 22/25 FP and 0 TP on sample8 only; prior 32-sample filter audit rejects this as a default due regressions |
| CAMI II HMP airskin sample8 FP diagnostic | `drop_prob_lt_0.7_and_breadth_lt_0.2` post-hoc | 8 | 0.953271 | 15.481815 pp | 0.982890 | removes all 25 FP but also 4 TP on sample8; useful adaptive-filter clue, not a promoted default |
| cross-panel call-filter audit with HMP omitted samples | MinCO current calls | 35 | mean F1 0.000000 delta | mean L1 0.000000 delta | Pearson 0.000000 delta | current comparison baseline across 4 panels plus HMP samples2/8/26 |
| cross-panel call-filter audit with HMP omitted samples | best mean-F1 filter `drop_prob_lt_0.7_and_breadth_lt_0.2` | 35 | +0.028620 delta | -3.542619 pp delta | NA | improves mean F1/L1 but worsens 10 samples and 2 panels; not a default |
| HMP airskin call-filter extended panel | `min_xny_ge_25` post-hoc | 27 | mean F1 0.919836 | L1 20.853452 pp | 0.980810 | HMP panel has 0 worsened samples and 16 improved, but global audit still rejects default promotion |
| CAMI II HMP airskin abundance variants | best simple fixed-call variant | 6,11 | 0.961058 | 22.780570 pp | 0.925614 | no simple HMP-only fixed-call rescaling beats current calibrated abundance after exact6 |
| CAMI III ToyGut C-direct profile gap | C `minco profile` direct | 3 | 0.284080 | 60.011872 pp | 0.875887 | too conservative for calibrated species default; documented as direct/domain path |
| CAMI III ToyGut C-direct profile gap | prior MinCO full-candidate offline gate | 3 | 0.555312 | NA | NA | much better than C profile but still below Sylph |
| CAMI III ToyGut C-direct profile gap | Sylph | 3 | 0.603485 | 27.789743 pp | 0.923617 | stronger than current MinCO outputs here |
| low-extra offline all29 | MinCO low-extra candidate | 29 | 0.728433 | 0.574735 | NA | small F1/L1 gain over guarded-tail all29 |
| near-split rescue all29 | best near-split candidate | 29 | 0.726088 | 0.580270 | NA | rejected; more FPs offset rescued TPs |
| loose CAMI3 split rescue crossval | readiness baseline | 26 | 0.709348 | 0.623300 | 0.939155 | baseline for the source-readmap-derived rescue |
| loose CAMI3 split rescue crossval | loose rescue | 26 | 0.702619 | 0.645890 | 0.937852 | rejected; CAMI3 source-readmap gain does not transfer |
| broad fixed-call abundance | max(split, unique) mean / zip AF^1 | 20 | 0.703609 | 0.607295 | 0.933795 | current abundance default |
| broad fixed-call abundance | split mean / zip AF^1 | 20 | 0.703609 | 0.609677 | 0.932528 | superseded abundance baseline |
| broad fixed-call abundance | split mean / zip AF^0.25 | 20 | 0.703609 | 0.686698 | 0.871859 | rejected; CAMI3 source-readmap abundance gain does not transfer |
| CAMI III source-readmap abundance | current calibrated abundance | 3 | 0.815957 | 62.732029 pp | 0.823915 | baseline in partial source-readmap namespace |
| CAMI III source-readmap abundance | split mean / zip AF^0.25 | 3 | 0.815957 | 61.280472 pp | 0.841450 | local improvement, not promoted after broad rejection |
| abundance error decomposition | Toy Mouse MinCO-Sylph delta | 3 | pooled F1 -0.066247 | L1 +6.888988 pp | NA | main gap is matched-TP mass allocation, not FP mass |
| abundance error decomposition | HMP airskin MinCO-Sylph delta | 24 | pooled F1 -0.116997 | L1 +17.181781 pp | NA | main gap is matched-TP mass allocation (+11.310500 pp), with extra predicted mass secondary (+4.921452 pp) |
| abundance error decomposition | HMP gastrooral MinCO-Sylph delta | 2 | pooled F1 -0.041393 | L1 +43.923862 pp | NA | main gap is matched-TP mass allocation |
| abundance error decomposition | CAMI3 source-readmap MinCO-Sylph delta | 3 | pooled F1 +0.072382 | L1 +30.099188 pp | NA | main abundance gap is missing truth mass despite higher MinCO F1 |
| abundance oracle bounds | fixed-call truth-renorm oracle overall | 32 | current calls retained | recoverable mean L1 22.825403 pp | NA | truth-aware allocation would beat Sylph on 2/4 panels only; diagnostic lower bound, not a deployable strategy |
| abundance oracle bounds | Toy Mouse fixed-call truth-renorm | 3 | current pooled F1 0.877647 | oracle L1 3.351852 pp vs Sylph 6.731523 pp | NA | allocation alone could close this panel's Sylph L1 gap if calls stayed fixed |
| abundance oracle bounds | HMP airskin fixed-call truth-renorm | 24 | current pooled F1 0.847987 | oracle L1 2.585540 pp vs Sylph 5.356745 pp | NA | allocation alone could close this panel's Sylph L1 gap if calls stayed fixed |
| abundance oracle bounds | HMP gastrooral fixed-call truth-renorm | 2 | current pooled F1 0.915663 | oracle L1 8.952874 pp vs Sylph 5.063129 pp | NA | fixed-call allocation is not enough; missing calls/normalization still matter |
| abundance oracle bounds | CAMI3 source-readmap fixed-call truth-renorm | 3 | current pooled F1 0.814570 | oracle L1 41.833504 pp vs Sylph 32.632841 pp | NA | fixed-call allocation is not enough; detected truth mass is only 79.083248% |
| missed-truth candidate audit | high-abundance missed truth visibility | 32 | visible 9/19 | emitted profiles with uncalled rows 9/32 | NA | diagnostic only; many HMP outputs are emitted-only, so raw-table follow-up is needed |
| missed-truth candidate audit | CAMI3 source-readmap emitted profile | 3 | high FN visible 8/8 | mean FN truth mass 20.916752 pp | NA | high-abundance CAMI3 misses are mostly present below the call gate |
| missed-truth candidate audit | HMP airskin emitted profile | 24 | high FN visible 1/10 | mean FN truth mass 1.292770 pp | NA | most high-abundance HMP misses are not visible in emitted profile rows |
| missed-truth candidate audit | HMP gastrooral emitted profile | 2 | high FN visible 0/1 | mean FN truth mass 4.476437 pp | NA | emitted profile rows do not expose the dominant missed species |
| HMP raw-table missed-truth audit | profile-absent high-FN raw visibility | 10 targets | raw visible 10/10 | raw absent 0 | NA | all high-abundance HMP misses absent from emitted profiles are present in raw unique/split evidence |
| HMP raw-table missed-truth audit | HMP airskin raw evidence | 9 targets | raw visible 9/9 | visible truth mass 19.711829 pp | mean ANI 0.933001 | mostly weak split/exact rows below raw default thresholds |
| HMP raw-table missed-truth audit | HMP gastrooral raw evidence | 1 target | raw visible 1/1 | visible truth mass 8.611111 pp | ANI 0.987029 | dominant Clostridium_F botulinum miss is present as a major raw row |
| cross-panel fixed-call abundance sweep | current calibrated abundance | 4 panels / 32 samples | unchanged | mean official L1 37.006346 pp | mean Pearson 0.918773 | retained as default |
| cross-panel fixed-call abundance sweep | best panel-mean diagnostic, genus-XnY realloc | 4 panels / 32 samples | unchanged | mean official L1 36.291143 pp | mean Pearson 0.920526 | not promoted: zero nontrivial variants are sample-safe and this row worsens one panel |
| adaptive abundance switch audit | best output-feature switch | 4 panels / 32 samples | unchanged | LOPO mean panel delta -0.737243 pp | NA | candidate only; worsens one held-out panel and remains overfit |
| cross-panel output call-filter sweep | best mean-F1 filter, drop P<0.7 and breadth<0.2 | 4 panels / 32 samples | mean delta F1 +0.027950 | mean L1 delta -3.631006 pp | NA | not promoted: worsens 8 samples and 2 panels; zero sample-safe or panel-safe filters |
| HMP airskin output call-filter diagnostic | same filter | 24 samples | F1 0.928581 | L1 18.130992 pp | 0.983506 | HMP improves, but the filter is not universal |
| adaptive call-filter switch audit | best sample-safe switch, min XnY>=25 when max XnY median>=253 | 4 panels / 32 samples | mean delta F1 +0.038640 | mean L1 delta -1.979829 pp | mean Pearson delta +0.001604 | candidate only: LOPO mean delta F1 +0.012057 with zero worsened panels, but needs raw-wrapper and independent holdout validation |
| adaptive call-filter switch with HMP omitted samples | best sample-safe switch, min breadth>=0.05 when max XnY median>=387.5 | 4 panels / 35 samples | mean delta F1 +0.043659 in-sample | mean L1 delta -2.127068 pp | LOPO mean F1 delta -0.046909 | rejected: leave-one-panel-out worsens two holdout panels; current call gate remains default |
| adaptive call-filter wrapper validation | `lopo-min-xny25` wrapper replay | 4 panels / 10 cached-table profiles | expected delta F1 +0.000926 | expected delta L1 -0.006534 pp | expected Pearson delta -0.000005 | call mask matches offline expected rows (counts delta 0, F1 delta 1.11e-16); adaptive switch remains off by default and needs independent holdout |
| adaptive call-filter external exact-split stress | `lopo-min-xny25` on plant/strain cached profiles | plant3-5, strain0-2 | mean delta F1 -0.030282 | switched 4 samples, removed 14 rows | baseline validation counts delta 0 | rejected as default: worsens both external diagnostic datasets, including strain mean F1 -0.059176 |
| supervised abundance calibrator | HGB log l2=1 blend 0.75 LOPO | 4 panels / 11 samples | unchanged | mean L1 delta -0.090593 pp | Pearson delta +0.000712 | candidate only; best mean-L1 model still worsens one panel and needs independent holdout |
| supervised abundance sample28 holdout | RF log leaf3 blend 0.75 trained without sample28 | sample28 | unchanged | L1 18.057211 pp; delta -0.295422 pp | Pearson delta +0.000346 | same-panel support only; Sylph remains much better at L1 1.803063 pp |
| supervised abundance external exact-split | RF log leaf3 blend 0.75 | plant3-5, strain0-2 | unchanged | plant delta +0.001128; strain delta -0.001264 | plant Pearson delta -0.000591; strain +0.001373 | rejected as default: one external diagnostic panel worsens |
| guarded feature allocator external exact-split | implemented `guarded-genus-hit-breadth-a002` | plant3-5, strain0-2 | unchanged | mean L1 delta -0.004134 pp | Pearson delta +0.002907 | stays experimental/off: both diagnostic panel means improve but 2/6 samples have small L1 regressions, max +0.000249 |
| guarded feature allocator release candidate | wrapper parity plus external exact-split | 38 profiles | unchanged | mean L1 delta -0.269171 pp | switched 36/38; improved 34; worsened 2 | do not promote: strict external sample safety fails on plant_holdout3 and plant_holdout4 despite zero F1 change |
| refined guarded feature allocator | output-only combined guard search | 38 profiles | unchanged | best strict rule mean L1 delta -0.262955 pp | `s_xny_median>=230.3`; cached gain preserved 0.976818; external worsened 0 | candidate: removes cached external regressions; current status is set by the score replay row below |
| refined guarded feature allocator parity | implemented `guarded-genus-hit-breadth-a002-xny230` | 38 profiles | unchanged | application mismatches 0 | applied total=33, selected=29, external=4 | wrapper switch validated as opt-in; remains off by default |
| refined guarded feature allocator score replay | implemented `guarded-genus-hit-breadth-a002-xny230` | 38 profiles | unchanged | mean L1 delta -0.262955 pp | switched 33/38; improved 33; worsened 0; guard estimate match max delta 3.11e-7 | strongest opt-in abundance candidate; still needs independent holdout before default |
| refined guarded feature allocator marine diagnostic | independent cached marine exact-split replay | 2 profiles | unchanged | mean L1 delta -0.000522 fraction units | switched 2/2; improved 2; worsened 0; min mapped B/A truth mass 74.120670% | supportive diagnostic only; truth transfer is nonrelease |
| refined allocator independent holdout inventory | release-grade cache audit | 4 release panels / 32 samples | NA | independent release-ready samples 0 | release samples overlap guard cache 32/32; HMP omitted IDs profile pairs 3/5 | explains why refined allocator remains opt-in |
| default release snapshot | selected candidate preset | 4 matched GTDB panels + 5 exact-split diagnostic datasets | candidate beats previous MinCO on all matched panel metrics | vs Sylph: F1 wins 1/4, L1 wins 0/4, Pearson wins 0/4 | runtime: 2.515x Sylph seconds, 0.191x Sylph RSS on 7 timed samples | selected default is stable versus prior MinCO, but not a broad Sylph-beating release claim; experimental allocators remain off |
| universal strategy release gate | code/docs/cached evidence gate | 17 gates | pass=15, fail=0, expected_gap=2 | release-grade panels=4 and abundance L1 losses=4/4 | docs, ANI, default, default-contract tests, allocator-off, abundance boundary, runtime boundaries, CAMI3 post-recovery runner, and next-evidence route contract pass | current default is supported with expected release gaps; goal remains active until more clean holdouts or stronger abundance evidence close the gap |
| universal strategy goal completion audit | active-goal state from cached decision TSVs | 9 checks | pass=6, expected_gap=3, fail=0 | goal_not_complete_expected_gaps_remain | current candidate default supported; release-holdout and abundance claims remain expected gaps; marine, HMP omitted, and CAMI3 samples3-5 routes are completed negative | keep current default for users; next useful work is a strategy change or a new clean same-namespace holdout before any final universal or broad external-baseline superiority claim |
| universal strategy next-evidence routes | expected-gap route audit | 5 routes | top route refined allocator candidate | completed routes marine truth upgrade, HMP omitted samples2/8/12/26/27, and CAMI3 source-readmap samples3-5 | ready without external restore: completed marine route, completed HMP route, completed CAMI3 route, and current default audit | do not treat additional local threshold sweeps as release evidence; keep default and refined allocator opt-in |
| holdout gap action plan | holdout component of release-gate expected gaps | 8 manifest panels | release=4, diagnostic=0, nonrelease=4 | top next actions are nonrelease/mixed-readiness records; marine is completed negative | HMP airskin release panel retained as 24-sample stress evidence; abundance blocker tracked separately | no local clean unfinished route remains; next useful evidence needs strategy change or a new clean same-namespace holdout |
| CAMI3 source-readmap scope audit | cached source rows + CAMI taxonomic profiles | samples0-2 | all-row mapped pct 79.832/79.059/84.914 | profile-scope mapped pct 93.827/94.554/96.870 | 324026 additional in-scope rows needed for samples0/1 to reach 95% before fallback | superseded by exact-binomial fallback/policy audit for samples0-2; still useful as provenance for the release-upgrade path |
| CAMI3 source-readmap extension cache audit | local file-presence audit | samples0-5 | accepted subset ready 3/3 | extension subset ready 0/3 at audit time | superseded by post-recovery scorer below | preserved as provenance for the pre-recovery blocker; do not use this older cache-only row as the current route state |
| CAMI3 source-profile extension diagnostic | local taxonomic-profile source rows | samples3-5 | mapped in-scope GTDB truth 62.729/77.005/67.404% | MinCO F1 wins 0/3 and L1 wins 0/3 versus Sylph | refined allocator guard does not pass | rejected as release substitute; keep pursuing per-read source-readmap truth |
| CAMI3 source-readmap recovery inputs | local file-presence and URL manifest audit | samples3-5 | remote archive URLs 3/3 | local readmaps 3/3 after stream extraction | rescoring artifacts ready 3/3; anonymous reads optional; recovered readmap files total 742M under `/mnt/new3T` | recovery completed without storing full archives; use post-recovery scorer result below |
| CAMI3 source-readmap post-recovery scorer | completed extension scorer | samples3-5 | MinCO pooled F1 0.709193 vs baseline 0.783582 | MinCO L1_union 43.530195 pp vs baseline 21.448926 pp | refined allocator delta 0.000000 pp; F1 wins 0/3 and L1 wins 0/3 | negative independent holdout for allocator promotion; move next route to HMP omitted-sample profile-pair recovery |
| marine source-readmap feasibility | sampled local read maps from marine archives | samples3-5 | sampled rows 6000000 | unresolved transfer rows sampled 314110 | direct source-sequence-to-GTDB matches 0 | local readmaps are useful but not sufficient; marine needs contig/OTU-to-assembly mapping, source FASTA provenance, source-specific truth mapping, or a different clean holdout |
| marine source-mapping local inventory | local file and archive-member audit | samples3-5 | setup `genome_to_id.tsv` and `metadata.tsv` present | pooled mapping remains missing | inventory alone does not promote the route | use setup-metadata truth-upgrade audit for the threshold decision |
| marine setup metadata truth upgrade | setup metadata plus local assembly summary | all 10 gold-profile samples | exact source-name unique 300/977; assembly-summary unique 472/977 | strict source-name-or-assembly rule min mapped truth 95.029253% | 10/10 samples pass the 95% mapped-truth threshold | truth mapping became sufficient; selected-default profile scoring is recorded in the next row |
| marine setup truth profile rescore | selected-default same-namespace profiles | samples3-5 | MinCO mean F1 0.781358 vs baseline 0.838716 | MinCO L1 41.299694 pp vs baseline 40.997305 pp | MinCO Pearson 0.886706 vs baseline 0.923667 | completed negative route; marine does not support default promotion |
| remaining holdout route options | cached route audit after marine profile rescore | 8 routes | local viable uncompleted route: none | locally exhausted routes include HMP omitted, marine selected-default, CAMI3 samples3-5, plant, strain, and mixed-readiness | top remaining route refined_allocator_independent_release_holdout | keep current default with expected gaps; new strategy or new clean holdout needed |
| completed-negative route failure mode | raw best-diff visibility audit | marine3-5 and CAMI3 samples3-5 | FN raw visible 336/338 | high-abundance FN raw visible 13/13 | HMP context high-abundance raw visible 10/10 | missed truth is mostly present before final emission; prioritize in-pass candidate retention/gate recovery |
| raw-side candidate retention sweep | truth-aware offline diagnostic | marine3-5 and CAMI3 samples3-5 | 59/78 simple rules strict-pass | best rule `raw_ani0.95_xny100_breadth0.05` mean route F1 delta +0.030664 | min route F1 delta +0.005061; sample worsens 0; added TP/FP 130/53 | promising implementation candidate, but diagnostic only and not a default change |
| cached exact-split diagnostic | current MinCO exact-split vs Sylph | 5 comparable diagnostic datasets | F1 wins 3/5 datasets | L1 wins 2/5 datasets | mixed | diagnostic only; not release-grade evidence |
| wrapper abundance blend validation | genus-XnY alpha 0.25 | 4 panels | unchanged | mean official L1 41.901839 pp | mean Pearson 0.900292 | experimental only; wrapper replay does not validate as an all-panel default |
| cross-domain edge-EM policy | adaptive edge-EM versus edge-marker baseline | 3 spot samples | mean delta F1 0.000000 | mean delta L1 -0.123827 pp | NA | narrow diagnostic win: improves 2/3 and worsens none versus its own edge-marker baseline |
| cross-domain edge-EM policy | adaptive edge-EM versus Sylph | 3 spot samples | mean delta F1 -0.051425 | mean delta L1 +1.111455 pp | NA | rejected as default: lower F1 than Sylph on all three spot samples |
| CAMI III source/ref ANI | MinCO `reported_ani` basis | n=164 | NA | MAE 0.007627 | Pearson 0.478258 | ZIP-AAF ANI is the useful MinCO continuous ANI signal here |
| CAMI III source/ref ANI | Sylph Adjusted_ANI | n=165 | NA | MAE 0.008038 | Pearson 0.462588 | close external baseline |
| Toy Mouse source/rep ANI | MinCO ZIP-AAF ANI | n=69 | NA | MAE 0.006771 | Pearson 0.492383 | better than raw/emitted MinCO ANI, but not best |
| Toy Mouse source/rep ANI | MinCO raw/emitted ANI | n=69 | NA | MAE 0.025913 | Pearson 0.423006 | rejected for continuous ANI reporting |
| Toy Mouse source/rep ANI | Sylph Adjusted_ANI | n=69 | NA | MAE 0.003862 | Pearson 0.802689 | stronger on this panel |

A separate call-filter sweep gives the same conclusion for F1 gates: no
output-only filter is sample-safe or panel-safe across all panels. The best
mean-F1 filter improves HMP by pruning false positives but worsens CAMI3, Toy
Mouse, or HMP gastrooral samples, so it is not a universal default.

Runtime and memory evidence from completed `/usr/bin/time -v` logs is recorded
in `results/runtime_memory_minco_vs_sylph.tsv`. On the comparable HMP airskin
sample6/11 logs, MinCO current strategy took 4:00 and 4:04 with 3.47-3.63 GiB
peak RSS, while Sylph sketch+profile took 1:27 and 1:31 with 18.55-18.57 GiB
peak RSS. The wrapper now auto-selects the initial evidence scheduler:
`-p1` uses same-stream unique sidecar generation, while `-p>=2` runs the
initial best-diff-unique and best-diff-split raw-read passes concurrently and
splits the `-p/--threads` budget between them. On CAMI II Toy Mouse sample6,
the p16 concurrent scheduler reduced current wrapper wall time from 6:33 to
4:10 at the same measured peak RSS of 3.74 GiB; the generated unique, split,
and exact-split MinCO tables have identical MD5s to the prior sequential run.
Sylph sketch+profile for the same sample took 1:51 and
18.81 GiB RSS. The same-release GTDB r232 chunked Sylph DB build now completed:
20 chunks, 24.097 GiB total output, 25:54 wall time from the driver log, and
2.616 GiB maximum per-chunk RSS. HMP r232 chunked Sylph profiling took 3:24.81
for sample6 and 3:22.43 for sample11, with 25.919 GiB peak RSS. Additional
same-release HMP airskin raw-read runs are mixed: MinCO is faster than Sylph
sketch+profile on samples22 and 13, roughly tied/slower on sample21, and slower
on samples0, 1, 3, 4, 5, 15, 16, 17, 18, 19, 20, 23, 24, and 25, while using about 3.27-4.04 GiB
RSS instead of Sylph's about 25.9 GiB peak RSS. The sample17 run took 3:31.93
wall and 3,839,524 kB RSS for MinCO versus 2:59.16 total sketch+profile wall
and 27,164,768 kB RSS peak for Sylph. The sample24 run took 3:30.53 wall and
3,613,232 kB RSS for MinCO versus 2:33.08 total sketch+profile wall and
27,169,636 kB RSS peak for Sylph. The sample16 run took 3:29.27 wall and
3,636,720 kB RSS for MinCO versus 2:55.53 total sketch+profile wall and
27,142,948 kB RSS peak for Sylph. The sample15 run took 3:30.25 wall and
3,691,312 kB RSS for MinCO versus 2:54.23 total sketch+profile wall and
27,147,624 kB peak RSS for Sylph. The sample4 run took 3:36.19 wall and
3,950,380 kB RSS for MinCO versus 2:25.11 total sketch+profile wall and
27,236,604 kB peak RSS for Sylph. The sample23 run took 3:33.17 wall and
3,732,784 kB RSS for MinCO versus 3:03.23 total sketch+profile wall and
27,169,836 kB peak RSS for Sylph. The sample20 run took 3:35.20 wall and
3,943,344 kB RSS for MinCO versus 3:02.95 total sketch+profile wall and
27,234,152 kB peak RSS for Sylph. The sample19 run took 3:37.65 wall and
4,044,720 kB RSS for MinCO versus 3:13.38 total sketch+profile wall and
27,232,320 kB peak RSS for Sylph.

An experimental same-stream unique sidecar was then added:
`minco ani --readwise-unique-out` writes the best-diff-unique profile table
while the split pass streams the reads. The current implementation fuses the
candidate traversal so best-diff-split and best-diff-unique evidence are
updated from the same candidate group. Smoke tests and the CAMI II Toy Mouse
sample6 benchmark showed identical MD5s for the unique, split, and exact-split
raw evidence tables versus the concurrent scheduler. It was not a speed win:
sample6 took 5:18.90 wall time and 4.69 GiB peak RSS, compared with 4:09.65
and 3.74 GiB for the concurrent dual-pass default. The likely reason is that
the sidecar path loses pass-level concurrency and keeps dual accumulators in
one process. On a bounded sample6 first50k-read p1 subset, however,
same-stream took 8.10s versus 10.68s for sequential unique+split, with
identical unique and split MD5s. Therefore the sidecar path is promoted only
for the auto `-p1` initial scheduler; `-p>=2` remains concurrent dual-pass.

Important decision checks:

| check | status | interpretation |
|---|---|---|
| default launcher forces `universal-auto-exact` | pass | no dataset-specific manual strategy is needed for the wrapper |
| default launcher help explains species boundary | pass | wrapper help separates calibrated species profiling from direct domain profiling |
| docs recommend the default launcher | pass | README and manual point calibrated users to `scripts/minco_profile_default.py` |
| default launcher discovers packaged sidecars | pass | packaged refdbs can provide taxmap/training sidecars so users do not need calibration-path flags |
| C profile boundary documented as conservative direct | pass | C help, README, and manual prevent treating `minco profile` as the calibrated default |
| mixed-panel F1 beats legacy probability gate | pass | supports F1-priority promotion over older MinCO gate |
| low-extra F1 beats guarded-tail all29 | pass | supports adding guarded low-extra split rescue |
| near-split rejected by F1 | pass | prevents lowering split ANI as a default rescue |
| loose CAMI3 source-readmap split rescue rejected by crossval F1 | pass | prevents promoting a sample/namespace-specific CAMI3 rescue |
| adaptive call-filter wrapper matches expected call mask | pass | experimental `lopo-min-xny25` wrapper replay has counts delta 0 and F1 delta 1.11e-16 versus offline expected rows |
| adaptive call-filter external exact-split rejected as default | pass | plant/strain cached exact-split stress worsens both external diagnostic datasets, so the filter stays experimental/off |
| zip-power 0.25 abundance rejected by broad fixed-call L1 | pass | prevents promoting a CAMI3-specific abundance formula |
| zip-power 0.25 improves CAMI III source-readmap L1 | pass | confirms the local diagnostic result, not a universal rule |
| edge-EM adaptive improves edge-marker L1 without F1 cost | pass | confirms the selected_group_species2 edge-EM pilot has a narrow abundance signal |
| edge-EM adaptive beats Sylph F1 on cross-domain spots | known gap | edge-EM has lower F1 than Sylph on CAMI3, marine, and strain spot samples |
| edge-EM default not promoted | pass | preserves F1 priority because edge-EM loses F1 versus current/Sylph baselines |
| reported ANI ZIP-AAF output documented | pass | wrapper, README, and manual expose `reported_ani` as diagnostic output |
| reported ANI beats raw/emitted MinCO ANI on Toy Mouse MAE | pass | supports avoiding saturated raw/emitted ANI for continuous reporting |
| reported ANI close to Sylph on CAMI III source/ref MAE | pass | supports ZIP-AAF ANI as useful diagnostic output |
| reported ANI beats Sylph on Toy Mouse source/rep MAE | known gap | prevents a broad continuous-ANI winning claim |
| CAMI III F1 beats Sylph | pass | supports F1-priority behavior on one holdout panel |
| CAMI III abundance L1 beats Sylph | known gap | Sylph is better for abundance here |
| CAMI III GTDB-transfer F1 beats Sylph | pass | supports MinCO in a comparable but partial GTDB-transfer namespace |
| CAMI III GTDB-transfer truth quality | known gap | mean mapped truth mass is 13.4968% of total abundance, so this is not release-grade |
| CAMI III source-readmap F1 beats Sylph | pass | supports MinCO in a stronger partial source-readmap GTDB namespace |
| CAMI III source-readmap abundance L1 beats Sylph | known gap | Sylph remains better for abundance in this view |
| CAMI III source-readmap truth quality | known gap | mean mapped read-row mass is 81.2684%, so this is stronger but still partial |
| HMP airskin source-abundance truth quality | pass | mean mapped truth abundance is 98.4661%, useful diagnostic counterexample evidence |
| HMP airskin source-abundance F1 beats Sylph | contradicts broad Sylph-beating claim | Sylph has slightly higher mean and pooled F1 |
| HMP airskin source-abundance L1 beats Sylph | known gap | Sylph is much better for abundance in this view |
| mouse5-7 F1 beats Sylph | contradicts broad Sylph-beating claim | use conservative default claim |
| mouse5-7 abundance L1 beats Sylph | known gap | Sylph is better for abundance here |
| marine domain recipe F1 beats Sylph | pass | useful F1 evidence, but it is an older S1000 domain recipe in local species-taxid scoring |
| marine domain recipe is universal-default evidence | known gap | prevents treating the marine domain recipe as release-grade universal evidence |
| marine GTDB-transfer F1 beats Sylph | pass | MinCO beats Sylph on samples0,3,4,5 after conservative unique-taxid GTDB transfer |
| marine GTDB-transfer abundance L1 beats Sylph | pass | MinCO has lower union L1 in this partial transfer namespace |
| marine GTDB-transfer truth quality | known gap | mean scored Bacteria/Archaea mapped truth mass is 74.4053%, below the 95% release-grade threshold |
| plant calibrated gate F1 beats Sylph slightly | pass | supports calibrated rescue, but only in bacteria-scope local scoring |
| plant current exact-split F1 beats Sylph slightly | pass | supports current gate competitiveness, but still non-release evidence |
| plant GTDB transfer truth quality | known gap | only 19.0835-22.2414% of source abundance maps to unique GTDB species, so plant cannot be counted as release-grade GTDB truth |
| strain current exact-split F1 beats Sylph | pass | supports current exact-split F1 on strainmadness, but only as non-release evidence |
| strain current exact-split L1 beats Sylph | known gap | Sylph is much better for abundance in the strainmadness local scorer |
| strain current exact-split is release-grade | known gap | cached outputs use train12/default-equivalent calibration that includes strainmadness rows, plus a local scorer and Sylph r226/default comparison |
| strain GTDB transfer truth quality | known gap | only 0.4775-2.9454% of source abundance maps to unique GTDB species, so strainmadness cannot be counted as release-grade GTDB truth |
| plant C-direct profile recall gap | documented direct-profile gap | C `minco profile` is much lower recall than calibrated plant evidence and is not the calibrated species default |
| HMP pilot Sylph F1 beats calibrated | contradicts broad Sylph-beating claim | HMP gastrooral is a strong counterexample |
| HMP pilot unique direct beats calibrated | known gap | suggests an automatic unique-anchor path is still needed for HMP-like samples |
| CAMI3 C-direct profile below prior offline gate | documented direct-profile gap | C `minco profile` does not expose the calibrated/full-candidate species logic |
| CAMI3 C-direct profile F1 beats Sylph | documented direct-profile gap | C direct profile is below Sylph here, but it is not the calibrated species default |
| all holdout panels are release-grade GTDB | known gap | 4 of 8 panels audit as release-grade; remaining panels still lack clean same-namespace truth and comparable selected-default profiles |

Release-readiness gate:

| gate | status | interpretation |
|---|---|---|
| single default entrypoint | pass | wrapper/docs provide one no-manual-strategy default |
| C profile boundary documented | pass | conservative C profile is explicitly separated from the calibrated wrapper |
| F1 priority beats previous MinCO gate | pass | supports the current MinCO default candidate |
| unsafe rescues rejected | pass | near-split and loose CAMI3 rescues are not promoted |
| abundance default not replaced by local candidate | pass | fixed-call variants, adaptive switch, supervised calibrator, sample28 holdout, external exact-split stress test, and cross-domain edge-EM do not provide a promotable default replacement |
| diagnostic ANI reporting available | pass | `reported_ani` is implemented and documented |
| extended non-release panels audited | pass | marine, plant legacy/current-exactsplit plus plant GTDB-transfer feasibility, strain current-exactsplit plus strain GTDB-transfer feasibility, HMP gastrooral, HMP airskin source-abundance, and CAMI3 C-direct rows are represented in decision checks |
| C direct profile gap documented | pass | C `minco profile` is explicitly separated from the calibrated species default |
| clean release-grade holdout bundle | fail | 4 of 8 panels audit as release-grade GTDB-species evidence; several panels remain partial or mixed-namespace |
| abundance beats Sylph on key panels | known gap | Sylph remains better on important abundance checks |
| broad Sylph-beating claim supported | fail | CAMI II Toy Mouse, HMP, and strainmadness abundance evidence contradict this claim |
| final decision | pre-release candidate | use as best current MinCO default candidate, not a broad release-grade universal claim |

2026-06-30 wrapper-output replay for the combined raw-retention surface:

| method | samples/panels | mean delta F1 | min delta F1 | sample regressions | panel regressions | mean delta L1 pp | TP delta | FP delta | decision |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `accession-current-or-ani93-xny650-br20` + `normalized-depth-alpha2` | 43 / 7 | 0.006464 | -0.003754 | 3 | 1 | -0.444718 | 30 | 12 | diagnostic only |
| `accession-current-or-ani93-xny650-br20` + `zero` | 43 / 7 | 0.006464 | -0.003754 | 3 | 1 | -0.018777 | 30 | 12 | diagnostic only |

The wrapper replay completed all 86 expected profiles. The combined switch
preserves the current candidate surface and adds the stricter raw-retention
surface, so it is a fairer wrapper implementation of the prior offline replay
than replacing the current surface. It improved mean F1 on the aggregate panel
and improved six of seven pooled panels, but all three marine samples regressed
slightly by adding 2-3 extra calls each and no true-call rescues. Therefore
this result does not justify a default or preset change. The result is useful
as implementation evidence for an in-wrapper candidate side channel, and it
narrows the next algorithmic direction to avoiding the marine false-call
pattern while preserving the true-call gains on HMP and CAMI3-like panels.

2026-06-30 max-called-species guard follow-up:

| method | samples/panels | mean delta F1 | min delta F1 | sample regressions | panel regressions | mean delta L1 pp | TP delta | FP delta | decision |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| combined surface active only when pre-surface called species count <=250 | 43 / 7 | 0.006671 | 0.000000 | 0 | 0 | -0.488058 | 30 | 5 | diagnostic candidate |

The guard is implemented as opt-in
`--candidate-surface-max-called-species`; default `0` disables it. The cached
threshold sweep showed thresholds 175-250 give the same strict-pass result.
Threshold 250 was selected for the implementation smoke because it blocks the
three high-complexity marine regressions while preserving every non-marine
candidate-surface gain in this replay. Marine samples3/4/5 were rerun with
`--candidate-surface-max-called-species 250`; all three emitted
`candidate_surface_guard_blocked=True`, had zero candidate-surface additions,
and completed in 0:17.01-0:17.62 wall time with max RSS 3,452,476 KB. This is
still not a default promotion because it is a cached replay plus targeted
implementation smoke, not a fresh independent holdout.

## Validation

- Focused pytest for `tests/test_minco_profile_calibrated_auto_exact.py`
  passed with 54 tests on 2026-06-29. The suite now includes a
  selected-default contract check that reads
  `results/current_default_strategy_manifest.tsv` and verifies the launcher
  constants still match the candidate preset, `universal-auto-exact`, the
  selected candidate rescue/surface switches, normalized-depth-alpha2
  candidate abundance, and off-by-default experimental allocators.
- `score_cami3_gtdb_taxid_transfer.py` generated CAMI III Toy Human Gut
  samples0-5 partial GTDB-transfer truth, quality, per-sample scores, and
  summary tables.
- `score_cami3_gtdb_source_readmap.py` generated CAMI III Toy Human Gut
  samples0-2 partial source-readmap GTDB truth, quality, per-sample scores, and
  summary tables.
- `score_hmp_gtdb_source_abundance.py` generated CAMI II HMP airskin samples6
  and 11 source-genome abundance GTDB truth, quality, per-sample scores, and
  summary tables.
- `summarize_plant_holdout_exactsplit.py` copied compact current exact-split
  plant holdout3-5 MinCO/Sylph score rows from cached `/tmp` files into
  note-local result TSVs. This is non-release evidence because the source files
  are temporary and the scorer is plant bacteria-scope.
- `audit_plant_gtdb_transfer_feasibility.py` generated the plant
  source-abundance GTDB-transfer audit. With CAMI `camitax.tsv` taxID plus
  numeric filename-taxID fallback, only 19.0835%, 20.3716%, and 22.2414% of
  source abundance maps to unique GTDB species for samples3-5.
- `summarize_strain_holdout_exactsplit.py` copied compact current exact-split
  strainmadness samples0-2 MinCO/Sylph score rows from cached `/tmp` files into
  note-local result TSVs. This is non-release evidence because the cached
  train12/default-equivalent calibration includes strainmadness rows and the
  scorer is local.
- `audit_strain_gtdb_transfer_feasibility.py` generated the strainmadness
  source-abundance GTDB-transfer audit from local CAMISIM metadata and coverage
  files. Only 2.9454%, 0.4775%, and 1.3037% of source abundance maps to unique
  GTDB species for samples0-2.
- The current-code HMP airskin refresh regenerated sample11 exact split from
  FASTQ, rebuilt sample11 with the current wrapper, rebuilt sample6 from the
  existing unique and block-mode split raw tables because sample6 FASTQ was no
  longer local, and rescored both samples under the
  `hmp_current_refresh_source_abundance` prefix.
- `cross_validate_cami3_loose_split_rescue.py` generated a 26-sample
  cross-validation for the loose CAMI3 source-readmap split rescue; the rule
  lowered mean F1 from 0.709348 to 0.702619 and worsened L1 from 0.623300 to
  0.645890.
- `search_fixed_call_abundance.py` was rerun after adding `zip_power=0.25`.
  The CAMI3 source-readmap-favored formula worsened broad-panel fixed-call L1
  from 0.609677 to 0.686698 and Pearson from 0.932528 to 0.871859.
- The same broad fixed-call sweep supports the promoted max split/unique
  ZIP-depth abundance formula: mean L1 improves from 0.609677 to 0.607295 and
  Pearson from 0.932528 to 0.933795 with unchanged calls.
- `summarize_edge_em_cross_domain.py` imported the 2026-06-24 cross-domain
  edge-EM spot comparison into `results/edge_em_cross_domain_methods.tsv` and
  `results/edge_em_cross_domain_policy_summary.tsv`. Adaptive edge-EM improves
  L1 versus its own edge-marker baseline on 2/3 spot samples without F1 loss,
  but has lower F1 than Sylph on all three spot samples.
- `scripts/minco_profile_calibrated.py` now reports `reported_ani` and
  `reported_ani_source`, backed by focused regression tests. This is a
  reporting-only change for ANI; call logic is unchanged.
- `validate_holdout_bundle.py` regenerated `results/holdout_bundle_audit.tsv`
  and `results/holdout_bundle_summary.tsv`; the current summary reports 4 of
  8 panels release-grade and 0 diagnostic-grade panels in the manifest.
- `audit_holdout_resources.py` generated
  `results/holdout_resource_audit.tsv` and
  `results/holdout_resource_audit_summary.tsv`. HMP airskin cached inputs,
  chunked GTDB r232 Sylph DB list, and r232 sample6/sample11 profile outputs
  are complete; `hmp_airskin_same_release_sylph_ready=true`.
- `preflight_gtdb232_sylph_db.py` generated
  `results/gtdb232_sylph_db_preflight.tsv`,
  `results/gtdb232_sylph_db_missing_paths.tsv`, and
  monolithic plus chunked GTDB r232 Sylph build/profile/rescore command scripts.
  The preflight reports 199,923/199,923 representative paths present,
  178.909 GiB compressed input, `disk_heuristic_pass=true`,
  `chunk_existing_syldb_count=20`, `chunked_output_complete=true`, and
  `hmp_sample_sketches_ready=true`.
- The chunked GTDB r232 Sylph DB build completed on 2026-06-27 with 20 chunk
  `.syldb` files totaling 24.097 GiB. Per-chunk peak RSS was at most
  2.616 GiB; the driver log ran from 10:39:01 to 11:04:55 UTC.
- HMP r232 chunked Sylph profiling completed for sample6 and sample11. The
  profiles have 18 and 33 output rows, respectively; peak RSS was about
  25.9 GiB for each profile.
- The same-release current-code HMP source-abundance scorer generated
  `results/hmp_current_refresh_r232_source_abundance_summary.tsv`. MinCO
  remained below Sylph on this panel: F1 0.961058 versus 0.970588 and L1
  22.780570 pp versus 5.133659 pp. The sample6 MinCO row now uses an exact
  split table after restoring the temporary FASTQ by streaming selected BAMs
  from the CAMI II archive.
- `audit_hmp_airskin_source_truth_candidates.py` generated
  `results/hmp_airskin_source_truth_candidate_audit.tsv` and
  `results/hmp_airskin_source_truth_candidate_summary.tsv`. Of 29 available
  HMP airskin truth tables, 27 have release-quality GTDB r232 truth transfer,
  and seventeen now have same-release r232 outputs. Samples1/3/4/13/15/16/17/18/21/24/25 are scored
  from the alternate `/mnt/new3T/minco_release_holdouts_20260628` work root.
  The audit ignores local FASTQ placeholders smaller than 1024 bytes.
- HMP airskin sample28 restore/profile results:
  selected-BAM stream restore took 6:41.46 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 1:39.98 wall and 3,506,992 kB RSS; Sylph
  r232 chunked profile took 2:07.01 wall and 27,164,692 kB RSS. Sample28
  source-abundance scoring gives MinCO F1 0.888889, L1 18.352634 pp, Pearson
  0.989034, versus Sylph F1 0.923077, L1 1.803063 pp, Pearson 0.999933.
- The expanded same-release HMP airskin panel (samples6,11,28) gives MinCO
  mean F1 0.937001, L1 21.304591 pp, Pearson 0.946754, versus Sylph mean F1
  0.954751, L1 4.023461 pp, Pearson 0.997693. This strengthens the
  release-readiness boundary: current MinCO remains a default candidate, not a
  broad Sylph-beating abundance strategy.
- HMP airskin sample22 restore/profile results:
  selected-BAM stream restore took 6:48.32 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:28.88 wall and 3,764,272 kB RSS; Sylph
  r232 sketch plus chunked profile took 3:59.34 wall and 27,163,848 kB peak RSS.
  Sample22 source-abundance scoring gives MinCO F1 0.644068, L1 20.648591 pp,
  Pearson 0.988592, versus Sylph F1 0.909091, L1 6.043692 pp, Pearson 0.997079.
  The corrected expanded four-sample same-release HMP airskin panel (samples6,11,22,28)
  gives MinCO mean F1 0.863768, L1 21.140591 pp, Pearson 0.957213, versus Sylph
  mean F1 0.943336, L1 4.528518 pp, Pearson 0.997539. This is a new
  release-grade counterexample for both F1 and abundance.
- HMP airskin sample5 restore/profile results:
  selected-BAM stream restore took 8:07.59 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:32.57 wall and 3,811,504 kB RSS; Sylph
  r232 sketch plus chunked profile took 2:52.33 wall and 27,170,992 kB peak RSS.
  Sample5 source-abundance scoring gives MinCO F1 0.961538, L1 10.337736 pp,
  Pearson 0.998382, versus Sylph F1 0.943396, L1 5.965327 pp, Pearson
  0.998278. The corrected expanded five-sample same-release HMP airskin panel
  (samples5,6,11,22,28) gives MinCO mean F1 0.883322, L1 18.980020 pp, Pearson
  0.965447, versus Sylph mean F1 0.943348, L1 4.815880 pp, Pearson 0.997687.
  Sample5 is a MinCO single-sample F1 win, but it does not overturn the HMP
  airskin abundance or five-sample F1 gap.
- HMP airskin sample0 restore/profile results:
  selected-BAM stream restore took 6:35.72 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:33.54 wall and 3,963,224 kB RSS; Sylph
  r232 sketch plus profile took 3:18.89 wall and 27,177,688 kB peak RSS.
  Sample0 source-abundance scoring gives MinCO F1 0.840580, L1 18.780365 pp,
  Pearson 0.989504, versus Sylph F1 0.903226, L1 2.510444 pp, Pearson
  0.999779. The expanded six-sample same-release HMP airskin panel
  (samples0,5,6,11,22,28) gives MinCO mean F1 0.876198, L1 18.946744 pp,
  Pearson 0.969457, versus Sylph mean F1 0.936661, L1 4.431641 pp, Pearson
  0.998036. Sample0 strengthens the HMP airskin counterexample.
- HMP airskin sample21 restore/profile results:
  selected-BAM stream restore took 6:45.19 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:32.88 wall and 4,001,584 kB RSS; Sylph
  r232 sketch plus profile took 3:30.07 wall and 27,181,760 kB peak RSS.
  Sample21 source-abundance scoring gives MinCO F1 0.880000, L1 34.744967 pp,
  Pearson 0.962690, versus Sylph F1 0.970588, L1 2.197983 pp, Pearson
  0.999889. The expanded seven-sample same-release HMP airskin panel
  (samples0,5,6,11,21,22,28) gives MinCO mean F1 0.876741, L1 21.203633 pp,
  Pearson 0.968490, versus Sylph mean F1 0.941508, L1 4.112547 pp, Pearson
  0.998300.
- HMP airskin sample18 restore/profile results:
  selected-BAM stream restore took 6:01.23 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:31.30 wall and 3,819,440 kB RSS; Sylph
  r232 sketch plus profile took 2:15.22 wall and 27,162,472 kB peak RSS.
  Sample18 source-abundance scoring gives MinCO F1 0.943396, L1 12.023035 pp,
  Pearson 0.994733, versus Sylph F1 1.000000, L1 1.882765 pp, Pearson
  0.999930. The corrected expanded eight-sample same-release HMP airskin panel
  (samples0,5,6,11,18,21,22,28) gives MinCO mean F1 0.885073, L1 20.056058 pp,
  Pearson 0.971770, versus Sylph mean F1 0.948819, L1 3.833824 pp, Pearson
  0.998504.
- HMP airskin sample13 restore/profile results:
  selected-BAM stream restore took 6:12.98 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 1:40.59 wall and 3,431,608 kB RSS; Sylph
  r232 sketch plus profile took 3:20.23 wall and 27,154,424 kB peak RSS.
  Sample13 source-abundance scoring gives MinCO F1 0.925000, L1 54.747984 pp,
  Pearson 0.910524, versus Sylph F1 0.987654, L1 8.344359 pp, Pearson
  0.997403. The corrected expanded nine-sample same-release HMP airskin panel
  (samples0,5,6,11,13,18,21,22,28) gives MinCO mean F1 0.889510, L1
  23.910717 pp, Pearson 0.964965, versus Sylph mean F1 0.953134, L1
  4.334995 pp, Pearson 0.998382.
- HMP airskin sample25 restore/profile results:
  selected-BAM stream restore took 5:38.23 wall and 11,136 kB RSS; MinCO
  current `universal-auto-exact` took 3:34.91 wall and 3,939,252 kB RSS; Sylph
  r232 sketch plus profile took 3:10.13 wall and 27,184,872 kB peak RSS.
  Sample25 source-abundance scoring gives MinCO F1 0.720000, L1 23.876830 pp,
  Pearson 0.984476, versus Sylph F1 0.969072, L1 4.836648 pp, Pearson
  0.998867. The corrected expanded ten-sample same-release HMP airskin panel
  (samples0,5,6,11,13,18,21,22,25,28) gives MinCO mean F1 0.872559, L1
  23.907328 pp, Pearson 0.966916, versus Sylph mean F1 0.954728, L1
  4.385160 pp, Pearson 0.998430.
- HMP airskin sample3 restore/profile results:
  selected-BAM stream restore took 5:27.53 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:30.24 wall and 3,688,364 kB RSS; Sylph
  r232 sketch plus profile took 2:23.51 wall and 27,160,752 kB peak RSS.
  Sample3 source-abundance scoring gives MinCO F1 0.957447, L1 12.670667 pp,
  Pearson 0.997025, versus Sylph F1 0.989247, L1 2.961929 pp, Pearson
  0.999926. The corrected expanded eleven-sample same-release HMP airskin panel
  (samples0,3,5,6,11,13,18,21,22,25,28) gives MinCO mean F1 0.880276, L1
  22.885814 pp, Pearson 0.969653, versus Sylph mean F1 0.957866, L1
  4.255775 pp, Pearson 0.998566.
- HMP airskin sample1 restore/profile results:
  selected-BAM stream restore took 5:41.29 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:33.32 wall and 3,873,456 kB RSS; Sylph
  r232 sketch plus profile took 3:10.56 wall and 27,163,132 kB peak RSS.
  Sample1 source-abundance scoring gives MinCO F1 0.680000, L1 18.105884 pp,
  Pearson 0.997394, versus Sylph F1 0.990291, L1 2.178406 pp, Pearson
  0.999907. The corrected expanded twelve-sample same-release HMP airskin panel
  (samples0,1,3,5,6,11,13,18,21,22,25,28) gives MinCO mean F1 0.863586, L1
  22.487486 pp, Pearson 0.971965, versus Sylph mean F1 0.960568, L1
  4.082661 pp, Pearson 0.998678.
- HMP airskin sample17 restore/profile results:
  selected-BAM stream restore took 5:38.61 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:31.93 wall and 3,839,524 kB RSS; Sylph
  r232 sketch plus profile took 2:59.16 wall and 27,164,768 kB peak RSS.
  Sample17 source-abundance scoring gives MinCO F1 0.773109, L1 26.525181 pp,
  Pearson 0.996010, versus Sylph F1 0.978723, L1 1.872577 pp, Pearson
  0.999967. The corrected expanded thirteen-sample same-release HMP airskin
  panel (samples0,1,3,5,6,11,13,17,18,21,22,25,28) gives MinCO mean F1
  0.856626, L1 22.798078 pp, Pearson 0.973815, versus Sylph mean F1 0.961965,
  L1 3.912655 pp, Pearson 0.998777.
- HMP airskin sample24 restore/profile results:
  selected-BAM stream restore took 5:44.30 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:30.53 wall and 3,613,232 kB RSS; Sylph
  r232 sketch plus profile took 2:33.08 wall and 27,169,636 kB peak RSS.
  Sample24 source-abundance scoring gives MinCO F1 0.909091, L1 13.776347 pp,
  Pearson 0.994696, versus Sylph F1 0.915254, L1 4.668799 pp, Pearson
  0.999353. The corrected expanded fourteen-sample same-release HMP airskin
  panel (samples0,1,3,5,6,11,13,17,18,21,22,24,25,28) gives MinCO mean F1
  0.860374, L1 22.153669 pp, Pearson 0.975306, versus Sylph mean F1 0.958628,
  L1 3.966665 pp, Pearson 0.998818.
- HMP airskin sample16 restore/profile results:
  selected-BAM stream restore took 8:11.14 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:29.27 wall and 3,636,720 kB RSS; Sylph
  r232 sketch plus profile took 2:55.53 wall and 27,142,948 kB peak RSS.
  Sample16 source-abundance scoring gives MinCO F1 1.000000, L1 6.914590 pp,
  Pearson 0.999835, versus Sylph F1 1.000000, L1 1.288589 pp, Pearson
  0.999960. The corrected expanded fifteen-sample same-release HMP airskin
  panel (samples0,1,3,5,6,11,13,16,17,18,21,22,24,25,28) gives MinCO mean F1
  0.869682, L1 21.137730 pp, Pearson 0.976941, versus Sylph mean F1 0.961386,
  L1 3.788127 pp, Pearson 0.998894.
- HMP airskin sample15 restore/profile results:
  selected-BAM stream restore took 5:27.43 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:30.25 wall and 3,691,312 kB RSS; Sylph
  r232 sketch plus profile took 2:54.23 wall and 27,147,624 kB peak RSS.
  Sample15 source-abundance scoring gives MinCO F1 0.913043, L1 8.896280 pp,
  Pearson 0.999525, versus Sylph F1 0.938776, L1 0.978464 pp, Pearson
  0.999984. The corrected expanded sixteen-sample same-release HMP airskin
  panel (samples0,1,3,5,6,11,13,15,16,17,18,21,22,24,25,28) gives MinCO mean
  F1 0.872392, L1 20.372639 pp, Pearson 0.978353, versus Sylph mean F1
  0.959973, L1 3.612523 pp, Pearson 0.998962.
- HMP airskin sample4 restore/profile results:
  selected-BAM stream restore took 5:35.00 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:36.19 wall and 3,950,380 kB RSS; Sylph
  r232 sketch plus profile took 2:25.11 wall and 27,236,604 kB peak RSS.
  Sample4 source-abundance scoring gives MinCO F1 0.839506, L1 28.726741 pp,
  Pearson 0.988378, versus Sylph F1 0.971963, L1 4.592948 pp, Pearson
  0.999771. The corrected expanded seventeen-sample same-release HMP airskin
  panel (samples0,1,3,4,5,6,11,13,15,16,17,18,21,22,24,25,28) gives MinCO mean
  F1 0.870458, L1 20.864057 pp, Pearson 0.978943, versus Sylph mean F1
  0.960679, L1 3.670195 pp, Pearson 0.999010.
- HMP airskin sample23 restore/profile results:
  selected-BAM stream restore took 5:32.89 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:33.17 wall and 3,732,784 kB RSS; Sylph
  r232 sketch plus profile took 3:03.23 wall and 27,169,836 kB peak RSS.
  Sample23 source-abundance scoring gives MinCO F1 0.898551, L1 28.419232 pp,
  Pearson 0.984447, versus Sylph F1 0.969697, L1 20.895611 pp, Pearson
  0.983133. The corrected expanded eighteen-sample same-release HMP airskin
  panel (samples0,1,3,4,5,6,11,13,15,16,17,18,21,22,23,24,25,28) gives MinCO
  mean F1 0.872019, L1 21.283789 pp, Pearson 0.979248, versus Sylph mean F1
  0.961180, L1 4.627162 pp, Pearson 0.998128.
- HMP airskin sample20 restore/profile results:
  selected-BAM stream restore took 5:41.98 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:35.20 wall and 3,943,344 kB RSS; Sylph
  r232 sketch plus profile took 3:02.95 wall and 27,234,152 kB peak RSS.
  Sample20 source-abundance scoring gives MinCO F1 0.821622, L1 25.952986 pp,
  Pearson 0.990482, versus Sylph F1 0.968153, L1 4.085099 pp, Pearson
  0.999816. The corrected expanded nineteen-sample same-release HMP airskin
  panel (samples0,1,3,4,5,6,11,13,15,16,17,18,20,21,22,23,24,25,28) gives
  MinCO mean F1 0.869366, L1 21.529536 pp, Pearson 0.979840, versus Sylph mean
  F1 0.961547, L1 4.598633 pp, Pearson 0.998217.
- HMP airskin sample14 restore/profile results:
  selected-BAM stream restore took 5:30.26 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:31.61 wall and 3,759,664 kB RSS; Sylph
  r232 sketch plus profile took 3:07.66 wall and 27,161,108 kB peak RSS.
  Sample14 source-abundance scoring gives MinCO F1 0.862069, L1 14.824168 pp,
  Pearson 0.998636, versus Sylph F1 0.949153, L1 4.072864 pp, Pearson
  0.999842. The corrected expanded 21-sample same-release HMP airskin
  panel (samples0,1,3,4,5,6,7,11,13,14,15,16,17,18,20,21,22,23,24,25,28)
  gives MinCO mean F1 0.867690, L1 21.499750 pp, Pearson 0.981102, versus
  Sylph mean F1 0.960908, L1 4.602505 pp, Pearson 0.998353.
- HMP airskin sample9 restore/profile results:
  selected-BAM stream restore took 5:37.94 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:33.90 wall and 3,856,592 kB RSS; Sylph
  r232 sketch plus profile took 3:03.06 wall and 27,231,128 kB peak RSS.
  Sample9 source-abundance scoring gives MinCO F1 0.929293, L1 31.526930 pp,
  Pearson 0.986667, versus Sylph F1 0.973822, L1 7.687184 pp, Pearson
  0.998098. The corrected expanded twenty-two-sample same-release HMP airskin
  panel (samples0,1,3,4,5,6,7,9,11,13,14,15,16,17,18,20,21,22,23,24,25,28)
  gives MinCO mean F1 0.870490, L1 21.955531 pp, Pearson 0.981355, versus
  Sylph mean F1 0.961495, L1 4.742718 pp, Pearson 0.998341.
- HMP airskin sample10 restore/profile results:
  selected-BAM stream restore took 5:51.91 wall and 11,264 kB RSS; MinCO
  current `universal-auto-exact` took 3:36.69 wall and 3,833,648 kB RSS; Sylph
  r232 sketch plus profile took 2:21.70 wall and 27,234,804 kB peak RSS.
  Sample10 source-abundance scoring gives MinCO F1 0.913043, L1 20.729208 pp,
  Pearson 0.988950, versus Sylph F1 0.941176, L1 11.172738 pp, Pearson
  0.993135. The corrected expanded 23-sample same-release HMP airskin panel
  (samples0,1,3,4,5,6,7,9,10,11,13,14,15,16,17,18,20,21,22,23,24,25,28)
  gives MinCO mean F1 0.872340, L1 21.902213 pp, Pearson 0.981685, versus
  Sylph mean F1 0.960611, L1 5.022284 pp, Pearson 0.998115.
- `diagnose_hmp_abundance_gap.py` generated fixed-call abundance variant
  sweeps for the original HMP r232 samples6/11 panel. After the exact6
  refresh, simple fixed-call rescaling no longer improves on current
  calibrated abundance:
  the current MinCO L1 remains 22.780570 pp, still far behind Sylph's
  5.133659 pp. The top baseline errors are
  species-mass allocation errors, not a global scale issue: sample6 assigns
  11.0124 pp to a false-positive `s__Bifidobacterium sp947289115` and
  underestimates `s__Prevotella intermedia` by 8.8666 pp; sample11 overweights
  `s__Moraxella bovoculi` by 7.6034 pp and underweights `s__Moraxella ovis` by
  7.1536 pp.
- `build_strategy_decision_summary.py` regenerated `summary.tsv`,
  `results/decision_checks.tsv`, and `results/objective_audit.tsv` from source
  TSVs, including marine, plant, strainmadness, HMP, and CAMI3 C-direct counterexample
  rows.
- `build_release_readiness.py` generated `results/release_readiness.tsv`; the
  final status is `pre_release_candidate`.
- Source wrapper inspection confirms `scripts/minco_profile_default.py`
  delegates to `scripts/minco_profile_calibrated.py --strategy
  universal-auto-exact`.
- `scripts/minco_profile_default.py --help` now advertises itself as the
  recommended no-manual-strategy entry point for calibrated species profiling
  and tells users to keep direct AMR/gene/virus/mixed-domain profiling on the C
  `minco profile` path.
- `scripts/minco_profile_default.py` now auto-discovers packaged
  `species_taxmap.tsv` and `joined_feature_training/` sidecars beside `--ref`
  when explicit flags and environment variables are absent.
- C profile wrapper inspection confirms `minco profile` is documented as the
  conservative direct readwise path, not the calibrated `universal-auto-exact`
  default candidate.
- Source wrapper inspection confirms the low-extra split rescue constants and
  output detail columns exist in `scripts/minco_profile_calibrated.py`.
- Experiment-folder file-size hygiene was checked. The oversized exploratory
  CAMI3 source-readmap split-rescue raw grid TSVs were moved to
  `/tmp/minco_universal_strategy_large_artifacts_20260627/`; compact summaries
  remain in the repo.

## Important Artifacts

See `artifacts.md` for source notes, result paths, and validation paths.

## Conclusion

The stable user-facing default should be `universal-auto-exact` with guarded
exact split, high-extra low-uAF probability-tail rescue, and low-extra split
rescue. This is the best known stable MinCO default under the requested
priority order because it prioritizes F1 and avoids rejected broad rescue rules.

The release-readiness status is `pre_release_candidate`, not `release_ready`.
Default-candidate gates pass, but broad release gates fail because the clean
GTDB holdout bundle is incomplete and Sylph remains stronger for abundance on
important panels.

The claim must stay conservative. Current evidence does not prove a universal
strategy that beats Sylph in most cases. Sylph remains stronger on CAMI II Toy
Mouse Gut samples5-7, HMP gastrooral, HMP airskin source-abundance, and
strainmadness abundance, and usually stronger for abundance L1/Pearson.
MinCO's defensible current claim is "best known stable MinCO default so far,
optimized for F1 first", not
"generally beats Sylph".

The extended non-release panels clarify what should happen next. Marine and
plant contain useful MinCO-supporting signals. Marine now has a conservative
GTDB taxid-transfer check where MinCO beats Sylph on F1 and union L1, but it is
still partial: only about 74.4% of scored Bacteria/Archaea truth mass maps
uniquely to GTDB species, and the MinCO rows are cached S1000 unique ZIP-AAF
domain-recipe outputs rather than the current universal default. Plant remains
bacteria-scope local scoring. The current exact-split plant check is
competitive with Sylph, but it still depends on cached temporary exact-split
outputs rather than a clean release holdout. A source-abundance GTDB transfer
audit confirms plant cannot be
promoted from local files: only 19-22% of sample abundance maps to unique GTDB
species, with most unmapped mass coming from RNODE/no-taxid sources or
genus/family-only CAMI taxids. Strainmadness current exact-split evidence has
the same non-release pattern in the opposite metric direction: MinCO has higher
mean F1 than Sylph, 0.628737 versus 0.527374, but Sylph has much better L1 and
Pearson, and only 0.4775-2.9454% of source abundance maps uniquely into GTDB
species. HMP gastrooral shows that unique direct evidence can beat the
calibrated train12 gate inside MinCO, while Sylph still wins overall. CAMI3
ToyGut shows that the C `minco profile` command is too conservative for
calibrated species calling and should remain the direct/domain profiler rather
than being treated as the calibrated species default.

The new CAMI III GTDB taxid-transfer diagnostic strengthens the statement that
MinCO's current default can outperform Sylph on F1 when predictions and truth
are compared in one GTDB-species namespace. It does not close the holdout gap:
only 13.4968% mean truth mass maps by unique NCBI taxid, so this is a partial
diagnostic rather than a release-grade clean GTDB truth set.

The CAMI III source-readmap diagnostic is stronger and more relevant: MinCO
beats Sylph on F1 after accession-aware transfer, 0.815957 versus 0.741521 on
samples0-2. However, Sylph remains better for abundance L1/Pearson, and only
81.2684% mean read-row mass maps to GTDB species, so this still cannot be
counted as a release-grade clean holdout panel.

The HMP airskin source-abundance panel is now a same-release GTDB r232
counterexample: 98.4661% mean source abundance maps into GTDB species, the
chunked r232 Sylph DB and HMP profiles are complete, and Sylph has slightly
higher F1 plus much better abundance than MinCO on samples6 and 11. This
promotes the panel from "diagnostic because Sylph was r226" to release-grade
external evidence against a broad Sylph-beating default claim. The earlier
sample6 non-exact caveat was removed by restoring a temporary 4.2 GiB FASTQ
through streamed CAMI II BAM extraction and rerunning exact split; both sample6
and sample11 now use exact split tables in the promoted current-code MinCO row.

The loose split rescue discovered from CAMI III source-readmap false negatives
should not be promoted. It improves that partial namespace, but cross-validation
on the broader 26-sample readiness panel reduces mean F1 and worsens abundance
L1.

The same applies to the CAMI3 source-readmap abundance formula
`s_mean / zip_af^0.25`: it improves the partial CAMI3 source-readmap L1 from
62.732029 pp to 61.280472 pp, but fails the broader 20-sample fixed-call
abundance panel. The promoted default abundance formula is now
`max(s_mean / s_zip_af, u_mean / u_zip_af)` for normal rows, with
`s_mean_depth` retained for tail-rescue-added rows.

The same-release HMP r232 abundance-gap diagnostic rejects another easy
abundance-only rescue. Power transforms, split/unique/ZIP-depth variants,
representative output-only quality weights, and genus-local reallocation using
XnY, breadth, and mutation/diff penalties on the fixed MinCO call set do not
close the gap. The best diagnostic variant improves L1 by only 3.0087 pp on
samples6 and 11, leaving a 29.4984 pp mean-L1 gap to Sylph. The dominant
failures are local species-mass allocation errors inside close/reference-
substituted groups, so the next abundance work should target true context/read-
level ambiguity modeling rather than a global scalar, per-row quality
correction, or genus-local post-hoc reweighting.

The runtime comparison also identifies an implementation gap. The current
`universal-auto-exact` default is slower than Sylph because it gets its
features from repeated raw-read scans: one best-diff-unique pass, one
best-diff-split pass, and sometimes a third exact split pass. The older faster
MinCO measurements used simpler direct/single-pass modes and are not the same
default. This is not a necessary algorithmic tradeoff: the calibrated gate does
need both unique and split evidence, but it should not need duplicated FASTQ
parsing or duplicated candidate lookup to get them.

The model-cache benchmark removes one avoidable Python-side cost. On Toy Mouse
sample6, cached table scoring is 15.85 s versus 39.97 s uncached. Combined
with `--strategy universal` and calls-only output, MinCO ran in 1:48.19 versus
Sylph sketch+profile at 1:50.98 on the same read file, with much lower peak
RSS. This does not change the main default decision because the same sample's
GTDB truth still favors Sylph for F1 and abundance accuracy, and the
F1-priority exact default remains slower. The Toy Mouse samples5-7 exact-skip
audit also shows why `universal` is not promoted as the F1-priority default:
exact split adds one true low-abundance species (`s__Lachnoanaerobaculum
orale`) in sample5.

The obvious shortcut of deriving the unique channel from split-pass columns is
unsafe. `results/one_scan_feature_feasibility.tsv` compares cached separate
best-diff-unique and best-diff-split tables. On HMP airskin sample11,
`unique.XnY_ctx` differs from `split.Raw_XnY_ctx` for 32,684 of 72,486 shared
refs (45.090086%); on sample6 it differs for 18,910 of 62,505 shared refs
(30.253580%). ANI, breadth, mean depth, and ZIP-AF also differ. Therefore the
speed fix must be a real one-pass dual accumulator that updates best-diff-
unique and best-diff-split evidence from the same candidate groups, not a
post-hoc remapping of the split output.

To support the next speed and abundance step, the existing internal
context-ambiguity edge dump is now exposed as experimental CLI output:
`minco ani --readwise-edge-out` and `minco profile --edge-out`. The edge sidecar
records per-read/per-context `qctx`, candidate `gid`, `diff`, `best_diff`,
candidate/selected counts, selection flag, and coverage increment. It forces
`--density-block-ctx 0`, is not consumed by the calibrated default, and is
intended as a diagnostic/EM substrate. The audit file
`results/read_context_ambiguity_output_audit.tsv` records why existing
`--track` and aggregate profile tables are insufficient for context-level EM:
they lose the full candidate reference set per context.

The new `--readwise-unique-out` sidecar is a correctness-preserving
single-thread speed fix for initial evidence generation. It proves the unique
table can be emitted from the same read stream and candidate group. The p16
benchmark also shows that same-stream initial evidence alone does not beat the
pass-level concurrency of the current p16 default on exact-mode samples.

For continuous ANI reporting, the wrapper should expose ZIP-AAF ANI as
`reported_ani`. This is better supported than the raw/emitted readwise ANI field,
which can saturate, but it remains diagnostic. It is not a universal
Sylph-beating ANI claim because Toy Mouse source-to-representative ANIm still
favors Sylph Adjusted_ANI.

The exact-split speed optimization now has a working sidecar implementation.
`minco ani --readwise-exact-split-out` writes exact per-read best-diff-split
evidence while the main split pass uses block mode, and
`scripts/minco_profile_calibrated.py --same-stream-exact-split` reuses that
sidecar when `universal-auto-exact` fires. On a CAMI II Toy Mouse first50k
subset, the exact sidecar matched standalone exact output by MD5 and block
output remained unchanged; block+exact sidecar took 8.42 s versus 10.79 s for
separate block+exact C passes. At wrapper level, forced exact first50k took
18.30 s with the sidecar versus 20.64 s with the legacy exact rerun. On full
Toy Mouse sample6, p16 sidecar mode reduced the exact default path from the
previous 4:09.65 to 2:28.51, but it still did not beat Sylph sketch+profile at
1:50.98 on the same sample. Accuracy was unchanged from exact MinCO on sample6
and remained below Sylph: MinCO F1 0.854749, L1 20.062956 pp, Pearson 0.973755;
Sylph F1 0.939086, L1 11.759075 pp, Pearson 0.984876.

The HMP gastrooral sample0 exact-needed run shows why the sidecar remains worth
keeping as a speed substrate. Current default auto-exact used a legacy third
FASTQ pass and took 4:00.22 with 3.7958 GiB RSS. The same raw-read run with
`--same-stream-exact-split` reused the sidecar, produced identical calls and
calibrated abundance to the exact rerun, and took 2:25.93 with 4.6302 GiB RSS.
Chunked r232 Sylph profile time for the same sample was 3:27.18 with 25.9178
GiB RSS. A block-only replay kept the same 34 called species but changed
abundance enough to worsen the two-sample HMP source-abundance L1 from
48.986991 pp to 68.950069 pp, so simply dropping exact for speed is not
acceptable on this panel. This supports a future automatic exact-sidecar policy,
but not unconditional sidecar default, because Toy Mouse sample6 shows eager
sidecar overhead when exact is not needed.

The sidecar-policy audit formalizes that boundary over six available instances.
The existing block-feature policy marks three exact candidates and three
avoid-exact rows, with no observed avoid-exact contradiction in this small set.
However, those block features are known only after the initial block/unique
passes, while eager sidecar must be requested before the split pass starts. Toy
Mouse sample6 remains the key negative control: eager sidecar costs +40.61 s,
+1.0685 GiB RSS, and +0.231379 L1 pp relative to the low-extra exact-skip
default. Therefore auto-sidecar is not promoted without either a validated
preflight exact-need predictor or a cheaper conditional exact sidecar.

I tested the most obvious cheap preflight: run the block/unique wrapper on the
first 50k or first 200k HMP gastrooral reads and use the same exact-trigger
features. This is not stable enough for default promotion. Without a support
floor, the prefix policy mismatches the full decision on sample6 at both 50k
and 200k by incorrectly calling it an exact candidate; with a support floor
(`joined_base_n >= 10` and `raw_unique95_n >= 10`), all four prefix runs are
insufficient. The preflight itself costs about 11-12 s and about 2.44-2.46 GiB
RSS at `-p4`, so it is not useful unless it makes a reliable decision. The
decision remains to reject prefix-preflight default and instead look for a
cheaper conditional exact-sidecar implementation or a stronger preflight signal.

A narrower exact-trigger optimization is now implemented in the default
wrapper. `--exact-split-low-extra-mode skip` is the default and suppresses the
exact split rerun when the block-mode low-extra split rescue already added
candidates; `allow` restores the older behavior. On CAMI II Toy Mouse sample6,
this skipped exact mode (`auto_exact_split_unavailable_reason =
block_low_extra_split_rescue_already_active`) and ran in 1:47.90 with 3.4335
GiB RSS, compared with Sylph sketch+profile at 1:50.98 and 18.8135 GiB. Calls
were unchanged relative to the exact sidecar path on sample6 (TP/FP/FN
153/4/48, F1 0.854749). L1/Pearson matched the block-mode abundance path
better than the exact sidecar path: 19.831577 pp / 0.974331 versus 20.062956
pp / 0.973755. This is a speed improvement, not a broad accuracy improvement:
Sylph remains better on sample6 for F1 and abundance. A forced one-process
p16 run that emitted unique, block split, and exact split sidecars together
took 3:55.94 and 6.0030 GiB, so that one-stream p16 path is not promoted.

The abundance-error decomposition clarifies why the current default is not yet
a Sylph accuracy replacement. On Toy Mouse, HMP airskin samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28,
and HMP gastrooral, MinCO already calls most true species, but assigns too much
or too little mass among true positives: the MinCO-Sylph matched-TP error
deltas are +5.481572 pp, +11.310500 pp, and +38.274541 pp, respectively. On CAMI3
source-readmap, MinCO has higher pooled F1 than Sylph but loses abundance
mainly through missing truth mass: +19.641662 pp versus Sylph. Therefore the
next accuracy work should target true read/context-level mass allocation and
high-abundance FN rescue, not simply another global abundance scalar.

The fixed-call cross-panel abundance sweep reinforces that conclusion. After
adding HMP airskin sample19, no nontrivial variant is strict-safe at the
individual-sample level across 32 evaluated samples. Six variants are
strict-safe at the panel-mean level, and the best panel-mean diagnostic,
`genus_realloc_xny`, reduced mean official L1 by
0.715203 pp but worsened one panel by 1.218914 pp; sample-level safety is still
zero.
Current calibrated abundance remains the stable default until a real
ambiguity-aware allocation model is validated. A simple output-feature adaptive
switch still has a leave-one-panel-out signal, improving three held-out panels
but worsening one with mean delta -0.737243 pp in the current cache. That
pattern is still treated as overfit, so current calibrated abundance remains
the default.

The first supervised abundance-calibrator audit found a stronger candidate than
the fixed formulas or simple switches: `lopo_hgb_log_l2_1_blend0.75` has the
best mean L1 delta at -0.090593 pp with mean Pearson delta +0.000712, but it
still worsens one held-out panel by 0.062803 pp. The stricter all-panel-safe
RF blend has only a smaller mean L1 gain. These effects are trained and
evaluated inside the existing decision panels, and no wrapper implementation has
yet been validated on an independent holdout, so supervised abundance remains
candidate-only and does not replace the current calibrated abundance default.

A stricter post-hoc sample28 holdout trained the same model family without HMP
airskin sample28. The predefined RF log leaf3 blend 0.75 improved sample28 L1
from 18.352634 pp to 18.057211 pp with unchanged F1 and Pearson delta
+0.000346. This means the candidate survives the newly restored HMP sample, but
the improvement is still small and Sylph's same-release r232 L1 remains much
lower at 1.803063 pp. The next promotion test must therefore be a new
independent dataset or newly profiled release-grade samples, not another
same-panel rescore.

The external exact-split stress test rejects the supervised RF abundance
candidate as a stable default. Applying the predefined RF log leaf3 model
without retuning to cached plant and strainmadness panels leaves F1 unchanged,
but none of the tested blends is strict-safe externally. Blend 0.25 already
worsens plant L1 by +0.000126 while improving strainmadness by -0.000141;
blend 0.75 has the best external mean delta, -0.000068, but still worsens plant
by +0.001128 and improves strainmadness by -0.001264. Because the requested
default prioritizes stable cross-panel behavior over mean-only gains, this
tradeoff keeps current calibrated abundance as the default.

A cached exact-split diagnostic now summarizes the older exact-split output
folder as a decision-boundary audit rather than treating those files as release
evidence. Across five comparable diagnostic datasets, current exact-split MinCO
wins F1 on three datasets and loses F1 on two; for abundance L1 it wins two and
loses three. This mixed result supports keeping those cached panels in the
non-release evidence bucket and does not alter the `universal-auto-exact`
default.

The existing cross-domain edge-EM pilot is now part of the decision evidence.
The adaptive selected_group_species2 policy has a real but narrow signal: versus
its own S2000 edge-marker beta=0 baseline it improves mean L1 by 0.123827 pp
across CAMI3 ToyGut, CAMI2 Marine, and CAMI2 Strain Madness spot samples without
changing F1. That is not enough for default promotion under the requested
priority order. Against Sylph, the same adaptive edge-EM policy has lower F1 on
all three spot samples and a worse mean L1 by 1.111455 pp. Against the current
non-edge MinCO rows it loses F1 on CAMI3 and marine. Therefore edge-EM remains a
diagnostic abundance substrate, not the selected default strategy.

I also refreshed the HMP gastrooral pilot with current table-mode
`universal-auto-exact` outputs. This matters because the older HMP note scored
a plain calibrated train12 table, not the current universal fallback logic.
Current MinCO improves the two-sample HMP gastrooral mean F1 from 0.779956 to
0.861646 and beats the unique-direct MinCO signal, but Sylph remains stronger:
mean F1 0.901587 and abundance L1 0.311628 versus current MinCO L1 0.766191.
This weakens the old "calibrated gate fails" wording, but it still remains a
counterexample to a broad Sylph-beating claim.

The stronger HMP gastrooral source-abundance scorer confirms the counterexample
in a cleaner namespace. Chunked GTDB r232 Sylph profiling from existing `.sylsp`
sketches took 3:27.18 for sample0 and 1:58.40 for sample6, both with about
25.9 GiB peak RSS. The raw-read current default rerun improves MinCO's GTDB
source-abundance L1 from the earlier table-mode 68.950069 pp to 48.986991 pp
with unchanged mean F1 0.919568. It still remains far behind Sylph r232 at mean
F1 0.955300 and L1 5.063129 pp. Runtime was 4:00.22 and 3.7958 GiB RSS for
sample0, where auto-exact fired and reran exact split, and 2:03.11 and 3.3268
GiB RSS for sample6, where auto-exact did not fire. This adds a third
release-grade panel to the holdout bundle and strengthens the conclusion that
abundance allocation, not only call gating, is the main remaining
universal-strategy gap.

The follow-up HMP sample0 sidecar run shows the exact-needed speed path can be
faster than Sylph without changing MinCO's output on that sample: 2:25.93 versus
Sylph's 3:27.18, with 4.6302 GiB versus 25.9178 GiB RSS. This does not change
the accuracy conclusion, because Sylph still has better two-sample F1 and much
better abundance. It also does not make sidecar the default yet, because the
existing Toy Mouse sample6 benchmark shows eager sidecar is slower when exact
is skipped.

The generated runtime-readiness gate makes the speed boundary explicit. Across
the three currently timed current-default cases, MinCO is faster than Sylph in
1/3 cases and uses less memory in 3/3 cases. Across two opt-in exact-sidecar
cases, the sidecar is faster than Sylph in 1/2 cases and uses less memory in
2/2 cases. Therefore the current default should be described as a low-memory
default candidate, not yet as a generally faster-than-Sylph default. The next
speed target remains a conditional exact sidecar or cheaper exact-evidence path
that keeps the HMP exact-needed speed win without paying the Toy Mouse
no-exact overhead.

I also audited the most tempting exact-mode shortcut: deriving the normal
density-block split evidence from exact per-read split rows. This is not safe.
Across HMP gastrooral sample0, Toy Mouse sample6, and the Toy Mouse first50k
sidecar smoke pair, 3/3 block/exact table pairs are non-interchangeable. The
largest shared-row normalized-abundance absolute delta is 0.871859, and HMP
sample0 has 118,118/137,162 shared refs with a changed normalized abundance
column. This matches the C semantics: density-block mode merges reads into a
unit before best-diff/split assignment, while exact sidecar assigns each read
separately. Therefore the next speed fix must either predict exact need before
the split pass or reduce the cost of preserving both block and exact evidence
channels; it should not replace block evidence with exact evidence.

I then audited a narrower exact-mode shortcut: run exact split only against
candidate/final references. Existing full exact tables show that only
0.140213-0.154453% of scoped exact rows belong to final called taxids on HMP
sample0 and Toy Mouse sample6, so post-hoc output filtering would be large.
That does not reduce scan cost, and it does not prove a safe candidate-only
rerun. Best-diff assignment depends on references that competed for the same
read/context, while the aggregate per-reference tables do not retain those
competitor groups. Therefore candidate-restricted exact is not promoted unless
the C layer records context/read-level competitor closure or performs the exact
side channel during the full candidate scan.

## Paper-Relevant Claim

Preliminary: A guarded `universal-auto-exact` MinCO wrapper is the current best
single default candidate for F1-priority species profiling, but cross-dataset
evidence is mixed and does not support a broad Sylph-beating claim.

## Caveats

- Scoring namespaces differ across source panels: GTDB species, NCBI/CAMI
  taxids, and source-aware diagnostics are not interchangeable.
- The CAMI III GTDB taxid-transfer diagnostic is conservative and partial:
  most truth taxids are ambiguous against GTDB metadata, and the mean mapped
  truth mass is only 13.4968% of total abundance.
- The CAMI III source-readmap diagnostic covers only samples0-2 because those
  are the locally available read-mapping files. It maps 81.2684% mean read-row
  mass; the remaining 15-21% per sample is excluded from normalized truth.
- The CAMI II marine GTDB taxid-transfer diagnostic is partial. It scores only
  samples0,3,4,5 with paired local MinCO/Sylph profiles, and the scored
  Bacteria/Archaea truth mapping is 71.3731-77.0927%, below the 95% release
  threshold. It supports MinCO in a marine-like namespace but cannot be used as
  release-grade universal-default evidence.
- The plant-associated panel is non-release for GTDB source-abundance claims:
  only 19.0835-22.2414% of source abundance maps uniquely into GTDB species from
  local CAMI taxid metadata.
- The strainmadness panel is non-release for GTDB source-abundance claims:
  only 0.4775-2.9454% of source abundance maps uniquely into GTDB species, and
  the cached exact-split check uses train12/default-equivalent calibration that
  includes strainmadness rows.
- The cached exact-split diagnostic reuses
  `/tmp/minco_exactsplit_universal_20260626` score files. It is useful for
  checking whether older cached panels agree with the current decision
  boundary, but it is not release-grade because it depends on temporary cached
  outputs and heterogeneous local scorers.
- The HMP airskin source-abundance panel is same-release GTDB r232 evidence
  for the external MinCO/Sylph comparison. Truth mapping is high but not
  perfect: sample11 has 3.0679% source abundance unmapped to GTDB r232.
- The HMP airskin sample1, sample3, sample4, sample7, sample9, sample10, sample13, sample14, sample15, sample16, sample17, sample18, sample19, sample20, sample21, sample23, sample24, and sample25 release-grade runs
  are complete and recorded as same-release GTDB r232 evidence under
  `/mnt/new3T/minco_release_holdouts_20260628`. The moving next-action script
  now targets sample26. It requires network access to the CAMI II archive and
  enough scratch headroom to stream selected BAMs, write the selected-source
  FASTQ, and keep temporary MinCO/Sylph outputs. The latest headroom audit
  reports both `/tmp` and the alternate work root
  `/mnt/new3T/minco_release_holdouts_20260628` as below the conservative floor,
  so cleanup is needed before the next network-backed run.
- The HMP gastrooral source-abundance panel is same-release chunked GTDB r232
  evidence with high mapped truth mass, 99.1736% and 99.2685%. It uses cached
  `.sylsp` sketches and newly generated chunked r232 Sylph profiles under
  `/tmp/cami2_hmp_pilot_20260625/run_sylph_r232`. The current MinCO evidence
  now includes raw-read wrapper reruns under
  `/tmp/minco_current_code_hmp_gastrooral_20260627`; those outputs are
  temporary and should be preserved or rerun before release packaging.
- The refreshed current-code HMP airskin MinCO result now has exact split for
  both samples6 and 11. Sample6 required streaming selected BAMs from
  `https://frl.publisso.de/data/frl:6425518/airskinurogenital/sample_6.tar.gz`
  through `samtools fastq`; the regenerated FASTQ is temporary under
  `/tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample6.nonzero.fastq.gz`.
- The GTDB r232 Sylph monolithic DB does not exist locally. The completed
  same-release artifact is the chunked DB list under
  `/mnt/new3T/sylph_db/gtdb-r232-c200-dbv1_chunks/`, with 20 `.syldb` chunks
  totaling 24.097 GiB.
- The r232 representative manifest paths include `GTDBr226_genomes` in the
  directory name. This is a local storage path naming issue in the r232
  manifest, not evidence that the chunked Sylph DB was built from the r226
  representative set.
- The loose CAMI III source-readmap split rescue is diagnostic-only. It is not
  part of the default because it fails cross-panel validation.
- The `zip_power=0.25` abundance alternative is diagnostic-only. It helps the
  partial CAMI III source-readmap namespace but fails broad fixed-call
  abundance validation.
- The HMP r232 abundance variant sweep is diagnostic-only. It shows that simple
  fixed-call rescaling does not close the Sylph abundance gap and should not
  replace the current broad fixed-call abundance default.
- The cross-panel fixed-call abundance variant sweep is diagnostic-only. It
  finds no sample-safe replacement for current calibrated abundance across 32
  evaluated Toy Mouse, CAMI3 source-readmap, HMP airskin, and HMP gastrooral
  samples; some
  panel-mean candidates improve all panels, but they still fail individual
  sample safety or leave-one-panel-out validation.
- The output-only call-filter sweep is diagnostic-only. It proves simple
  post-hoc pruning is not enough for a universal F1 default; future F1 work
  needs richer ambiguity/context evidence or dataset-independent exact evidence,
  not a single threshold.
- The adaptive call-filter switch is also diagnostic-only. Its
  leave-one-panel-out result does not worsen any current held-out panel, but the rule family and
  candidate filters were discovered from the same current panels. It needs
  wrapper-realistic raw-output validation and a new independent holdout before
  any default change.
- The cross-domain edge-EM pilot is diagnostic-only. It improves abundance only
  relative to its own S2000 edge-marker baseline and fails the F1-first promotion
  rule against Sylph/current non-edge baselines.
- `--abundance-genus-xny-blend-alpha` is an experimental abundance-only wrapper
  option with default 0. Alpha 0.25 was tested and rejected as the default
  because it is not an all-panel improvement.
- Several source panels reuse cached tables generated before the latest
  low-extra rescue; the low-extra rule was verified offline on all29 and
  concretely on mouse samples5-7, not by rerunning every raw-read panel.
- The current no-manual-strategy species default is the Python calibrated
  wrapper. The C `minco profile` subcommand remains conservative and lower
  recall on plant and CAMI3 ToyGut, but that is now documented as the
  direct/domain path rather than the calibrated species default.
- The speed-priority cached `--strategy universal` mode is only a measured
  speed option, not a new accuracy default. On Toy Mouse sample6 it is faster
  than Sylph by 2.79 s and uses far less memory, but its F1 and abundance
  accuracy remain below Sylph on that same sample.
- The abundance-error decomposition is diagnostic. It reproduces the official
  score TSVs to numerical roundoff, but it does not itself improve calls or
  abundance. It points to matched-TP mass allocation on Toy Mouse/HMP and
  missing truth mass on CAMI3 as the next targets.
- The fixed-call abundance oracle is also diagnostic and truth-aware. It proves
  that better allocation could close the Toy Mouse and HMP airskin L1 gaps, but
  it also proves that fixed-call allocation alone is not enough on HMP
  gastrooral or CAMI3 source-readmap.
- The missed-truth candidate audit is limited to emitted profile rows. It shows
  that CAMI3 high-abundance misses are usually visible below the gate, but many
  HMP outputs do not expose uncalled candidate rows, so HMP absence must be
  followed up against raw tables or an in-pass candidate side channel before
  designing a rescue rule.
- The HMP raw-table missed-truth audit is also truth-aware and diagnostic. It
  shows that all high-abundance HMP misses absent from emitted profiles are
  present in raw evidence, but it does not yet define a safe output-only or
  raw-side-channel rescue rule.
- The HMP raw-side-channel rescue threshold sweep used the corrected GTDB
  namespace mapping for current calls and cached one best raw row per GTDB
  species under `/tmp/minco_hmp_raw_candidate_cache`. Across HMP airskin
  samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28
  plus HMP gastrooral samples0/6, the best HMP-only rule was
  `raw_rescue_ani0.9_xny100_br0.05_af0.7_zero_rescue_mass`, improving mean
  F1 by 0.007757 with TP +13, FP +1, FN -13, and no per-sample F1 decrease.
  It leaves abundance unchanged by assigning zero rescued mass. Raw-depth
  abundance was rejected by this diagnostic because it severely degraded L1
  and Pearson. The corrected cached sweep took 2:33.55 wall time and
  4,143,500 kB max RSS.
- The cross-panel zero-mass candidate rescue sweep combined HMP raw-cache
  candidates with emitted Toy Mouse and CAMI3 candidate rows. Baseline
  validation against official score TSVs matched with max absolute delta
  7.673861546209082e-13. Across Toy Mouse samples5-7, HMP airskin
  samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28,
  HMP gastrooral samples0/6, and CAMI3 source-readmap samples0-2, the best
  rule was `cross_rescue_ani0.9_xny100_br0.01_af0.7_zero_mass`. It improved
  mean F1 by 0.007319 with TP +19, FP +1, FN -19, no sample-level F1 decrease,
  and no panel-level F1 decrease. L1 and Pearson were unchanged because rescued
  calls carried zero abundance mass. This is a promising F1 candidate, but it
  is not a default until implemented/validated in the wrapper with a real
  raw-side-channel source and a nonzero abundance/reporting policy. The sweep
  took 3:27.49 wall time and 5,250,044 kB max RSS.
- The rescued-species detail audit for that rule found 20 rescued species calls:
  19 true positives and 1 false positive. The rescued true positives sum to
  57.219220% truth abundance across sample-specific truth profiles. By panel:
  Toy Mouse rescued 1/1 TP with 0.775595% truth mass; CAMI3 rescued 5/5 TP with
  38.222911% truth mass; HMP airskin rescued 11/12 TP/FP with 9.272691% truth
  mass; HMP gastrooral rescued 2/2 TP with 8.948024% truth mass. The only FP was
  `s__Haemophilus influenzae_I` in HMP airskin sample19. This supports
  prioritizing a real abundance policy for rescued calls rather than leaving
  them at zero mass.
- `scripts/minco_profile_calibrated.py` now has an off-by-default experimental
  `--candidate-rescue-switch emitted-ani90-xny100-br01-af70`. It adds uncalled
  emitted-profile candidates with max(split,unique) ZIP-AAF ANI >=0.90,
  XnY >=100, breadth >=0.01, and real AF >=0.70; raw MinCO ANI is deliberately
  not used because it can saturate in readwise profiles. Rescued rows receive
  zero abundance mass. On the 10 cached-table wrapper replay profiles
  (Toy Mouse 5/6/7, HMP airskin 6/11, HMP gastrooral 0/6, CAMI3 0/1/2), the
  switch rescued 7 rows with 0 zero-mass violations. Mean wrapper delta versus
  current was F1 +0.006471, L1 -0.061420 pp, and Pearson +0.001163. The replay
  matched the offline cross-panel rule on Toy Mouse, HMP airskin 6/11, and
  CAMI3, but missed two HMP gastrooral raw-cache rescues; max count delta versus
  the offline rule was 1 and max F1 delta was 0.0153453. Therefore this validates
  the implementation path but not a new default. A real default would need an
  in-pass raw-candidate side channel or equivalent, plus an abundance/reporting
  policy for rescued calls. Replay runtime was 1:57.36 with 1,984,092 kB max
  RSS; scoring runtime was 1:47.14 with 5,393,784 kB max RSS.
- The raw-side candidate visibility audit explains why emitted-profile replay
  cannot fully reproduce the raw-cache rescue. Under the selected raw-side rule
  (ANI >=0.90, XnY >=100, breadth >=0.01, real AF >=0.70), there are 14 HMP
  raw-side candidates not already called by the current profile: 13 true
  positives and 1 false positive. Seven are NCBI-taxid rows already called but
  collapsed to another GTDB species/best accession, five have no emitted profile
  row for the raw accession's NCBI species taxid, and only two are ordinary
  profile-row uncalled/other cases. The `Clostridium_F botulinum` HMP
  gastrooral sample0 miss is a representative example: raw accession
  `GCF_001276985.1` and called accession `GCF_000827935.1` share NCBI taxid
  1491 but map to different GTDB species. This makes an accession-level or
  GTDB-species-level raw-side candidate surface the next implementation target;
  another threshold over emitted taxid rows cannot recover most of these cases.
  The audit took 1:02.63 wall time and 4,014,288 kB max RSS.
- A focused accession-level profile-row replay on HMP gastrooral samples0/6
  validates that target mechanism. The experimental postprocessor appended two
  zero-mass raw-side rows to the current profiles, one for
  `s__Clostridium_F botulinum` and one for
  `s__Flavobacterium psychrophilum`. The postprocessed profiles matched the
  offline raw-cache rule exactly for TP/FP/FN and F1 (max count delta 0; max F1
  delta 1.11e-16), with zero zero-mass violations. Pooled HMP gastrooral F1 in
  this focused replay was 0.928571 (TP=78, FP=9, FN=3), while L1/Pearson stayed
  unchanged because rescued rows had zero abundance. This proves that an
  accession-level output surface can represent the raw-cache rescue calls, but
  it is still not a default: it is a focused postprocessor, not an integrated
  raw-side candidate-generation path, and it still has no nonzero abundance
  policy. Runtime was 0:44.07 with 4,010,936 kB max RSS.
- The all-cached-HMP accession-level profile-row replay extends that validation
  to 26 HMP source-abundance profiles: airskin samples
  0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28 and
  gastrooral samples0/6. The postprocessor appended 14 zero-mass raw-side rows
  with zero zero-mass violations and matched the offline raw-cache rule exactly
  for TP/FP/FN counts and F1 (max count delta 0; max F1 delta 1.11e-16). Mean
  F1 delta versus current MinCO was +0.007757. HMP airskin pooled F1 was
  0.852821 (TP=1043, FP=308, FN=52); HMP gastrooral pooled F1 was 0.928571
  (TP=78, FP=9, FN=3). L1/Pearson are unchanged by construction because
  rescued rows carry zero abundance mass. This is stronger evidence that an
  accession-level raw-side output surface can reproduce the raw-cache rescue
  calls, but it remains experimental and not default because it is still a
  postprocessed cached-HMP replay with no integrated raw-side candidate
  generation and no nonzero abundance policy. Rerun runtime was 1:10.19 with
  4,071,784 kB max RSS.
- The wrapper now has an integrated, off-by-default candidate-surface switch:
  `--candidate-surface-switch accession-ani90-xny100-br01-af70`. On focused
  HMP gastrooral cached-table replay, the switch needed the regular
  model/taxid map for model features plus `--candidate-surface-taxmap` for
  accession-level GTDB labels. Without that second label map, the replay added
  zero rows because both target accessions collapsed under the normal taxid
  labels. With the GTDB label map, the wrapper added exactly two zero-mass rows
  (`s__Clostridium_F botulinum` and `s__Flavobacterium psychrophilum`) and
  reproduced the focused postprocessor result: pooled F1=0.928571 (TP=78,
  FP=9, FN=3), mean L1=48.986991 pp, mean Pearson=0.882277. Sample0 cached
  replay took 0:16.41 and 3,385,392 kB RSS; sample6 took 0:14.23 and
  3,078,136 kB RSS. This is the first in-wrapper proof of the accession-level
  surface, but it remains experimental/not default because it is focused and a
  nonzero abundance policy is not validated.
- The all-HMP integrated wrapper replay is now complete for 26 cached profiles
  (HMP airskin samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28
  plus HMP gastrooral samples0/6). With `--candidate-surface-taxmap`, the
  wrapper added 14 zero-mass rows, had zero zero-mass violations, and matched
  both the postprocessor and offline raw-cache rule exactly for TP/FP/FN counts
  with max F1 delta 1.11e-16. HMP airskin pooled F1 was 0.852821 (TP=1043,
  FP=308, FN=52), and HMP gastrooral pooled F1 was 0.928571 (TP=78, FP=9,
  FN=3). L1/Pearson stayed unchanged by construction because rescued rows
  carry zero abundance mass. The cached replay took 430.64 s wall time with
  3,821,056 kB max RSS; `/tmp/minco_candidate_surface_integrated_hmp_20260629`
  was 2.8 MiB after the run. This validates the integrated surface mechanics
  across cached HMP panels, but it remains experimental/not default because it
  is F1-only until a nonzero rescued-row abundance policy is validated.
- A first HMP-only nonzero abundance sweep for integrated candidate-surface
  rows tested 114 candidate-only mass policies without writing profile copies.
  The best tested policy was `normalized_depth_alpha2`, which assigns appended
  rows `2 * max(s,u Normalized_abundance_depth_max)` before profile
  normalization.
  It improved mean per-profile L1 by 0.310173 percentage points versus the
  zero-mass surface with no per-sample L1 increase among the 26 cached HMP
  profiles (mean added mass before normalization 0.018568). The effect is positive but
  small and remains diagnostic/not default because it is HMP-only and does not
  address the broader abundance gap. Runtime was 49.78 s with 3,805,916 kB max
  RSS.
- The `normalized_depth_alpha2` abundance policy is now implemented in the
  wrapper as `--candidate-abundance-policy normalized-depth-alpha2`, still off
  by default (`zero`). The policy treats `2 * max(s,u Normalized_abundance_depth)`
  as normalized-scale candidate mass and converts it to raw scale using the
  base called raw mass before final output normalization. On the same 26 cached
  HMP integrated candidate-surface profiles, wrapper output matched the
  fixed-call policy sweep with count deltas 0, max L1 delta 3.38e-14, and max
  Pearson delta 4.44e-16; raw-mass formula audit max delta was
  8.88e-16. HMP airskin mean L1/Pearson improved from 22.538526/0.980924
  (zero-mass surface) to 22.421594/0.981091; HMP gastrooral improved from
  48.986991/0.882277 to 46.357942/0.896281. The replay took 432.55 s wall time
  with 3,801,932 kB max RSS. This validates HMP wrapper mechanics for the
  nonzero policy, but it remains experimental/not default because full
  cross-panel wrapper validation and a raw-side candidate channel beyond HMP
  are still incomplete.
- The same nonzero-mass idea was then stress-tested across the full cached
  cross-panel selected-candidate set: CAMI II Toy Mouse samples5/6/7, HMP
  airskin samples0/1/3/4/5/6/7/9/10/11/13/14/15/16/17/18/19/20/21/22/23/24/25/28,
  HMP gastrooral samples0/6, and CAMI3 source-readmap samples0/1/2. The call
  set was fixed to the selected zero-mass rule
  `ANI>=0.90; XnY>=100; breadth>=0.01; real_af>=0.70`. The best mean-L1
  policy, `raw_abundance_scaled_alpha2_cap0.05`, improved mean per-profile L1
  by 1.079573 pp but worsened one sample and one panel, so it is unsafe. The
  best sample-safe policy was again `normalized_depth_alpha2`, improving mean
  per-profile official L1 by 0.743320 pp with no per-sample L1 increase and no
  panel-mean L1 increase. This was the first cross-panel support for a nonzero
  rescued-row abundance policy and motivated the wrapper validation below.
  Runtime was 1:56.76 with 5,318,964 kB max RSS.
- The cross-panel wrapper validation is now complete for the same 32 cached
  profiles. Toy Mouse and CAMI3 were replayed through the emitted-profile
  candidate-rescue switch with `--candidate-abundance-policy
  normalized-depth-alpha2`; HMP reused the integrated accession-level
  candidate-surface abundance profiles. The wrapper call set matched the
  fixed-call policy exactly (counts delta 0; max F1 delta 1.11e-16), but
  abundance differed from the older fixed-call diagnostic because that
  diagnostic used previously emitted base profiles while the wrapper replay
  regenerated base abundances from cached tables. Against a matched wrapper
  zero-mass baseline, `normalized-depth-alpha2` improved mean official L1 by
  0.722425 pp with no sample-level L1 regressions and improved mean Pearson by
  0.004448. The wrapper added 20 candidate rows and all 20 received nonzero
  mass; the raw-mass formula audit passed with max delta 3.55e-15. Panel
  official L1/Pearson were: Toy Mouse 13.530191/0.988216, CAMI3
  57.747234/0.864321, HMP airskin 22.421594/0.981091, and HMP gastrooral
  46.357942/0.896281. Runtime was 4:54.54 with 5,479,568 kB max RSS. This
  validates the candidate-abundance wrapper mechanics across the cached
  cross-panel set and supports the selected candidate preset, but broader
  release gaps remain and this is not a broad Sylph-beating claim.
- A direct candidate-default comparison now places the wrapper candidate
  strategy against current MinCO and Sylph on the same four panel/sample sets.
  The candidate strategy dominates current MinCO on all four panels under the
  priority order used here: no panel has lower pooled F1, higher official L1,
  or lower official Pearson. Mean deltas versus current MinCO are
  +0.007199 pooled F1, -1.992106 official L1 percentage points, and +0.013704
  Pearson. Per panel, the candidate improved current MinCO pooled F1/L1/Pearson
  on Toy Mouse by +0.001319/-0.237646/+0.000239, CAMI3 by
  +0.009733/-4.984795/+0.040406, HMP airskin by
  +0.004834/-0.116933/+0.000168, and HMP gastrooral by
  +0.012909/-2.629049/+0.014004. Against Sylph, it only wins F1 on 1/4 panels
  and wins L1/Pearson on 0/4, so this is evidence to promote the strategy for
  MinCO default study, not evidence for a broad Sylph-beating release claim.
- `scripts/minco_profile_default.py` now uses `--profile-preset candidate` as
  the selected default. The preset keeps `universal-auto-exact`, adds the
  validated candidate rescue/surface switches, uses
  `normalized-depth-alpha2` candidate abundance, and discovers
  `candidate_surface_taxmap.tsv`, `accession_species_taxmap.tsv`, or
  `accession_taxmap.tsv` sidecars beside `--ref` unless
  `MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP` is set. Use
  `--profile-preset current` or `MINCO_PROFILE_PRESET=current` to reproduce the
  previous universal-auto-exact default without candidate additions. This
  resolves the wrapper-ergonomics part of the candidate path, but it does not
  change the broader release-readiness conclusion.
- The one-flag candidate preset was replayed on the same 32 cached table-mode
  profiles used for the cross-panel candidate-abundance validation. It matched
  the prior explicit-wrapper candidate strategy in calls exactly
  (`counts=0`, max F1 delta 1.11e-16), added the same 20 nonzero candidate
  rows, and preserved the candidate raw-mass formula with max delta 3.55e-15.
  The only score differences were small abundance shifts in HMP airskin
  (max L1 delta 0.14 percentage points; max Pearson delta 6.88e-05) because
  the preset uses `universal-auto-exact` while the earlier explicit HMP replay
  used `universal`; the HMP airskin panel L1 improved slightly from 22.421594
  to 22.413096. The called-row replay took 9:36.79 with 5,122,588 kB max RSS.
  An earlier report-all replay filled `/tmp` on HMP airskin sample25, so broad
  preset validation should use called-row output unless full candidate tables
  are specifically required.
- `results/current_default_strategy_manifest.tsv` and the tracked
  `CURRENT_DEFAULT_STRATEGY.md` are the compact recall records for the selected
  default. They are generated from the readiness and candidate audit TSVs by
  `write_current_default_strategy_manifest.py` and record the default entry
  point, `candidate` preset, `universal-auto-exact` base strategy, candidate
  rescue/surface switches, `normalized-depth-alpha2` candidate abundance,
  previous-default compatibility path, evidence deltas, ANI reporting rule, and
  `pre_release_candidate` claim boundary.
- `results/abundance_allocation_gap_manifest.tsv` and the tracked
  `ABUNDANCE_ALLOCATION_GAP.md` are the compact recall records for the
  remaining abundance blocker. They are generated from cached abundance
  decomposition, fixed-call allocator sweep, and candidate-vs-baseline TSVs by
  `write_abundance_allocation_gap_manifest.py`. The selected default still
  loses abundance L1 to Sylph on all four cached comparison panels; the main
  gap is matched-call mass allocation on most panels. The later selected-call
  oracle narrows call-recovery-first behavior to one cached HMP airskin sample.
  The best cached fixed-call allocator signal is
  `blend_current_genus_realloc_xny_a0.25`, but it remains an unpromoted
  follow-up until wrapper/raw-output validation and independent holdout checks
  close.
- `audit_candidate_preset_genus_xny_blend.py` tested that follow-up directly
  on the 32 selected-default candidate-preset profiles by post-processing only
  the base called-row abundance mass and preserving candidate-row mass. Calls
  were unchanged, mean L1 improved by 0.780908 percentage points, and all four
  panel means improved, but 6/32 individual samples regressed with a maximum
  L1 increase of 0.993420 percentage points. Therefore
  `--abundance-genus-xny-blend-alpha 0.25` is still not promoted into the
  selected default. The abundance target remains a sample-safe matched-call
  mass allocator, not a panel-mean-only allocator.
- `audit_candidate_preset_genus_xny_blend_guard.py` then tested whether an
  output-only guard could rescue the same blend. It searched 3,722
  single-feature and simple two-feature guards using MinCO output features and
  deterministic blend-perturbation features, with truth used only for scoring.
  There were 394 in-panel sample-safe guards and 3,509 panel-safe guards, but
  leave-one-panel-out selection still produced held-out sample regressions in
  two holdout panels (`max_sample_worse=0.220479` percentage points). The guard
  remains diagnostic and is not promoted.
- `audit_candidate_preset_genus_xny_alpha_sweep.py` tested whether a smaller
  or capped genus-XnY redistribution would be stable enough for the selected
  default. Six variants were tested on the same 32 candidate-preset profiles:
  alpha 0.02, 0.05, 0.10, 0.15, alpha 0.25 with a per-genus movement cap, and
  alpha 0.25 uncapped. Every tested variant improved all four panel mean L1
  values, but none was sample-safe. Even alpha 0.02 regressed 4/32 samples
  (`max_sample_worse=0.044407` percentage points). This closes the simple
  genus-XnY redistribution family as a default candidate for now.
- `write_abundance_strategy_decision_matrix.py` consolidates the abundance
  strategy branches into `ABUNDANCE_STRATEGY_DECISION.md` and
  `results/abundance_strategy_decision_matrix.tsv`. The matrix records the
  current decision: keep the selected candidate preset, keep
  `normalized-depth-alpha2` only as the candidate-row abundance component, and
  do not promote additional fixed-call formula sweeps, adaptive output
  switches, supervised calibrators, edge-level redistribution, or genus-XnY
  blend variants until they pass sample-safe external checks.
- `audit_candidate_callset_oracle_feasibility.py` then asks whether the
  selected candidate call set itself is the limiting factor for abundance. It
  reuses the official GTDB mapping helpers and validates selected-profile
  rescoring with max absolute delta `4.54747350886e-13`. The truth-aware
  selected-call oracle can close the selected-candidate L1 gap versus Sylph on
  31/32 cached samples; only HMP airskin sample20 remains call-recovery-first
  (`selected_oracle_minus_sylph_L1=+3.492771` percentage points). This shifts
  the next abundance target toward a sample-safe matched-call mass allocator,
  with targeted call recovery for the one insufficient-call case.
- `sweep_selected_call_mass_transforms.py` tests whether that allocation
  target can be captured by simple fixed-call transforms on the selected
  candidate profiles. It preserves candidate-added row mass and sweeps small
  base-row power, uniform, global XnY/probability, and within-genus XnY
  transforms. The baseline validates with max absolute delta
  `1.14930287509e-12`, but none of 14 variants is sample-safe. The best ranked
  variant, `base_genus_xny_a0.10`, improves mean L1 by 0.382481 percentage
  points but regresses 4/32 samples. This rejects simple unguarded base-row
  transforms as the default allocator.
- `audit_selected_call_mass_transform_guard.py` then tested whether
  output-derived guards can apply those selected-call mass transforms only on
  stable samples. Cached F1/L1/Pearson score columns are excluded from guard
  rule features. The corrected audit scored 5,659 nontrivial guard rules in
  memory, saved the top 1,000 scored rules, and tested the top 512 by
  leave-one-panel-out. It found 822 in-panel sample-safe guards and 3,600
  panel-safe guards, but leave-one-panel-out still had one holdout mean
  regression and two holdouts with sample regressions
  (`max_sample_worse=0.825999` percentage points). Guarded selected-call mass
  transforms therefore remain diagnostic and are not promoted into the
  default.
- `sweep_selected_call_feature_allocators.py` widened the selected-call
  allocation search to use abundance-like MinCO output features as target mass:
  normalized/relative depth, mean-depth-by-breadth, hit-depth-by-breadth,
  ZIP-AF, probability, real-AF, and XnY variants, both globally and within
  genus. It tested 264 fixed-call variants while preserving candidate-added row
  mass. The best ranked variant,
  `genus_raw_hit_breadth_a0.02`, improved mean L1 by 0.318875 percentage
  points and all four panel means, but it still regressed 1/32 samples
  (`max_sample_worse=0.040714` percentage points). Forty-six variants were
  panel-safe, but none was strict sample-safe, so this allocator family remains
  diagnostic and is not promoted into the default.
- `audit_selected_call_feature_allocator_guard.py` then tested whether
  output-derived guards can apply the feature allocators only on stable
  samples. Cached F1/L1/Pearson score columns are excluded from guard rule
  features. It scored 5,687 nontrivial guard rules, saved the top 1,000, and
  tested the top 512 by leave-one-panel-out. The best sample-safe rule was
  `genus_raw_hit_breadth_a0.02` with
  `base_mass_multi_genus_frac>=0.375015`, improving mean L1 by 0.312324
  percentage points while switching 28/32 cached samples. Leave-one-panel-out
  selected output-derived rules with no held-out sample regressions
  (`mean_delta=-0.146115`, `max_sample_worse=0`). This is the first
  allocation branch to pass the cached LOPO screen, but it is not promoted
  until implemented in the wrapper and validated on independent holdouts.
- `--abundance-feature-allocator-switch guarded-genus-hit-breadth-a002` was
  added as an explicit experimental wrapper option. It keeps the default `off`,
  uses the no-leak output guard
  `base_mass_multi_genus_frac>=0.375015`, and applies the
  `genus_raw_hit_breadth_a0.02` allocator only to base called rows while
  preserving candidate-added row mass. The allocator uses accession/taxmap
  species labels for genus grouping when an accession-level taxmap is available,
  falling back to output `species_name`.
- `validate_feature_allocator_wrapper_parity.py` validates that implemented
  switch on the 32 cached selected-candidate profiles without rerunning raw
  profiling. With accession/taxmap genus labels, the wrapper-function replay
  switches 30/32 samples, improves mean L1 by 0.318866 percentage points, and
  has zero sample regressions (`max_worse=0`). This confirms the implementation
  matches the intended sample-safe cached signal better than a display-name
  implementation, but it still remains experimental until independent holdouts
  pass.
- `validate_refined_feature_allocator_score_replay.py` validates the refined
  `guarded-genus-hit-breadth-a002-xny230` switch through the same fixed-call
  selected-profile and external diagnostic scorers. Baseline validation passes
  (`counts=0`, `F1=1.11e-16`, `L1=1.15e-12`, `Pearson=5.55e-16`). The refined
  replay switches 33/38 profiles, improves 33, worsens 0, keeps F1 unchanged,
  and gives mean L1 delta -0.262955 percentage points. The replay matches the
  rounded guard estimate within `3.10619365229e-07`, so this is now the
  strongest opt-in abundance candidate, but it remains off by default until
  independent holdout evidence passes.
- `validate_refined_feature_allocator_marine_diagnostic.py` adds an
  independent cached marine exact-split stress replay on samples 0 and 2, which
  were not part of the 38-profile refined-guard selection cache. It uses the
  conservative GTDB taxid-transfer marine truth, so it is not release-grade
  evidence (`min_truth_mass_mapped_pct_bacteria_archaea=74.120670`). Within
  that diagnostic namespace, the refined allocator switches 2/2 samples,
  improves 2, worsens 0, leaves F1 unchanged, and improves mean L1 by
  0.000521634 fraction units.
- `audit_refined_allocator_independent_holdout_inventory.py` checks the
  release-grade cache coverage for this allocator. All four release-grade
  manifest panels, totaling 32 samples, already overlap the refined
  guard-selection cache. The omitted HMP airskin sample IDs 2, 8, 12, 26, and
  27 have truth files but no cached MinCO/Sylph profile pairs. Therefore the
  refined allocator has independent diagnostic support but still lacks unused
  release-grade holdout evidence.
- The experimental `--same-stream-exact-split` path is correct and improves
  exact-mode runtime, but it is not the no-expertise default. It made HMP
  gastrooral sample0 faster than Sylph with identical MinCO output, but it adds
  work when the exact gate does not fire and increased Toy Mouse sample6 RSS to
  4.5020 GiB. The new default low-extra exact skip makes Toy Mouse sample6
  faster than Sylph, but this is a sample-level runtime win with lower
  F1/abundance accuracy than Sylph.
- The runtime-readiness audit is small and timing-path specific. It supports a
  consistent memory advantage in the timed cases, but not a general default
  speed claim: current default is faster than Sylph in only 1/3 formal timed
  cases. The later HMP airskin runs are mixed: sample22 and sample13 are MinCO
  speed wins versus Sylph sketch+profile, sample21 is roughly tied/slower, and
  samples0, 1, 3, 4, 5, 7, 9, 10, 14, 15, 16, 17, 18, 19, 20, 23, 24, and 25 are slower than Sylph sketch+profile. All keep the
  memory advantage, and Sylph remains more accurate on most of those HMP
  samples.
- The block/exact split semantics audit rejects an optimization shortcut, not
  exact sidecar itself. Exact sidecar can replace a later exact rerun, but exact
  per-read split rows cannot replace the normal density-block split rows used
  for block-mode feature generation.
- The candidate-restricted exact audit rejects another speed shortcut. Although
  post-hoc filtering keeps only 0.140-0.154% of scoped exact rows for final
  called taxids in the two audited cases, current aggregate tables do not prove
  that rerunning exact against only those candidates would preserve best-diff
  assignments.
- Abundance remains weaker than Sylph on important panels.
- Continuous ANI reporting remains diagnostic. `reported_ani` exposes the best
  supported MinCO continuous ANI signal, but no universal MinCO ANI reporter is
  promoted as a Sylph-beating default.
- The CAMI3 source-readmap release-upgrade path is now sharper. The cached
  scope audit showed profile-scope mapping is close to the 95% release target
  for samples 0 and 1. A new cached resolver-candidate audit reviewed the top
  20 in-scope unmapped sources per sample and found that an exact
  GTDB-binomial metadata fallback would resolve 1,732,020 reviewed rows. The
  projected profile-scope mapping after that candidate fallback is sample0
  95.941%, sample1 96.806%, and sample2 99.908%. This is strong enough to run
  a rescored truth-transfer audit, but it is not promoted as release evidence
  until that audit proves the scores.
- The CAMI3 source-readmap exact-binomial fallback rescore used only cached
  source-level tables and cached profiler outputs. It adds 1,732,020 truth
  rows and reaches the 95% profile-scope coverage threshold for samples 0, 1,
  and 2. Against that adjusted truth, MinCO has higher F1 than Sylph
  (`0.853058` versus `0.762024` mean F1), but Sylph remains much better for
  abundance (`31.103234` versus MinCO `61.884168` mean L1 percentage points).
  This strengthens CAMI3 as a potential F1 holdout and reinforces the
  remaining abundance-allocation gap; it still does not justify changing the
  selected default.
- The stricter selected-default CAMI3 rescore uses the cached `candidate`
  preset profiles instead of the older CAMI3 profile set. Under the accepted
  exact-binomial fallback truth policy, the selected default has mean F1
  `0.862183` versus Sylph `0.762024`, but mean L1 is still worse
  (`55.919344` versus Sylph `31.103234`) and mean Pearson is worse
  (`0.866593` versus Sylph `0.947203`).
- The exact-binomial fallback truth policy passes the cached audit: 51/51
  fallback rows are unique GTDB-binomial assignments; no source-priority
  conflicts remain after allowing only same-species unmapped remainders; and
  samples 0, 1, and 2 all exceed the 95% profile-scope coverage threshold.
  The CAMI3 source-readmap row was therefore promoted from diagnostic to
  release-grade evidence in `holdout_bundle_manifest.tsv`.
- After regenerating the holdout and release gates, the clean matched snapshot
  has four release-grade panels and zero diagnostic matched panels. The broad
  claim boundary is unchanged: candidate default still beats previous MinCO on
  the matched panels without regressions, but versus Sylph it wins F1 on only
  1/4 panels and wins neither L1 nor Pearson on any matched panel.
- The abundance release-blocker audit now makes that secondary-priority gap
  explicit in the release gate: selected candidate has 0/4 F1/L1/Pearson
  regressions versus previous MinCO, but loses abundance L1 and Pearson on 4/4
  matched panels versus the external baseline. The cached oracle says
  allocation-only could close 31/32 selected-candidate samples, while 1/32
  first needs call recovery; this keeps abundance as an expected gap rather
  than a default-blocking failure.
- The refined guarded abundance allocator remains opt-in. Its CAMI3
  source-readmap samples 3, 4, and 5 validation route is now completed and
  negative for promotion. The selected-default table-mode replay profiles were
  generated from cached raw MinCO tables with exit status 0 for all three
  samples, 7.99-8.67 seconds wall time per sample, and max RSS 1,601,816 KB.
  The later readmap-only recovery restored all 3 required read maps without
  storing full archives, and the post-recovery scorer found that candidate
  MinCO wins F1 on 0/3 and L1 on 0/3 samples versus the external baseline.
  The refined allocator has no effect on this route (`mean_L1_delta_pp=0`).
- A local CAMI3 samples3-5 source-profile fallback was tested as a possible
  shortcut, using taxonomic-profile strain/source rows plus deterministic GTDB
  transfer. It is not adequate: mapped in-scope GTDB truth is only 62.728832%,
  77.004747%, and 67.403900% for samples 3, 4, and 5. In that diagnostic
  namespace, Sylph beats candidate MinCO on all three samples for both F1 and
  L1, and the refined allocator does not apply because its guard does not pass.
  Therefore this path should not replace per-read source-readmap truth.
- A source-readmap recovery-input audit made the CAMI3 samples3-5 blocker
  concrete, and the blocker has now been removed. The existing download
  manifest contains URLs for all three sample read archives. The recovery
  command file streams each archive URL through `tar` and extracts only
  `reads_mapping.tsv.gz`, avoiding full archive storage and FASTQ extraction.
  The recovered external files are
  `/mnt/new3T/minco_cami3_toygut_extra_20260621/sample_{3,4,5}_reads/reads_mapping.tsv.gz`,
  gzip-validated, and total 742M.
- The post-recovery extension scorer now records the completed independent
  extension result. For samples3-5, the external baseline has mean F1
  `0.784196`, pooled F1 `0.783582`, mean L1_union `21.448926` pp, and
  Pearson_union `0.981150`. Candidate MinCO has mean F1 `0.706093`, pooled F1
  `0.709193`, mean L1_union `43.530195` pp, and Pearson_union `0.922077`.
  The refined allocator row matches candidate MinCO exactly on this route, so
  it remains opt-in and should not be promoted from this evidence.
- `tests/test_cami3_extension_after_recovery.py` now protects this route with
  focused regression tests for the clean missing-readmap blocker, accepted
  exact-binomial fallback on an unmapped source, and rejection when fallback
  conflicts with a prior WGS/unique-taxid source mapping.
- The CAMI II marine exact-binomial transfer audit improves truth mapping but
  does not clear the release threshold. Across all 10 gold-profile samples, the
  fallback rescues 59.034500 abundance-percentage-points of species truth mass,
  but the minimum mapped Bacteria/Archaea truth mass is only 91.005806%, below
  the 95% threshold. The scored cached marine subset remains below threshold as
  well, so marine stays nonrelease and should not be promoted without stronger
  truth mapping plus selected-default profiles in the same namespace.
- HMP airskin sample26 was restored and scored from the generated
  next-sample command after freeing scratch space. The stream-to-FASTQ step
  took 5:48.61 wall time with 11 MB peak RSS and produced a 4.2 GB temporary
  FASTQ. The current MinCO default took 3:34.08 wall time with 3.83 GB peak
  RSS. The external baseline took 0:50.90 for sketching plus 3:00.89 for
  profiling, with 27.23 GB peak RSS during profiling. On sample26, the
  baseline scored F1=0.926471 and L1_union=20.141523 pp; MinCO current default
  scored F1=0.882353 and L1_union=31.580382 pp. On the 25-sample HMP aggregate
  including sample26, the baseline scored mean_F1=0.959554 and
  mean_L1_union=5.948136 pp; MinCO current default scored mean_F1=0.870095 and
  mean_L1_union=22.900201 pp. This is another negative abundance/F1 holdout
  for a broad external-baseline superiority claim, but it provides one of the five
  omitted profile pairs needed for the independent HMP route.
- The route audit bookkeeping was corrected after sample26: the HMP omitted
  inventory now reads truth=5/5, MinCO profile=1/5, baseline profile=1/5, and
  ready pairs=1/5. The next generated action is sample2, with `/tmp` free
  12.348 GiB against a 6.531 GiB estimated requirement and a 12 GiB
  conservative floor. The release gate remains pass=15, fail=0,
  expected_gap=2, so the current default stays supported but the active goal is
  not complete.
- HMP airskin sample2 was then restored and scored with the same generated
  next-sample command. The stream-to-FASTQ step took 5:38.11 wall time with
  11 MB peak RSS and produced a 4.2 GB temporary FASTQ. The current MinCO
  default took 3:33.10 wall time with 3.69 GB peak RSS. The external baseline
  took 0:50.07 for sketching plus 3:02.62 for profiling, with 27.17 GB peak
  RSS during profiling. On sample2, the baseline scored F1=0.892857 and
  L1_union=6.585788 pp; MinCO current default scored F1=0.872727 and
  L1_union=35.247753 pp. On the 25-sample HMP aggregate including sample2,
  the baseline scored mean_F1=0.958210 and mean_L1_union=5.405907 pp; MinCO
  current default scored mean_F1=0.869710 and mean_L1_union=23.046896 pp.
  This is another negative holdout for a broad external-baseline superiority
  claim. The omitted-profile-pair route now has 2/5 ready pairs.
- Before sample8, two approved completed FASTQ caches were removed from `/tmp`:
  `/tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample0.nonzero.fastq.gz`
  and `/tmp/cami2_hmp_unseen_transfer_20260626/reads/airskinurogenital_sample5.nonzero.fastq.gz`.
  This freed about 8.278 GiB of scratch data that can be regenerated from the
  recorded commands and left the completed profile/score artifacts in place.
- HMP airskin sample8 was then restored and scored with the same generated
  next-sample command. The stream-to-FASTQ step took 6:09.88 wall time with
  11 MB peak RSS and produced a 4.2 GB temporary FASTQ. The current MinCO
  default took 3:36.38 wall time with 3.90 GB peak RSS. The external baseline
  took 0:50.18 for sketching plus 3:01.05 for profiling, with 27.18 GB peak
  RSS during profiling. On sample8, the baseline scored F1=0.982456 and
  L1_union=7.967710 pp; MinCO current default scored F1=0.808824 and
  L1_union=19.710972 pp. On the 25-sample HMP aggregate including sample8,
  the baseline scored mean_F1=0.961794 and mean_L1_union=5.461184 pp; MinCO
  current default scored mean_F1=0.867154 and mean_L1_union=22.425424 pp.
  This run is high recall for MinCO but has 25 extra calls, so it further
  rejects a broad external-baseline superiority claim. The omitted-profile-pair
  route now has 3/5 ready pairs.
- The sample8 extra-call diagnostic shows the extra calls are mainly
  low-support retained rows, not rescue rows: FP median split XnY is 4 and
  median split breadth is 0.004, while TP medians are split XnY 575 and split
  breadth 0.473972. On sample8 alone, `min_breadth_ge_0.05` would improve F1
  from 0.808824 to 0.964912 by removing 22 FP and no TP; the existing
  32-sample cross-panel call-filter audit still rejects that simple filter as
  a default replacement because it has sample regressions. Treat this as
  evidence for a future adaptive/conditional filter, not a current default
  change.
- Re-running the cross-panel call-filter sweep with the restored HMP omitted
  samples2/8/26 strengthens the same conclusion. On the expanded 27-sample HMP
  panel, `min_xny_ge_25` improves or preserves F1 for every HMP sample and
  raises mean F1 from 0.867923 to 0.919836. Across all 35 samples and four
  panels, however, no tested output-only filter is sample-safe; the best
  mean-F1 filter worsens 10 samples and 2 panels. Keep the current call gate
  as default and treat the HMP signal as input for an adaptive trigger rather
  than a global filter.
- The expanded adaptive call-filter switch audit is also negative for default
  promotion. In-sample, the best sample-safe switch improves mean F1 by
  0.043659 and mean L1 by 2.127068 pp, but leave-one-panel-out validation has
  mean F1 delta -0.046909 and worsens two held-out panels. This supersedes the
  earlier in-panel adaptive-switch optimism for default selection; keep
  adaptive call-filter switches experimental/off.
- The next-action planner was corrected after sample8 so it derives completed
  HMP samples from the candidate audit and falls back to the explicit omitted
  sample route instead of defaulting to sample0. The generated action now
  targets sample27. Sample27 has truth_mapped_pct=94.722468, below the 95%
  release-grade threshold, so this is independent route-completion evidence
  rather than a clean release-grade upgrade. Scratch headroom is still blocked:
  `/tmp` has 11.734 GiB free versus a 6.531 GiB estimated requirement and a
  12 GiB conservative floor; the alternate root has 10.217 GiB free and also
  fails the floor. Do not start another recovery run until cleanup is approved
  or another work root is available.
- After cleanup approval, old generated FASTQ caches were removed rather than
  profile/score artifacts. Two completed `/tmp` read caches were removed first
  to unblock sample27, and 18 generated HMP FASTQ caches were then removed
  under `/mnt/new3T/minco_release_holdouts_20260628`. `/mnt/new3T` free space
  increased from roughly 11 GiB to 86 GiB, and no `*.fastq*` or `*.fq*` files
  remain under that holdout workspace. These read caches are regenerable from
  the recorded stream commands.
- HMP airskin sample27 was restored and scored after the cleanup. The
  stream-to-FASTQ step took 6:04.58 wall time with 11 MB peak RSS. The current
  MinCO default took 3:38.49 wall time with 4.04 GB peak RSS. The external
  baseline took 0:51.17 for sketching plus 3:08.02 for profiling, with
  27.23 GB peak RSS during profiling. On sample27, the baseline scored
  F1=0.955056 and L1_union=11.662543 pp; MinCO current default scored
  F1=0.723214 and L1_union=35.635571 pp. The 28-sample aggregate after
  sample27 had baseline mean_F1=0.957830 and mean_L1_union=6.247123 pp versus
  MinCO mean_F1=0.862755 and mean_L1_union=23.682118 pp.
- HMP airskin sample12 was then restored and scored. The stream-to-FASTQ step
  took 5:44.84 wall time with 11 MB peak RSS. The current MinCO default took
  3:33.48 wall time with 3.69 GB peak RSS. The external baseline took 0:49.86
  for sketching plus 3:00.76 for profiling, with 27.15 GB peak RSS during
  profiling. On sample12, the baseline scored F1=0.964706 and L1_union=
  13.765698 pp; MinCO current default scored F1=0.901099 and L1_union=
  27.438082 pp. The 29-sample HMP aggregate after sample12 had baseline
  mean_F1=0.958067 and mean_L1_union=6.506384 pp versus MinCO mean_F1=
  0.864077 and mean_L1_union=23.811634 pp.
- The route audit now counts the HMP omitted sample route
  samples2/8/12/26/27 as completed and scored negative, rather than treating a
  nonexistent sample -1 as a scratch-headroom blocker. Completed routes are now
  HMP omitted and CAMI3 source-readmap samples3-5; the top unresolved route is
  marine truth upgrade. The release gate passes with expected gaps only, and
  the active goal remains incomplete.
- The marine truth-upgrade route was tested from cached metadata with a weaker
  diagnostic fallback: within an ambiguous NCBI taxid, assign only when a
  single GTDB candidate has the same species epithet. This rescues 51 rows and
  3.146300 abundance-percentage-points in total, but still leaves only 2/10
  marine samples above the 95% mapped Bacteria/Archaea truth threshold; the
  minimum mapped fraction rises only from 91.005806% to 91.292431%. An unsafe
  upper bound that assigns every ambiguous taxid would pass all samples, so the
  blocker is not missing GTDB metadata alone; it is unresolved split ambiguity.
  The route audit now records `truth_upgrade_needs_external_truth_source` for
  marine rather than treating it as an untested local threshold problem.
- A marine source-readmap feasibility audit then sampled 6,000,000 readmap rows
  from local marine archives for samples3-5. The readmaps cover 314,110 sampled
  unresolved transfer rows, but direct source-sequence-to-GTDB matches are zero
  (`sequence_accession_rows=0;core_rows=0`). This sharpens the marine route
  decision to `truth_upgrade_needs_source_mapping`: the next useful step is
  contig/OTU-to-assembly metadata, source FASTA provenance, source-specific
  truth mapping, or a different clean holdout, not another local threshold
  sweep.
- A local source-mapping inventory audit then checked the marine workspace and
  the three local marine read archives. It found no local
  `genome_to_id.tsv`, `metadata.tsv`, or `gsa_pooled_mapping.tsv.gz` crosswalk
  for marine, and all three read archives expose only `reads_mapping.tsv.gz`
  among relevant members. Sibling CAMI2 plant/strain workspaces contain
  crosswalk files, but OTU IDs are reused across panels and cannot resolve the
  marine source IDs.
- Superseding 2026-06-30 audit: the marine setup archive metadata was extracted
  to `/tmp/cami2_marine_samples3_5_20260625/setup_metadata`, giving local
  `genome_to_id.tsv` and `metadata.tsv` for 977 setup genomes. Exact normalized
  source-name matching to GTDB r232 resolves 300 setup genomes uniquely, 8
  ambiguously, and leaves 669 without an exact source-name match. The local
  RefSeq assembly-summary table adds strain/isolate evidence: 472 setup
  genomes resolve uniquely, 6 ambiguously, and 499 remain unmatched under that
  source. A release-style strict source-name-or-assembly rule reaches minimum
  mapped in-scope truth 95.029253% with 10/10 samples passing the 95%
  threshold. That route is now completed with selected-default same-namespace
  profile scoring.
- The selected-default marine profile rescore ran samples3/4/5 with exact split
  in about 4:17 wall time per sample and 4.15-4.20 GB peak RSS. It is a
  completed negative route: MinCO mean F1 0.781358, L1 41.299694 pp, Pearson
  0.886706, versus the external baseline mean F1 0.838716, L1 40.997305 pp,
  Pearson 0.923667.
- A remaining-holdout route audit now finds no local viable uncompleted route.
  Plant and strain cached transfers are not close to release-grade truth
  coverage (plant minimum mapped GTDB truth 19.083472%; strain minimum
  0.477515%). HMP omitted, marine selected-default scoring, and CAMI3
  samples3-5 are completed negative. Therefore the next useful action is a
  strategy change or a new clean same-namespace holdout, not another local
  threshold sweep.
- A completed-negative route failure-mode audit found that the issue is not
  mostly missing raw evidence. For marine and CAMI3 samples3-5, 336/338 missed
  truth species are visible in best-diff raw rows, and all high-abundance
  misses at >=1% truth abundance are raw-visible (13/13). The existing HMP
  raw-table context has the same pattern for high-abundance misses (10/10).
  A limited raw-side retention sweep then found 59/78 simple rules with no
  route or sample F1 regression on the two newly scored routes. The best
  diagnostic rule, `raw_ani0.95_xny100_breadth0.05`, gives mean route pooled
  F1 delta +0.030664 and minimum route delta +0.005061, adding 130 TP and 53
  FP. This is a candidate implementation direction, not a default change.
- The broader cached raw-side retention replay extends that diagnostic to 43
  samples across 7 panels. It found 39/78 strict-pass rules with no sample or
  panel F1 regression. The best strict-pass rule,
  `raw_ani0.93_xny650_breadth0.20`, gives mean sample F1 delta +0.007355,
  minimum sample delta 0, no panel regressions, and adds 54 TP / 11 FP. Because
  rescued rows are assigned zero abundance in this replay, the result supports
  an in-wrapper raw-candidate retention implementation target but still does
  not justify a default change by itself.
- `scripts/minco_profile_calibrated.py` now has an opt-in
  `--candidate-surface-switch accession-ani93-xny650-br20` mode that implements
  the broad replay rule shape on raw accession rows: ANI >= 0.93, XnY >= 650,
  Ref_breadth >= 0.20, and no real-AF gate. The default candidate preset still
  uses `accession-ani90-xny100-br01-af70`. A one-sample cached table-mode smoke
  on Toy Mouse sample5 with zero candidate mass passed in 13.69 s wall time and
  2,815,608 KB peak RSS, adding 1 zero-mass surface row. This is implementation
  progress toward the next validation pass, not a new default decision.
- The full wrapper-output replay for
  `accession-current-or-ani93-xny650-br20` plus
  `--candidate-surface-max-called-species 250` is now complete across the same
  43-sample / 7-panel cached replay. Both zero-mass and
  `normalized-depth-alpha2` candidate abundance policies strict-pass with no
  sample or panel F1 regression versus the current default. The better
  abundance policy is `normalized-depth-alpha2`: mean sample F1 delta
  +0.006671, minimum sample delta 0, panel regressions 0, total TP/FP delta
  +30/+5, mean L1 delta -0.488058 pp, and mean Pearson delta +0.004763. The
  guard blocked candidate-surface additions in 6 emitted profiles, including
  the marine samples that regressed in the unguarded replay. This is the
  strongest opt-in candidate so far, but it remains validation-only because it
  is still the same cached panel family and not an independent holdout.
- An independent same-namespace HMP gastrooral sample1 holdout was then added
  from `https://frl.publisso.de/data/frl:6425518/gastrooral/sample_1.tar.gz`.
  The streamed nonzero-source read file is 4.2G; extraction took 5:30.77 with
  11264 KB peak RSS and `gzip -t` passed. GTDB source-abundance truth maps
  98.150316% of source mass (39 GTDB species), so this is release-grade for
  the source-abundance comparison. The raw default command, without a
  GTDB-accession candidate-surface sidecar beside the reference, scored F1
  0.886076 (TP/FP/FN 35/5/4), L1 18.244592 pp, Pearson 0.985778 in 3:38.76
  wall / 4000660 KB RSS. Replaying the selected current surface
  (`accession-ani90-xny100-br01-af70`) with an explicit GTDB accession-level
  candidate-surface taxmap rescued `s__Clostridium_F botulinum` without adding
  an FP, scoring F1 0.900000 (36/5/3), L1 17.455403 pp, Pearson 0.986396.
  The guarded `accession-current-or-ani93-xny650-br20` candidate with
  normalized-depth-alpha2 scored F1 0.888889 (36/6/3), L1 18.564655 pp,
  Pearson 0.985790; zero candidate mass kept the same F1 and had L1
  18.585816 pp. Sylph r232 scored F1 0.938272 (38/4/1), L1 8.099610 pp,
  Pearson 0.998654; Sylph sketch+profile took 4:18.02 wall and 27186576 KB
  peak RSS. The stricter no-realAF branch rescued the same true species but
  added `s__Bacteroides sp947646015` as an FP, so this holdout is negative for
  promoting that guarded branch to default.
- `scripts/minco_profile_calibrated.py` now disables candidate-surface rows
  unless `--candidate-surface-taxmap` is supplied or the normal `--taxmap`
  contains GTDB-style `s__` species labels. This prevents a non-GTDB taxmap
  from being silently used as an accession-level GTDB surface map. The practical
  deployment requirement is unchanged: downloadable MinCO GTDB reference
  bundles need to include a `candidate_surface_taxmap.tsv`/
  `accession_species_taxmap.tsv` sidecar, or users must set
  `MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP`.
- The GTDB r232 S1000 reference bundle under
  `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
  now has `candidate_surface_taxmap.tsv` copied from
  `/tmp/cami2_toymouse_current_default_20260626/gtdb_r232_species.taxmap.tsv`
  (388M, 1,348,479 rows, sha256
  `5b461a4f807600d0fee7691c76550b9d7646b04c23534e72ae1f7fea848746f8`).
  A plain `scripts/minco_profile_default.py` table-mode replay with `-r` set
  to that reference and no explicit `--candidate-surface-taxmap` discovered the
  sidecar and matched the explicit accession-surface profile in calls,
  abundance, ANI, and the single added `s__Clostridium_F botulinum` row.
- `scripts/minco_profile` is now the user-facing calibrated species-profile
  launcher. It is a thin executable shim over `scripts/minco_profile_default.py`
  and keeps the same selected default: `universal-auto-exact`,
  `--profile-preset candidate`, accession-level candidate-surface sidecar
  discovery, and normalized-depth-alpha2 candidate abundance. The old
  `scripts/minco_profile_default.py` command remains the compatibility
  implementation path. README, the user manual, wrapper help, smoke tests, and
  the current-default manifest now point users to `scripts/minco_profile` as
  the simplest no-manual-strategy command.
- `scripts/minco_profile --check-ref` now performs a wrapper-level preflight.
  It validates the reference sketch, reads input when supplied, species taxmap,
  model cache or training sidecar, candidate-surface taxmap for the selected
  candidate preset, MinCO binary, and model-cache metadata against the resolved
  `scope`, `train_pool`, and `filter_training_scope`.
- The same local GTDB r232 S1000 bundle is now self-contained for the selected
  default: it contains `species_taxmap.tsv`, `candidate_surface_taxmap.tsv`,
  `minco_profile_rf_hgb.train12.unfiltered.joblib`, and
  `minco_profile_defaults.tsv`. The first no-extra-flags profile attempt found
  a real packaging defect: the model cache was trained for `scope='bacteria'`
  while the parser default was `scope='all'`. The fix is a whitelisted
  packaged-default sidecar with `scope	bacteria`; the wrapper injects only
  allowed domain defaults from that file, and explicit user options still
  override it.
- A real no-extra-flags packaged smoke completed on
  `/tmp/toymouse_sample0_50000_reads.fq.gz` using only `scripts/minco_profile
  -r REF --reads READS -p16 -o OUT`. It wrote
  `/tmp/minco_packaged_default_smoke_20260630/toymouse50k.default.v2.tsv` with
  6 called rows, `calibrated_abundance` sum 1.0, `scope=bacteria`, and
  `profile_strategy=universal-auto-exact`. Runtime was 0:22.49 wall with
  2,754,736 KB max RSS.
- The repository is dirty, so this note records a stage decision before a clean
  release commit.

## Achievable Release Goal

The original research goal, a stable universal MinCO strategy that broadly
beats Sylph and closes the abundance-release gap, is still incomplete and
should remain a research target. The achievable release goal is narrower:
make the current MinCO profiler easy for a non-expert user to run, with one
default command, reference preflight, packaged sidecar defaults, documented
claim boundaries, and cached cross-panel evidence supporting the chosen MinCO
default against previous MinCO defaults.

`validate_achievable_release_goal.py` records this scope explicitly. On
2026-06-30 it passed 9/9 checks and wrote
`results/achievable_release_goal.tsv` plus `ACHIEVABLE_RELEASE_GOAL.md`.
The decision is `achievable_release_goal_supported`. The same gate also records
that the original broad goal remains separate:
`do_not_mark_broad_goal_complete`, because the release-grade holdout and broad
abundance claims remain expected gaps in
`results/universal_strategy_release_gate.tsv`.

## Next Experiment

Before calling the goal complete, expand the clean release-grade holdout bundle
beyond the current same-release/GTDB-species panels. The current action plan has
four release panels and four nonrelease panels. The next evidence target is not
another local threshold sweep; it is stronger same-namespace truth and
selected-default profile coverage for the nonrelease panels.

CAMI3 source-readmap samples0-2 are now release-grade under the accepted
exact-binomial fallback policy, and the samples3-5 extension route has now
been recovered and scored. That route is negative for promoting the refined
allocator and for any broad abundance claim: the external baseline wins F1 and
L1 on all three extension samples, and the refined allocator has zero measured
effect. The local source-profile fallback was also tested and rejected as a
release substitute because its in-scope GTDB truth mapping remains below
threshold. The HMP omitted-sample profile-pair recovery path is complete for
samples2/8/12/26/27 and is negative for a broad default-promotion claim. The
marine setup+assembly truth route now clears the mapped-truth threshold and
has selected-default same-namespace profile scoring, but that scoring is also
negative for promotion. A separate remaining-route audit now rules out the
current local plant, strain, mixed-readiness, HMP-omitted, marine, and CAMI3
extension paths as immediate release-gap closures. The newest failure-mode
audit changes the next algorithmic target: implement or emulate an in-wrapper
raw-candidate retention side channel and validate it across the existing
release panels before searching more output-only thresholds. The first wrapper
rule to test is now `raw_ani0.93_xny650_breadth0.20`, because it is the best
strict-pass rule in the broader 43-sample replay. Keep
`raw_ani0.95_xny650_breadth0.05` and `raw_ani0.93_xny650_breadth0.30` as
lower-FP fallbacks if implementation-time scoring shows too many extra calls
or abundance side effects. The first opt-in implementation exists as
`--candidate-surface-switch accession-ani93-xny650-br20`. The follow-up
wrapper-output replay using combined mode
`accession-current-or-ani93-xny650-br20` is complete across the same 43-sample
panel. The first guard candidate, `--candidate-surface-max-called-species 250`,
blocks the unguarded marine regressions in targeted implementation smokes and
strict-passes the full wrapper-output replay across all 43 cached samples. The
new HMP gastrooral sample1 independent holdout changes the deployment lesson:
the selected current surface helps only when a GTDB accession-level
candidate-surface taxmap is available, while the guarded no-realAF branch adds
an FP and worsens L1. The local S1000 reference sidecar has now been installed
and discovery-validated, so the default usability issue is fixed for this
local reference bundle. Therefore keep the guarded raw-retention surface opt-in
and do not change the default preset thresholds. The next required experiment
is not another promotion replay of this exact guarded candidate; it is a new
strategy that can recover the remaining missed truth without extra
FPs/abundance drift, or additional independent same-namespace holdouts showing
the selected current surface remains robust when the sidecar is packaged.
Treat the current strainmadness files as non-release until a stronger truth
mapping is found; and test abundance algorithms that combine better
read/context-level allocation with high-abundance missed-call recovery beyond
the rejected selected_group_species2 edge-EM policy rather than aggregate
per-reference totals. Before tuning another output-only rescue threshold, add
or emulate an in-pass raw-candidate side channel and score candidate-retention
rules against both HMP and CAMI3. The emitted-profile wrapper replay validates
the zero-mass rescue mechanics but does not expose the HMP gastrooral raw-cache
rescues; the taxid-collapse audit shows the side channel must preserve
accession-level or GTDB-species-level candidates rather than only taxid-level
aggregates. The final decision table should keep reporting F1, L1/Pearson, ANI
diagnostics, runtime, and memory for MinCO versus Sylph. For speed, the next
deeper implementation target is deciding whether an
adaptive fast/exact trigger can skip exact mode without losing the panels where
exact split helps. The current sidecar-policy audit narrows that target: the
next implementation should test a stronger preflight exact-need predictor or a
cheaper conditional exact sidecar, because full block-mode features identify
exact candidates only after it is too late to decide eager sidecar for that
same pass, and first50k/first200k prefix preflight was too sparse on HMP
gastrooral. Candidate-only exact rerun should not be implemented from aggregate
profile outputs alone; it would need context/read-level competitor closure or
an in-pass full-candidate exact side channel.
