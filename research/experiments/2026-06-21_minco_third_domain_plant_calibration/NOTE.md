# minco Third-Domain Plant Calibration

Date: 2026-06-21

## Question

After the marine + toy gut calibration, test whether the minco readwise call model generalizes to a third CAMI domain and whether a deployable, no-dataset-label strategy can still match or beat Sylph.

## Data

- Existing six samples: CAMI marine short-read samples 0-2 and CAMI III toy human gut samples 0-2.
- Added third domain: CAMI II plant-associated short-read samples 0-2.
- Plant source URL: `https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/`
- Plant setup/gold profiles: `/mnt/new3T/minco_cami2_plant_20260621/simulation_short_read/`
- Plant reads:
  - `/mnt/new3T/minco_cami2_plant_20260621/sample_0/.../anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_plant_20260621/sample_1/.../anonymous_reads.fq.gz`
  - `/mnt/new3T/minco_cami2_plant_20260621/sample_2/.../anonymous_reads.fq.gz`

Plant profiles include bacteria, viruses, eukaryotes, and unidentified mass. This benchmark uses `scope=bacteria` for plant, matching the GTDB/Sylph bacterial reference comparison.

Plant bacterial gold counts:

| sample | all species | bacterial species | bacterial percentage |
| --- | ---: | ---: | ---: |
| plant0 | 86 | 45 | 25.8807 |
| plant1 | 87 | 47 | 30.6912 |
| plant2 | 78 | 37 | 24.7759 |

## Methods

Reference databases:

- minco S1000 GTDBr232 + RefSeq virus DB: `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno`
- Sylph GTDB r226 DB: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`

minco plant readwise outputs used both assignment modes:

- `--readwise-assign best-diff-unique`
- `--readwise-assign best-diff-split`
- `--readwise-ani zip-aaf`
- `-m0 -f0 -n0 -t0` to keep unfiltered candidate tables for model calibration.

Calibration design:

- Leave-one-sample-out across 9 samples: marine0-2, toy0-2, plant0-2.
- Removed dataset identity features (`is_marine`, `is_toy_gut`, `is_bacteria_scope`) from the learned model because they are not deployable signals.
- Added deployable rule-evidence features:
  - `u_direct_call`, `u_relaxed_call`
  - `s_direct_call`, `s_relaxed_call`
- Evaluated individual RF/HGB models and an averaged probability ensemble `loso_rf_hgb_avg`.

## Accuracy Results

Final no-leak, rule-feature, ensemble result:

| method | samples | mean precision | mean recall | mean F1 | mean FP | mean FN |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| minco `loso_rf_hgb_avg` | 9 | 0.7333 | 0.7001 | 0.6976 | 30.67 | 33.00 |
| minco `loso_rf` | 9 | 0.7756 | 0.6424 | 0.6771 | 23.00 | 38.44 |
| Sylph default | 9 | 0.7190 | 0.6646 | 0.6763 | 32.00 | 36.44 |
| minco unique direct | 9 | 0.7196 | 0.5546 | 0.5980 | 23.44 | 46.11 |

Dataset-level F1:

| dataset | minco `rf_hgb_avg` | Sylph default | minco unique direct |
| --- | ---: | ---: | ---: |
| marine | 0.8518 | 0.8323 | 0.8612 |
| plant | 0.6143 | 0.5799 | 0.5013 |
| toy gut | 0.6267 | 0.6167 | 0.4314 |

The ensemble improves over Sylph on mean F1 across the 9-sample panel, and it improves both plant and toy gut relative to strict minco direct calls. Unique direct remains slightly better on marine alone.

## Timing

Plant minco readwise scan timings, one pass per assignment mode:

| sample/pass | wall time | max RSS |
| --- | ---: | ---: |
| plant0 unique | 1:56.23 | 3.62 GB |
| plant0 split | 1:55.84 | 4.07 GB |
| plant1 unique | 1:53.71 | 3.65 GB |
| plant1 split | 1:55.55 | 4.11 GB |
| plant2 unique | 1:53.84 | 3.61 GB |
| plant2 split | 1:55.70 | 4.06 GB |

Sylph plant timings:

| sample | sketch wall/RSS | profile wall/RSS |
| --- | ---: | ---: |
| plant0 | 0:54.65 / 0.99 GB | 2:27.76 / 19.53 GB |
| plant1 | 0:57.64 / 0.99 GB | 1:35.08 / 19.53 GB |
| plant2 | 0:57.27 / 0.99 GB | 1:26.56 / 19.53 GB |

## Interpretation

The earlier six-sample “beats Sylph” result was partly inflated by dataset identity features. After removing those non-deployable signals, RF alone nearly tied Sylph. Adding strict/relaxed rule-evidence features and averaging RF+HGB probabilities produced a cleaner improvement: F1 0.6976 vs Sylph 0.6763.

The main improvement mechanism is not a new ANI correction. It is call-level calibration that combines breadth, depth dispersion, XnY support, min align fraction, and rule-evidence features to recover low/medium abundance true positives while avoiding the large FP explosion from relaxed direct thresholds.

## Caveats

- Plant scoring is bacteria-only because the current comparison reference is GTDB/Sylph bacterial-focused.
- Sylph uses GTDB r226 while minco uses GTDBr232 plus RefSeq virus annotations; exact taxonomic naming/release differences can affect TP/FP/FN.
- This is still a 9-sample LOSO panel, not a production-trained classifier.

## Artifacts

Small archived files are in this directory. Large generated tables remain under:

`/mnt/new3T/minco_cami2_plant_20260621/calibration_marine_toy_plant_20260621/model_noleak_rules_ensemble/`
