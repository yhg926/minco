# External CAMI Strain Holdout

Date: 2026-06-21

## Question

Does the current minco readwise species-call calibration generalize to a new CAMI sample type that was not used in the marine/toy-gut/plant calibration panel?

## Dataset

CAMI II strain-madness short-read samples 0-2 were used as an external holdout. CAMI describes strain-madness as simulated short- and long-read shotgun metagenomes with large strain-level variation.

Download source:

- CAMI datasets page: https://cami-challenge.org/datasets/
- Strain-madness download root: https://frl.publisso.de/data/frl:6425521/strain/

Local data root:

- `/mnt/new3T/minco_cami2_strain_20260621`

The strain gold profiles have blank `@SampleID:` headers. Species-level scoring used the species taxid set from `taxonomic_profile_0.txt`, `taxonomic_profile_1.txt`, and `taxonomic_profile_2.txt`; each contains 20 species taxids. Three gold species were absent from the current GTDB/taxmap mapping: `1235803`, `481743`, and `97137`.

## Methods

Training panel:

- `research/experiments/2026-06-21_minco_multisample_call_calibration/manifest_marine_toy0_2_plant0_2.tsv`
- 9 samples: CAMI marine 0-2, toy gut 0-2, plant 0-2.

External test panel:

- `manifest_strain0_2.tsv`
- 3 samples: CAMI II strain-madness short-read 0-2.

minco settings:

- Reference: GTDBr232 S1000 annotated minco sketch.
- `minco ani -p16 --qraw FASTQ --query-density ref --abundance-est depth --readwise-profile-only --readwise-ani zip-aaf -m0 -f0 -n0 -t0`
- Assignment modes tested: `best-diff-unique` and `best-diff-split`.

Sylph baseline:

- `/home/ubuntu/yihuiguang/bin/sylph`
- Database: `/mnt/new3T/sylph_db/gtdb-r226-c200-dbv1.syldb`

Modeling:

- Train only on the 9-sample panel.
- Tune model thresholds only on the 9-sample training panel.
- Test once on strain0-2.
- Features intentionally exclude dataset identity.
- Models compared: logistic regression, RF, HGB, RF/HGB probability average.

## Results

External strain holdout, mean over 3 samples:

| Method | Mean precision | Mean recall | Mean F1 | Mean FP | Mean FN | Mean FP+FN |
|---|---:|---:|---:|---:|---:|---:|
| Sylph default | 0.5504 | 0.5167 | 0.5274 | 9.00 | 9.67 | 18.67 |
| train9 RF/HGB average | 0.5433 | 0.4167 | 0.4656 | 7.33 | 11.67 | 19.00 |
| minco split direct | 0.6007 | 0.3833 | 0.4632 | 5.33 | 12.33 | 17.67 |
| train9 HGB | 0.4977 | 0.4000 | 0.4415 | 8.00 | 12.00 | 20.00 |
| train9 RF | 0.5180 | 0.3333 | 0.4011 | 6.33 | 13.33 | 19.67 |
| minco unique direct | 0.5985 | 0.2667 | 0.3671 | 3.67 | 14.67 | 18.33 |
| minco unique relaxed | 0.2644 | 0.5167 | 0.3475 | 29.33 | 9.67 | 39.00 |
| train9 logistic regression | 0.2583 | 0.3167 | 0.2824 | 18.33 | 13.67 | 32.00 |

Per-sample metrics are in `summary.tsv`.

Candidate availability diagnostic:

| Sample | Gold species | Gold in taxmap | Gold in minco joined candidates | All bacterial candidate taxa |
|---|---:|---:|---:|---:|
| strain0 | 20 | 17 | 17 | 18158 |
| strain1 | 20 | 17 | 17 | 19873 |
| strain2 | 20 | 17 | 17 | 17962 |

Threshold oracle diagnostic, not a valid benchmark:

- RF/HGB average improves only from mean F1 0.4656 to 0.4784 when the threshold is tuned directly on the strain holdout.
- RF improves to 0.4759, HGB to 0.4627.

This means the external-test loss is not mostly a threshold calibration problem; the feature ranking/model itself is not generalizing well enough to strain-madness.

## Follow-Up Simple Rules

After the external test, I checked whether a simpler deployable rule could be a
better path than RF/HGB. The useful rule was:

- call if either unique or split direct evidence fires;
- then optionally suppress generic species labels whose names contain markers
  such as `uncultured` or `bacterium`.

This is recorded as `u_or_s_direct_validtax` in `simple_filter_diagnostics.tsv`.

On the external strain holdout:

| Method | Mean precision | Mean recall | Mean F1 | Mean FP | Mean FN | Mean FP+FN |
|---|---:|---:|---:|---:|---:|---:|
| Sylph default | 0.5504 | 0.5167 | 0.5274 | 9.00 | 9.67 | 18.67 |
| Sylph with same generic-name filter | 0.5966 | 0.5167 | 0.5493 | 7.33 | 9.67 | 17.00 |
| minco `u_or_s_direct_validtax` | 0.6677 | 0.4000 | 0.4966 | 4.00 | 12.00 | 16.00 |
| minco RF/HGB valid taxon | 0.5957 | 0.4167 | 0.4876 | 5.67 | 11.67 | 17.33 |
| minco `s_direct_validtax` | 0.6597 | 0.3833 | 0.4826 | 4.00 | 12.33 | 16.33 |

The filter helps minco precision and FP+FN, but it does not beat Sylph by F1.
On the nine-sample training panel, raw `u_or_s_direct` remained better than
`u_or_s_direct_validtax` by F1, so this should remain an optional/reporting
filter rather than a default.

I also tested genus-capped relaxed rescue rules: start from lower-ANI or relaxed
split evidence, then keep only the top K species per genus by support,
abundance, breadth, or ANI. The training-preferred rule was still raw
`u_or_s_direct_validtax`; genus-capped rescue did not improve the strain
holdout. The missed gold taxa mostly have substantial split support but ANI
below the current direct threshold, while the remaining false positives are
nearby species within gold genera/families. That points to reference/taxonomy
resolution rather than a simple threshold-only fix.

## Speed And Memory

minco on each 2.0 GB compressed strain FASTQ:

- Unique assignment: 47.8-48.8 s wall time, 3.13-3.16 GB peak RSS.
- Split assignment: 48.8-49.4 s wall time, 3.41-3.47 GB peak RSS.

Sylph:

- Sketch: about 21.6-21.7 s wall time, 0.50 GB peak RSS.
- Profile: 37.1-70.0 s wall time, 19.28 GB peak RSS.

## Conclusion

The current multi-sample learned model does not yet generalize better than Sylph on a new strain-heavy CAMI sample type. The RF/HGB average is slightly better than single RF or HGB on this holdout, but the gain is small and it still does not beat Sylph. The simpler `minco_split_direct` rule is nearly tied with the RF/HGB average by F1 and has the best FP+FN, so the current ensemble complexity is not justified as a default.

The useful signal is that minco candidates already contain 17 of 20 gold species for each strain sample. The next improvement should focus on robust candidate filtering/ranking across domains, not sketch recovery.

Follow-up diagnostics strengthen the same conclusion: a simple `u_or_s_direct`
rule with generic-name suppression is more attractive than the current RF/HGB
ensemble, but the remaining gap to Sylph is recall and species-resolution in
strain-heavy communities.

## Next Steps

1. Keep strain-madness as a locked external holdout while developing filters on marine/toy/plant.
2. Prefer a simple deployable model or rule until external validation shows a consistent gain.
3. Add more non-training CAMI domains before claiming broad performance, especially clinical pathogen and human-body-site datasets if their gold profiles can be mapped cleanly.
4. Investigate features for the strain FNs and FPs, with special attention to breadth/depth thresholds and strain-level shared-context ambiguity.
