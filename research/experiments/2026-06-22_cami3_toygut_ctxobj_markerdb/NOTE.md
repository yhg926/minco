# CAMI3 ToyGut Ctx+Obj Markerdb Comparison

Date: 2026-06-22

## Question

Compare the new S2000 ctx+obj markerdb against the current S2000 ctx-only markerdb and Sylph on the CAMI3 ToyGut samples.

## Inputs

- Workspace: `/home/ubuntu/yihuiguang/tools/KSSD3mini`
- Binary: `bin/minco`, version `minco 0.1`
- Binary sha256: `40af778fb7b9a122e2db513a27d652ef52caeffb7196bba4ea1f936258c65f2b`
- Ctx+obj markerdb: `/tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker`
- Ctx-only baseline outputs: `/tmp/cami3_toygut_current_minco_vs_sylph_20260622`
- Samples: CAMI3 ToyGut sample0, sample1, sample2.
- Truth: bacterial species rows from CAMI taxonomic profiles, scored by NCBI species taxid.

Only bacterial species were scored, matching the earlier CAMI3 ToyGut comparison. Viral, plasmid, and eukaryotic truth rows were excluded because the references tested here are GTDB references.

## Methods

The ctx+obj MinCO profiles used the same readwise command family as the current ctx-only baseline:

```bash
bin/minco ani -p16 -r /tmp/gtdb232_s2000_ctxobj_marker_20260622/sketch_T_S2000_aaf003_dedup_ctxobjmarker \
  --qraw <reads.fq.gz> \
  --query-density ref --abundance-est depth --readwise-profile-only \
  --readwise-assign best-diff-split --readwise-ani naive \
  --readwise-ctx-filter product-topfrac-median --readwise-fake-threshold 0.25 \
  -m0 -f0 -n0 -t0
```

Scored methods:

- `minco_ctx_only_current`: current S2000 ctx-only markerdb active gate.
- `minco_ctxobj_active`: ctx+obj markerdb with active gate, `Reliable_ztp_af >= 0.40`.
- `minco_ctxobj_relaxed`: ctx+obj markerdb with Toy Mouse relaxed gate, `Reliable_ztp_af >= 0.35` and `Reliable_Ref_breadth >= 0.03`.
- `sylph`: existing Sylph GTDB r226 c200 profiles from the previous comparison.

The common active gate components were:

```text
XnY_ctx >= 15
ANI_naive_calc > 0.95
active delta rule unchanged:
  if Reliable_Ref_hit_mean_depth > 3 and Reliable_Ref_hit_depth_variance / Reliable_Ref_hit_mean_depth > 50:
      require ANI_naive_calc - ANI_from_Reliable_ztp_af < 0.03
```

## Presence/Absence Results

Per-sample bacterial species scores:

| sample | method | truth | pred | TP | FP | FN | precision | recall | F1 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | minco_ctx_only_current | 118 | 61 | 49 | 12 | 69 | 0.803279 | 0.415254 | 0.547486 |
| 0 | minco_ctxobj_active | 118 | 68 | 53 | 15 | 65 | 0.779412 | 0.449153 | 0.569892 |
| 0 | minco_ctxobj_relaxed | 118 | 61 | 47 | 14 | 71 | 0.770492 | 0.398305 | 0.525140 |
| 0 | sylph | 118 | 70 | 57 | 13 | 61 | 0.814286 | 0.483051 | 0.606383 |
| 1 | minco_ctx_only_current | 125 | 70 | 55 | 15 | 70 | 0.785714 | 0.440000 | 0.564103 |
| 1 | minco_ctxobj_active | 125 | 81 | 65 | 16 | 60 | 0.802469 | 0.520000 | 0.631068 |
| 1 | minco_ctxobj_relaxed | 125 | 71 | 56 | 15 | 69 | 0.788732 | 0.448000 | 0.571429 |
| 1 | sylph | 125 | 77 | 63 | 14 | 62 | 0.818182 | 0.504000 | 0.623762 |
| 2 | minco_ctx_only_current | 127 | 57 | 51 | 6 | 76 | 0.894737 | 0.401575 | 0.554348 |
| 2 | minco_ctxobj_active | 127 | 69 | 59 | 10 | 68 | 0.855072 | 0.464567 | 0.602041 |
| 2 | minco_ctxobj_relaxed | 127 | 60 | 47 | 13 | 80 | 0.783333 | 0.370079 | 0.502674 |
| 2 | sylph | 127 | 66 | 56 | 10 | 71 | 0.848485 | 0.440945 | 0.580311 |

Mean over three samples:

| method | precision | recall | F1 | mean TP | mean FP | mean FN |
|---|---:|---:|---:|---:|---:|---:|
| minco_ctx_only_current | 0.827910 | 0.418943 | 0.555312 | 51.67 | 11.00 | 71.67 |
| minco_ctxobj_active | 0.812318 | 0.477906 | 0.601000 | 59.00 | 13.67 | 64.33 |
| minco_ctxobj_relaxed | 0.780853 | 0.405461 | 0.533081 | 50.00 | 14.00 | 73.33 |
| sylph | 0.826984 | 0.475999 | 0.603485 | 58.67 | 12.33 | 64.67 |

Totals over three samples:

| method | TP | FP | FN | predicted species | truth species |
|---|---:|---:|---:|---:|---:|
| minco_ctx_only_current | 155 | 33 | 215 | 188 | 370 |
| minco_ctxobj_active | 177 | 41 | 193 | 218 | 370 |
| minco_ctxobj_relaxed | 150 | 42 | 220 | 192 | 370 |
| sylph | 176 | 37 | 194 | 213 | 370 |

## Abundance

Mean abundance correlation over the three samples, with predictions renormalized to sum to 1:

| method | Pearson | Spearman | MAE percentage points | L1 percentage points |
|---|---:|---:|---:|---:|
| minco_ctx_only_current | 0.884011 | 0.545709 | 0.417232 | 51.237138 |
| minco_ctxobj_active | 0.894380 | 0.571113 | 0.379254 | 46.613987 |
| minco_ctxobj_relaxed | 0.896439 | 0.643682 | 0.372436 | 45.756950 |
| sylph | 0.923617 | 0.668254 | 0.225281 | 27.789743 |

Without renormalizing predictions, Sylph also had much better recovered abundance mass on bacterial truth (`mean_pred_sum_on_truth=0.758908`) than ctx-only MinCO (`0.384135`) or ctx+obj active MinCO (`0.302476`).

## Runtime

Ctx+obj MinCO runtime:

| sample | elapsed | max RSS KB | exit |
|---:|---:|---:|---:|
| 0 | 2:11.33 | 5944032 | 0 |
| 1 | 2:14.92 | 6007500 | 0 |
| 2 | 2:14.86 | 5890656 | 0 |

This is similar wall time to ctx-only MinCO but uses more memory. Previous ctx-only RSS was about `4.1G`; ctx+obj was about `5.9-6.0G`.

## Conclusion

Ctx+obj active is clearly better than the current ctx-only MinCO gate on CAMI3 ToyGut:

```text
mean F1: 0.555312 -> 0.601000
total TP/FN: 155/215 -> 177/193
```

It is nearly tied with Sylph by mean F1, but still slightly lower:

```text
ctx+obj active mean F1 = 0.601000
Sylph mean F1          = 0.603485
```

Ctx+obj active has one more total TP and one fewer total FN than Sylph (`177 TP, 193 FN` vs `176 TP, 194 FN`), but also has more FP (`41` vs `37`), which keeps its mean precision and mean F1 slightly below Sylph.

The Toy Mouse relaxed rule does not transfer to CAMI3. It improves renormalized abundance metrics among MinCO variants but hurts presence/absence calling (`mean F1=0.533081`).

## Caveats

- This is CAMI bacterial NCBI species-taxid scoring, not GTDB species scoring.
- The ctx+obj markerdb was built from the same S2000 deduplicated reference set used in the Toy Mouse experiment.
- The active gate was tuned on CAMI2 Toy Mouse sample0, so this CAMI3 improvement is encouraging but still needs broader validation.
- Sylph uses GTDB r226 c200 while MinCO uses the local GTDB r232-derived reference/metadata workflow; accession mapping uses the same metadata helper as the previous CAMI3 comparison.
