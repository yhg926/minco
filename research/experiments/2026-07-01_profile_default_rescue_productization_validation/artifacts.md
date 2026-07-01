# Artifacts

- `run_default_entrypoint_validation.py`: runs `scripts/minco_profile` default and opt-out commands.
- `run_raw_default_smoke.py`: builds a tiny temporary raw-input fixture and runs the default launcher end to end.
- `results/default_profile_runtime.tsv`: exact commands, elapsed time, and peak RSS.
- `results/default_profile_scores.tsv`: per-sample GTDB scores.
- `results/default_profile_panel_summary.tsv`: panel-level metric rollup.
- `results/default_profile_overall_summary.tsv`: overall rollup copied to `summary.tsv`.
- `results/default_profile_added.tsv`: rows added by the default rescue.
- `results/raw_default_smoke_summary.tsv`: raw-input smoke pass/fail, default switches, and runtime.
- `/tmp/minco_profile_default_rescue_productization_validation_20260701`: temporary large profile outputs.
- `/tmp/minco_default_raw_smoke_20260701`: temporary raw-input smoke fixture and profile output.
