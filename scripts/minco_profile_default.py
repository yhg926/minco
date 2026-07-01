#!/usr/bin/env python3
"""Default MinCO species-profile launcher.

This is the no-manual-strategy entry point for the calibrated profile default.
It delegates to ``minco_profile_calibrated.py`` and fills stable
deployment paths from packaged sidecars or environment variables when the
caller does not provide them explicitly:

  MINCO_PROFILE_TAXMAP
  MINCO_PROFILE_TRAIN_FEATURES
  MINCO_PROFILE_MODEL_CACHE
  MINCO_PROFILE_MINCO
  MINCO_PROFILE_PRESET
  MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP

The selected strategy remains explicit in the delegated argv so logs and
process tables show which default was used.
"""

from __future__ import annotations

import os
import shlex
import sys
from pathlib import Path
from typing import Mapping, Optional, Sequence

try:
    from scripts import minco_profile_calibrated as calibrated
except ModuleNotFoundError:  # direct execution from scripts/
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import minco_profile_calibrated as calibrated  # type: ignore


DEFAULT_STRATEGY = "universal-auto-exact"
CURRENT_PRESET = "current"
CANDIDATE_PRESET = "candidate"
DEFAULT_PRESET = CANDIDATE_PRESET
ENV_TAXMAP = "MINCO_PROFILE_TAXMAP"
ENV_TRAIN_FEATURES = "MINCO_PROFILE_TRAIN_FEATURES"
ENV_MODEL_CACHE = "MINCO_PROFILE_MODEL_CACHE"
ENV_MINCO = "MINCO_PROFILE_MINCO"
ENV_PROFILE_PRESET = "MINCO_PROFILE_PRESET"
ENV_CANDIDATE_SURFACE_TAXMAP = "MINCO_PROFILE_CANDIDATE_SURFACE_TAXMAP"
WRAPPER_PROFILE_PRESET_OPT = "--profile-preset"
WRAPPER_NO_PROFILE_RESCUE_OPT = "--no-profile-rescue"
WRAPPER_PREFLIGHT_OPTS = ("--check-ref", "--preflight")
VALID_PROFILE_PRESETS = {CURRENT_PRESET, CANDIDATE_PRESET}
MODEL_CACHE_NAMES = [
    "minco_profile_rf_hgb.train12.unfiltered.joblib",
    "minco_profile_rf_hgb.train12.joblib",
    "minco_profile_rf_hgb.joblib",
]
PROFILE_DEFAULT_FILENAMES = [
    "minco_profile_defaults.tsv",
    "profile_defaults.tsv",
    "minco_profile.defaults.tsv",
]
PROFILE_DEFAULT_OPTIONS = {
    "scope": "--scope",
    "train_pool": "--train-pool",
}
PROFILE_DEFAULT_FLAGS = {
    "filter_training_scope": "--filter-training-scope",
}


def default_help_banner() -> str:
    return f"""\
MinCO calibrated species-profile default

This launcher is the recommended no-manual-strategy entry point for calibrated
species profiling. It delegates to scripts/minco_profile_calibrated.py with
--strategy {DEFAULT_STRATEGY}.

Typical command:
  scripts/minco_profile -r ref.minco --reads reads.fastq.gz -p16 -o profile.tsv

Preflight only:
  scripts/minco_profile --check-ref -r ref.minco --reads reads.fastq.gz

Compatibility command:
  scripts/minco_profile_default.py -r ref.minco --reads reads.fastq.gz -p16 -o profile.tsv

Required inputs:
  -r/--ref REF.minco and --reads READS.fastq.gz, or precomputed
  --unique-table/--split-table inputs
  -o/--out output.tsv
  --taxmap TAXMAP, {ENV_TAXMAP}, or a taxmap sidecar beside --ref
  --model-cache MODEL.joblib / {ENV_MODEL_CACHE}, or
  --train-features DIR / --train-table TSV, {ENV_TRAIN_FEATURES}, or a
  joined_feature_training sidecar beside --ref
  Packaged references may also include minco_profile_defaults.tsv for
  domain-specific defaults such as scope=bacteria.

Optional deployment default:
  --minco BIN, or {ENV_MINCO}

Optional wrapper preset:
  {WRAPPER_PROFILE_PRESET_OPT} current|candidate, or {ENV_PROFILE_PRESET}
  candidate is the selected default. It enables strict split-evidence profile
  rescue, candidate-surface calls, abundance-only reliability guards, plus
  normalized-depth candidate abundance. current reproduces the previous
  calibrated default without those candidate additions. Use
  {WRAPPER_NO_PROFILE_RESCUE_OPT} to keep the candidate preset but disable
  profile-rescue additions. The candidate preset can discover
  candidate_surface_taxmap.tsv beside --ref, or use
  {ENV_CANDIDATE_SURFACE_TAXMAP}.

For direct AMR/gene/virus/mixed-domain profiling and read tracking, use the C
subcommand `minco profile`. For calibrated species profiling, use this launcher.
"""


def has_option(argv: Sequence[str], long_names: Sequence[str]) -> bool:
    for token in argv:
        for name in long_names:
            if token == name or token.startswith(f"{name}="):
                return True
    return False


def help_requested(argv: Sequence[str]) -> bool:
    return any(token in {"-h", "--help"} for token in argv)


def option_value(
    argv: Sequence[str],
    long_names: Sequence[str],
    short_names: Sequence[str] = (),
) -> Optional[str]:
    for i, token in enumerate(argv):
        for name in long_names:
            if token == name and i + 1 < len(argv):
                return argv[i + 1]
            if token.startswith(f"{name}="):
                return token.split("=", 1)[1]
        if token in short_names and i + 1 < len(argv):
            return argv[i + 1]
    return None


def option_values(
    argv: Sequence[str],
    long_names: Sequence[str],
    short_names: Sequence[str] = (),
) -> list[str]:
    values: list[str] = []
    for i, token in enumerate(argv):
        for name in long_names:
            if token == name and i + 1 < len(argv):
                values.append(argv[i + 1])
            elif token.startswith(f"{name}="):
                values.append(token.split("=", 1)[1])
        if token in short_names and i + 1 < len(argv):
            values.append(argv[i + 1])
    return values


def extract_wrapper_option(argv: Sequence[str], long_name: str) -> tuple[Optional[str], list[str]]:
    value: Optional[str] = None
    cleaned: list[str] = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == long_name:
            if i + 1 >= len(argv):
                raise SystemExit(f"{long_name} requires a value")
            value = argv[i + 1]
            i += 2
            continue
        if token.startswith(f"{long_name}="):
            value = token.split("=", 1)[1]
            i += 1
            continue
        cleaned.append(token)
        i += 1
    return value, cleaned


def extract_wrapper_flags(argv: Sequence[str], flags: Sequence[str]) -> tuple[bool, list[str]]:
    found = False
    cleaned: list[str] = []
    for token in argv:
        if token in flags:
            found = True
            continue
        cleaned.append(token)
    return found, cleaned


def first_existing_file(paths: Sequence[Path]) -> Optional[Path]:
    seen: set[Path] = set()
    for path in paths:
        normalized = path.expanduser()
        if normalized in seen:
            continue
        seen.add(normalized)
        if normalized.is_file():
            return normalized
    return None


def first_training_dir(paths: Sequence[Path]) -> Optional[Path]:
    seen: set[Path] = set()
    for path in paths:
        normalized = path.expanduser()
        if normalized in seen:
            continue
        seen.add(normalized)
        if normalized.is_dir() and (normalized / "train.joined_features.tsv").is_file():
            return normalized
    return None


def first_model_cache(paths: Sequence[Path]) -> Optional[Path]:
    seen: set[Path] = set()
    for path in paths:
        normalized = path.expanduser()
        if normalized in seen:
            continue
        seen.add(normalized)
        if normalized.is_file():
            return normalized
    return None


def unique_glob(root: Path, patterns: Sequence[str], require_train: bool = False) -> list[Path]:
    if not root.expanduser().is_dir():
        return []
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(sorted(root.expanduser().glob(pattern)))
    if require_train:
        matches = [path for path in matches if (path / "train.joined_features.tsv").is_file()]
    else:
        matches = [path for path in matches if path.is_file()]
    return matches if len(matches) == 1 else []


def discover_taxmap(argv: Sequence[str]) -> Optional[Path]:
    ref_value = option_value(argv, ["--ref"], ["-r"])
    if not ref_value:
        return None
    ref = Path(ref_value)
    parent = ref.parent if str(ref.parent) else Path(".")
    stem = ref.stem
    name = ref.name
    candidates: list[Path] = []
    if ref.expanduser().is_dir():
        candidates.extend(
            [
                ref / "species_taxmap.tsv",
                ref / "ref.species_taxmap.tsv",
                ref / "taxmap.tsv",
                ref / f"{stem}.species_taxmap.tsv",
                ref / f"{stem}.taxmap.tsv",
            ]
        )
        candidates.extend(unique_glob(ref, ["*.species_taxmap.tsv", "*.taxmap.tsv"]))
    candidates.extend(
        [
            parent / f"{name}.species_taxmap.tsv",
            parent / f"{stem}.species_taxmap.tsv",
            parent / f"{name}.taxmap.tsv",
            parent / f"{stem}.taxmap.tsv",
            parent / "species_taxmap.tsv",
            parent / "ref.species_taxmap.tsv",
            parent / "taxmap.tsv",
        ]
    )
    candidates.extend(unique_glob(parent, ["*.species_taxmap.tsv", "*.taxmap.tsv"]))
    return first_existing_file(candidates)


def discover_train_features(argv: Sequence[str]) -> Optional[Path]:
    ref_value = option_value(argv, ["--ref"], ["-r"])
    if not ref_value:
        return None
    ref = Path(ref_value)
    parent = ref.parent if str(ref.parent) else Path(".")
    stem = ref.stem
    name = ref.name
    dir_names = [
        "joined_feature_training",
        "minco_profile_training",
        "profile_training",
        "training",
        "calibration",
    ]
    candidates: list[Path] = []
    if ref.expanduser().is_dir():
        candidates.append(ref)
        candidates.extend(ref / dirname for dirname in dir_names)
        candidates.extend(unique_glob(ref, ["*joined_feature_training", "*profile_training"], require_train=True))
    candidates.extend(
        [
            parent / f"{name}.joined_feature_training",
            parent / f"{stem}.joined_feature_training",
            parent / f"{name}.profile_training",
            parent / f"{stem}.profile_training",
        ]
    )
    candidates.extend(parent / dirname for dirname in dir_names)
    candidates.append(parent)
    candidates.extend(unique_glob(parent, ["*joined_feature_training", "*profile_training"], require_train=True))
    return first_training_dir(candidates)


def discover_model_cache(argv: Sequence[str]) -> Optional[Path]:
    ref_value = option_value(argv, ["--ref"], ["-r"])
    if not ref_value:
        return None
    ref = Path(ref_value)
    parent = ref.parent if str(ref.parent) else Path(".")
    stem = ref.stem
    name = ref.name
    candidates: list[Path] = []
    training = discover_train_features(argv)
    if training is not None:
        candidates.extend(training / cache_name for cache_name in MODEL_CACHE_NAMES)
    if ref.expanduser().is_dir():
        candidates.extend(ref / cache_name for cache_name in MODEL_CACHE_NAMES)
        candidates.extend(ref / "joined_feature_training" / cache_name for cache_name in MODEL_CACHE_NAMES)
        candidates.extend(ref / "minco_profile_training" / cache_name for cache_name in MODEL_CACHE_NAMES)
    candidates.extend(parent / f"{name}.{cache_name}" for cache_name in MODEL_CACHE_NAMES)
    candidates.extend(parent / f"{stem}.{cache_name}" for cache_name in MODEL_CACHE_NAMES)
    candidates.extend(parent / cache_name for cache_name in MODEL_CACHE_NAMES)
    return first_model_cache(candidates)


def discover_candidate_surface_taxmap(argv: Sequence[str]) -> Optional[Path]:
    ref_value = option_value(argv, ["--ref"], ["-r"])
    if not ref_value:
        return None
    ref = Path(ref_value)
    parent = ref.parent if str(ref.parent) else Path(".")
    stem = ref.stem
    name = ref.name
    candidate_names = [
        "candidate_surface_taxmap.tsv",
        "accession_species_taxmap.tsv",
        "accession_taxmap.tsv",
        "ref.candidate_surface_taxmap.tsv",
        "ref.accession_species_taxmap.tsv",
    ]
    candidates: list[Path] = []
    if ref.expanduser().is_dir():
        candidates.extend(ref / sidecar for sidecar in candidate_names)
        candidates.extend(
            unique_glob(
                ref,
                [
                    "*.candidate_surface_taxmap.tsv",
                    "*.accession_species_taxmap.tsv",
                    "*.accession_taxmap.tsv",
                ],
            )
        )
    candidates.extend(
        [
            parent / f"{name}.candidate_surface_taxmap.tsv",
            parent / f"{stem}.candidate_surface_taxmap.tsv",
            parent / f"{name}.accession_species_taxmap.tsv",
            parent / f"{stem}.accession_species_taxmap.tsv",
        ]
    )
    candidates.extend(parent / sidecar for sidecar in candidate_names)
    candidates.extend(
        unique_glob(
            parent,
            [
                "*.candidate_surface_taxmap.tsv",
                "*.accession_species_taxmap.tsv",
                "*.accession_taxmap.tsv",
            ],
        )
    )
    return first_existing_file(candidates)


def discover_profile_defaults_path(argv: Sequence[str]) -> Optional[Path]:
    ref_value = option_value(argv, ["--ref"], ["-r"])
    if not ref_value:
        return None
    ref = Path(ref_value)
    parent = ref.parent if str(ref.parent) else Path(".")
    candidates: list[Path] = []
    if ref.expanduser().is_dir():
        candidates.extend(ref / name for name in PROFILE_DEFAULT_FILENAMES)
    candidates.extend(parent / name for name in PROFILE_DEFAULT_FILENAMES)
    return first_existing_file(candidates)


def read_profile_defaults(path: Path) -> dict[str, str]:
    defaults: dict[str, str] = {}
    with path.expanduser().open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 2:
                continue
            key = fields[0].strip().lstrip("-").replace("-", "_")
            value = fields[1].strip()
            if key == "option" and value == "value":
                continue
            if key in PROFILE_DEFAULT_OPTIONS or key in PROFILE_DEFAULT_FLAGS:
                defaults[key] = value
    return defaults


def apply_packaged_profile_defaults(
    args: Sequence[str],
) -> tuple[list[str], Optional[Path], dict[str, str]]:
    out = list(args)
    defaults_path = discover_profile_defaults_path(out)
    if defaults_path is None:
        return out, None, {}
    defaults = read_profile_defaults(defaults_path)
    applied: dict[str, str] = {}

    for key, option in PROFILE_DEFAULT_OPTIONS.items():
        value = defaults.get(key, "")
        if value and not has_option(out, [option]):
            out.extend([option, value])
            applied[key] = value

    for key, option in PROFILE_DEFAULT_FLAGS.items():
        value = defaults.get(key, "").strip().lower()
        enabled = value in {"1", "true", "yes", "on"}
        if enabled and not has_option(out, [option]):
            out.append(option)
            applied[key] = "true"

    return out, defaults_path, applied


def preflight_path_status(
    rows: list[tuple[str, str, str, str, str]],
    errors: list[str],
    name: str,
    source: str,
    value: str,
    kind: str,
    required: bool = True,
) -> bool:
    if not value:
        status = "missing" if required else "optional"
        rows.append((name, status, source, "", "not supplied"))
        if required:
            errors.append(f"{name} is required")
        return not required

    path = Path(value).expanduser()
    if kind == "dir":
        ok = path.is_dir()
        detail = "directory exists" if ok else "directory not found"
    elif kind == "training-dir":
        ok = path.is_dir() and (path / "train.joined_features.tsv").is_file()
        detail = "train.joined_features.tsv found" if ok else "training directory missing train.joined_features.tsv"
    elif kind == "executable":
        ok = path.is_file() and os.access(path, os.X_OK)
        detail = "executable" if ok else "not executable or not found"
    else:
        ok = path.is_file()
        detail = "file exists" if ok else "file not found"

    rows.append((name, "pass" if ok else "fail", source, str(path), detail))
    if required and not ok:
        errors.append(f"{name} {detail}: {path}")
    return ok


def source_value(
    args: Sequence[str],
    option_names: Sequence[str],
    env: Mapping[str, str],
    env_name: str,
    discover_value: Optional[Path],
    short_names: Sequence[str] = (),
) -> tuple[str, str]:
    explicit = option_value(args, option_names, short_names)
    if explicit:
        return "explicit", explicit
    env_value = env.get(env_name, "")
    if env_value:
        return f"env:{env_name}", env_value
    if discover_value is not None:
        return "sidecar", str(discover_value)
    return "missing", ""


def run_preflight(argv: Sequence[str], env: Optional[Mapping[str, str]] = None) -> int:
    env_map = os.environ if env is None else env
    rows: list[tuple[str, str, str, str, str]] = []
    errors: list[str] = []

    try:
        preset_arg, args = extract_wrapper_option(argv, WRAPPER_PROFILE_PRESET_OPT)
        no_profile_rescue, args = extract_wrapper_flags(args, (WRAPPER_NO_PROFILE_RESCUE_OPT,))
    except SystemExit as exc:
        print("MinCO profile preflight")
        print("status\tfail")
        print(f"error\t{exc}")
        return 1

    preset = preset_arg or env_map.get(ENV_PROFILE_PRESET, "") or DEFAULT_PRESET
    if preset not in VALID_PROFILE_PRESETS:
        errors.append(f"{WRAPPER_PROFILE_PRESET_OPT} must be one of: {', '.join(sorted(VALID_PROFILE_PRESETS))}")
    rows.append(("profile_preset", "pass" if preset in VALID_PROFILE_PRESETS else "fail", "wrapper", preset, "selected wrapper preset"))
    rows.append(
        (
            "profile_rescue",
            "pass",
            "wrapper",
            "off" if no_profile_rescue else "default",
            f"{WRAPPER_NO_PROFILE_RESCUE_OPT} opt-out",
        )
    )

    explicit_scope = has_option(args, ["--scope"])
    explicit_train_pool = has_option(args, ["--train-pool"])
    args, defaults_path, applied_defaults = apply_packaged_profile_defaults(args)
    if defaults_path is not None:
        rows.append(("profile_defaults", "pass", "sidecar", str(defaults_path), "packaged profile defaults"))
    for key, value in sorted(applied_defaults.items()):
        rows.append((f"default_{key}", "pass", "sidecar", value, "applied packaged default"))

    strategy = option_value(args, ["--strategy"]) or DEFAULT_STRATEGY
    rows.append(("strategy", "pass", "explicit" if has_option(args, ["--strategy"]) else "default", strategy, "delegated strategy"))
    scope = option_value(args, ["--scope"]) or "all"
    train_pool = option_value(args, ["--train-pool"]) or "train12"
    filter_training_scope = has_option(args, ["--filter-training-scope"])
    scope_source = "explicit" if explicit_scope else "sidecar" if "scope" in applied_defaults else "default"
    train_pool_source = "explicit" if explicit_train_pool else "sidecar" if "train_pool" in applied_defaults else "default"
    rows.append(("scope", "pass", scope_source, scope, "taxonomic output scope"))
    rows.append(("train_pool", "pass", train_pool_source, train_pool, "calibration training pool"))
    rows.append(
        (
            "filter_training_scope",
            "pass",
            "explicit" if filter_training_scope else "default",
            "true" if filter_training_scope else "false",
            "model-cache metadata setting",
        )
    )

    ref = option_value(args, ["--ref"], ["-r"]) or ""
    preflight_path_status(rows, errors, "ref", "explicit" if ref else "missing", ref, "dir")

    reads = option_value(args, ["--reads", "--qraw"]) or ""
    if reads:
        preflight_path_status(rows, errors, "reads", "explicit", reads, "file")

    taxmap_source, taxmap = source_value(
        args,
        ["--taxmap"],
        env_map,
        ENV_TAXMAP,
        discover_taxmap(args),
    )
    preflight_path_status(rows, errors, "species_taxmap", taxmap_source, taxmap, "file")

    model_source, model_cache = source_value(
        args,
        ["--model-cache"],
        env_map,
        ENV_MODEL_CACHE,
        discover_model_cache(args),
    )
    if model_cache:
        if preflight_path_status(rows, errors, "model_cache", model_source, model_cache, "file"):
            try:
                calibrated.load_model_cache(Path(model_cache), train_pool, scope, filter_training_scope)
            except SystemExit as exc:
                rows.append(("model_cache_metadata", "fail", model_source, model_cache, str(exc)))
                errors.append(str(exc))
            except Exception as exc:
                rows.append(("model_cache_metadata", "fail", model_source, model_cache, repr(exc)))
                errors.append(f"could not validate model cache metadata: {exc!r}")
            else:
                rows.append(("model_cache_metadata", "pass", model_source, model_cache, "metadata matches scope/train-pool"))
    else:
        train_source, train_features = source_value(
            args,
            ["--train-features"],
            env_map,
            ENV_TRAIN_FEATURES,
            discover_train_features(args),
        )
        train_tables = option_values(args, ["--train-table"])
        if train_features:
            preflight_path_status(rows, errors, "train_features", train_source, train_features, "training-dir")
        elif train_tables:
            ok_tables = True
            for index, table in enumerate(train_tables, start=1):
                ok_tables &= preflight_path_status(
                    rows,
                    errors,
                    f"train_table_{index}",
                    "explicit",
                    table,
                    "file",
                )
            if not ok_tables:
                errors.append("one or more --train-table inputs are missing")
        else:
            rows.append(("model_or_training", "fail", "missing", "", "model cache or training sidecar required"))
            errors.append(
                f"model cache or training features are required; provide --model-cache, set {ENV_MODEL_CACHE}, "
                f"or package joined_feature_training/"
            )

    surface_switch = option_value(args, ["--candidate-surface-switch"])
    rescue_switch = option_value(args, ["--candidate-rescue-switch"])
    if rescue_switch is None:
        if preset == CANDIDATE_PRESET:
            rescue_switch = (
                calibrated.CANDIDATE_RESCUE_SWITCH_OFF
                if no_profile_rescue
                else calibrated.CANDIDATE_RESCUE_SWITCH_SPLIT_P002_X300_ANI95_AF60_BR025_D1_TOP1
            )
        else:
            rescue_switch = calibrated.CANDIDATE_RESCUE_SWITCH_OFF
    rows.append(
        (
            "candidate_rescue_switch",
            "pass",
            "explicit" if has_option(args, ["--candidate-rescue-switch"]) else "preset",
            rescue_switch,
            "profile rescue mode",
        )
    )
    if surface_switch is None:
        surface_switch = (
            calibrated.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
            if preset == CANDIDATE_PRESET
            else calibrated.CANDIDATE_SURFACE_SWITCH_OFF
        )
    rows.append(("candidate_surface_switch", "pass", "explicit" if has_option(args, ["--candidate-surface-switch"]) else "preset", surface_switch, "candidate surface mode"))
    if surface_switch != calibrated.CANDIDATE_SURFACE_SWITCH_OFF:
        surface_source, surface_taxmap = source_value(
            args,
            ["--candidate-surface-taxmap"],
            env_map,
            ENV_CANDIDATE_SURFACE_TAXMAP,
            discover_candidate_surface_taxmap(args),
        )
        preflight_path_status(
            rows,
            errors,
            "candidate_surface_taxmap",
            surface_source,
            surface_taxmap,
            "file",
        )

    minco_value = option_value(args, ["--minco"]) or env_map.get(ENV_MINCO, "")
    minco_source = "explicit" if has_option(args, ["--minco"]) else f"env:{ENV_MINCO}" if minco_value else "default"
    if not minco_value:
        minco_value = str(Path(__file__).resolve().parents[1] / "bin/minco")
    preflight_path_status(rows, errors, "minco_binary", minco_source, minco_value, "executable")

    delegated_args: list[str] = []
    if not errors:
        try:
            delegated_args = build_calibrated_argv(argv, env_map)
        except SystemExit as exc:
            errors.append(str(exc))

    print("MinCO profile preflight")
    print(f"status\t{'pass' if not errors else 'fail'}")
    print("check\tstatus\tsource\tpath_or_value\tdetail")
    for row in rows:
        print("\t".join(row))
    if delegated_args:
        print("delegated_argv\t" + " ".join(shlex.quote(part) for part in delegated_args))
    for error in errors:
        print(f"error\t{error}")
    return 0 if not errors else 1


def apply_profile_preset(
    args: list[str],
    preset: str,
    env: Mapping[str, str],
    no_profile_rescue: bool = False,
) -> list[str]:
    if preset not in VALID_PROFILE_PRESETS:
        raise SystemExit(
            f"{WRAPPER_PROFILE_PRESET_OPT} must be one of: {', '.join(sorted(VALID_PROFILE_PRESETS))}"
        )
    if preset == CURRENT_PRESET:
        return args

    if not has_option(args, ["--candidate-rescue-switch"]):
        rescue_switch = (
            calibrated.CANDIDATE_RESCUE_SWITCH_OFF
            if no_profile_rescue
            else calibrated.CANDIDATE_RESCUE_SWITCH_SPLIT_P002_X300_ANI95_AF60_BR025_D1_TOP1
        )
        args.extend(
            [
                "--candidate-rescue-switch",
                rescue_switch,
            ]
        )
    if not has_option(args, ["--candidate-surface-switch"]):
        args.extend(
            [
                "--candidate-surface-switch",
                calibrated.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
            ]
        )
    if not has_option(args, ["--candidate-abundance-policy"]):
        args.extend(
            [
                "--candidate-abundance-policy",
                calibrated.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
            ]
        )
    if not has_option(args, ["--abundance-ani-floor"]):
        args.extend(["--abundance-ani-floor", "0.90"])
    if not has_option(args, ["--abundance-sparse-depth-cap"]):
        args.extend(
            [
                "--abundance-sparse-depth-cap",
                calibrated.ABUNDANCE_SPARSE_DEPTH_CAP_SWITCH_POISSON_BREADTH,
            ]
        )
    if not has_option(args, ["--abundance-sparse-breadth-max"]):
        args.extend(["--abundance-sparse-breadth-max", "0.15"])
    if not has_option(args, ["--abundance-sparse-depth-ratio-min"]):
        args.extend(["--abundance-sparse-depth-ratio-min", "200"])
    surface_switch = option_value(args, ["--candidate-surface-switch"]) or ""
    surface_enabled = surface_switch != calibrated.CANDIDATE_SURFACE_SWITCH_OFF
    if surface_enabled and not has_option(args, ["--candidate-surface-taxmap"]):
        surface_taxmap = env.get(ENV_CANDIDATE_SURFACE_TAXMAP, "")
        if not surface_taxmap:
            discovered = discover_candidate_surface_taxmap(args)
            surface_taxmap = str(discovered) if discovered else ""
        if surface_taxmap:
            args.extend(["--candidate-surface-taxmap", surface_taxmap])
    return args


def build_calibrated_argv(
    argv: Optional[Sequence[str]] = None,
    env: Optional[Mapping[str, str]] = None,
) -> list[str]:
    args = list(sys.argv[1:] if argv is None else argv)
    preset_arg, args = extract_wrapper_option(args, WRAPPER_PROFILE_PRESET_OPT)
    no_profile_rescue, args = extract_wrapper_flags(args, (WRAPPER_NO_PROFILE_RESCUE_OPT,))
    if help_requested(args):
        return args

    env_map = os.environ if env is None else env
    preset = preset_arg or env_map.get(ENV_PROFILE_PRESET, "") or DEFAULT_PRESET
    if preset not in VALID_PROFILE_PRESETS:
        raise SystemExit(
            f"{WRAPPER_PROFILE_PRESET_OPT} must be one of: {', '.join(sorted(VALID_PROFILE_PRESETS))}"
        )

    args, _, _ = apply_packaged_profile_defaults(args)

    if not has_option(args, ["--strategy"]):
        args.extend(["--strategy", DEFAULT_STRATEGY])

    if not has_option(args, ["--taxmap"]):
        taxmap = env_map.get(ENV_TAXMAP, "")
        if not taxmap:
            discovered = discover_taxmap(args)
            taxmap = str(discovered) if discovered else ""
        if not taxmap:
            raise SystemExit(
                f"pass --taxmap, set {ENV_TAXMAP}, or place a species_taxmap.tsv sidecar beside --ref"
            )
        args.extend(["--taxmap", taxmap])

    if not has_option(args, ["--model-cache"]):
        model_cache = env_map.get(ENV_MODEL_CACHE, "")
        if not model_cache:
            discovered_model = discover_model_cache(args)
            model_cache = str(discovered_model) if discovered_model else ""
        if model_cache:
            args.extend(["--model-cache", model_cache])

    if not has_option(args, ["--train-features", "--train-table"]) and not has_option(args, ["--model-cache"]):
        train_features = env_map.get(ENV_TRAIN_FEATURES, "")
        if not train_features:
            discovered = discover_train_features(args)
            train_features = str(discovered) if discovered else ""
        if not train_features:
            raise SystemExit(
                f"pass --model-cache, --train-features/--train-table, set {ENV_MODEL_CACHE} "
                f"or {ENV_TRAIN_FEATURES}, or place a model cache or joined_feature_training "
                "sidecar beside --ref"
            )
        args.extend(["--train-features", train_features])

    if not has_option(args, ["--minco"]):
        minco = env_map.get(ENV_MINCO, "")
        if minco:
            args.extend(["--minco", minco])

    args = apply_profile_preset(args, preset, env_map, no_profile_rescue)
    return args


def main(argv: Optional[Sequence[str]] = None, env: Optional[Mapping[str, str]] = None) -> int:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    preflight, raw_args = extract_wrapper_flags(raw_args, WRAPPER_PREFLIGHT_OPTS)
    if preflight:
        return run_preflight(raw_args, env)
    args = build_calibrated_argv(raw_args, env)
    if help_requested(args):
        print(default_help_banner(), file=sys.stderr)
    return calibrated.main(args)


if __name__ == "__main__":
    raise SystemExit(main())
