#!/usr/bin/env python3
"""Regression checks for the calibrated profile auto-exact wrapper path."""

from __future__ import annotations

import contextlib
import io
import tempfile
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import scripts.minco_profile_calibrated as wrapper  # noqa: E402
import scripts.minco_profile_default as default_wrapper  # noqa: E402


MINCO_COLS = [
    "Ref",
    "ANI",
    "XnY_ctx",
    "N_diff_obj",
    "N_diff_obj_section",
    "N_mut2_ctx",
    "Real_min_align_fraction",
    "Ref_breadth",
    "Ref_mean_depth",
    "Ref_hit_mean_depth",
    "Ref_depth_variance",
    "Ref_depth_cv",
    "Normalized_abundance_depth",
    "Ref_zip_af",
    "Ref_zip_aaf_ani",
]


@pytest.fixture
def work(tmp_path: Path) -> Path:
    return tmp_path


def accession(i: int) -> str:
    return f"GCF_{i:09d}.1"


def write_taxmap(path: Path, n: int = 25) -> None:
    with path.open("w") as fh:
        fh.write("ref_key\tTAXID\tRANK\tTAXPATH\tTAXPATHSN\n")
        for i in range(1, n + 1):
            acc = accession(i)
            taxid = str(1000 + i)
            pathsn = f"d__Bacteria|p__Test|c__Test|o__Test|f__Test|g__Minco|s__fixture_{i}"
            fh.write(f"{acc}\t{taxid}\tspecies\t{pathsn}\t{pathsn}\n")


def write_ncbi_style_taxmap(path: Path, n: int = 25) -> None:
    with path.open("w") as fh:
        fh.write("ref_key\tTAXID\tRANK\tTAXPATH\tTAXPATHSN\n")
        for i in range(1, n + 1):
            acc = accession(i)
            taxid = str(1000 + i)
            pathsn = f"root|Bacteria|Fixture species {i}"
            fh.write(f"{acc}\t{taxid}\tspecies\t1|2|{taxid}\t{pathsn}\n")


def minco_rows(extra_mass: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    unique = []
    split = []
    for i in range(1, 21):
        ref = f"/ref/{accession(i)}.fna"
        unique.append([ref, 0.99, 20, 0, 0, 0, 0.60, 0.80, 2.0, 2.0, 0.1, 0.1, 1.0, 0.80, 0.98])
        split.append([ref, 0.99, 25, 0, 0, 0, 0.60, 0.80, 2.0, 2.0, 0.1, 0.1, 1.0, 0.80, 0.98])
    ref = f"/ref/{accession(21)}.fna"
    unique.append([ref, 0.90, 5, 0, 0, 0, 0.01, 0.01, 0.2, 0.2, 0.1, 0.1, 0.0, 0.01, 0.90])
    split.append([ref, 0.99, 50, 0, 0, 0, 0.60, 0.50, 2.0, 2.0, 0.1, 0.1, extra_mass, 0.50, 0.98])
    return pd.DataFrame(unique, columns=MINCO_COLS), pd.DataFrame(split, columns=MINCO_COLS)


def minco_rows_low_extra() -> tuple[pd.DataFrame, pd.DataFrame]:
    unique = []
    split = []
    for i in range(1, 21):
        ref = f"/ref/{accession(i)}.fna"
        unique.append([ref, 0.99, 20, 0, 0, 0, 0.20, 0.20, 2.0, 2.0, 0.1, 0.1, 1.0, 0.50, 0.97])
        split.append([ref, 0.99, 25, 0, 0, 0, 0.20, 0.20, 2.0, 2.0, 0.1, 0.1, 1.0, 0.50, 0.97])
    ref = f"/ref/{accession(21)}.fna"
    unique.append([ref, 0.90, 5, 0, 0, 0, 0.01, 0.01, 0.2, 0.2, 0.1, 0.1, 0.0, 0.01, 0.90])
    split.append([ref, 0.99, 80, 0, 0, 0, 0.60, 0.30, 1.0, 1.0, 0.1, 0.1, 0.001, 0.30, 0.96])
    return pd.DataFrame(unique, columns=MINCO_COLS), pd.DataFrame(split, columns=MINCO_COLS)


def write_minco_tables(unique_path: Path, split_path: Path, extra_mass: float) -> None:
    unique, split = minco_rows(extra_mass)
    unique.to_csv(unique_path, sep="\t", index=False)
    split.to_csv(split_path, sep="\t", index=False)


def patch_model_hooks() -> tuple[object, object, object]:
    old_load = wrapper.load_training_features
    old_fit = wrapper.fit_rf_hgb
    old_predict = wrapper.predict_rf_hgb
    wrapper.load_training_features = lambda *args, **kwargs: pd.DataFrame({"label": [0, 1]})
    wrapper.fit_rf_hgb = lambda train: [object(), object()]
    wrapper.predict_rf_hgb = lambda models, features: np.full(len(features), 0.9)
    return old_load, old_fit, old_predict


def restore_model_hooks(old_hooks: tuple[object, object, object]) -> None:
    wrapper.load_training_features, wrapper.fit_rf_hgb, wrapper.predict_rf_hgb = old_hooks


@pytest.fixture(autouse=True)
def patched_model_hooks():
    old_hooks = patch_model_hooks()
    try:
        yield
    finally:
        restore_model_hooks(old_hooks)


def read_output(path: Path) -> pd.DataFrame:
    assert path.exists() and path.stat().st_size > 0, path
    return pd.read_csv(path, sep="\t")


def assert_bool_column(df: pd.DataFrame, col: str, expected: bool) -> None:
    values = {str(v).lower() for v in df[col].unique()}
    assert values == {str(expected).lower()}, (col, values, expected)


def adaptive_fixture(base_af: float) -> tuple[pd.DataFrame, dict[str, str]]:
    rows: list[dict[str, object]] = []
    names: dict[str, str] = {}
    for i in range(1, 6):
        taxid = str(1000 + i)
        names[taxid] = f"Fixture species_{i}"
        rows.append(
            {
                "taxid": taxid,
                "calibrated_probability": 0.9,
                "u_XnY_ctx_max": 20,
                "u_ANI_max": 0.99,
                "u_Real_min_align_fraction_max": base_af,
                "s_XnY_ctx_max": 25,
                "s_ANI_max": 0.99,
                "s_Real_min_align_fraction_max": base_af,
                "s_Ref_breadth_max": base_af,
                "s_Ref_mean_depth_max": 2.0,
                "s_Ref_zip_af_max": 0.50,
                "s_Ref_hit_mean_depth_max": 2.0,
                "s_Normalized_abundance_depth_max": 1.0,
            }
        )
    for i in range(6, 16):
        taxid = str(1000 + i)
        names[taxid] = f"Fixture species_{i}"
        rows.append(
            {
                "taxid": taxid,
                "calibrated_probability": 0.9,
                "u_XnY_ctx_max": 5,
                "u_ANI_max": 0.90,
                "u_Real_min_align_fraction_max": 0.01,
                "s_XnY_ctx_max": 60,
                "s_ANI_max": 0.99,
                "s_Real_min_align_fraction_max": 0.60,
                "s_Ref_breadth_max": 0.50,
                "s_Ref_mean_depth_max": 0.2,
                "s_Ref_zip_af_max": 0.01,
                "s_Ref_hit_mean_depth_max": 2.0,
                "s_Normalized_abundance_depth_max": 1.0,
            }
        )
    return pd.DataFrame(rows), names


def low_extra_fixture(base_af: float) -> tuple[pd.DataFrame, dict[str, str]]:
    rows: list[dict[str, object]] = []
    names: dict[str, str] = {}
    for i in range(1, 6):
        taxid = str(2000 + i)
        names[taxid] = f"Fixturebase species_{i}"
        rows.append(
            {
                "taxid": taxid,
                "calibrated_probability": 0.9,
                "u_XnY_ctx_max": 20,
                "u_ANI_max": 0.99,
                "u_Real_min_align_fraction_max": base_af,
                "s_XnY_ctx_max": 25,
                "s_ANI_max": 0.99,
                "s_Real_min_align_fraction_max": base_af,
                "s_Ref_breadth_max": base_af,
                "s_Ref_mean_depth_max": 2.0,
                "s_Ref_zip_af_max": 0.50,
                "s_Ref_hit_mean_depth_max": 2.0,
                "s_Normalized_abundance_depth_max": 1.0,
            }
        )
    for i in range(6, 9):
        taxid = str(2000 + i)
        names[taxid] = f"Fixtureextra{i} species"
        rows.append(
            {
                "taxid": taxid,
                "calibrated_probability": 0.25,
                "u_XnY_ctx_max": 5,
                "u_ANI_max": 0.90,
                "u_Real_min_align_fraction_max": 0.01,
                "s_XnY_ctx_max": 100,
                "s_ANI_max": 0.96,
                "s_Real_min_align_fraction_max": 0.10,
                "s_Ref_breadth_max": 0.10,
                "s_Ref_mean_depth_max": 1.0,
                "s_Ref_zip_af_max": 0.30,
                "s_Ref_hit_mean_depth_max": 1.5,
                "s_Normalized_abundance_depth_max": 0.001,
            }
        )
    return pd.DataFrame(rows), names


def test_high_extra_low_unique_af_uses_probability() -> None:
    features, names = adaptive_fixture(base_af=0.20)
    mask, details = wrapper.adaptive_sub95_mask(features, names, 0.35)
    assert details["adaptive_mode"] == "probability_high_extra_low_uaf_tail025"
    assert int(mask.sum()) == len(features)
    assert details["joined_base_median_uaf"] < 0.35
    assert details["tail_rescue_added_n"] > 0
    assert details["tail_rescue_probability_threshold"] == 0.25
    assert details["tail_rescue_median_uaf_threshold"] == 0.35
    raw = wrapper.panel_abundance_raw(features)
    tail_added = features["tail_rescue_added"].astype(bool).to_numpy()
    base_rows = wrapper.numeric(features, "u_XnY_ctx_max").to_numpy(dtype=float) >= 10.0
    assert np.allclose(raw[tail_added], 0.2)
    assert np.allclose(raw[base_rows], 4.0)


def test_high_extra_high_unique_af_keeps_rescue() -> None:
    features, names = adaptive_fixture(base_af=0.60)
    mask, details = wrapper.adaptive_sub95_mask(features, names, 0.35)
    assert details["adaptive_mode"] == "rescue"
    assert int(mask.sum()) < len(features)
    assert details["joined_base_median_uaf"] >= 0.35
    assert details["tail_rescue_added_n"] == 0


def test_low_extra_low_unique_af_adds_guarded_split_rescue() -> None:
    features, names = low_extra_fixture(base_af=0.20)
    mask, details = wrapper.adaptive_sub95_mask(features, names, 0.35)
    assert details["adaptive_mode"] == "base_low_extra_split_rescue"
    assert details["joined_base_median_uaf"] < 0.45
    assert details["low_extra_split_rescue_added_n"] == 3
    assert details["low_extra_split_rescue_probability_threshold"] == 0.20
    assert details["low_extra_split_rescue_median_uaf_threshold"] == 0.45
    assert details["low_extra_split_rescue_topn_per_genus"] == 1
    assert int(mask.sum()) == 8
    assert int(features["low_extra_split_rescue_added"].astype(bool).sum()) == 3


def test_low_extra_high_unique_af_keeps_base_only() -> None:
    features, names = low_extra_fixture(base_af=0.60)
    mask, details = wrapper.adaptive_sub95_mask(features, names, 0.35)
    assert details["adaptive_mode"] == "base"
    assert details["joined_base_median_uaf"] >= 0.45
    assert details["low_extra_split_rescue_added_n"] == 0
    assert int(mask.sum()) == 5


def test_reported_ani_prefers_split_zip_aaf_unless_raw_unique_fallback() -> None:
    features = pd.DataFrame(
        {
            "s_Ref_zip_aaf_ani_max": [0.97, 0.96, 0.0],
            "u_Ref_zip_aaf_ani_max": [0.99, 0.98, 0.95],
            "raw_unique_fallback": [False, True, False],
        }
    )
    ani, source = wrapper.reported_ani_values(features)
    assert np.allclose(ani, [0.97, 0.98, 0.95])
    assert list(source) == [
        "s_Ref_zip_aaf_ani_max",
        "u_Ref_zip_aaf_ani_max",
        "u_Ref_zip_aaf_ani_max",
    ]


def test_panel_abundance_uses_max_split_unique_zip_depth_and_tail_exception() -> None:
    features = pd.DataFrame(
        {
            "s_Ref_mean_depth_max": [2.0, 1.0, 10.0],
            "s_Ref_zip_af_max": [0.5, 0.5, 0.2],
            "u_Ref_mean_depth_max": [1.0, 3.0, 50.0],
            "u_Ref_zip_af_max": [0.5, 0.5, 0.1],
            "tail_rescue_added": [False, False, True],
        }
    )
    raw = wrapper.panel_abundance_raw(features)
    assert np.allclose(raw, [4.0, 6.0, 10.0])


def test_genus_xny_abundance_blend_is_disabled_at_alpha_zero() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Alpha two"],
            "s_XnY_ctx_max": [1000, 100],
        }
    )
    raw = np.array([10.0, 10.0])
    blended = wrapper.genus_xny_blended_abundance_raw(
        features,
        np.array([True, True]),
        raw,
        alpha=0.0,
    )
    assert np.allclose(blended, raw)


def test_genus_xny_abundance_blend_reallocates_only_within_called_genus() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Alpha two", "Beta one", "Alpha three"],
            "s_XnY_ctx_max": [1000, 100, 500, 1000],
        }
    )
    raw = np.array([10.0, 10.0, 7.0, 30.0])
    blended = wrapper.genus_xny_blended_abundance_raw(
        features,
        np.array([True, True, True, False]),
        raw,
        alpha=0.25,
    )
    expected_alpha_realloc = np.array([18.1818181818, 1.8181818182])
    expected = raw.copy()
    expected[:2] = 0.75 * raw[:2] + 0.25 * expected_alpha_realloc
    assert np.allclose(blended, expected)
    assert np.isclose(blended[:2].sum(), raw[:2].sum())
    assert blended[2] == raw[2]
    assert blended[3] == raw[3]


def test_guarded_feature_allocator_default_off_preserves_raw_abundance() -> None:
    features = pd.DataFrame({"species_name": ["Alpha one", "Alpha two"]})
    raw = np.array([10.0, 5.0])
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        features,
        np.array([True, True]),
        raw,
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF,
    )
    assert np.allclose(adjusted, raw)
    assert details["abundance_feature_allocator_applied"] is False
    assert details["abundance_feature_allocator_switch"] == wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF


def test_guarded_feature_allocator_requires_multi_genus_mass_guard() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Beta one"],
            "s_Ref_hit_mean_depth_max": [20.0, 1.0],
            "s_Ref_breadth_max": [1.0, 1.0],
            "u_Ref_hit_mean_depth_max": [0.0, 0.0],
            "u_Ref_breadth_max": [0.0, 0.0],
        }
    )
    raw = np.array([10.0, 5.0])
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        features,
        np.array([True, True]),
        raw,
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002,
    )
    assert np.allclose(adjusted, raw)
    assert details["abundance_feature_allocator_guard_passed"] is False
    assert details["abundance_feature_allocator_applied"] is False


def test_guarded_feature_allocator_reallocates_base_rows_and_preserves_candidates() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Alpha two", "Beta one", "Beta two", "Gamma one"],
            "candidate_surface_added": [False, False, False, True, False],
            "candidate_rescue_added": [False, False, False, False, False],
            "s_Ref_hit_mean_depth_max": [10.0, 1.0, 1.0, 100.0, 1.0],
            "s_Ref_breadth_max": [1.0, 1.0, 1.0, 1.0, 1.0],
            "u_Ref_hit_mean_depth_max": [0.0, 0.0, 0.0, 0.0, 0.0],
            "u_Ref_breadth_max": [0.0, 0.0, 0.0, 0.0, 0.0],
        }
    )
    raw = np.array([10.0, 10.0, 5.0, 5.0, 2.0])
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        features,
        np.array([True, True, True, True, True]),
        raw,
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002,
    )
    alpha_target = np.array([18.1818181818, 1.8181818182])
    expected = raw.copy()
    expected[:2] = 0.98 * raw[:2] + 0.02 * alpha_target
    assert np.allclose(adjusted, expected)
    assert adjusted[3] == raw[3]
    assert adjusted.sum() == pytest.approx(raw.sum())
    assert details["abundance_feature_allocator_guard_passed"] is True
    assert details["abundance_feature_allocator_applied"] is True
    assert details["abundance_feature_allocator_adjusted_rows_n"] == 2
    assert details["abundance_feature_allocator_base_mass_multi_genus_frac"] == pytest.approx(20.0 / 27.0)


def test_refined_feature_allocator_requires_split_xny_median_guard() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Alpha two", "Beta one"],
            "candidate_surface_added": [False, False, False],
            "candidate_rescue_added": [False, False, False],
            "s_XnY_ctx_max": [120.0, 140.0, 1000.0],
            "s_Ref_hit_mean_depth_max": [10.0, 1.0, 1.0],
            "s_Ref_breadth_max": [1.0, 1.0, 1.0],
            "u_Ref_hit_mean_depth_max": [0.0, 0.0, 0.0],
            "u_Ref_breadth_max": [0.0, 0.0, 0.0],
        }
    )
    raw = np.array([10.0, 10.0, 1.0])
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        features,
        np.array([True, True, True]),
        raw,
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230,
    )
    assert np.allclose(adjusted, raw)
    assert details["abundance_feature_allocator_guard_passed"] is False
    assert details["abundance_feature_allocator_applied"] is False
    assert details["abundance_feature_allocator_s_xny_median"] == pytest.approx(140.0)

    features["s_XnY_ctx_max"] = [260.0, 300.0, 1000.0]
    adjusted, details = wrapper.guarded_feature_allocator_abundance_raw(
        features,
        np.array([True, True, True]),
        raw,
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230,
    )
    assert details["abundance_feature_allocator_guard_passed"] is True
    assert details["abundance_feature_allocator_applied"] is True
    assert details["abundance_feature_allocator_s_xny_median"] == pytest.approx(300.0)


def test_feature_allocator_genus_labels_can_use_accession_taxmap() -> None:
    features = pd.DataFrame(
        {
            "species_name": ["Alpha one", "Alpha two"],
            "s_best_accession": ["GCF_000000001.1", "GCF_000000002.1"],
        }
    )
    taxmap = {
        "GCF_000000001.1": {"species_name": "s__Alpha_A one"},
        "GCF_000000002.1": {"species_name": "s__Alpha_B two"},
    }
    raw = np.array([10.0, 5.0])
    without_taxmap = wrapper.base_multi_genus_mass_fraction(
        features,
        np.array([True, True]),
        raw,
    )
    with_taxmap = wrapper.base_multi_genus_mass_fraction(
        features,
        np.array([True, True]),
        raw,
        taxmap,
    )
    assert without_taxmap == pytest.approx(1.0)
    assert with_taxmap == pytest.approx(0.0)


def test_adaptive_call_filter_switch_default_off_preserves_calls() -> None:
    features = pd.DataFrame(
        {
            "s_XnY_ctx_max": [500.0, 10.0],
            "u_XnY_ctx_max": [500.0, 10.0],
        }
    )
    call = np.array([True, True])
    filtered, details = wrapper.apply_adaptive_call_filter_switch(
        features,
        call,
        wrapper.ADAPTIVE_CALL_FILTER_SWITCH_OFF,
    )
    assert filtered.tolist() == [True, True]
    assert details["adaptive_call_filter_applied"] is False
    assert details["adaptive_call_filter_removed_n"] == 0


def test_adaptive_call_filter_switch_prunes_only_when_sample_gate_fires() -> None:
    features = pd.DataFrame(
        {
            "s_XnY_ctx_max": [500.0, 500.0, 10.0],
            "u_XnY_ctx_max": [500.0, 20.0, 10.0],
        }
    )
    call = np.array([True, True, False])
    filtered, details = wrapper.apply_adaptive_call_filter_switch(
        features,
        call,
        wrapper.ADAPTIVE_CALL_FILTER_SWITCH_LOPO_MIN_XNY25,
    )
    assert filtered.tolist() == [True, False, False]
    assert details["adaptive_call_filter_applied"] is True
    assert details["adaptive_call_filter_removed_n"] == 1
    assert details["adaptive_call_filter_input_call_n"] == 2
    assert details["adaptive_call_filter_output_call_n"] == 1


def test_adaptive_call_filter_switch_does_not_prune_low_signal_sample() -> None:
    features = pd.DataFrame(
        {
            "s_XnY_ctx_max": [200.0, 200.0],
            "u_XnY_ctx_max": [200.0, 20.0],
        }
    )
    call = np.array([True, True])
    filtered, details = wrapper.apply_adaptive_call_filter_switch(
        features,
        call,
        wrapper.ADAPTIVE_CALL_FILTER_SWITCH_LOPO_MIN_XNY25,
    )
    assert filtered.tolist() == [True, True]
    assert details["adaptive_call_filter_applied"] is False
    assert details["adaptive_call_filter_removed_n"] == 0


def test_candidate_rescue_switch_default_off_preserves_calls_and_abundance() -> None:
    features = pd.DataFrame(
        {
            "s_Ref_zip_aaf_ani_max": [0.99, 0.99],
            "u_Ref_zip_aaf_ani_max": [0.99, 0.99],
            "s_XnY_ctx_max": [200.0, 200.0],
            "u_XnY_ctx_max": [200.0, 200.0],
            "s_Ref_breadth_max": [0.5, 0.5],
            "u_Ref_breadth_max": [0.5, 0.5],
            "s_Real_min_align_fraction_max": [0.8, 0.8],
            "u_Real_min_align_fraction_max": [0.8, 0.8],
        }
    )
    call = np.array([True, False])
    abundance = np.array([10.0, 5.0])
    rescued, rescued_abundance, details = wrapper.apply_candidate_rescue_switch(
        features,
        call,
        abundance,
        wrapper.CANDIDATE_RESCUE_SWITCH_OFF,
    )
    assert rescued.tolist() == [True, False]
    assert np.allclose(rescued_abundance, abundance)
    assert features["candidate_rescue_added"].tolist() == [False, False]
    assert details["candidate_rescue_applied"] is False
    assert details["candidate_rescue_added_n"] == 0


def test_candidate_rescue_switch_adds_strong_uncalled_candidates_with_zero_mass() -> None:
    features = pd.DataFrame(
        {
            "s_Ref_zip_aaf_ani_max": [0.99, 0.91, 0.99],
            "u_Ref_zip_aaf_ani_max": [0.99, 0.89, 0.99],
            "s_XnY_ctx_max": [200.0, 150.0, 150.0],
            "u_XnY_ctx_max": [200.0, 50.0, 150.0],
            "s_Ref_breadth_max": [0.5, 0.02, 0.02],
            "u_Ref_breadth_max": [0.5, 0.00, 0.02],
            "s_Real_min_align_fraction_max": [0.8, 0.75, 0.69],
            "u_Real_min_align_fraction_max": [0.8, 0.10, 0.69],
        }
    )
    call = np.array([True, False, False])
    abundance = np.array([10.0, 5.0, 7.0])
    rescued, rescued_abundance, details = wrapper.apply_candidate_rescue_switch(
        features,
        call,
        abundance,
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70,
    )
    assert rescued.tolist() == [True, True, False]
    assert features["candidate_rescue_added"].tolist() == [False, True, False]
    assert np.allclose(rescued_abundance, [10.0, 0.0, 7.0])
    assert details["candidate_rescue_applied"] is True
    assert details["candidate_rescue_added_n"] == 1
    assert details["candidate_rescue_input_call_n"] == 1
    assert details["candidate_rescue_output_call_n"] == 2
    assert details["candidate_rescue_zero_mass"] is True


def test_candidate_rescue_switch_normalized_depth_policy_assigns_candidate_mass() -> None:
    features = pd.DataFrame(
        {
            "s_Ref_zip_aaf_ani_max": [0.99, 0.91],
            "u_Ref_zip_aaf_ani_max": [0.99, 0.89],
            "s_XnY_ctx_max": [200.0, 150.0],
            "u_XnY_ctx_max": [200.0, 50.0],
            "s_Ref_breadth_max": [0.5, 0.02],
            "u_Ref_breadth_max": [0.5, 0.00],
            "s_Real_min_align_fraction_max": [0.8, 0.75],
            "u_Real_min_align_fraction_max": [0.8, 0.10],
            "s_Normalized_abundance_depth_max": [1.0, 0.20],
            "u_Normalized_abundance_depth_max": [1.0, 0.35],
        }
    )
    rescued, rescued_abundance, details = wrapper.apply_candidate_rescue_switch(
        features,
        np.array([True, False]),
        np.array([10.0, 5.0]),
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70,
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    )
    assert rescued.tolist() == [True, True]
    assert np.allclose(rescued_abundance, [10.0, 0.0])
    assert np.allclose(features["candidate_abundance_norm_mass"], [0.0, 0.70])
    assert details["candidate_rescue_zero_mass"] is False
    assert details["candidate_abundance_policy"] == wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    assert details["candidate_abundance_policy_alpha"] == 2.0


def test_candidate_rescue_switch_uses_zip_ani_not_saturated_raw_ani() -> None:
    features = pd.DataFrame(
        {
            "s_ANI_max": [1.0],
            "u_ANI_max": [1.0],
            "s_Ref_zip_aaf_ani_max": [0.89],
            "u_Ref_zip_aaf_ani_max": [0.89],
            "s_XnY_ctx_max": [500.0],
            "u_XnY_ctx_max": [500.0],
            "s_Ref_breadth_max": [0.5],
            "u_Ref_breadth_max": [0.5],
            "s_Real_min_align_fraction_max": [0.9],
            "u_Real_min_align_fraction_max": [0.9],
        }
    )
    rescued, rescued_abundance, details = wrapper.apply_candidate_rescue_switch(
        features,
        np.array([False]),
        np.array([3.0]),
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70,
    )
    assert rescued.tolist() == [False]
    assert np.allclose(rescued_abundance, [3.0])
    assert details["candidate_rescue_added_n"] == 0
    assert details["candidate_rescue_ani_source_rule"] == "max(s_Ref_zip_aaf_ani_max,u_Ref_zip_aaf_ani_max)"


def accession_surface_taxmap() -> dict[str, dict[str, str]]:
    return {
        "GCF_000000101.1": {
            "taxid": "9001",
            "rank": "species",
            "taxpath": "",
            "taxpathsn": "d__Bacteria|g__Fixture|s__Fixture alpha_A",
            "species_name": "s__Fixture alpha_A",
        },
        "GCF_000000102.1": {
            "taxid": "9001",
            "rank": "species",
            "taxpath": "",
            "taxpathsn": "d__Bacteria|g__Fixture|s__Fixture alpha_B",
            "species_name": "s__Fixture alpha_B",
        },
    }


def accession_surface_raw_rows() -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = {
        "Ref": ["/ref/GCF_000000102.1.fna"],
        "accession": ["GCF_000000102.1"],
        "taxid": ["9001"],
        "species_name": ["s__Fixture alpha_B"],
        "ANI": [0.97],
        "XnY_ctx": [150],
        "Real_min_align_fraction": [0.80],
        "Ref_breadth": [0.15],
        "Ref_mean_depth": [2.0],
        "Ref_hit_mean_depth": [2.0],
        "Ref_depth_cv": [0.2],
        "Normalized_abundance_depth": [0.1],
        "Ref_zip_af": [0.3],
        "Ref_zip_aaf_ani": [0.96],
    }
    return pd.DataFrame(), pd.DataFrame(cols)


def test_accession_candidate_surface_default_off_preserves_empty_surface() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})
    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_OFF,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
    )
    assert surface.empty
    assert details["candidate_surface_switch"] == wrapper.CANDIDATE_SURFACE_SWITCH_OFF
    assert details["candidate_surface_added_n"] == 0
    assert details["candidate_surface_zero_mass"] is False


def test_accession_candidate_surface_adds_strong_uncalled_same_taxid_species() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})
    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
    )
    assert len(surface) == 1
    assert surface.iloc[0]["candidate_surface_species_name"] == "s__Fixture alpha_B"
    assert surface.iloc[0]["taxid"] == "9001"
    assert details["candidate_surface_added_n"] == 1
    assert details["candidate_surface_zero_mass"] is True


def test_accession_candidate_surface_cross_panel_mode_uses_strict_xny_breadth_without_real_af_gate() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    split_rows.loc[0, "ANI"] = 0.94
    split_rows.loc[0, "XnY_ctx"] = 700
    split_rows.loc[0, "Ref_breadth"] = 0.25
    split_rows.loc[0, "Real_min_align_fraction"] = 0.05
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})

    old_surface, old_details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
    )
    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI93_XNY650_BR20,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
    )

    assert old_surface.empty
    assert old_details["candidate_surface_real_af_min"] == 0.70
    assert len(surface) == 1
    assert surface.iloc[0]["candidate_surface_species_name"] == "s__Fixture alpha_B"
    assert details["candidate_surface_ani_min"] == 0.93
    assert details["candidate_surface_xny_min"] == 650.0
    assert details["candidate_surface_breadth_min"] == 0.20
    assert details["candidate_surface_real_af_min"] == 0.0


def test_accession_candidate_surface_combined_mode_preserves_current_and_adds_strict_rows() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    strict = split_rows.iloc[0].copy()
    strict["Ref"] = "/ref/GCF_000000103.1.fna"
    strict["accession"] = "GCF_000000103.1"
    strict["species_name"] = "s__Fixture alpha_C"
    strict["ANI"] = 0.94
    strict["XnY_ctx"] = 700
    strict["Ref_breadth"] = 0.25
    strict["Real_min_align_fraction"] = 0.05
    split_rows = pd.concat([split_rows, strict.to_frame().T], ignore_index=True)
    taxmap = accession_surface_taxmap()
    taxmap["GCF_000000103.1"] = {
        "taxid": "9001",
        "rank": "species",
        "taxpath": "",
        "taxpathsn": "d__Bacteria|g__Fixture|s__Fixture alpha_C",
        "species_name": "s__Fixture alpha_C",
    }
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})

    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20,
        taxmap,
        {"9001": "s__Fixture alpha_A"},
    )

    assert set(surface["candidate_surface_species_name"]) == {
        "s__Fixture alpha_B",
        "s__Fixture alpha_C",
    }
    assert details["candidate_surface_switch"] == (
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_CURRENT_OR_ANI93_XNY650_BR20
    )
    assert "passing either current surface rule" in details["candidate_surface_rule"]


def test_accession_candidate_surface_max_called_species_guard_blocks_surface() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    taxmap = accession_surface_taxmap()
    taxmap["GCF_000000104.1"] = {
        "taxid": "9002",
        "rank": "species",
        "taxpath": "",
        "taxpathsn": "d__Bacteria|g__Fixture|s__Fixture beta_A",
        "species_name": "s__Fixture beta_A",
    }
    features = pd.DataFrame(
        {
            "taxid": ["9001", "9002"],
            "s_best_accession": ["GCF_000000101.1", "GCF_000000104.1"],
        }
    )

    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True, True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        taxmap,
        {"9001": "s__Fixture alpha_A", "9002": "s__Fixture beta_A"},
        max_called_species=1,
    )

    assert surface.empty
    assert details["candidate_surface_guard_blocked"] is True
    assert details["candidate_surface_called_species_n"] == 2
    assert details["candidate_surface_max_called_species"] == 1
    assert "called species count > 1" in details["candidate_surface_rule"]


def test_accession_candidate_surface_appended_rows_are_called_zero_mass() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})
    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
    )
    base = pd.DataFrame(
        {
            "taxid": ["9001"],
            "species_name": ["s__Fixture alpha_A"],
            "calibrated_call": [True],
            "calibrated_probability": [0.9],
            "calibrated_abundance": [1.0],
            "calibrated_abundance_raw": [5.0],
            "reported_ani": [0.98],
            "reported_ani_source": ["s_Ref_zip_aaf_ani_max"],
            **details,
        }
    )
    out = wrapper.append_accession_candidate_surface_rows(base, surface)
    added = out.loc[out["candidate_surface_added"].astype(bool)]
    assert len(added) == 1
    assert bool(added.iloc[0]["calibrated_call"]) is True
    assert added.iloc[0]["calibrated_abundance"] == 0.0
    assert added.iloc[0]["calibrated_abundance_raw"] == 0.0
    assert added.iloc[0]["s_best_accession"] == "GCF_000000102.1"


def test_accession_candidate_surface_appended_rows_can_use_normalized_depth_mass() -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    features = pd.DataFrame({"taxid": ["9001"], "s_best_accession": ["GCF_000000101.1"]})
    surface, details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    )
    base = pd.DataFrame(
        {
            "taxid": ["9001"],
            "species_name": ["s__Fixture alpha_A"],
            "calibrated_call": [True],
            "calibrated_probability": [0.9],
            "calibrated_abundance": [1.0],
            "calibrated_abundance_raw": [5.0],
            "reported_ani": [0.98],
            "reported_ani_source": ["s_Ref_zip_aaf_ani_max"],
            **details,
        }
    )
    out = wrapper.append_accession_candidate_surface_rows(
        base,
        surface,
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    )
    added = out.loc[out["candidate_surface_added"].astype(bool)]
    assert len(added) == 1
    assert details["candidate_surface_zero_mass"] is False
    assert added.iloc[0]["calibrated_abundance_raw"] == pytest.approx(0.2)
    assert bool(added.iloc[0]["candidate_surface_zero_mass"]) is False
    assert added.iloc[0]["candidate_abundance_policy"] == wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2


def test_write_output_normalizes_candidate_surface_mass(work: Path) -> None:
    unique_rows, split_rows = accession_surface_raw_rows()
    features = pd.DataFrame(
        {
            "taxid": ["9001"],
            "calibrated_probability": [0.9],
            "s_best_accession": ["GCF_000000101.1"],
            "s_Ref_zip_aaf_ani_max": [0.98],
            "u_Ref_zip_aaf_ani_max": [0.0],
        }
    )
    surface, surface_details = wrapper.build_accession_candidate_surface(
        unique_rows,
        split_rows,
        features,
        np.array([True]),
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70,
        accession_surface_taxmap(),
        {"9001": "s__Fixture alpha_A"},
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    )
    out = work / "profile.tsv"
    strategy_details = {
        **surface_details,
        "candidate_abundance_policy": wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
        "candidate_abundance_policy_alpha": 2.0,
        "candidate_abundance_policy_rule": wrapper.candidate_abundance_policy_rule(
            wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
        ),
    }
    wrapper.write_output(
        out,
        features,
        0.35,
        "train12",
        "all",
        {"9001": "s__Fixture alpha_A"},
        False,
        "universal-auto-exact",
        np.array([True]),
        np.array([5.0]),
        strategy_details,
        candidate_surface_rows=surface,
        candidate_abundance_policy=wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2,
    )
    result = pd.read_csv(out, sep="\t")
    added = result.loc[result["candidate_surface_added"].astype(str) == "True"].iloc[0]
    base = result.loc[result["species_name"] == "s__Fixture alpha_A"].iloc[0]
    assert added["candidate_abundance_norm_mass"] == pytest.approx(0.2)
    assert added["calibrated_abundance_raw"] == pytest.approx(1.0)
    assert added["calibrated_abundance"] == pytest.approx(1.0 / 6.0)
    assert base["calibrated_abundance"] == pytest.approx(5.0 / 6.0)
    assert result["calibrated_abundance"].sum() == pytest.approx(1.0)


def test_initial_pass_thread_budget_is_split_without_exceeding_total() -> None:
    assert wrapper.split_initial_pass_threads(1) == (1, 1)
    assert wrapper.split_initial_pass_threads(2) == (1, 1)
    assert wrapper.split_initial_pass_threads(3) == (1, 2)
    assert wrapper.split_initial_pass_threads(16) == (8, 8)


def test_minco_pass_command_can_request_exact_split_sidecar(work: Path) -> None:
    exact = work / "exact.tsv"
    cmd = wrapper.minco_pass_command(
        "minco",
        work / "ref.minco",
        work / "reads.fq",
        work / "split.tsv",
        "best-diff-split",
        4,
        "",
        exact_split_sidecar_out=exact,
    )
    assert "--readwise-exact-split-out" in cmd
    assert cmd[cmd.index("--readwise-exact-split-out") + 1] == str(exact)


def test_same_stream_initial_pass_uses_unique_sidecar(work: Path) -> None:
    calls: list[tuple[str, int, str, str, object]] = []
    old_run = wrapper.run_minco_pass

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
    ):
        assert density_block_ctx is None
        calls.append(
            (
                assign_mode,
                threads,
                Path(out_path).name,
                Path(unique_sidecar_out).name if unique_sidecar_out else "",
                density_block_ctx,
            )
        )

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": True,
            "parallel_readwise_passes": False,
            "minco": "minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 7,
            "pipecmd": "",
        },
    )()
    wrapper.run_minco_pass = fake_run
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_pass = old_run
    assert calls == [
        (
            "best-diff-split",
            7,
            "minco.best_diff_split.unfiltered.tsv",
            "minco.best_diff_unique.unfiltered.tsv",
            None,
        ),
    ]


def test_same_stream_initial_pass_can_write_exact_split_sidecar(work: Path) -> None:
    calls: list[tuple[str, str, str]] = []
    old_run = wrapper.run_minco_pass

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        _threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
        exact_split_sidecar_out=None,
    ):
        assert density_block_ctx is None
        calls.append(
            (
                assign_mode,
                Path(unique_sidecar_out).name if unique_sidecar_out else "",
                Path(exact_split_sidecar_out).name if exact_split_sidecar_out else "",
            )
        )

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": True,
            "parallel_readwise_passes": False,
            "minco": "minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 7,
            "pipecmd": "",
        },
    )()
    wrapper.run_minco_pass = fake_run
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
            work / "minco.best_diff_split.exact.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_pass = old_run
    assert calls == [
        (
            "best-diff-split",
            "minco.best_diff_unique.unfiltered.tsv",
            "minco.best_diff_split.exact.unfiltered.tsv",
        )
    ]


def test_auto_initial_pass_uses_same_stream_at_one_thread(work: Path) -> None:
    calls: list[tuple[str, int, str, str]] = []
    old_run = wrapper.run_minco_pass

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
    ):
        assert density_block_ctx is None
        calls.append(
            (
                assign_mode,
                threads,
                Path(out_path).name,
                Path(unique_sidecar_out).name if unique_sidecar_out else "",
            )
        )

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": None,
            "parallel_readwise_passes": True,
            "minco": "minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 1,
            "pipecmd": "",
        },
    )()
    wrapper.run_minco_pass = fake_run
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_pass = old_run
    assert calls == [
        (
            "best-diff-split",
            1,
            "minco.best_diff_split.unfiltered.tsv",
            "minco.best_diff_unique.unfiltered.tsv",
        )
    ]


def test_legacy_sequential_initial_passes_use_requested_threads(work: Path) -> None:
    calls: list[tuple[str, int, str]] = []
    old_run = wrapper.run_minco_pass

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
    ):
        assert density_block_ctx is None
        assert unique_sidecar_out is None
        calls.append((assign_mode, threads, Path(out_path).name))

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": None,
            "parallel_readwise_passes": False,
            "minco": "minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 7,
            "pipecmd": "",
        },
    )()
    wrapper.run_minco_pass = fake_run
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_pass = old_run
    assert calls == [
        ("best-diff-unique", 7, "minco.best_diff_unique.unfiltered.tsv"),
        ("best-diff-split", 7, "minco.best_diff_split.unfiltered.tsv"),
    ]


def test_parallel_initial_passes_delegate_to_parallel_helper(work: Path) -> None:
    calls: list[tuple[str, str, str, str, int, str]] = []
    old_parallel = wrapper.run_minco_passes_parallel

    def fake_parallel(minco, ref, reads, unique_out, split_out, threads, pipecmd):
        calls.append(
            (
                minco,
                str(ref),
                str(reads),
                f"{Path(unique_out).name}|{Path(split_out).name}",
                threads,
                pipecmd,
            )
        )

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": None,
            "parallel_readwise_passes": True,
            "minco": "bin/minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 16,
            "pipecmd": "zcat",
        },
    )()
    wrapper.run_minco_passes_parallel = fake_parallel
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_passes_parallel = old_parallel
    assert calls == [
        (
            "bin/minco",
            str(work / "ref.minco"),
            str(work / "reads.fq"),
            "minco.best_diff_unique.unfiltered.tsv|minco.best_diff_split.unfiltered.tsv",
            16,
            "zcat",
        )
    ]


def test_parallel_initial_passes_can_forward_exact_sidecar(work: Path) -> None:
    calls: list[str] = []
    old_parallel = wrapper.run_minco_passes_parallel

    def fake_parallel(
        _minco,
        _ref,
        _reads,
        _unique_out,
        _split_out,
        _threads,
        _pipecmd,
        exact_split_sidecar_out=None,
    ):
        calls.append(Path(exact_split_sidecar_out).name if exact_split_sidecar_out else "")

    args = type(
        "Args",
        (),
        {
            "same_stream_readwise_passes": None,
            "parallel_readwise_passes": True,
            "minco": "bin/minco",
            "ref": work / "ref.minco",
            "reads": work / "reads.fq",
            "threads": 16,
            "pipecmd": "",
        },
    )()
    wrapper.run_minco_passes_parallel = fake_parallel
    try:
        wrapper.run_initial_minco_passes(
            args,
            work / "minco.best_diff_unique.unfiltered.tsv",
            work / "minco.best_diff_split.unfiltered.tsv",
            work / "minco.best_diff_split.exact.unfiltered.tsv",
        )
    finally:
        wrapper.run_minco_passes_parallel = old_parallel
    assert calls == ["minco.best_diff_split.exact.unfiltered.tsv"]


def test_exact_split_guard_rejects_sample4_like_case() -> None:
    assert not wrapper.exact_split_guard_passes(
        {
            "joined_base_median_uaf": 0.3565,
            "raw_unique_to_base_ratio": 0.74,
        }
    )


def test_exact_split_guard_accepts_sparse_or_strong_unique_cases() -> None:
    assert wrapper.exact_split_guard_passes(
        {
            "joined_base_median_uaf": 0.34,
            "raw_unique_to_base_ratio": 0.50,
        }
    )
    assert wrapper.exact_split_guard_passes(
        {
            "joined_base_median_uaf": 0.60,
            "raw_unique_to_base_ratio": 0.94,
        }
    )


def test_table_mode_trigger_falls_back(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    unique = work / "unique.tsv"
    split = work / "split.tsv"
    out = work / "table_mode.tsv"
    write_taxmap(taxmap)
    write_minco_tables(unique, split, extra_mass=1.0)
    rc = wrapper.main(
        [
            "--unique-table",
            str(unique),
            "--split-table",
            str(split),
            "--taxmap",
            str(taxmap),
            "--train-table",
            str(work / "unused.tsv"),
            "--strategy",
            "universal-auto-exact",
            "--report-all",
            "-o",
            str(out),
        ]
    )
    assert rc == 0
    df = read_output(out)
    assert_bool_column(df, "auto_exact_split_requested", True)
    assert_bool_column(df, "auto_exact_split_used", False)
    assert set(df["auto_exact_split_unavailable_reason"]) == {"missing_ref_or_reads"}
    assert set(df["profile_strategy"]) == {"universal-auto-exact"}
    assert {"u_best_accession", "s_best_accession", "s_best_ref"}.issubset(df.columns)
    assert {"reported_ani", "reported_ani_source"}.issubset(df.columns)
    row = df.loc[df["taxid"].astype(str) == "1021"].iloc[0]
    assert row["s_best_accession"] == accession(21)
    assert np.isclose(float(row["reported_ani"]), 0.98)
    assert row["reported_ani_source"] == "s_Ref_zip_aaf_ani_max"


def test_default_strategy_matches_explicit_autoexact(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    unique = work / "unique.tsv"
    split = work / "split.tsv"
    default_out = work / "default.tsv"
    explicit_out = work / "explicit.tsv"
    write_taxmap(taxmap)
    write_minco_tables(unique, split, extra_mass=1.0)
    base_args = [
        "--unique-table",
        str(unique),
        "--split-table",
        str(split),
        "--taxmap",
        str(taxmap),
        "--train-table",
        str(work / "unused.tsv"),
        "--report-all",
    ]
    rc = wrapper.main(base_args + ["-o", str(default_out)])
    assert rc == 0
    rc = wrapper.main(base_args + ["--strategy", "universal-auto-exact", "-o", str(explicit_out)])
    assert rc == 0

    default_df = read_output(default_out)
    explicit_df = read_output(explicit_out)
    assert set(default_df["profile_strategy"]) == {"universal-auto-exact"}
    pd.testing.assert_frame_equal(default_df, explicit_df, check_dtype=False)


def test_default_launcher_candidate_preset_runs_table_mode(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    unique = work / "unique.tsv"
    split = work / "split.tsv"
    out = work / "candidate.tsv"
    write_taxmap(taxmap)
    write_minco_tables(unique, split, extra_mass=1.0)

    rc = default_wrapper.main(
        [
            "--profile-preset",
            "candidate",
            "--unique-table",
            str(unique),
            "--split-table",
            str(split),
            "--taxmap",
            str(taxmap),
            "--train-table",
            str(work / "unused.tsv"),
            "--report-all",
            "-o",
            str(out),
        ],
        {},
    )
    assert rc == 0
    df = read_output(out)
    assert set(df["profile_strategy"]) == {"universal-auto-exact"}
    assert set(df["candidate_rescue_switch"]) == {
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70
    }
    assert set(df["candidate_surface_switch"]) == {
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
    }
    assert set(df["candidate_abundance_policy"]) == {
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    }


def test_candidate_surface_disables_without_gtdb_label_taxmap(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "ncbi_taxmap.tsv"
    unique = work / "unique.tsv"
    split = work / "split.tsv"
    out = work / "candidate.tsv"
    write_ncbi_style_taxmap(taxmap)
    write_minco_tables(unique, split, extra_mass=1.0)

    rc = default_wrapper.main(
        [
            "--profile-preset",
            "candidate",
            "--unique-table",
            str(unique),
            "--split-table",
            str(split),
            "--taxmap",
            str(taxmap),
            "--train-table",
            str(work / "unused.tsv"),
            "--report-all",
            "-o",
            str(out),
        ],
        {},
    )
    assert rc == 0
    df = read_output(out)
    assert set(df["candidate_surface_switch"]) == {wrapper.CANDIDATE_SURFACE_SWITCH_OFF}
    assert set(df["candidate_surface_taxmap"].fillna("")) == {""}


def test_raw_mode_trigger_runs_exact(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    reads = work / "reads.fq"
    ref = work / "ref.minco"
    out = work / "raw_exact.tsv"
    run_dir = work / "raw_exact_work"
    run_dir.mkdir()
    reads.write_text("@r1\nACGT\n+\n!!!!\n")
    ref.mkdir()
    write_taxmap(taxmap)
    calls: list[tuple[str, int | None, str]] = []
    old_run = wrapper.run_minco_pass
    old_initial = wrapper.run_initial_minco_passes

    def fake_initial(_args, unique_out, split_out):
        calls.append(("best-diff-unique", None, Path(unique_out).name))
        calls.append(("best-diff-split", None, Path(split_out).name))
        unique, split = minco_rows(extra_mass=1.0)
        unique.to_csv(unique_out, sep="\t", index=False)
        split.to_csv(split_out, sep="\t", index=False)

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        _threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
    ):
        calls.append((assign_mode, density_block_ctx, Path(out_path).name))
        if assign_mode == "best-diff-unique":
            unique, _split = minco_rows(extra_mass=1.0)
            unique.to_csv(out_path, sep="\t", index=False)
        else:
            _unique, split = minco_rows(extra_mass=1.0)
            split.to_csv(out_path, sep="\t", index=False)

    wrapper.run_initial_minco_passes = fake_initial
    wrapper.run_minco_pass = fake_run
    try:
        rc = wrapper.main(
            [
                "-r",
                str(ref),
                "--reads",
                str(reads),
                "--taxmap",
                str(taxmap),
                "--train-table",
                str(work / "unused.tsv"),
                "--strategy",
                "universal-auto-exact",
                "--workdir",
                str(run_dir),
                "--report-all",
                "-o",
                str(out),
            ]
        )
    finally:
        wrapper.run_initial_minco_passes = old_initial
        wrapper.run_minco_pass = old_run
    assert rc == 0
    assert ("best-diff-split", 0, "minco.best_diff_split.exact.unfiltered.tsv") in calls
    df = read_output(out)
    assert_bool_column(df, "auto_exact_split_requested", True)
    assert_bool_column(df, "auto_exact_split_used", True)
    assert set(df["auto_exact_split_path"]) == {str(run_dir / "minco.best_diff_split.exact.unfiltered.tsv")}


def test_raw_mode_low_extra_rescue_skips_exact_by_default(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    reads = work / "reads.fq"
    ref = work / "ref.minco"
    out = work / "raw_low_extra_skip.tsv"
    run_dir = work / "raw_low_extra_skip_work"
    run_dir.mkdir()
    reads.write_text("@r1\nACGT\n+\n!!!!\n")
    ref.mkdir()
    write_taxmap(taxmap)
    calls: list[tuple[str, int | None, str]] = []
    old_run = wrapper.run_minco_pass
    old_initial = wrapper.run_initial_minco_passes

    def fake_initial(_args, unique_out, split_out):
        calls.append(("best-diff-unique", None, Path(unique_out).name))
        calls.append(("best-diff-split", None, Path(split_out).name))
        unique, split = minco_rows_low_extra()
        unique.to_csv(unique_out, sep="\t", index=False)
        split.to_csv(split_out, sep="\t", index=False)

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        _threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
        exact_split_sidecar_out=None,
    ):
        calls.append((assign_mode, density_block_ctx, Path(out_path).name))
        unique, split = minco_rows_low_extra()
        (unique if assign_mode == "best-diff-unique" else split).to_csv(out_path, sep="\t", index=False)

    wrapper.run_initial_minco_passes = fake_initial
    wrapper.run_minco_pass = fake_run
    try:
        rc = wrapper.main(
            [
                "-r",
                str(ref),
                "--reads",
                str(reads),
                "--taxmap",
                str(taxmap),
                "--train-table",
                str(work / "unused.tsv"),
                "--strategy",
                "universal-auto-exact",
                "--workdir",
                str(run_dir),
                "--report-all",
                "-o",
                str(out),
            ]
        )
    finally:
        wrapper.run_initial_minco_passes = old_initial
        wrapper.run_minco_pass = old_run
    assert rc == 0
    assert all(call[1] != 0 for call in calls), calls
    df = read_output(out)
    assert_bool_column(df, "auto_exact_split_requested", False)
    assert_bool_column(df, "auto_exact_split_used", False)
    assert set(df["auto_exact_split_unavailable_reason"]) == {"block_low_extra_split_rescue_already_active"}
    assert set(df["auto_exact_split_low_extra_mode"]) == {"skip"}
    assert set(df["auto_exact_split_low_extra_added_n"]) == {1}
    assert_bool_column(df, "auto_exact_split_low_extra_gate_passed", False)


def test_raw_mode_low_extra_rescue_allow_runs_exact(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    reads = work / "reads.fq"
    ref = work / "ref.minco"
    out = work / "raw_low_extra_allow.tsv"
    run_dir = work / "raw_low_extra_allow_work"
    run_dir.mkdir()
    reads.write_text("@r1\nACGT\n+\n!!!!\n")
    ref.mkdir()
    write_taxmap(taxmap)
    calls: list[tuple[str, int | None, str]] = []
    old_run = wrapper.run_minco_pass
    old_initial = wrapper.run_initial_minco_passes

    def fake_initial(_args, unique_out, split_out):
        calls.append(("best-diff-unique", None, Path(unique_out).name))
        calls.append(("best-diff-split", None, Path(split_out).name))
        unique, split = minco_rows_low_extra()
        unique.to_csv(unique_out, sep="\t", index=False)
        split.to_csv(split_out, sep="\t", index=False)

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        _threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
        exact_split_sidecar_out=None,
    ):
        calls.append((assign_mode, density_block_ctx, Path(out_path).name))
        unique, split = minco_rows_low_extra()
        (unique if assign_mode == "best-diff-unique" else split).to_csv(out_path, sep="\t", index=False)

    wrapper.run_initial_minco_passes = fake_initial
    wrapper.run_minco_pass = fake_run
    try:
        rc = wrapper.main(
            [
                "-r",
                str(ref),
                "--reads",
                str(reads),
                "--taxmap",
                str(taxmap),
                "--train-table",
                str(work / "unused.tsv"),
                "--strategy",
                "universal-auto-exact",
                "--exact-split-low-extra-mode",
                "allow",
                "--workdir",
                str(run_dir),
                "--report-all",
                "-o",
                str(out),
            ]
        )
    finally:
        wrapper.run_initial_minco_passes = old_initial
        wrapper.run_minco_pass = old_run
    assert rc == 0
    assert ("best-diff-split", 0, "minco.best_diff_split.exact.unfiltered.tsv") in calls
    df = read_output(out)
    assert_bool_column(df, "auto_exact_split_requested", True)
    assert_bool_column(df, "auto_exact_split_used", True)
    assert set(df["auto_exact_split_low_extra_mode"]) == {"allow"}
    assert_bool_column(df, "auto_exact_split_low_extra_gate_passed", True)


def test_raw_mode_trigger_off_skips_exact(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    taxmap = work / "taxmap.tsv"
    reads = work / "reads.fq"
    ref = work / "ref.minco"
    out = work / "raw_no_exact.tsv"
    run_dir = work / "raw_no_exact_work"
    run_dir.mkdir()
    reads.write_text("@r1\nACGT\n+\n!!!!\n")
    ref.mkdir()
    write_taxmap(taxmap)
    calls: list[tuple[str, int | None, str]] = []
    old_run = wrapper.run_minco_pass
    old_initial = wrapper.run_initial_minco_passes

    def fake_initial(_args, unique_out, split_out):
        calls.append(("best-diff-unique", None, Path(unique_out).name))
        calls.append(("best-diff-split", None, Path(split_out).name))
        unique, split = minco_rows(extra_mass=3.0)
        unique.to_csv(unique_out, sep="\t", index=False)
        split.to_csv(split_out, sep="\t", index=False)

    def fake_run(
        _minco,
        _ref,
        _reads,
        out_path,
        assign_mode,
        _threads,
        _pipecmd,
        density_block_ctx=None,
        unique_sidecar_out=None,
    ):
        calls.append((assign_mode, density_block_ctx, Path(out_path).name))
        if assign_mode == "best-diff-unique":
            unique, _split = minco_rows(extra_mass=3.0)
            unique.to_csv(out_path, sep="\t", index=False)
        else:
            _unique, split = minco_rows(extra_mass=3.0)
            split.to_csv(out_path, sep="\t", index=False)

    wrapper.run_initial_minco_passes = fake_initial
    wrapper.run_minco_pass = fake_run
    try:
        rc = wrapper.main(
            [
                "-r",
                str(ref),
                "--reads",
                str(reads),
                "--taxmap",
                str(taxmap),
                "--train-table",
                str(work / "unused.tsv"),
                "--strategy",
                "universal-auto-exact",
                "--workdir",
                str(run_dir),
                "--report-all",
                "-o",
                str(out),
            ]
        )
    finally:
        wrapper.run_initial_minco_passes = old_initial
        wrapper.run_minco_pass = old_run
    assert rc == 0
    assert all(call[1] != 0 for call in calls), calls
    assert not (run_dir / "minco.best_diff_split.exact.unfiltered.tsv").exists()
    df = read_output(out)
    assert_bool_column(df, "auto_exact_split_requested", False)
    assert_bool_column(df, "auto_exact_split_used", False)


def test_model_cache_roundtrip_and_metadata_validation(work: Path) -> None:
    cache = work / "minco_profile_rf_hgb.train12.unfiltered.joblib"
    models = ["rf-model", "hgb-model"]
    wrapper.write_model_cache(cache, models, "train12", "all", False)

    loaded = wrapper.load_model_cache(cache, "train12", "all", False)
    assert list(loaded) == models

    with pytest.raises(SystemExit, match="train_pool"):
        wrapper.load_model_cache(cache, "train9", "all", False)
    with pytest.raises(SystemExit, match="scope"):
        wrapper.load_model_cache(cache, "train12", "bacteria", False)
    with pytest.raises(SystemExit, match="filter_training_scope"):
        wrapper.load_model_cache(cache, "train12", "all", True)


def test_load_or_fit_models_prefers_discovered_train_feature_cache(work: Path) -> None:
    training = work / "joined_feature_training"
    training.mkdir(parents=True)
    cache = training / "minco_profile_rf_hgb.train12.unfiltered.joblib"
    wrapper.write_model_cache(cache, ["rf-cache", "hgb-cache"], "train12", "all", False)
    args = type(
        "Args",
        (),
        {
            "model_cache": None,
            "train_features": training,
            "train_table": [],
            "train_pool": "train12",
            "scope": "all",
            "filter_training_scope": False,
            "write_model_cache": None,
        },
    )()

    loaded = wrapper.load_or_fit_models(args, {})
    assert list(loaded) == ["rf-cache", "hgb-cache"]


def test_default_launcher_injects_universal_autoexact_and_env_defaults() -> None:
    args = default_wrapper.build_calibrated_argv(
        ["-r", "ref.minco", "--reads", "reads.fq.gz", "-o", "out.tsv"],
        {
            default_wrapper.ENV_TAXMAP: "ref.species.taxmap.tsv",
            default_wrapper.ENV_TRAIN_FEATURES: "joined_training",
            default_wrapper.ENV_MINCO: "bin/minco",
        },
    )
    assert args[:5] == ["-r", "ref.minco", "--reads", "reads.fq.gz", "-o"]
    assert args[5] == "out.tsv"
    assert args.count("--strategy") == 1
    assert args[args.index("--strategy") + 1] == "universal-auto-exact"
    assert args[args.index("--taxmap") + 1] == "ref.species.taxmap.tsv"
    assert args[args.index("--train-features") + 1] == "joined_training"
    assert args[args.index("--minco") + 1] == "bin/minco"
    assert args[args.index("--candidate-rescue-switch") + 1] == (
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70
    )
    assert args[args.index("--candidate-surface-switch") + 1] == (
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
    )
    assert args[args.index("--candidate-abundance-policy") + 1] == (
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    )


def test_default_launcher_constants_match_current_strategy_manifest() -> None:
    manifest_path = (
        ROOT
        / "research/experiments/2026-06-27_universal_strategy_decision/results/current_default_strategy_manifest.tsv"
    )
    manifest = pd.read_csv(manifest_path, sep="\t").set_index("item")["value"].to_dict()

    assert manifest["selected_entrypoint"] == "scripts/minco_profile"
    assert manifest["selected_profile_preset"] == default_wrapper.DEFAULT_PRESET
    assert default_wrapper.DEFAULT_PRESET == default_wrapper.CANDIDATE_PRESET
    assert manifest["base_strategy"] == default_wrapper.DEFAULT_STRATEGY
    assert manifest["candidate_rescue_switch"] == (
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70
    )
    assert manifest["candidate_surface_switch"] == (
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
    )
    assert manifest["candidate_abundance_policy"] == (
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    )
    assert wrapper.DEFAULT_ABUNDANCE_GENUS_XNY_BLEND_ALPHA == 0.0
    assert wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_OFF == "off"
    assert (
        wrapper.ABUNDANCE_FEATURE_ALLOCATOR_SWITCH_GUARDED_GENUS_HIT_BREADTH_A002_XNY230
        == "guarded-genus-hit-breadth-a002-xny230"
    )
    assert wrapper.ADAPTIVE_CALL_FILTER_SWITCH_OFF == "off"


def test_default_launcher_current_preset_reproduces_previous_default() -> None:
    args = default_wrapper.build_calibrated_argv(
        [
            "--profile-preset",
            "current",
            "-r",
            "ref.minco",
            "--reads",
            "reads.fq.gz",
            "-o",
            "out.tsv",
        ],
        {
            default_wrapper.ENV_TAXMAP: "ref.species.taxmap.tsv",
            default_wrapper.ENV_TRAIN_FEATURES: "joined_training",
        },
    )
    assert "--profile-preset" not in args
    assert args[args.index("--strategy") + 1] == "universal-auto-exact"
    assert "--candidate-rescue-switch" not in args
    assert "--candidate-surface-switch" not in args
    assert "--candidate-abundance-policy" not in args


def test_default_launcher_discovers_packaged_ref_sidecars(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = work / "ref.species_taxmap.tsv"
    training = work / "joined_feature_training"
    training.mkdir()
    write_taxmap(taxmap)
    (training / "train.joined_features.tsv").write_text("label\n1\n")

    args = default_wrapper.build_calibrated_argv(
        ["-r", str(ref), "--reads", "reads.fq.gz", "-o", "out.tsv"],
        {},
    )
    assert args[args.index("--strategy") + 1] == "universal-auto-exact"
    assert args[args.index("--taxmap") + 1] == str(taxmap)
    assert args[args.index("--train-features") + 1] == str(training)
    assert args[args.index("--candidate-rescue-switch") + 1] == (
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70
    )


def test_default_launcher_discovers_packaged_model_cache_without_training(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = work / "ref.species_taxmap.tsv"
    cache = ref / "minco_profile_rf_hgb.train12.unfiltered.joblib"
    write_taxmap(taxmap)
    cache.write_text("placeholder\n")

    args = default_wrapper.build_calibrated_argv(
        ["-r", str(ref), "--reads", "reads.fq.gz", "-o", "out.tsv"],
        {},
    )
    assert args[args.index("--strategy") + 1] == "universal-auto-exact"
    assert args[args.index("--taxmap") + 1] == str(taxmap)
    assert args[args.index("--model-cache") + 1] == str(cache)
    assert args[args.index("--candidate-abundance-policy") + 1] == (
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    )
    assert "--train-features" not in args


def test_default_launcher_applies_packaged_scope_default(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = ref / "species_taxmap.tsv"
    cache = ref / "minco_profile_rf_hgb.train12.unfiltered.joblib"
    defaults = ref / "minco_profile_defaults.tsv"
    write_taxmap(taxmap)
    cache.write_text("placeholder\n")
    defaults.write_text("option\tvalue\nscope\tbacteria\n")

    args = default_wrapper.build_calibrated_argv(
        ["-r", str(ref), "--reads", "reads.fq.gz", "-o", "out.tsv"],
        {},
    )
    assert args[args.index("--scope") + 1] == "bacteria"
    assert args[args.index("--model-cache") + 1] == str(cache)


def test_default_launcher_explicit_scope_overrides_packaged_default(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    write_taxmap(ref / "species_taxmap.tsv")
    (ref / "minco_profile_rf_hgb.train12.unfiltered.joblib").write_text("placeholder\n")
    (ref / "minco_profile_defaults.tsv").write_text("option\tvalue\nscope\tbacteria\n")

    args = default_wrapper.build_calibrated_argv(
        ["-r", str(ref), "--reads", "reads.fq.gz", "--scope", "all", "-o", "out.tsv"],
        {},
    )
    assert args.count("--scope") == 1
    assert args[args.index("--scope") + 1] == "all"


def test_default_launcher_candidate_preset_injects_candidate_switches_and_sidecar(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = work / "ref.species_taxmap.tsv"
    surface_taxmap = work / "ref.candidate_surface_taxmap.tsv"
    training = work / "joined_feature_training"
    training.mkdir()
    write_taxmap(taxmap)
    write_taxmap(surface_taxmap)
    (training / "train.joined_features.tsv").write_text("label\n1\n")

    args = default_wrapper.build_calibrated_argv(
        [
            "--profile-preset",
            "candidate",
            "-r",
            str(ref),
            "--reads",
            "reads.fq.gz",
            "-o",
            "out.tsv",
        ],
        {},
    )
    assert "--profile-preset" not in args
    assert args[args.index("--strategy") + 1] == "universal-auto-exact"
    assert args[args.index("--candidate-rescue-switch") + 1] == (
        wrapper.CANDIDATE_RESCUE_SWITCH_EMITTED_ANI90_XNY100_BR01_AF70
    )
    assert args[args.index("--candidate-surface-switch") + 1] == (
        wrapper.CANDIDATE_SURFACE_SWITCH_ACCESSION_ANI90_XNY100_BR01_AF70
    )
    assert args[args.index("--candidate-abundance-policy") + 1] == (
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    )
    assert args[args.index("--candidate-surface-taxmap") + 1] == str(surface_taxmap)


def test_default_launcher_candidate_preset_can_use_env_surface_taxmap(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = work / "ref.species_taxmap.tsv"
    training = work / "joined_feature_training"
    training.mkdir()
    write_taxmap(taxmap)
    (training / "train.joined_features.tsv").write_text("label\n1\n")

    args = default_wrapper.build_calibrated_argv(
        ["-r", str(ref), "--reads", "reads.fq.gz", "-o", "out.tsv"],
        {
            default_wrapper.ENV_PROFILE_PRESET: "candidate",
            default_wrapper.ENV_CANDIDATE_SURFACE_TAXMAP: "env.surface.taxmap.tsv",
        },
    )
    assert args[args.index("--candidate-surface-taxmap") + 1] == "env.surface.taxmap.tsv"
    assert args[args.index("--candidate-abundance-policy") + 1] == (
        wrapper.CANDIDATE_ABUNDANCE_POLICY_NORMALIZED_DEPTH_ALPHA2
    )


def test_default_launcher_candidate_preset_respects_surface_off_with_sidecar(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    taxmap = work / "ref.species_taxmap.tsv"
    surface_taxmap = work / "candidate_surface_taxmap.tsv"
    training = work / "joined_feature_training"
    training.mkdir()
    write_taxmap(taxmap)
    write_taxmap(surface_taxmap)
    (training / "train.joined_features.tsv").write_text("label\n1\n")

    args = default_wrapper.build_calibrated_argv(
        [
            "--profile-preset",
            "candidate",
            "-r",
            str(ref),
            "--reads",
            "reads.fq.gz",
            "--candidate-surface-switch",
            "off",
            "-o",
            "out.tsv",
        ],
        {},
    )
    assert args[args.index("--candidate-surface-switch") + 1] == "off"
    assert "--candidate-surface-taxmap" not in args


def test_default_launcher_preflight_passes_packaged_ref(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    reads = work / "reads.fq.gz"
    minco = work / "minco"
    write_taxmap(ref / "species_taxmap.tsv")
    write_taxmap(ref / "candidate_surface_taxmap.tsv")
    (ref / "minco_profile_defaults.tsv").write_text("option\tvalue\nscope\tbacteria\n")
    wrapper.write_model_cache(
        ref / "minco_profile_rf_hgb.train12.unfiltered.joblib",
        ["rf-cache", "hgb-cache"],
        "train12",
        "bacteria",
        False,
    )
    reads.write_text("placeholder\n")
    minco.write_text("#!/usr/bin/env sh\nexit 0\n")
    minco.chmod(0o755)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = default_wrapper.run_preflight(
            ["-r", str(ref), "--reads", str(reads), "--minco", str(minco)],
            {},
        )

    text = out.getvalue()
    assert rc == 0
    assert "status\tpass" in text
    assert "species_taxmap\tpass\tsidecar" in text
    assert "model_cache\tpass\tsidecar" in text
    assert "model_cache_metadata\tpass\tsidecar" in text
    assert "candidate_surface_taxmap\tpass\tsidecar" in text
    assert "default_scope\tpass\tsidecar\tbacteria" in text
    assert "--candidate-surface-taxmap" in text
    assert "--scope bacteria" in text


def test_default_launcher_preflight_fails_missing_candidate_surface(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    reads = work / "reads.fq.gz"
    minco = work / "minco"
    write_taxmap(ref / "species_taxmap.tsv")
    wrapper.write_model_cache(
        ref / "minco_profile_rf_hgb.train12.unfiltered.joblib",
        ["rf-cache", "hgb-cache"],
        "train12",
        "all",
        False,
    )
    reads.write_text("placeholder\n")
    minco.write_text("#!/usr/bin/env sh\nexit 0\n")
    minco.chmod(0o755)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = default_wrapper.run_preflight(
            ["-r", str(ref), "--reads", str(reads), "--minco", str(minco)],
            {},
        )

    text = out.getvalue()
    assert rc == 1
    assert "status\tfail" in text
    assert "candidate_surface_taxmap\tmissing" in text
    assert "error\tcandidate_surface_taxmap is required" in text


def test_default_launcher_preflight_current_preset_does_not_require_surface(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "ref.minco"
    ref.mkdir()
    minco = work / "minco"
    write_taxmap(ref / "species_taxmap.tsv")
    wrapper.write_model_cache(
        ref / "minco_profile_rf_hgb.train12.unfiltered.joblib",
        ["rf-cache", "hgb-cache"],
        "train12",
        "all",
        False,
    )
    minco.write_text("#!/usr/bin/env sh\nexit 0\n")
    minco.chmod(0o755)

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = default_wrapper.run_preflight(
            ["--profile-preset", "current", "-r", str(ref), "--minco", str(minco)],
            {},
        )

    text = out.getvalue()
    assert rc == 0
    assert "profile_preset\tpass\twrapper\tcurrent" in text
    assert "candidate_surface_switch\tpass\tpreset\toff" in text
    assert "candidate_surface_taxmap" not in text


def test_default_launcher_preserves_explicit_advanced_options() -> None:
    args = default_wrapper.build_calibrated_argv(
        [
            "--profile-preset",
            "candidate",
            "--strategy",
            "probability",
            "--taxmap",
            "explicit.taxmap.tsv",
            "--train-table",
            "explicit.training.tsv",
            "--candidate-rescue-switch",
            "off",
            "--candidate-surface-switch",
            "off",
            "--candidate-abundance-policy",
            "zero",
            "--minco",
            "explicit-minco",
            "-o",
            "out.tsv",
        ],
        {
            default_wrapper.ENV_TAXMAP: "env.taxmap.tsv",
            default_wrapper.ENV_TRAIN_FEATURES: "env.training",
            default_wrapper.ENV_MINCO: "env-minco",
        },
    )
    assert args.count("--strategy") == 1
    assert args[args.index("--strategy") + 1] == "probability"
    assert args[args.index("--taxmap") + 1] == "explicit.taxmap.tsv"
    assert "--train-features" not in args
    assert args[args.index("--train-table") + 1] == "explicit.training.tsv"
    assert args[args.index("--minco") + 1] == "explicit-minco"
    assert args[args.index("--candidate-rescue-switch") + 1] == "off"
    assert args[args.index("--candidate-surface-switch") + 1] == "off"
    assert args[args.index("--candidate-abundance-policy") + 1] == "zero"
    assert "--profile-preset" not in args


def test_default_launcher_requires_taxmap_or_env_unless_help() -> None:
    assert default_wrapper.build_calibrated_argv(["--help"], {}) == ["--help"]
    assert default_wrapper.build_calibrated_argv(["--profile-preset", "candidate", "--help"], {}) == ["--help"]
    with pytest.raises(SystemExit, match=default_wrapper.ENV_TAXMAP):
        default_wrapper.build_calibrated_argv(["-o", "out.tsv"], {})
    with pytest.raises(SystemExit, match="profile-preset"):
        default_wrapper.build_calibrated_argv(["--profile-preset", "unknown", "-o", "out.tsv"], {})


def test_default_launcher_help_banner_names_strategy_and_boundary() -> None:
    text = default_wrapper.default_help_banner()
    assert "recommended no-manual-strategy entry point" in text
    assert "calibrated species profiling" in text
    assert f"--strategy {default_wrapper.DEFAULT_STRATEGY}" in text
    assert "scripts/minco_profile -r ref.minco" in text
    assert "scripts/minco_profile_default.py -r ref.minco" in text
    assert "scripts/minco_profile --check-ref -r ref.minco" in text
    assert "candidate is the selected default" in text
    assert "current reproduces the previous" in text
    assert default_wrapper.ENV_TAXMAP in text
    assert default_wrapper.ENV_TRAIN_FEATURES in text
    assert default_wrapper.ENV_MODEL_CACHE in text
    assert default_wrapper.ENV_MINCO in text
    assert default_wrapper.ENV_PROFILE_PRESET in text
    assert default_wrapper.ENV_CANDIDATE_SURFACE_TAXMAP in text
    assert "--profile-preset current|candidate" in text
    assert "sidecar beside --ref" in text
    assert "AMR/gene/virus/mixed-domain profiling" in text
    assert "minco profile" in text


def main() -> int:
    old_hooks = patch_model_hooks()
    try:
        with tempfile.TemporaryDirectory(prefix="minco_auto_exact_test_") as tmp:
            base = Path(tmp)
            test_table_mode_trigger_falls_back(base / "table")
            test_default_strategy_matches_explicit_autoexact(base / "default")
            test_default_launcher_candidate_preset_runs_table_mode(base / "candidate_table")
            (base / "raw_exact").mkdir()
            test_raw_mode_trigger_runs_exact(base / "raw_exact")
            (base / "raw_no_exact").mkdir()
            test_raw_mode_trigger_off_skips_exact(base / "raw_no_exact")
            test_model_cache_roundtrip_and_metadata_validation(base / "model_cache")
            test_load_or_fit_models_prefers_discovered_train_feature_cache(base / "model_cache_discover")
            test_high_extra_low_unique_af_uses_probability()
            test_high_extra_high_unique_af_keeps_rescue()
            test_low_extra_low_unique_af_adds_guarded_split_rescue()
            test_low_extra_high_unique_af_keeps_base_only()
            test_reported_ani_prefers_split_zip_aaf_unless_raw_unique_fallback()
            test_initial_pass_thread_budget_is_split_without_exceeding_total()
            test_same_stream_initial_pass_uses_unique_sidecar(base / "same_stream_pass")
            test_auto_initial_pass_uses_same_stream_at_one_thread(base / "auto_p1_pass")
            test_legacy_sequential_initial_passes_use_requested_threads(base / "sequential_passes")
            test_parallel_initial_passes_delegate_to_parallel_helper(base / "parallel_passes")
            test_exact_split_guard_rejects_sample4_like_case()
            test_exact_split_guard_accepts_sparse_or_strong_unique_cases()
            test_default_launcher_injects_universal_autoexact_and_env_defaults()
            test_default_launcher_constants_match_current_strategy_manifest()
            test_default_launcher_current_preset_reproduces_previous_default()
            test_default_launcher_discovers_packaged_ref_sidecars(base / "sidecars")
            test_default_launcher_discovers_packaged_model_cache_without_training(base / "cache_sidecar")
            test_default_launcher_applies_packaged_scope_default(base / "packaged_scope")
            test_default_launcher_explicit_scope_overrides_packaged_default(base / "packaged_scope_override")
            test_default_launcher_candidate_preset_injects_candidate_switches_and_sidecar(base / "candidate_sidecars")
            test_default_launcher_candidate_preset_can_use_env_surface_taxmap(base / "candidate_env")
            test_default_launcher_candidate_preset_respects_surface_off_with_sidecar(base / "candidate_surface_off")
            test_default_launcher_preflight_passes_packaged_ref(base / "preflight_pass")
            test_default_launcher_preflight_fails_missing_candidate_surface(base / "preflight_missing_surface")
            test_default_launcher_preflight_current_preset_does_not_require_surface(base / "preflight_current")
            test_default_launcher_preserves_explicit_advanced_options()
            test_default_launcher_requires_taxmap_or_env_unless_help()
            test_default_launcher_help_banner_names_strategy_and_boundary()
            test_panel_abundance_uses_max_split_unique_zip_depth_and_tail_exception()
    finally:
        restore_model_hooks(old_hooks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
