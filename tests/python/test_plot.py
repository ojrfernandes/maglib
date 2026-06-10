"""
Tests for maglib.plot — plot_footprint and plot_manifold.

Rendering is not tested (matplotlib's responsibility). Tests cover:
  - source dispatch (ndarray, object with output_data, file path)
  - input validation (TypeError / ValueError for bad arguments)
  - _apply_cut helper
  - smoke tests: function runs without error for valid inputs
"""

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pytest

from maglib import plot_footprint, plot_manifold
from maglib.plot import _apply_cut

# Suppress the "cannot show figure" warning emitted by the Agg backend.
pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_fp_data(n_phi=4, n_y=3):
    """Synthetic (n_phi*n_y, 6) footprint array on a horizontal plate."""
    phi   = np.tile(np.linspace(0.0, 270.0, n_phi), n_y)
    R     = np.repeat(np.linspace(1.4, 1.6, n_y), n_phi)
    Z     = np.zeros(n_phi * n_y)
    CL    = np.linspace(1.0, 20.0, n_phi * n_y)
    psi   = np.linspace(0.1, 0.8, n_phi * n_y)
    turns = np.linspace(2.0, 15.0, n_phi * n_y)
    return np.column_stack([R, Z, phi, CL, psi, turns])


def _make_mf_segs():
    seg0 = np.array([[0.50, -0.20], [0.51, -0.21], [0.52, -0.22], [0.53, -0.23]])
    seg1 = np.array([[0.60, -0.10], [0.61, -0.11], [0.62, -0.12]])
    return [seg0, seg1]


class _FakeFootprint:
    def __init__(self, data):
        self.output_data = data


class _FakeManifold:
    def __init__(self, segs):
        self.output_data = segs


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


@pytest.fixture
def fp_data():
    return _make_fp_data()


@pytest.fixture
def mf_segs():
    return _make_mf_segs()


@pytest.fixture
def tmp_fp_file(fp_data):
    with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as f:
        path = f.name
    np.savetxt(path, fp_data)
    yield path
    os.unlink(path)


@pytest.fixture
def tmp_mf_files(mf_segs):
    paths = []
    for seg in mf_segs:
        with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as f:
            path = f.name
        np.savetxt(path, seg)
        paths.append(path)
    yield paths
    for p in paths:
        os.unlink(p)


# ── _apply_cut ────────────────────────────────────────────────────────────────

def test_apply_cut_none_returns_full_array(fp_data):
    assert _apply_cut(fp_data, None) is fp_data


def test_apply_cut_zero_returns_full_array(fp_data):
    assert _apply_cut(fp_data, 0) is fp_data


def test_apply_cut_positive_slices_from_start():
    arr = np.arange(10).reshape(5, 2)
    result = _apply_cut(arr, 3)
    assert result.shape == (3, 2)
    np.testing.assert_array_equal(result, arr[:3])


def test_apply_cut_negative_slices_from_end():
    arr = np.arange(10).reshape(5, 2)
    result = _apply_cut(arr, -2)
    assert result.shape == (3, 2)
    np.testing.assert_array_equal(result, arr[:-2])


# ── plot_footprint — source dispatch ─────────────────────────────────────────

def test_plot_footprint_accepts_ndarray(fp_data):
    plot_footprint(fp_data, which_plot="cl")


def test_plot_footprint_accepts_object_with_output_data(fp_data):
    plot_footprint(_FakeFootprint(fp_data), which_plot="cl")


def test_plot_footprint_accepts_file_path(tmp_fp_file):
    plot_footprint(tmp_fp_file, which_plot="cl")


# ── plot_footprint — input validation ────────────────────────────────────────

def test_plot_footprint_invalid_source_type():
    with pytest.raises(TypeError):
        plot_footprint(42, which_plot="cl")


def test_plot_footprint_wrong_ndim():
    with pytest.raises(ValueError):
        plot_footprint(np.zeros((10, 4)), which_plot="cl")  # only 4 cols, needs 6


def test_plot_footprint_invalid_which_plot(fp_data):
    with pytest.raises(ValueError):
        plot_footprint(fp_data, which_plot="invalid")


def test_plot_footprint_invalid_xaxis(fp_data):
    with pytest.raises(ValueError):
        plot_footprint(fp_data, which_plot="cl", xaxis="invalid")


def test_plot_footprint_vmin_ge_vmax(fp_data):
    with pytest.raises(ValueError):
        plot_footprint(fp_data, which_plot="au", v_min=0.8, v_max=0.2)


def test_plot_footprint_au_raises_when_all_values_identical():
    # CL=2, psi=0.5 for all → 1/(CL*psi)=1 everywhere → lo==hi → ValueError
    n = 4 * 3
    phi   = np.tile(np.linspace(0.0, 270.0, 4), 3)
    R     = np.repeat(np.linspace(1.4, 1.6, 3), 4)
    Z     = np.zeros(n)
    CL    = np.full(n, 2.0)
    psi   = np.full(n, 0.5)
    turns = np.full(n, 5.0)
    data  = np.column_stack([R, Z, phi, CL, psi, turns])
    with pytest.raises(ValueError):
        plot_footprint(data, which_plot="au")


# ── plot_footprint — smoke tests ──────────────────────────────────────────────

@pytest.mark.parametrize("which_plot", ["cl", "psi", "turns", "au"])
def test_plot_footprint_each_subplot(fp_data, which_plot):
    plot_footprint(fp_data, which_plot=which_plot)


def test_plot_footprint_psi_cap(fp_data):
    plot_footprint(fp_data, which_plot="cl", psi_cap=True)


def test_plot_footprint_turn_cap(fp_data):
    plot_footprint(fp_data, which_plot="turns", turn_cap=(2, 10))


def test_plot_footprint_vertical_plate():
    # Constant R → vertical plate detection path
    n_phi, n_y = 4, 3
    phi   = np.tile(np.linspace(0.0, 270.0, n_phi), n_y)
    R     = np.full(n_phi * n_y, 1.5)
    Z     = np.repeat(np.linspace(-0.3, -0.1, n_y), n_phi)
    CL    = np.linspace(1.0, 20.0, n_phi * n_y)
    psi   = np.linspace(0.1, 0.8, n_phi * n_y)
    turns = np.linspace(2.0, 15.0, n_phi * n_y)
    plot_footprint(np.column_stack([R, Z, phi, CL, psi, turns]), which_plot="cl")


# ── plot_manifold — source dispatch ──────────────────────────────────────────

def test_plot_manifold_accepts_list_of_arrays(mf_segs):
    plot_manifold(mf_segs)


def test_plot_manifold_accepts_manifold_object(mf_segs):
    plot_manifold(_FakeManifold(mf_segs))


def test_plot_manifold_accepts_single_file(tmp_mf_files):
    plot_manifold(tmp_mf_files[0])


def test_plot_manifold_accepts_list_of_files(tmp_mf_files):
    plot_manifold(tmp_mf_files)


# ── plot_manifold — input validation ─────────────────────────────────────────

def test_plot_manifold_invalid_source_type():
    with pytest.raises(TypeError):
        plot_manifold(42)


def test_plot_manifold_empty_list():
    with pytest.raises(ValueError):
        plot_manifold([])


def test_plot_manifold_invalid_linestyle(mf_segs):
    with pytest.raises(ValueError):
        plot_manifold(mf_segs, linestyle="dotted")


def test_plot_manifold_cuts_length_mismatch(mf_segs):
    with pytest.raises(ValueError):
        plot_manifold(mf_segs, cuts=[10])  # 2 segments, 1 cut


# ── plot_manifold — smoke tests ───────────────────────────────────────────────

@pytest.mark.parametrize("linestyle", ["line", "scatter", "both"])
def test_plot_manifold_each_linestyle(mf_segs, linestyle):
    plot_manifold(mf_segs, linestyle=linestyle)


@pytest.mark.parametrize("linestyle", ["line", "scatter", "both"])
def test_plot_manifold_object_each_linestyle(mf_segs, linestyle):
    plot_manifold(_FakeManifold(mf_segs), linestyle=linestyle)


def test_plot_manifold_with_cuts(mf_segs):
    plot_manifold(mf_segs, cuts=[2, 0])


def test_plot_manifold_with_labels(mf_segs):
    plot_manifold(mf_segs, labels=["Stable", "Unstable"])
