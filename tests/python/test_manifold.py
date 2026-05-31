"""
Python tests for Manifold.

Reference values mirror those in tests/test_mfgen.cpp.
"""

import math
import os
import tempfile

import numpy as np
import pytest
import maglib
from maglib import M3DC1Source, Maglit, Manifold

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
C1_H5    = os.path.join(DATA_DIR, "C1.h5")

TOL = 1e-6   # matches C++ reference-value tolerance


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def tracer_and_source():
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    return src, t


@pytest.fixture(scope="module")
def manifold_at_xpoint(tracer_and_source):
    """Unstable manifold with confirmed X-point."""
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.find_x_point(0.497999, -0.218603)
    return mf


@pytest.fixture(scope="module")
def primary_seg(manifold_at_xpoint):
    return manifold_at_xpoint.primary_segment(10)


# ── Construction and configuration ───────────────────────────────────────────

def test_construction(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    assert mf is not None


def test_configure(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(1e-6, 1e-8, 1e-14, 50, 1e-14, 100)


# ── find_x_point ─────────────────────────────────────────────────────────────

def test_find_x_point_returns_true(manifold_at_xpoint):
    _, t = pytest.importorskip("maglib"), None
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    mf = Manifold(t, phi=0.0, stability=0)
    result = mf.find_x_point(0.497999, -0.218603)
    assert result is True


def test_find_x_point_reference_values(manifold_at_xpoint):
    xp = manifold_at_xpoint.x_point
    assert abs(xp[0] - 0.4979691771716279) < TOL   # R
    assert abs(xp[1] - (-0.2185980054447758)) < TOL  # Z


# ── x_point property ──────────────────────────────────────────────────────────

def test_x_point_is_numpy_array(manifold_at_xpoint):
    xp = manifold_at_xpoint.x_point
    assert isinstance(xp, np.ndarray)
    assert xp.shape == (2,)
    assert xp.dtype == np.float64


def test_x_point_setter(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    mf.x_point = np.array([0.498, -0.219])
    xp = mf.x_point
    assert abs(xp[0] - 0.498) < 1e-15
    assert abs(xp[1] - (-0.219)) < 1e-15


def test_x_point_setter_wrong_size_raises(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    with pytest.raises((ValueError, RuntimeError)):
        mf.x_point = np.array([0.498])


# ── primary_segment ───────────────────────────────────────────────────────────

def test_primary_segment_shape(primary_seg):
    assert isinstance(primary_seg, np.ndarray)
    assert primary_seg.shape == (11, 2)   # n_intervals=10 → 11 points
    assert primary_seg.dtype == np.float64


def test_primary_segment_first_point(primary_seg):
    assert abs(primary_seg[0, 0] - 0.4979701320801332) < TOL
    assert abs(primary_seg[0, 1] - (-0.2185977085445460)) < TOL


def test_primary_segment_last_point(primary_seg):
    assert abs(primary_seg[-1, 0] - 0.4979855389268979) < TOL
    assert abs(primary_seg[-1, 1] - (-0.2185929033960707)) < TOL


# ── new_segment (interpolant method) ─────────────────────────────────────────

def test_new_segment_returns_tuple_of_arrays(manifold_at_xpoint, primary_seg):
    result = manifold_at_xpoint.new_segment(primary_seg, l_lim=0.005, theta_lim=20.0)
    assert isinstance(result, tuple) and len(result) == 2
    prev_out, new_seg = result
    assert isinstance(prev_out, np.ndarray) and prev_out.ndim == 2 and prev_out.shape[1] == 2
    assert isinstance(new_seg,  np.ndarray) and new_seg.ndim  == 2 and new_seg.shape[1]  == 2


def test_new_segment_reference_values(manifold_at_xpoint, primary_seg):
    """Three iterations of the interpolant method; spot-check new_seg_3."""
    seg = primary_seg.copy()
    _, seg1 = manifold_at_xpoint.new_segment(seg,  l_lim=0.005, theta_lim=20.0)
    _, seg2 = manifold_at_xpoint.new_segment(seg1, l_lim=0.005, theta_lim=20.0)
    _, seg3 = manifold_at_xpoint.new_segment(seg2, l_lim=0.005, theta_lim=20.0)

    assert abs(seg3[0,  0] - 0.5012464706502761) < TOL
    assert abs(seg3[0,  1] - (-0.2174148574827388)) < TOL
    assert abs(seg3[4,  0] - 0.5030591235329525) < TOL
    assert abs(seg3[4,  1] - (-0.2168648161175224)) < TOL
    assert abs(seg3[9,  0] - 0.5133430437264513) < TOL
    assert abs(seg3[9,  1] - (-0.2130464750286543)) < TOL
    assert abs(seg3[14, 0] - 0.5298324820702656) < TOL
    assert abs(seg3[14, 1] - (-0.208634750788582)) < TOL
    assert abs(seg3[19, 0] - 0.5429352654656256) < TOL
    assert abs(seg3[19, 1] - (-0.199387939879125)) < TOL
    assert abs(seg3[24, 0] - 0.5594125465315225) < TOL
    assert abs(seg3[24, 1] - (-0.192432899160264)) < TOL


# ── new_segment (exact-map method) ───────────────────────────────────────────

def test_new_segment_from_primary_reference_values(manifold_at_xpoint, primary_seg):
    """Exact-map method with n_seg=3; spot-check the output."""
    _, seg3 = manifold_at_xpoint.new_segment(primary_seg.copy(),
                                              n_seg=3, l_lim=0.005, theta_lim=20.0)

    assert abs(seg3[0,  0] - 0.5037194321479069) < TOL
    assert abs(seg3[0,  1] - (-0.2167180039510279)) < TOL
    assert abs(seg3[4,  0] - 0.5177907149090401) < TOL
    assert abs(seg3[4,  1] - (-0.2125646567853695)) < TOL
    assert abs(seg3[9,  0] - 0.5328067328393506) < TOL
    assert abs(seg3[9,  1] - (-0.2069275324563348)) < TOL
    assert abs(seg3[14, 0] - 0.5448185439171523) < TOL
    assert abs(seg3[14, 1] - (-0.1979494329010073)) < TOL


# ── input validation ──────────────────────────────────────────────────────────

def test_new_segment_wrong_shape_raises(manifold_at_xpoint):
    bad = np.zeros((5, 3))
    with pytest.raises((ValueError, RuntimeError)):
        manifold_at_xpoint.new_segment(bad, l_lim=0.005, theta_lim=20.0)


# ── output_data ───────────────────────────────────────────────────────────────

def test_output_data_empty_before_computation(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    assert mf.output_data == []


def test_output_data_accumulates(manifold_at_xpoint, primary_seg):
    # manifold_at_xpoint fixture already called find_x_point; primary_seg called
    # primary_segment once — outputData has at least 1 entry from the fixture.
    # Build a fresh manifold so we control the count exactly.
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.find_x_point(0.497999, -0.218603)

    assert mf.output_data == []

    seg0 = mf.primary_segment(10)
    assert len(mf.output_data) == 1
    assert mf.output_data[0].shape == (11, 2)

    mf.new_segment(seg0, l_lim=0.005, theta_lim=20.0)
    assert len(mf.output_data) == 2

    mf.new_segment(seg0, n_seg=2, l_lim=0.005, theta_lim=20.0)
    assert len(mf.output_data) == 3


def test_output_data_types(manifold_at_xpoint, primary_seg):
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.find_x_point(0.497999, -0.218603)
    seg0 = mf.primary_segment(10)
    mf.new_segment(seg0, l_lim=0.005, theta_lim=20.0)

    for arr in mf.output_data:
        assert isinstance(arr, np.ndarray)
        assert arr.ndim == 2
        assert arr.shape[1] == 2
        assert arr.dtype == np.float64


# ── save ──────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def manifold_with_output():
    """Manifold with 2 segments (primary + 1 new) for save tests."""
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.find_x_point(0.497999, -0.218603)
    seg0 = mf.primary_segment(10)
    mf.new_segment(seg0, l_lim=0.005, theta_lim=20.0)
    return mf


def test_save_dat(manifold_with_output):
    with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as tmp:
        path = tmp.name
    manifold_with_output.save(path)
    with open(path) as f:
        lines = f.readlines()
    assert lines[0].startswith("#")           # header
    assert lines[1].startswith("0 ")          # first data row: segment 0
    assert any(l.startswith("1 ") for l in lines)  # segment 1 present
    os.unlink(path)


def test_save_txt(manifold_with_output):
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tmp:
        path = tmp.name
    manifold_with_output.save(path)
    with open(path) as f:
        header = f.readline()
    assert header.startswith("#")
    os.unlink(path)


def test_save_csv(manifold_with_output):
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        path = tmp.name
    manifold_with_output.save(path)
    with open(path) as f:
        lines = f.readlines()
    assert lines[0].strip() == "seg,R,Z"
    assert lines[1].startswith("0,")
    assert any(l.startswith("1,") for l in lines)
    os.unlink(path)


def test_save_npz(manifold_with_output):
    with tempfile.NamedTemporaryFile(suffix=".npz", delete=False) as tmp:
        path = tmp.name
    manifold_with_output.save(path)
    data = np.load(path)
    assert "seg_0" in data
    assert "seg_1" in data
    assert data["seg_0"].shape == (11, 2)
    assert data["seg_0"].dtype == np.float64
    os.unlink(path)


def test_save_unsupported_raises(manifold_with_output):
    with pytest.raises((ValueError, RuntimeError)):
        manifold_with_output.save("/tmp/manifold_test.xyz")


# ── run() ─────────────────────────────────────────────────────────────────────

def _fresh_manifold():
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(1e-2, 1e-6, 1e-2)
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.find_x_point(0.497999, -0.218603)
    return mf


def test_run_output_data_size_interpolant():
    mf = _fresh_manifold()
    mf.run(n_intervals=10, n_segments=4, method=1, l_lim=0.005, theta_lim=20.0)
    assert len(mf.output_data) == 4
    assert mf.output_data[0].shape == (11, 2)  # primary: n_intervals+1 points


def test_run_output_data_size_exact_map():
    mf = _fresh_manifold()
    mf.run(n_intervals=10, n_segments=3, method=0, l_lim=0.005, theta_lim=20.0)
    assert len(mf.output_data) == 3


def test_run_interpolant_matches_manual_loop():
    """run(method=1) should give the same segments as calling new_segment() manually."""
    mf_run = _fresh_manifold()
    mf_run.run(n_intervals=10, n_segments=3, method=1, l_lim=0.005, theta_lim=20.0)

    mf_manual = _fresh_manifold()
    seg = mf_manual.primary_segment(10)
    _, seg = mf_manual.new_segment(seg, l_lim=0.005, theta_lim=20.0)
    mf_manual.new_segment(seg, l_lim=0.005, theta_lim=20.0)

    for i in range(3):
        np.testing.assert_array_equal(mf_run.output_data[i], mf_manual.output_data[i])


def test_run_exact_map_matches_manual():
    """run(method=0) should give the same segments as calling new_segment(n_seg=i) manually."""
    mf_run = _fresh_manifold()
    mf_run.run(n_intervals=10, n_segments=3, method=0, l_lim=0.005, theta_lim=20.0)

    mf_manual = _fresh_manifold()
    primary = mf_manual.primary_segment(10)
    mf_manual.new_segment(primary, n_seg=1, l_lim=0.005, theta_lim=20.0)
    mf_manual.new_segment(primary, n_seg=2, l_lim=0.005, theta_lim=20.0)

    for i in range(3):
        np.testing.assert_array_equal(mf_run.output_data[i], mf_manual.output_data[i])
