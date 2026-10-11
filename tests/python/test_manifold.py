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
THETA_LIM = math.radians(20.0)   # turning-angle refinement threshold, rad

# Fixed inputs for the new_segment reference-value tests; identical to the
# hard-coded inputs in tests/test_mfgen.cpp. Using them (instead of the output
# of primary_segment()) keeps the references independent of how primarySegment
# parametrizes its interior points: a 1e-8 change there is amplified to ~1e-5
# by three map iterations.
HARD_XPOINT = np.array([0.4979691771716279, -0.2185980054447758])   # shape=(2,), [R, Z], units=m
HARD_PRIMARY = np.array([                                           # shape=(11,2), [R, Z], units=m
    [0.4979701320801332, -0.218597708544546],
    [0.4979716727648097, -0.2185972280296985],
    [0.4979732134494862, -0.2185967475148509],
    [0.4979747541341626, -0.2185962670000034],
    [0.4979762948188391, -0.2185957864851559],
    [0.4979778355035156, -0.2185953059703084],
    [0.497979376188192,  -0.2185948254554608],
    [0.4979809168728685, -0.2185943449406133],
    [0.497982457557545,  -0.2185938644257657],
    [0.4979839982422214, -0.2185933839109182],
    [0.4979855389268979, -0.2185929033960707],
])


def _fixed_input_manifold(t):
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.x_point = HARD_XPOINT
    return mf


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


def test_phi_is_in_radians(tracer_and_source):
    """Manifold(phi=) is in radians: phi=2*pi is the same section as phi=0."""
    _, t = tracer_and_source
    xp = []
    for phi in (0.0, 2 * math.pi):
        mf = Manifold(t, phi=phi, stability=0)
        mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                     max_iter=50, precision_limit=1e-14, max_insertions=100)
        assert mf.find_x_point(0.497999, -0.218603)
        xp.append(np.array(mf.x_point))   # shape=(2,), [R, Z], units=m
    assert np.allclose(xp[0], xp[1], atol=1e-8)


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
    result = manifold_at_xpoint.new_segment(primary_seg, l_lim=0.005, theta_lim=THETA_LIM)
    assert isinstance(result, tuple) and len(result) == 2
    prev_out, new_seg = result
    assert isinstance(prev_out, np.ndarray) and prev_out.ndim == 2 and prev_out.shape[1] == 2
    assert isinstance(new_seg,  np.ndarray) and new_seg.ndim  == 2 and new_seg.shape[1]  == 2


def test_new_segment_reference_values(tracer_and_source):
    """Three iterations of the interpolant method; spot-check new_seg_3."""
    _, t = tracer_and_source
    mf = _fixed_input_manifold(t)
    _, seg1 = mf.new_segment(HARD_PRIMARY.copy(), l_lim=0.005, theta_lim=THETA_LIM)
    _, seg2 = mf.new_segment(seg1, l_lim=0.005, theta_lim=THETA_LIM)
    _, seg3 = mf.new_segment(seg2, l_lim=0.005, theta_lim=THETA_LIM)

    # Reference values mirror Manifold_NewSegment_FromPrevious in tests/test_mfgen.cpp
    # (signed-turning-angle interpolant; indices 9-34 differ ~1e-5 m from the old unsigned one)
    assert abs(seg3[0 , 0] - 0.5037194236606823) < TOL
    assert abs(seg3[0 , 1] - (-0.2167180059957065)) < TOL
    assert abs(seg3[4 , 0] - 0.5177907016811655) < TOL
    assert abs(seg3[4 , 1] - (-0.21256465989197604)) < TOL
    assert abs(seg3[9 , 0] - 0.5329413313161456) < TOL
    assert abs(seg3[9 , 1] - (-0.20684384922194668)) < TOL
    assert abs(seg3[14, 0] - 0.544888620002773) < TOL
    assert abs(seg3[14, 1] - (-0.19789861738018513)) < TOL
    assert abs(seg3[19, 0] - 0.5623249055389825) < TOL
    assert abs(seg3[19, 1] - (-0.19188669144018536)) < TOL
    assert abs(seg3[24, 0] - 0.5773761419825298) < TOL
    assert abs(seg3[24, 1] - (-0.18797022952026268)) < TOL
    assert abs(seg3[29, 0] - 0.5936666484608587) < TOL
    assert abs(seg3[29, 1] - (-0.1813653812146642)) < TOL
    assert abs(seg3[34, 0] - 0.6134021821972827) < TOL
    assert abs(seg3[34, 1] - (-0.17094516309055374)) < TOL


# ── new_segment (exact-map method) ───────────────────────────────────────────

def test_new_segment_from_primary_reference_values(tracer_and_source):
    """Exact-map method with n_seg=3; spot-check the output."""
    _, t = tracer_and_source
    mf = _fixed_input_manifold(t)
    _, seg3 = mf.new_segment(HARD_PRIMARY.copy(),
                             n_seg=3, l_lim=0.005, theta_lim=THETA_LIM)

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
        manifold_at_xpoint.new_segment(bad, l_lim=0.005, theta_lim=THETA_LIM)


def test_new_segment_failed_map_raises(manifold_at_xpoint):
    # shape=(3,2), units=m; far outside the source domain, so tracing fails
    bad = np.array([[50.0, 50.0], [50.1, 50.0], [50.2, 50.0]])
    with pytest.raises(RuntimeError, match="map evaluation failed"):
        manifold_at_xpoint.new_segment(bad, l_lim=0.005, theta_lim=THETA_LIM)


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

    mf.new_segment(seg0, l_lim=0.005, theta_lim=THETA_LIM)
    assert len(mf.output_data) == 2

    mf.new_segment(seg0, n_seg=2, l_lim=0.005, theta_lim=THETA_LIM)
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
    mf.new_segment(seg0, l_lim=0.005, theta_lim=THETA_LIM)

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
    mf.new_segment(seg0, l_lim=0.005, theta_lim=THETA_LIM)
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
    mf.run(n_intervals=10, n_segments=4, method=1, l_lim=0.005, theta_lim=THETA_LIM)
    assert len(mf.output_data) == 4
    assert mf.output_data[0].shape == (11, 2)  # primary: n_intervals+1 points


def test_run_output_data_size_exact_map():
    mf = _fresh_manifold()
    mf.run(n_intervals=10, n_segments=3, method=0, l_lim=0.005, theta_lim=THETA_LIM)
    assert len(mf.output_data) == 3


def test_run_interpolant_matches_manual_loop():
    """run(method=1) should give the same segments as calling new_segment() manually."""
    mf_run = _fresh_manifold()
    mf_run.run(n_intervals=10, n_segments=3, method=1, l_lim=0.005, theta_lim=THETA_LIM)

    mf_manual = _fresh_manifold()
    seg = mf_manual.primary_segment(10)
    _, seg = mf_manual.new_segment(seg, l_lim=0.005, theta_lim=THETA_LIM)
    mf_manual.new_segment(seg, l_lim=0.005, theta_lim=THETA_LIM)

    for i in range(3):
        np.testing.assert_array_equal(mf_run.output_data[i], mf_manual.output_data[i])


def test_run_exact_map_matches_manual():
    """run(method=0) should give the same segments as calling new_segment(n_seg=i) manually."""
    mf_run = _fresh_manifold()
    mf_run.run(n_intervals=10, n_segments=3, method=0, l_lim=0.005, theta_lim=THETA_LIM)

    mf_manual = _fresh_manifold()
    primary = mf_manual.primary_segment(10)
    mf_manual.new_segment(primary, n_seg=1, l_lim=0.005, theta_lim=THETA_LIM)
    mf_manual.new_segment(primary, n_seg=2, l_lim=0.005, theta_lim=THETA_LIM)

    for i in range(3):
        np.testing.assert_array_equal(mf_run.output_data[i], mf_manual.output_data[i])


# ── read_xpoint_from_hdf5 ─────────────────────────────────────────────────────

from maglib import read_xpoint_from_hdf5

# Reference values: xnull[0] and znull[0] as stored in the test C1.h5 (float32)
_XN1 = 0.49799895
_ZN1 = -0.21860312


def test_read_xpoint_auto_single_null():
    R, Z = read_xpoint_from_hdf5(C1_H5)
    assert abs(R - _XN1) < 1e-5
    assert abs(Z - _ZN1) < 1e-5


def test_read_xpoint_null1():
    R, Z = read_xpoint_from_hdf5(C1_H5, null=1)
    assert abs(R - _XN1) < 1e-5
    assert abs(Z - _ZN1) < 1e-5


def test_read_xpoint_null2_raises_for_single_null():
    with pytest.raises(ValueError, match="single-null"):
        read_xpoint_from_hdf5(C1_H5, null=2)


# ── set_branch ────────────────────────────────────────────────────────────────

def _manifold_at_xpoint(tracer_and_source):
    _, t = tracer_and_source
    mf = Manifold(t, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)
    mf.x_point = np.array([0.4979691771716279, -0.2185980054447758])
    return mf


def test_set_branch_valid(tracer_and_source):
    mf = _manifold_at_xpoint(tracer_and_source)
    mf.set_branch(1, 1)
    mf.set_branch(-1, 1)
    mf.set_branch(1, -1)
    mf.set_branch(-1, -1)


def test_set_branch_invalid_raises(tracer_and_source):
    mf = _manifold_at_xpoint(tracer_and_source)
    with pytest.raises((ValueError, RuntimeError)):
        mf.set_branch(0, 1)
    with pytest.raises((ValueError, RuntimeError)):
        mf.set_branch(1, 2)


def test_set_branch_full_flip_mirrors_pivot(tracer_and_source):
    xp = np.array([0.4979691771716279, -0.2185980054447758])

    mf_def = _manifold_at_xpoint(tracer_and_source)
    seg_def = mf_def.primary_segment(2)

    mf_flip = _manifold_at_xpoint(tracer_and_source)
    mf_flip.set_branch(-1, -1)
    seg_flip = mf_flip.primary_segment(2)

    # First points must be mirror images through xPoint.
    np.testing.assert_allclose(seg_flip[0], 2 * xp - seg_def[0], atol=1e-10)


# ── segment continuity ────────────────────────────────────────────────────────

@pytest.mark.parametrize("method", [0, 1])
def test_run_segments_are_connected(method):
    """Last point of segment k must coincide (within l_lim) with first point of segment k+1."""
    l_lim = 0.005  # m
    mf = _fresh_manifold()
    mf.run(n_intervals=10, n_segments=5, method=method, l_lim=l_lim, theta_lim=THETA_LIM)
    seg = mf.output_data  # each shape=(npts,2), cols=[R,Z], units=m
    for k in range(len(seg) - 1):
        gap = np.linalg.norm(seg[k][-1] - seg[k + 1][0])
        assert gap <= l_lim, f"gap {gap:.4g} m between segments {k} and {k + 1}"
