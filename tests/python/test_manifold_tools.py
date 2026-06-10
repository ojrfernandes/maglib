"""
Tests for maglib.manifold_tools — simplify and to_linestrings.

All tests require shapely; the entire module is skipped when it is absent.
Tests verify geometric properties of RDP (point counts, subset guarantee,
endpoint preservation) and LineString fidelity, not visual output.
"""

import numpy as np
import pytest

pytest.importorskip("shapely", reason="shapely not installed")

from maglib import simplify, to_linestrings
from maglib.manifold_tools import _extract_segments


# ── Helpers ───────────────────────────────────────────────────────────────────

class _FakeManifold:
    def __init__(self, segs):
        self.output_data = segs


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def straight_line():
    """5 collinear points along Z=0. Any positive tolerance removes interior points."""
    return np.array([[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.5, 0.0], [2.0, 0.0]])


@pytest.fixture
def bent_line():
    """3 points with a 0.1 m perpendicular deviation at the midpoint."""
    return np.array([[0.0, 0.0], [1.0, 0.1], [2.0, 0.0]])


@pytest.fixture
def two_segs(straight_line, bent_line):
    return [straight_line, bent_line]


# ── _extract_segments ─────────────────────────────────────────────────────────

def test_extract_from_list_returns_same_object(two_segs):
    assert _extract_segments(two_segs) is two_segs


def test_extract_from_manifold_object(two_segs):
    assert _extract_segments(_FakeManifold(two_segs)) is two_segs


def test_extract_invalid_type_raises():
    with pytest.raises(TypeError):
        _extract_segments("not_valid")


# ── simplify — geometric properties ──────────────────────────────────────────

def test_simplify_collinear_reduces_to_endpoints(straight_line):
    # All interior points lie on the line (distance = 0 < any positive tolerance)
    result = simplify([straight_line], tolerance=0.001)
    assert len(result[0]) == 2
    np.testing.assert_allclose(result[0][0], straight_line[0])
    np.testing.assert_allclose(result[0][-1], straight_line[-1])


def test_simplify_preserves_endpoints_under_large_tolerance(bent_line):
    result = simplify([bent_line], tolerance=10.0)
    np.testing.assert_allclose(result[0][0], bent_line[0])
    np.testing.assert_allclose(result[0][-1], bent_line[-1])


def test_simplify_output_is_subset_of_input(straight_line):
    """RDP only removes existing points — no new coordinates are introduced."""
    result = simplify([straight_line], tolerance=0.001)
    input_pts = set(map(tuple, straight_line))
    for pt in result[0]:
        assert tuple(pt) in input_pts


def test_simplify_keeps_bend_above_tolerance(bent_line):
    # Perpendicular deviation is 0.1 m; tolerance 0.05 m → midpoint must be kept
    result = simplify([bent_line], tolerance=0.05)
    assert len(result[0]) == 3


def test_simplify_removes_bend_below_tolerance(bent_line):
    # Perpendicular deviation is 0.1 m; tolerance 0.2 m → midpoint is removed
    result = simplify([bent_line], tolerance=0.2)
    assert len(result[0]) == 2


def test_simplify_never_increases_point_count(two_segs):
    for tol in (1e-4, 1e-2, 0.5):
        result = simplify(two_segs, tol)
        for orig, simp in zip(two_segs, result):
            assert len(simp) <= len(orig)


def test_simplify_returns_ndarray_with_correct_shape(straight_line):
    result = simplify([straight_line], tolerance=0.001)
    assert isinstance(result[0], np.ndarray)
    assert result[0].ndim == 2
    assert result[0].shape[1] == 2


def test_simplify_single_point_segment():
    seg = np.array([[0.5, 0.3]])
    result = simplify([seg], tolerance=0.01)
    assert len(result[0]) == 1
    np.testing.assert_allclose(result[0], seg)


def test_simplify_multiple_segments_preserves_count(two_segs):
    result = simplify(two_segs, tolerance=0.01)
    assert len(result) == len(two_segs)


def test_simplify_from_manifold_object(two_segs):
    result = simplify(_FakeManifold(two_segs), tolerance=0.01)
    assert len(result) == 2
    assert all(isinstance(s, np.ndarray) for s in result)


# ── to_linestrings ────────────────────────────────────────────────────────────

def test_to_linestrings_returns_list(two_segs):
    assert isinstance(to_linestrings(two_segs), list)


def test_to_linestrings_count_matches_input(two_segs):
    assert len(to_linestrings(two_segs)) == len(two_segs)


def test_to_linestrings_type():
    from shapely.geometry import LineString
    result = to_linestrings([np.array([[0.0, 0.0], [1.0, 1.0]])])
    assert isinstance(result[0], LineString)


def test_to_linestrings_coordinates_match_input(straight_line):
    result = to_linestrings([straight_line])
    coords = np.array(result[0].coords)
    np.testing.assert_allclose(coords, straight_line)


def test_to_linestrings_no_interpolation(bent_line):
    """LineString must store exactly the input points, no extras added."""
    result = to_linestrings([bent_line])
    coords = np.array(result[0].coords)
    assert len(coords) == len(bent_line)
    np.testing.assert_allclose(coords, bent_line)


def test_to_linestrings_from_manifold_object(two_segs):
    result = to_linestrings(_FakeManifold(two_segs))
    assert len(result) == len(two_segs)
