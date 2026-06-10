"""
Tests for maglib.manifold_tools — simplify, to_linestrings, and lobe_map.

All tests require shapely; the entire module is skipped when it is absent.
Tests verify geometric properties of RDP (point counts, subset guarantee,
endpoint preservation), LineString fidelity, and lobe_map correctness
(area, perimeter, h_param, midpoint, angle, output shape).
"""

import numpy as np
import pytest

pytest.importorskip("shapely", reason="shapely not installed")

from maglib import simplify, to_linestrings, lobe_map
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


# ── lobe_map fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def rectangular_lobe():
    """
    Equilibrium: horizontal line with 5 points from R=0 to R=4 at Z=0.
    Perturbed: three-sided path that rises above Z=0 between R=1 and R=3.

    Intersection points at (1, 0) and (3, 0) — one lobe.
    Lobe polygon vertices: (1,0), (2,0), (3,0), (3,1), (1,1).
    Area  = 2 × 1 = 2.0 m²
    Perimeter = 1 + 1 + 1 + 2 + 1 = 6.0 m  (closing edge 1→0 not needed — shapely closes)

    Wait — shapely Polygon.length includes the closing edge.
    Vertices (1,0)→(2,0)→(3,0)→(3,1)→(1,1)→close→(1,0)
    Edges: 1 + 1 + 1 + 2 + 1 = 6 (closing edge has length sqrt((1-1)²+(0-1)²)=1)
    Total perimeter = 6.0 m.
    """
    eq = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]])
    ptb = np.array([[1.0, -1.0], [1.0, 0.0], [1.0, 1.0], [3.0, 1.0], [3.0, 0.0], [3.0, -1.0]])
    mag_axis = np.array([2.0, -2.0])
    x_point = np.array([1.0, -1.0])
    return eq, ptb, mag_axis, x_point


@pytest.fixture
def two_lobe_setup():
    """
    Two lobes: one above Z=0 (R=1–3) and one below (R=3–5).

    Equilibrium: R=0..6, Z=0 with 7 points.
    Perturbed: goes above between R=1–3, below between R=3–5.
    Intersections at R=1, 3, 5 → two lobes.
    """
    eq = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0],
                   [4.0, 0.0], [5.0, 0.0], [6.0, 0.0]])
    ptb = np.array([[1.0, -0.5], [1.0, 0.0], [1.0, 1.0], [3.0, 1.0], [3.0, 0.0],
                    [3.0, -1.0], [5.0, -1.0], [5.0, 0.0], [5.0, 0.5]])
    mag_axis = np.array([3.0, -3.0])
    x_point = np.array([1.0, -0.5])
    return eq, ptb, mag_axis, x_point


# ── lobe_map — output shape and types ─────────────────────────────────────────

def test_lobe_map_returns_ndarray(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    assert isinstance(result, np.ndarray)


def test_lobe_map_shape_one_lobe(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    assert result.shape == (1, 6)


def test_lobe_map_shape_two_lobes(two_lobe_setup):
    eq, ptb, mag_axis, x_point = two_lobe_setup
    result = lobe_map(eq, ptb, mag_axis, x_point)
    assert result.shape == (2, 6)


def test_lobe_map_no_intersection_returns_empty():
    eq = np.array([[0.0, 1.0], [4.0, 1.0]])   # horizontal line at Z=1
    ptb = np.array([[0.0, -1.0], [4.0, -1.0]])  # horizontal line at Z=-1
    mag_axis = np.array([2.0, 0.0])
    result = lobe_map(eq, ptb, mag_axis)
    assert result.shape == (0, 6)


def test_lobe_map_single_intersection_returns_empty():
    eq = np.array([[0.0, 0.0], [4.0, 0.0]])
    ptb = np.array([[0.0, -1.0], [2.0, 0.0], [4.0, -1.0]])  # touches Z=0 once
    mag_axis = np.array([2.0, -2.0])
    result = lobe_map(eq, ptb, mag_axis)
    assert result.shape == (0, 6)


# ── lobe_map — geometric correctness ─────────────────────────────────────────

def test_lobe_map_area(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    # 2×1 rectangle
    np.testing.assert_allclose(result[0, 4], 2.0, rtol=1e-10)


def test_lobe_map_perimeter(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    # edges: 1+1+1+2+1 = 6.0 (see fixture docstring)
    np.testing.assert_allclose(result[0, 3], 6.0, rtol=1e-10)


def test_lobe_map_h_param(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    # h = 2*area / base = 2*2 / 2 = 2.0
    np.testing.assert_allclose(result[0, 5], 2.0, rtol=1e-10)


def test_lobe_map_midpoint(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    # eq sub-curve: (1,0),(2,0),(3,0); midpoint = mean of (1,0),(2,0) = (1.5, 0.0)
    np.testing.assert_allclose(result[0, 0], 1.5, rtol=1e-10)
    np.testing.assert_allclose(result[0, 1], 0.0, atol=1e-14)


def test_lobe_map_angle_is_finite_and_in_range(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    angle = result[0, 2]
    assert np.isfinite(angle)
    assert 0.0 <= angle <= np.pi


def test_lobe_map_all_positive(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result = lobe_map(eq, ptb, mag_axis, x_point)
    # area, perimeter, h_param must be positive for a non-degenerate lobe
    assert result[0, 3] > 0.0   # perimeter
    assert result[0, 4] > 0.0   # area
    assert result[0, 5] > 0.0   # h_param


# ── lobe_map — x_point default ────────────────────────────────────────────────

def test_lobe_map_x_point_default_uses_perturbed_first(rectangular_lobe):
    eq, ptb, mag_axis, x_point = rectangular_lobe
    result_explicit = lobe_map(eq, ptb, mag_axis, x_point)
    result_default = lobe_map(eq, ptb, mag_axis)  # x_point defaults to ptb[0]
    np.testing.assert_allclose(result_explicit, result_default)
