"""
Python tests for the Collider vessel-wall boundary tester.

Mirrors the C++ GoogleTest suite in tests/test_collider.cpp.
Synthetic shapes are created in a temporary directory; the TCABR wall file
is loaded from tests/data/.
"""

import math
import os
import tempfile

import numpy as np
import pytest
import maglib
from maglib import Collider

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def square_file(tmp_path):
    """Unit square (0,0)–(1,1), counter-clockwise."""
    p = tmp_path / "square.txt"
    p.write_text("0.0 0.0\n1.0 0.0\n1.0 1.0\n0.0 1.0\n")
    return str(p)


@pytest.fixture
def octagon_file(tmp_path):
    """Approximate circle: regular octagon centred at (0.5,0.5), radius 0.4."""
    cx, cz, r = 0.5, 0.5, 0.4
    lines = "\n".join(
        f"{cx + r * math.cos(2 * math.pi * i / 8):.15f} "
        f"{cz + r * math.sin(2 * math.pi * i / 8):.15f}"
        for i in range(8)
    )
    p = tmp_path / "octagon.txt"
    p.write_text(lines + "\n")
    return str(p)


@pytest.fixture
def invalid_file(tmp_path):
    """Only two vertices — too few for a polygon."""
    p = tmp_path / "invalid.txt"
    p.write_text("0.0 0.0\n1.0 0.0\n")
    return str(p)


@pytest.fixture
def reversed_square_file(tmp_path):
    """Unit square with clockwise (reversed) winding order."""
    p = tmp_path / "reversed.txt"
    p.write_text("0.0 1.0\n1.0 1.0\n1.0 0.0\n0.0 0.0\n")
    return str(p)


@pytest.fixture
def tcabr_file():
    return os.path.join(DATA_DIR, "tcabr_first_wall.txt")


# ── Construction ──────────────────────────────────────────────────────────────

def test_default_construction():
    c = Collider()
    assert not c.is_loaded()


# ── Loading ───────────────────────────────────────────────────────────────────

def test_load_square(square_file):
    c = Collider()
    assert c.load_shape(square_file)
    assert c.is_loaded()
    assert len(c.get_vertices()) == 4


def test_load_octagon(octagon_file):
    c = Collider()
    assert c.load_shape(octagon_file)
    assert c.is_loaded()
    assert len(c.get_vertices()) == 8


def test_load_invalid(invalid_file):
    c = Collider()
    assert not c.load_shape(invalid_file)


def test_load_nonexistent():
    c = Collider()
    assert not c.load_shape("/nonexistent/path.txt")


def test_load_tcabr(tcabr_file):
    c = Collider()
    assert c.load_shape(tcabr_file)
    assert c.is_loaded()
    assert len(c.get_vertices()) == 68


# ── Inside/outside tests ──────────────────────────────────────────────────────

def test_square_inside(square_file):
    c = Collider()
    c.load_shape(square_file)
    assert c.inside(0.5, 0.5)   # centre
    assert not c.inside(1.5, 0.5)  # clearly outside
    assert c.inside(0.0, 0.5)   # on left edge


def test_octagon_inside(octagon_file):
    c = Collider()
    c.load_shape(octagon_file)
    assert c.inside(0.5, 0.5)   # centre
    assert not c.inside(0.1, 0.1)  # outside


def test_inside_without_loading():
    c = Collider()
    assert not c.inside(0.5, 0.5)


def test_reversed_winding_square(reversed_square_file):
    c = Collider()
    assert c.load_shape(reversed_square_file)
    assert c.inside(0.5, 0.5)
    assert not c.inside(1.5, 0.5)


# ── get_vertices shape ────────────────────────────────────────────────────────

def test_get_vertices_returns_numpy_array(square_file):
    c = Collider()
    c.load_shape(square_file)
    verts = c.get_vertices()
    assert isinstance(verts, np.ndarray)
    assert verts.shape == (4, 2)
    assert verts.dtype == np.float64


def test_get_vertices_values(square_file):
    c = Collider()
    c.load_shape(square_file)
    verts = c.get_vertices()
    expected = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    np.testing.assert_array_equal(verts, expected)
