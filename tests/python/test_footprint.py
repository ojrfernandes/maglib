"""
Python tests for Footprint.

Covers serial and parallel execution, output_data shape and values,
and all four save formats (.dat, .txt, .csv, .npy, .npz).
"""

import math
import os

import numpy as np
import pytest
import maglib
from maglib import M3DC1Source, Maglit, Footprint

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
C1_H5    = os.path.join(DATA_DIR, "C1.h5")
WALL_TXT = os.path.join(DATA_DIR, "tcabr_first_wall.txt")

TOL = 1e-3   # matches C++ reference-value tolerance


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_tracer():
    src = M3DC1Source(C1_H5, 1)
    t   = Maglit(src)
    t.configure(0.01, 1e-6, 0.1)
    t.set_monitor(WALL_TXT)
    return src, t   # caller must keep src alive


@pytest.fixture(scope="module")
def tracer_and_source():
    src, t = _make_tracer()
    return src, t


@pytest.fixture(scope="module")
def wall_fp(tracer_and_source):
    """Pre-run 2x2 wall footprint (stable manifold, 50 turns)."""
    _, t = tracer_and_source
    fp = Footprint(0, 0.435, -0.239, 0.435, -0.232, 2, 2, 50)
    fp.run([t])
    return fp


# ── Construction ──────────────────────────────────────────────────────────────

def test_construction():
    fp = Footprint(0, 0.435, -0.239, 0.435, -0.232, 2, 2, 50)
    assert fp is not None


def test_output_data_shape_before_run():
    fp = Footprint(0, 0.435, -0.239, 0.435, -0.232, 3, 5, 50)
    data = fp.output_data
    assert isinstance(data, np.ndarray)
    assert data.shape == (15, 6)   # nRZ * nPhi = 3 * 5
    assert data.dtype == np.float64


# ── Serial execution ──────────────────────────────────────────────────────────

def test_serial_run_shape(wall_fp):
    assert wall_fp.output_data.shape == (4, 6)   # 2 * 2


def test_serial_run_initial_positions(wall_fp):
    """R0, Z0, phi0 must match the grid definition."""
    data = wall_fp.output_data
    # row 0: i=0, j=0 → phi=0, (R,Z) = grid_R1, grid_Z1
    assert abs(data[0, 0] - 0.435) < TOL
    assert abs(data[0, 1] - (-0.239)) < TOL
    assert abs(data[0, 2] - 0.0) < TOL
    # row 2: i=1, j=0 → phi=pi
    assert abs(data[2, 2] - math.pi) < TOL


def test_serial_run_reference_values(wall_fp):
    """Mirror the C++ Footprint_RunGrid_Simple2x2_Wall reference values."""
    data = wall_fp.output_data
    assert abs(data[0, 3] - 4.699479)  < TOL   # connection length
    assert abs(data[1, 3] - 10.802152) < TOL
    assert abs(data[0, 4] - 0.995964)  < TOL   # psiMin
    assert abs(data[1, 4] - 1.001697)  < TOL
    assert data[0, 5] == 1                       # turn count
    assert data[1, 5] == 3


# ── Parallel execution ────────────────────────────────────────────────────────

def test_parallel_matches_serial(tracer_and_source):
    """2-thread result must be bit-for-bit identical to serial."""
    _, t1 = tracer_and_source
    src2, t2 = _make_tracer()

    fp_serial   = Footprint(0, 0.435, -0.239, 0.435, -0.232, 2, 2, 50)
    fp_parallel = Footprint(0, 0.435, -0.239, 0.435, -0.232, 2, 2, 50)

    fp_serial.run([t1])
    fp_parallel.run([t1, t2])

    np.testing.assert_array_equal(fp_serial.output_data, fp_parallel.output_data)


# ── Save formats ──────────────────────────────────────────────────────────────

def _round_trip_text(wall_fp, tmp_path, suffix, delimiter):
    """Save to text format and verify the header and data round-trip."""
    path = str(tmp_path / f"out{suffix}")
    wall_fp.save(path)
    assert os.path.exists(path)

    with open(path) as f:
        lines = f.readlines()

    # header line present
    assert len(lines) == 5   # 1 header + 4 data rows

    if delimiter == ",":
        assert lines[0].strip() == "R0,Z0,phi0,length,psiMin,turn"
        loaded = np.loadtxt(path, delimiter=delimiter, skiprows=1)
    else:
        assert lines[0].startswith("#")
        loaded = np.loadtxt(path, delimiter=delimiter, comments="#")

    np.testing.assert_allclose(loaded, wall_fp.output_data, rtol=1e-14)


def test_save_dat(wall_fp, tmp_path):
    _round_trip_text(wall_fp, tmp_path, ".dat", " ")


def test_save_txt(wall_fp, tmp_path):
    _round_trip_text(wall_fp, tmp_path, ".txt", " ")


def test_save_csv(wall_fp, tmp_path):
    _round_trip_text(wall_fp, tmp_path, ".csv", ",")


def test_save_npy(wall_fp, tmp_path):
    path = str(tmp_path / "out.npy")
    wall_fp.save(path)
    assert os.path.exists(path)
    loaded = np.load(path)
    np.testing.assert_array_equal(loaded, wall_fp.output_data)


def test_save_npz(wall_fp, tmp_path):
    path = str(tmp_path / "out.npz")
    wall_fp.save(path)
    assert os.path.exists(path)
    loaded = np.load(path)["data"]
    np.testing.assert_array_equal(loaded, wall_fp.output_data)


def test_save_unknown_extension_raises(wall_fp, tmp_path):
    with pytest.raises(RuntimeError):
        wall_fp.save(str(tmp_path / "out.xyz"))
