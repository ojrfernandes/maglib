"""
Python tests for M3DC1Source and Maglit.

Mirrors the C++ GoogleTest suite in tests/test_maglit.cpp. Requires
tests/data/C1.h5 and tests/data/tcabr_first_wall.txt.
"""

import math
import os

import numpy as np
import pytest
import maglib
from maglib import M3DC1Source, Maglit, SodeStatus

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
C1_H5    = os.path.join(DATA_DIR, "C1.h5")
WALL_TXT = os.path.join(DATA_DIR, "tcabr_first_wall.txt")

TOL = 1e-6


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def source_vac():
    src = M3DC1Source(C1_H5, -1)
    assert src.is_valid(), "vacuum M3DC1Source failed to open"
    return src


@pytest.fixture(scope="module")
def source_resp():
    src = M3DC1Source(C1_H5, 1)
    assert src.is_valid(), "response M3DC1Source failed to open"
    return src


@pytest.fixture(scope="module")
def tracer_vac(source_vac):
    t = Maglit(source_vac)
    t.configure(0.01, 1e-6, 0.1)
    return t


@pytest.fixture(scope="module")
def tracer_resp(source_resp):
    t = Maglit(source_resp)
    t.configure(0.01, 1e-6, 0.1)
    return t


# ── M3DC1Source ───────────────────────────────────────────────────────────────

def test_m3dc1_construction_valid(source_vac):
    assert source_vac.is_valid()


def test_m3dc1_construction_invalid():
    src = M3DC1Source("/nonexistent/path.h5", -1)
    assert not src.is_valid()


def test_m3dc1_two_timeslices(source_vac, source_resp):
    assert source_vac.is_valid()
    assert source_resp.is_valid()


# ── Maglit construction ───────────────────────────────────────────────────────

def test_maglit_construction(tracer_vac):
    assert tracer_vac is not None


# ── Magnetic field evaluation ─────────────────────────────────────────────────

def test_calc_mag_field(tracer_vac):
    ok, B = tracer_vac.calc_mag_field(R=0.7, phi=0.0, Z=0.0)
    assert ok
    assert isinstance(B, np.ndarray)
    assert B.shape == (3,)
    assert abs(B[0] - (-0.00401528)) < TOL   # B_R
    assert abs(B[1] - (-0.95980778)) < TOL   # B_phi
    assert abs(B[2] - (-0.14235661)) < TOL   # B_Z


def test_calc_mag_field_nonzero(tracer_vac):
    ok, B = tracer_vac.calc_mag_field(0.62, 0.0, 0.0)
    assert ok
    assert np.all(np.isfinite(B))
    assert np.linalg.norm(B) > 0.0


# ── Flux evaluation ───────────────────────────────────────────────────────────

def test_psi_eval(tracer_vac):
    ok, psi = tracer_vac.psi_eval(R=0.7, phi=0.0, Z=0.0)
    assert ok
    assert isinstance(psi, float)
    assert abs(psi - (-0.00988012)) < TOL


def test_psin_eval(tracer_vac):
    ok, psin = tracer_vac.psin_eval(R=0.7, phi=0.0, Z=0.0)
    assert ok
    assert isinstance(psin, float)
    assert abs(psin - 0.328436) < TOL


# ── Step integration ──────────────────────────────────────────────────────────

def _integrate(tracer, R, Z, phi, phi_max, dir=0):
    """Drive step() in a loop until a terminal status is reached."""
    tracer.reset()
    status = SodeStatus.CONTINUE_GOOD
    while status in (SodeStatus.CONTINUE_GOOD, SodeStatus.CONTINUE_BAD):
        R, Z, phi, status = tracer.step(R, Z, phi, phi_max, dir)
    return R, Z, phi, status


def test_step_returns_tuple(tracer_vac):
    tracer_vac.reset()
    result = tracer_vac.step(0.7, 0.0, 0.0, math.pi / 2)
    assert len(result) == 4
    R, Z, phi, status = result
    assert isinstance(R, float)
    assert isinstance(Z, float)
    assert isinstance(phi, float)


def test_step_reaches_phi_max(tracer_vac):
    phi_max = math.pi / 2
    R, Z, phi, status = _integrate(tracer_vac, 0.7, 0.0, 0.0, phi_max)
    assert status == SodeStatus.SUCCESS_TIME
    assert abs(phi - phi_max) < TOL


def test_inverse_map_roundtrip(tracer_vac):
    """Forward half-turn then backward half-turn must return to start."""
    R0, Z0, phi0 = 0.7, 0.0, 0.0
    phi_half = math.pi / 2

    try:
        tracer_vac.inverse_map(False)
        R1, Z1, phi1, status = _integrate(tracer_vac, R0, Z0, phi0, phi_half)
        assert status == SodeStatus.SUCCESS_TIME

        tracer_vac.inverse_map(True)
        phi_back = phi1 + phi_half
        R2, Z2, phi2, status = _integrate(tracer_vac, R1, Z1, phi1, phi_back)
        assert status == SodeStatus.SUCCESS_TIME
        assert abs(R2 - R0) < TOL
        assert abs(Z2 - Z0) < TOL
    finally:
        tracer_vac.inverse_map(False)  # restore default even on failure


# ── Boundary monitor ──────────────────────────────────────────────────────────

def test_set_monitor(tracer_resp):
    tracer_resp.set_monitor(WALL_TXT)


def test_monitor_stops_at_wall(tracer_resp):
    """A field line starting inside the vessel should eventually hit the wall."""
    tracer_resp.set_monitor(WALL_TXT)
    phi_max = 10000 * 2 * math.pi

    last_R = last_Z = 0.0
    R, Z, phi = 0.5, -0.225, 0.0
    tracer_resp.reset()
    tracer_resp.inverse_map(False)

    status = SodeStatus.CONTINUE_GOOD
    while status in (SodeStatus.CONTINUE_GOOD, SodeStatus.CONTINUE_BAD):
        last_R, last_Z = R, Z
        R, Z, phi, status = tracer_resp.step(R, Z, phi, phi_max, -1)

    assert status == SodeStatus.SUCCESS_MONITOR
    assert tracer_resp.boundary.inside(last_R, last_Z)
    assert not tracer_resp.boundary.inside(R, Z)
