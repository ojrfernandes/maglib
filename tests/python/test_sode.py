"""
Python tests for the Sode ODE solver.

Mirror the C++ GoogleTest suite in tests/test_sode.cpp so that correctness is
verified at both the C++ and Python binding layers.
"""

import math
import numpy as np
import pytest
import maglib
from maglib import Sode, SodeMethod, SodeStatus


TOL = 1e-6  # floating-point comparison tolerance


# ── Helper ────────────────────────────────────────────────────────────────────

def make_sho_solver(omega=1.0):
    """Cash-Karp solver pre-configured for a simple harmonic oscillator."""
    solver = Sode(SodeMethod.RK56_CK, dim=2)
    solver.configure_steps(1e-3, 1e-6, 0.1)
    solver.configure_tol(1e-10, 1e-12, 1e-12, 0.9)
    solver.set_system(lambda x, t: np.array([x[1], -(omega**2) * x[0]]))
    return solver


# ── Construction and configuration ───────────────────────────────────────────

def test_construction_all_methods():
    for method in (SodeMethod.RK56_FB, SodeMethod.RK56_CK, SodeMethod.RK78_DP):
        s = Sode(method, dim=3)
        s.configure_steps(0.01, 1e-6, 0.1)
        s.configure_tol(1e-8, 1e-12, 1e-12, 0.9)


def test_set_system_required():
    solver = Sode(SodeMethod.RK56_CK, dim=1)
    solver.configure_steps(1e-3, 1e-6, 0.1)
    with pytest.raises(RuntimeError, match="set_system"):
        solver.integrate(np.array([1.0]), t0=0.0, t_end=1.0)


# ── Linear ODE: dx/dt = -x  →  x(t) = exp(-t) ───────────────────────────────

def test_linear_ode():
    solver = Sode(SodeMethod.RK56_CK, dim=1)
    solver.configure_steps(1e-3, 1e-6, 0.1)
    solver.configure_tol(1e-10, 1e-12, 1e-12, 0.9)
    solver.set_system(lambda x, t: np.array([-x[0]]))

    x, t, status = solver.integrate(np.array([1.0]), t0=0.0, t_end=1.0)

    assert status == SodeStatus.SUCCESS_TIME
    assert abs(t - 1.0) < TOL
    assert abs(x[0] - math.exp(-1.0)) < TOL


# ── Simple harmonic oscillator ────────────────────────────────────────────────

def test_sho_returns_to_initial_position():
    """After one full period the SHO must return to (x=1, v=0)."""
    solver = make_sho_solver(omega=1.0)
    x, t, status = solver.integrate(np.array([1.0, 0.0]), t0=0.0, t_end=2 * math.pi)

    assert status == SodeStatus.SUCCESS_TIME
    assert abs(x[0] - 1.0) < TOL
    assert abs(x[1] - 0.0) < TOL


def test_sho_energy_conservation():
    """Mechanical energy must be conserved over one full period."""
    omega = 1.0
    x0 = np.array([1.0, 0.0])
    E0 = 0.5 * (x0[1] ** 2 + omega**2 * x0[0] ** 2)

    solver = make_sho_solver(omega=omega)
    x, _, status = solver.integrate(x0, t0=0.0, t_end=2 * math.pi)

    E1 = 0.5 * (x[1] ** 2 + omega**2 * x[0] ** 2)
    assert status == SodeStatus.SUCCESS_TIME
    assert abs(E1 - E0) < TOL


# ── Event detection ───────────────────────────────────────────────────────────

def test_event_detection_zero_crossing():
    """SHO starting at x=1 crosses zero at t=π/2."""
    solver = make_sho_solver(omega=1.0)
    solver.set_monitor(lambda x, t: x[0] <= 0.0)

    # Start with monitor=False (x[0]=1>0); stop on False→True (dir=1)
    x, t, status = solver.integrate(np.array([1.0, 0.0]), t0=0.0, t_end=10.0, monitor_dir=1)

    assert status == SodeStatus.SUCCESS_MONITOR
    assert abs(x[0]) < TOL
    assert abs(t - math.pi / 2) < TOL


def test_event_detection_no_trigger():
    """With monitor_dir=0 the monitor is ignored and integration runs to t_end."""
    solver = make_sho_solver(omega=1.0)
    solver.set_monitor(lambda x, t: x[0] <= 0.0)

    x, t, status = solver.integrate(np.array([1.0, 0.0]), t0=0.0, t_end=2 * math.pi,
                                    monitor_dir=0)

    assert status == SodeStatus.SUCCESS_TIME
    assert abs(t - 2 * math.pi) < TOL


# ── Solver consistency across methods ────────────────────────────────────────

def test_solver_consistency_lorenz():
    """All three RK methods must agree on the Lorenz attractor at t=1."""
    sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0
    x0 = np.array([0.1, 0.1, 0.1])

    def lorenz(x, t):
        return np.array([
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])

    results = []
    for method in (SodeMethod.RK56_FB, SodeMethod.RK56_CK, SodeMethod.RK78_DP):
        s = Sode(method, dim=3)
        s.configure_steps(1e-3, 1e-6, 0.1)
        s.configure_tol(1e-10, 1e-12, 1e-12, 0.9)
        s.set_system(lorenz)
        x, _, _ = s.integrate(x0.copy(), t0=0.0, t_end=1.0)
        results.append(x)

    for i in range(3):
        assert abs(results[0][i] - results[1][i]) < TOL
        assert abs(results[0][i] - results[2][i]) < TOL


def test_lorenz_stays_bounded():
    """Lorenz attractor is bounded; the solution must not blow up."""
    sigma, rho, beta = 10.0, 28.0, 8.0 / 3.0

    solver = Sode(SodeMethod.RK56_CK, dim=3)
    solver.configure_steps(1e-3, 1e-6, 0.1)
    solver.configure_tol(1e-10, 1e-12, 1e-12, 0.9)
    solver.set_system(lambda x, t: np.array([
        sigma * (x[1] - x[0]),
        x[0] * (rho - x[2]) - x[1],
        x[0] * x[1] - beta * x[2],
    ]))

    x, _, status = solver.integrate(np.array([0.1, 0.1, 0.1]), t0=0.0, t_end=10.0)

    assert status == SodeStatus.SUCCESS_TIME
    assert all(abs(xi) < 50.0 for xi in x)


# ── Step-by-step interface ────────────────────────────────────────────────────

def test_manual_step_loop():
    """step() in a manual loop must give the same result as integrate()."""
    solver_a = make_sho_solver()
    solver_b = make_sho_solver()

    x0 = np.array([1.0, 0.0])
    t_end = math.pi

    x_int, t_int, _ = solver_a.integrate(x0.copy(), t0=0.0, t_end=t_end)

    # Manual loop with step()
    solver_b.reset()
    x_step = x0.copy()
    t_step = 0.0
    status = SodeStatus.CONTINUE_GOOD
    while status in (SodeStatus.CONTINUE_GOOD, SodeStatus.CONTINUE_BAD):
        x_step, t_step, status = solver_b.step(x_step, t_step, t_end)

    assert abs(t_step - t_int) < TOL
    assert abs(x_step[0] - x_int[0]) < TOL
    assert abs(x_step[1] - x_int[1]) < TOL


# ── Output shapes ─────────────────────────────────────────────────────────────

def test_output_is_numpy_array():
    solver = make_sho_solver()
    x, t, status = solver.integrate(np.array([1.0, 0.0]), t0=0.0, t_end=1.0)
    assert isinstance(x, np.ndarray)
    assert x.shape == (2,)
    assert isinstance(t, float)
