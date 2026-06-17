"""
Lorenz attractor — demonstrates the Sode step() loop for trajectory accumulation.

Usage
-----

The Lorenz system:
    dx/dt = σ(y - x)
    dy/dt = x(ρ - z) - y
    dz/dt = xy - βz

with σ=10, ρ=28, β=8/3 exhibits deterministic chaos. The attractor is bounded
but sensitive to initial conditions (Lyapunov exponent λ₁ ≈ 0.906).
"""

import numpy as np
import matplotlib.pyplot as plt

from maglib import Sode, SodeMethod, SodeStatus

# ── Lorenz parameters ─────────────────────────────────────────────────────────

SIGMA = 10.0
RHO   = 28.0
BETA  = 8.0 / 3.0


def lorenz(r, t):
    # r: [x, y, z], shape (3,); t unused (autonomous system)
    return np.array([
        SIGMA * (r[1] - r[0]),
        r[0] * (RHO - r[2]) - r[1],
        r[0] * r[1] - BETA * r[2],
    ])


# ── Solver setup ──────────────────────────────────────────────────────────────

def make_solver():
    solver = Sode(SodeMethod.RK78_DP, dim=3)
    solver.configure_steps(h_init=1e-3, h_min=1e-9, h_max=0.05)
    solver.configure_tol(tol_sol=1e-10, tol_mon=1e-12, tol_end=1e-12, damp=0.9)
    solver.set_system(lorenz)
    return solver


# ── Trajectory integration ────────────────────────────────────────────────────

def integrate_trajectory(x0: np.ndarray, t_end: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Integrate the Lorenz system from x0 to t_end, collecting all adaptive steps.

    Returns
    -------
    traj : ndarray, shape (N, 3)
        State vectors at each accepted step.
    times : ndarray, shape (N,)
        Corresponding times.
    """
    solver = make_solver()
    solver.reset()

    xs = [x0.copy()]
    ts = [0.0]
    x  = x0.copy()
    t  = 0.0
    last_report = 0.0

    while True:
        x, t, status = solver.step(x, t, t_end)
        xs.append(x)
        ts.append(t)

        if t - last_report >= 0.1 * t_end:
            print(f"  t = {t:.1f} / {t_end:.1f}", flush=True)
            last_report = t

        if status == SodeStatus.SUCCESS_TIME:
            break
        if status == SodeStatus.FAILED:
            raise RuntimeError(f"Sode step failed at t={t:.4f} (step size below h_min).")
        if status not in (SodeStatus.CONTINUE_GOOD, SodeStatus.CONTINUE_BAD):
            raise RuntimeError(f"Unexpected status: {status}")

    return np.array(xs), np.array(ts)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    T_END = 50.0

    # Two trajectories from nearby initial conditions — illustrates sensitivity
    x0_a = np.array([0.1,  0.1,  0.1])
    x0_b = np.array([0.1 + 1e-8, 0.1, 0.1])

    print(f"Integrating Lorenz attractor to t={T_END} ...")
    traj_a, times_a = integrate_trajectory(x0_a, T_END)
    traj_b, times_b = integrate_trajectory(x0_b, T_END)
    print(f"  Trajectory A: {len(times_a)} steps")
    print(f"  Trajectory B: {len(times_b)} steps")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    projections = [
        (0, 1, "x", "y"),
        (0, 2, "x", "z"),
        (1, 2, "y", "z"),
    ]
    for ax, (i, j, xl, yl) in zip(axes.flat, projections):
        ax.plot(traj_a[:, i], traj_a[:, j], lw=0.3, alpha=0.8)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.set_title(f"{xl}-{yl} projection")

    # Divergence of the two trajectories
    ax_div = axes[1, 1]
    n = min(len(times_a), len(times_b))
    ax_div.semilogy(times_a[:n], np.linalg.norm(traj_a[:n] - traj_b[:n], axis=1) + 1e-16,
                    lw=1.0, color="tab:red", label=r"$\|\delta x(t)\|$")
    ax_div.set_xlabel("t")
    ax_div.set_ylabel(r"$\|\delta x\|$")
    ax_div.set_title("Sensitivity to initial conditions")
    ax_div.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
