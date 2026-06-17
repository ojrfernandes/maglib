"""
Simple harmonic oscillator — demonstrates event detection with set_monitor().

Usage
-----

System: q'' + ω²q = 0  →  state [q, p] with dq/dt = p, dp/dt = -ω²q
Exact solution from (q₀, p₀) = (1, 0):  q(t) = cos(ωt),  p(t) = -ω sin(ωt)

Demonstrated here:
  1. Full trajectory via step() loop — compare with analytical solution.
  2. Event detection: locate q=0 zero-crossings using set_monitor().
  3. Period estimation from successive crossing times.
"""

import math
import numpy as np
import matplotlib.pyplot as plt

from maglib import Sode, SodeMethod, SodeStatus

OMEGA = 1.0          # angular frequency (rad / time unit)
T_END = 4 * math.pi  # two full periods


# ── System function ───────────────────────────────────────────────────────────

def sho_rhs(x, t):
    # x: [q, p]
    return np.array([x[1], -(OMEGA**2) * x[0]])


# ── Solver factory ────────────────────────────────────────────────────────────

def make_solver():
    solver = Sode(SodeMethod.RK56_CK, dim=2)
    solver.configure_steps(h_init=1e-3, h_min=1e-9, h_max=0.1)
    solver.configure_tol(tol_sol=1e-10, tol_mon=1e-12, tol_end=1e-12, damp=0.9)
    solver.set_system(sho_rhs)
    return solver


# ── Part 1: Full trajectory ───────────────────────────────────────────────────

def build_trajectory(x0: np.ndarray, t_end: float) -> tuple[np.ndarray, np.ndarray]:
    """Collect state at every accepted adaptive step."""
    solver = make_solver()
    solver.reset()

    xs = [x0.copy()]
    ts = [0.0]
    x  = x0.copy()
    t  = 0.0

    while True:
        x, t, status = solver.step(x, t, t_end)
        xs.append(x)
        ts.append(t)

        if status == SodeStatus.SUCCESS_TIME:
            break
        if status == SodeStatus.FAILED:
            raise RuntimeError(f"Sode step failed at t={t:.4f}")
        if status not in (SodeStatus.CONTINUE_GOOD, SodeStatus.CONTINUE_BAD):
            raise RuntimeError(f"Unexpected status: {status}")

    return np.array(xs), np.array(ts)


# ── Part 2: Event detection — locate q = 0 zero-crossings ────────────────────

def find_zero_crossings(n: int = 4) -> list[float]:
    """
    Find the times t_k where q crosses zero from positive to negative.
    With (q₀, p₀) = (1, 0) and ω=1 these occur at t = π/2 + k·π, k=0,1,...

    Strategy: after each crossing q just turned negative (monitor True).
    Reset; integrate with monitor_dir=-1 (True→False) to find the next
    positive-to-negative transition pair and so on.
    """
    monitor = lambda x, t: x[0] <= 0.0   # True when q ≤ 0

    crossings = []
    x = np.array([1.0, 0.0])
    t = 0.0

    for _ in range(n):
        solver = make_solver()
        solver.set_monitor(monitor)
        solver.reset()

        # Integrate until q goes from positive (False) to ≤ 0 (True)
        x_c, t_c, status = solver.integrate(x.copy(), t0=t, t_end=t + 10.0,
                                             monitor_dir=1)
        if status != SodeStatus.SUCCESS_MONITOR:
            print(f"Warning: expected SUCCESS_MONITOR, got {status}")
            break
        crossings.append(t_c)

        # Advance just past the crossing so the monitor starts False again
        # (q continues past zero and becomes negative; we need it to come back)
        solver2 = make_solver()
        solver2.set_monitor(monitor)
        solver2.reset()
        x_next, t_next, status2 = solver2.integrate(x_c.copy(), t0=t_c,
                                                     t_end=t_c + 10.0,
                                                     monitor_dir=-1)
        if status2 != SodeStatus.SUCCESS_MONITOR:
            break
        x = x_next
        t = t_next

    return crossings


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    x0 = np.array([1.0, 0.0])

    # --- Trajectory ---
    print("Building full trajectory ...")
    traj, times = build_trajectory(x0, T_END)
    print(f"  {len(times)} adaptive steps over t ∈ [0, {T_END:.2f}]")

    q_num = traj[:, 0]
    p_num = traj[:, 1]
    q_ana = np.cos(OMEGA * times)
    p_ana = -OMEGA * np.sin(OMEGA * times)

    max_err_q = np.max(np.abs(q_num - q_ana))
    max_err_p = np.max(np.abs(p_num - p_ana))
    print(f"  max |q_num - q_ana| = {max_err_q:.2e}")
    print(f"  max |p_num - p_ana| = {max_err_p:.2e}")

    # Energy: E = (p² + ω²q²)/2 — should be constant at 0.5
    E = 0.5 * (p_num**2 + OMEGA**2 * q_num**2)
    print(f"  Energy drift: {np.max(np.abs(E - 0.5)):.2e}  (initial E = 0.5)")

    # --- Zero crossings ---
    print("\nFinding zero-crossings (q: + → -) ...")
    crossings = find_zero_crossings(n=4)
    analytical = [math.pi / 2 + k * math.pi for k in range(len(crossings))]
    for k, (tc, ta) in enumerate(zip(crossings, analytical)):
        print(f"  crossing {k}: t_num = {tc:.8f}  t_ana = {ta:.8f}  "
              f"err = {abs(tc - ta):.2e}")

    # Period estimates from consecutive crossings separated by 2π
    if len(crossings) >= 3:
        periods = [crossings[i + 2] - crossings[i] for i in range(len(crossings) - 2)]
        print(f"\n  Period estimates (from every-other crossing): "
              f"{[f'{p:.6f}' for p in periods]}")
        print(f"  Analytical period = {2 * math.pi:.6f}")

    # --- Plots ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Phase portrait (q-p plane): should be a circle of radius 1
    ax = axes[0]
    ax.plot(q_num, p_num, lw=1.0, label="numerical")
    ax.plot(q_ana, p_ana, lw=1.5, ls="--", alpha=0.6, label="analytical")
    for tc in crossings:
        ax.axvline(0, color="gray", lw=0.5, ls=":")
    ax.set_xlabel("q (displacement)")
    ax.set_ylabel("p (momentum)")
    ax.set_title("SHO phase portrait")
    ax.set_aspect("equal")
    ax.legend()

    # q(t): numerical vs analytical
    ax = axes[1]
    ax.plot(times, q_num, lw=1.0, label="numerical")
    ax.plot(times, q_ana, lw=1.5, ls="--", alpha=0.6, label=r"$\cos(\omega t)$")
    for tc in crossings:
        ax.axvline(tc, color="tab:red", lw=0.8, ls=":", alpha=0.7)
    ax.set_xlabel("t")
    ax.set_ylabel("q")
    ax.set_title("q(t): numerical vs analytical  (red lines: zero-crossings)")
    ax.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
