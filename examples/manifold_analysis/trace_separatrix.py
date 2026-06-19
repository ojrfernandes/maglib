"""
Equilibrium separatrix tracing example.

Calls maglib.trace_separatrix() with the equilibrium field (timeslice=-1) to
produce a closed (N, 2) curve that can be loaded by lobe_analysis.py.

Usage
-----
    PYTHONPATH=<maglib_root>/python python examples/manifold_analysis/trace_separatrix.py

Outputs (written to the current working directory)
------
    separatrix.npz  — single-key archive (key "seg_0")
    separatrix.dat  — space-separated text, 3 columns: seg R Z
"""

from pathlib import Path
import numpy as np
import maglib

# ── Data paths ────────────────────────────────────────────────────────────────

_REPO    = Path(__file__).resolve().parents[2]
DATA_DIR = _REPO / "tests" / "data"

HDF5_PATH = str(DATA_DIR / "C1.h5")
WALL_PATH = str(DATA_DIR / "tcabr_first_wall.txt")

# ── Parameters ────────────────────────────────────────────────────────────────

TIMESLICE = -1    # equilibrium field (no perturbation)
STABILITY = 0     # 0 = stable manifold; try 1 if curve goes the wrong way
PHI       = 0.0   # Poincaré section toroidal angle (rad)

N_SEGMENTS  = 30  # total segments including the primary; adjust until curve closes
N_INTERVALS = 9

L_LIM     = 0.005
THETA_LIM = 20.0

DPHI_INIT = 1e-2
DPHI_MIN  = 1e-6
DPHI_MAX  = 1e-2

R_XPOINT_GUESS = 0.497999
Z_XPOINT_GUESS = -0.218603

# ── Run ───────────────────────────────────────────────────────────────────────

def main():
    print(f"Loading M3DC1Source: {HDF5_PATH}  (timeslice {TIMESLICE})")
    source = maglib.M3DC1Source(HDF5_PATH, TIMESLICE)
    if not source.is_valid():
        raise RuntimeError(f"Failed to load M3DC1Source from {HDF5_PATH}")

    tracer = maglib.Maglit(source)
    tracer.configure(DPHI_INIT, DPHI_MIN, DPHI_MAX)
    tracer.set_monitor(WALL_PATH)

    print(f"\nTracing equilibrium separatrix "
          f"(stability={STABILITY}, {N_SEGMENTS} segments) ...")
    closed, x_point = maglib.trace_separatrix(
        tracer, phi=PHI,
        r_xpoint=R_XPOINT_GUESS, z_xpoint=Z_XPOINT_GUESS,
        stability=STABILITY,
        n_intervals=N_INTERVALS, n_segments=N_SEGMENTS,
        l_lim=L_LIM, theta_lim=THETA_LIM,
        verbose=True,
    )
    print(f"\nX-point: R = {x_point[0]:.6f} m,  Z = {x_point[1]:.6f} m")
    print(f"Separatrix: {len(closed)-1} traced + 1 closure = {len(closed)} points")

    np.savez_compressed("separatrix.npz", seg_0=closed)

    with open("separatrix.dat", "w") as f:
        f.write("# seg                  R                   Z\n")
        for row in closed:
            f.write(f"0  {row[0]:.16f}  {row[1]:.16f}\n")

    print("Saved: separatrix.npz  separatrix.dat")


if __name__ == "__main__":
    main()
