"""
Magnetic footprint plotting — demonstrates all input modes of maglib.plot_footprint().

Usage
-----
    # First generate the output files:
    PYTHONPATH=<maglib_root>/python python examples/fpgen/run.py

    # Then plot:
    PYTHONPATH=<maglib_root>/python python examples/fpgen/plot.py

plot_footprint() accepts three source types:

  1. File path   — load from .dat / .csv / .npy / .npz on disk
  2. numpy array — (N, 6) array held in memory
  3. Footprint object — only available immediately after run.py in the same session

All four sub-plots are shown: "cl" (connection length), "psi" (psi_N min),
"turns" (toroidal turns), and "au" (combined proxy). The "all" shortcut
produces all four in sequence.
"""

from pathlib import Path
import numpy as np
import maglib

# ── Case 3: Footprint object (requires a live computation) ───────────────────
# Only available immediately after running the integration in the same session.
# Uncomment this block and comment out cases 1/2 if you want to avoid
# loading from disk.
#
#   from run import build_tracer, Footprint, NRZ, NPHI, MAX_TURNS
#   from run import MANIFOLD, GRID_R1, GRID_Z1, GRID_R2, GRID_Z2, DPHI_INIT, DPHI_MIN, DPHI_MAX, N_THREADS
#   pairs   = [build_tracer() for _ in range(N_THREADS)]
#   tracers = [t for _, t in pairs]
#   fp = Footprint(MANIFOLD, GRID_R1, GRID_Z1, GRID_R2, GRID_Z2, NRZ, NPHI, MAX_TURNS)
#   fp.run(tracers)
#   maglib.plot_footprint(fp, which_plot="all")


# ── Case 1: file path ─────────────────────────────────────────────────────────
# Load directly from the .dat file saved by run.py. Produces all four sub-plots.

def plot_from_file():
    print("Case 1: all sub-plots from footprint.dat")
    maglib.plot_footprint("footprint.dat", which_plot="all")


# ── Case 2: numpy array ────────────────────────────────────────────────────────
# Load the array manually, then pass it to plot_footprint. Useful when you want
# to inspect or filter the data before plotting.

def plot_from_array():
    print("Case 2: connection-length sub-plot from numpy array")
    # .dat format: space-separated, first line is a comment
    data = np.loadtxt("footprint.dat", comments="#")
    maglib.plot_footprint(data, which_plot="cl")


def plot_from_npz():
    print("Case 3: psi sub-plot from .npz file")
    arr = np.load("footprint.npz")["data"]
    maglib.plot_footprint(arr, which_plot="psi")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    plot_from_file()
    plot_from_array()
    plot_from_npz()


if __name__ == "__main__":
    main()
