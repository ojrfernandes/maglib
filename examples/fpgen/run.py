"""
Magnetic footprint computation — runs the integration and saves outputs.

Usage
-----

Traces a 2-D grid of field lines from a divertor target segment through nPhi
toroidal angles. One Maglit tracer per CPU thread is created for parallel
execution. Outputs are saved to .dat and .npz for use with plot.py.

Grid here mirrors the fpgen_input.txt reference case: a horizontal segment
at Z = -0.24 m spanning R = 0.50 – 0.55 m (TCABR lower divertor target).
"""

import os
from pathlib import Path

from maglib import M3DC1Source, Maglit, Footprint

# ── Data paths ────────────────────────────────────────────────────────────────

_REPO    = Path(__file__).resolve().parents[2]
DATA_DIR = _REPO / "tests" / "data"

HDF5_PATH = str(DATA_DIR / "C1.h5")
WALL_PATH = str(DATA_DIR / "tcabr_first_wall.txt")

# ── Parameters (mirror fpgen_input.txt defaults) ──────────────────────────────

TIMESLICE  = 1        # 0 = vacuum, 1 = plasma response
MANIFOLD   = 1        # 0 = stable (backward map), 1 = unstable (forward map)

GRID_R1    = 0.50     # divertor target first endpoint (m)
GRID_Z1    = -0.24
GRID_R2    = 0.55     # divertor target second endpoint (m)
GRID_Z2    = -0.24

NRZ        = 20       # grid points along the target segment
NPHI       = 20       # toroidal starting angles (uniformly spaced over 2π)
MAX_TURNS  = 50       # field line lost after this many toroidal turns

DPHI_INIT  = 1e-2    # initial step size
DPHI_MIN   = 1e-6    # minimum step size
DPHI_MAX   = 1e-2    # maximum step size

N_THREADS  = min(4, os.cpu_count() or 1)


# ── Helpers ───────────────────────────────────────────────────────────────────

def build_tracer() -> tuple:
    """Create one (M3DC1Source, Maglit) pair. Both must be kept alive together."""
    src = M3DC1Source(HDF5_PATH, TIMESLICE)
    if not src.is_valid():
        raise RuntimeError(f"Failed to load M3DC1Source from {HDF5_PATH}")
    t = Maglit(src)
    t.configure(DPHI_INIT, DPHI_MIN, DPHI_MAX)
    t.set_monitor(WALL_PATH)
    return src, t


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print(f"Loading {N_THREADS} M3DC1Source(s): {HDF5_PATH}  (timeslice {TIMESLICE})")

    # Each thread needs its own source + tracer pair. The source must remain
    # alive as long as the tracer is in use, so we keep all pairs in a list.
    pairs   = [build_tracer() for _ in range(N_THREADS)]
    tracers = [t for _, t in pairs]

    label = "unstable" if MANIFOLD else "stable"
    print(f"\nGrid: {label} manifold, R=[{GRID_R1}, {GRID_R2}] m, Z={GRID_Z1} m")
    print(f"nRZ={NRZ}, nPhi={NPHI}, max_turns={MAX_TURNS}, threads={N_THREADS}")
    print("Running footprint integration... (this may take a moment)")

    fp = Footprint(MANIFOLD, GRID_R1, GRID_Z1, GRID_R2, GRID_Z2, NRZ, NPHI, MAX_TURNS)
    fp.run(tracers)

    data = fp.output_data
    print(f"Done. Output shape: {data.shape}  (nRZ={NRZ} × nPhi={NPHI} = {NRZ*NPHI} rows)")

    fp.save("footprint.dat")
    fp.save("footprint.npz")
    print("Saved: footprint.dat  footprint.npz")


if __name__ == "__main__":
    main()
