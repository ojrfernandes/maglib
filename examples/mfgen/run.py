"""
Invariant manifold computation — runs the integration and saves outputs.

Usage
-----

Computes both the unstable (stability=0) and stable (stability=1) manifolds
of the X-point and saves them to .npz and .dat files for use with plot.py.

Data paths default to tests/data/ relative to the repo root.
Override DATA_DIR at the top of the script to point to a different dataset.

The number of mapped segments (N_SEGMENTS) is set to 6 for faster computing
upon testing. Resolving the full topology of the manifolds may require resolving
more segments depending on the value set to EPSILON.
"""

import math
from pathlib import Path

from maglib import M3DC1Source, SuperpositionSource, Maglit, Manifold

# ── Data paths ────────────────────────────────────────────────────────────────

_REPO = Path(__file__).resolve().parents[2]
DATA_DIR = _REPO / "tests" / "data"

WALL_PATH = str(DATA_DIR / "tcabr_first_wall.txt")

# ── Field sources (mirror the [M3DC1 SOURCE] section of the input file) ──────
#
# Single source (nsources = 1): M3DC1Source; phase/amplitude are not applied.
# Multiple sources (nsources > 1): SuperpositionSource,
#     B(R,φ,Z) = B_eq(R,φ,Z) + Σ_i A_i · [B_i(R, φ−δ_i, Z) − B_eq(R,φ,Z)]
# The equilibrium B_eq is loaded automatically (timeslice = -1) from the first
# component's file. Use timeslice 0 (vacuum) or 1 (plasma response) per component.
#   path      : M3DC1 HDF5 file
#   timeslice : 0 = vacuum, 1 = full single-fluid response
#   phase     : toroidal phase shift δ_i (radians)
#   amplitude : linear scale factor A_i (dimensionless, may be negative)

SOURCES = [
    dict(path=str(DATA_DIR / "C1.h5"), timeslice=1, phase=0.0, amplitude=1.0),
]

# --- Three-coil superposition example (replace paths, then uncomment) ---
# SOURCES = [
#     dict(path="/path/to/IM_C1.h5", timeslice=1, phase=0.0, amplitude=1.0),
#     dict(path="/path/to/IL_C1.h5", timeslice=1, phase=math.radians(-100.0), amplitude=1.0),
#     dict(path="/path/to/IU_C1.h5", timeslice=1, phase=math.radians(80.0), amplitude=1.0),
# ]

# ── Tracing parameters (mirror mfgen_input.txt defaults) ──────────────────────

PHI            = 0.0      # Poincaré section toroidal angle (rad)

DPHI_INIT      = 1e-2     # initial stepsize for the integration
DPHI_MIN       = 1e-6     # minimal stepsize for the integration
DPHI_MAX       = 1e-2     # maximum stepsize for the integration

EPSILON        = 1e-8     # pivot distance from X-point along eigenvector
H_DERIV        = 1e-8     # step for numerical Jacobian
TOL_NEWTON     = 1e-14    # Newton convergence tolerance
MAX_ITER       = 50       # Maximum iteration for Newton method
PRECISION      = 1e-14    # Tolerance for floating point comparison
MAX_INSERTIONS = 50       # Maximum number of inserted points during refinement

N_INTERVALS    = 9        # points on primary segment = N_INTERVALS + 1
N_SEGMENTS     = 6        # total segments (including primary)
L_LIM          = 0.005    # arc-length refinement threshold (m)
THETA_LIM      = math.radians(20.0)  # turning-angle refinement threshold (rad)

# Initial X-point guess (from TCABR equilibrium metadata); the X-point is that of
# the equilibrium shared by all superposed components
R_XPOINT_GUESS = 0.497999
Z_XPOINT_GUESS = -0.218603


# ── Helpers ───────────────────────────────────────────────────────────────────

def build_source():
    """Create a FieldSource from SOURCES. Call once per thread (Fusion-IO is not thread-safe)."""
    print(f"Loading {len(SOURCES)} field source component(s):")
    for s in SOURCES:
        print(f"  {s['path']}  ts={s['timeslice']}  phase={s['phase']} rad  amp={s['amplitude']}")
    for s in SOURCES:
        if not Path(s["path"]).is_file():
            raise FileNotFoundError(s["path"])
    if len(SOURCES) == 1:
        src = M3DC1Source(SOURCES[0]["path"], SOURCES[0]["timeslice"])
    else:
        src = SuperpositionSource()
        for s in SOURCES:
            src.add_component(s["path"], s["timeslice"], s["phase"], s["amplitude"])
    if not src.is_valid():
        raise RuntimeError("Failed to load field source(s): " + ", ".join(s["path"] for s in SOURCES))
    return src


def build_tracer(source) -> Maglit:
    tracer = Maglit(source)
    tracer.configure(DPHI_INIT, DPHI_MIN, DPHI_MAX)
    tracer.set_monitor(WALL_PATH)
    return tracer


def grow_manifold(tracer: Maglit, stability: int) -> Manifold:
    label = "stable" if stability == 0 else "unstable"
    print(f"\n{'─'*60}")
    print(f"Computing {label} manifold (stability={stability}) ...")

    mf = Manifold(tracer, phi=PHI, stability=stability)
    mf.configure(
        epsilon         = EPSILON,
        h               = H_DERIV,
        tol             = TOL_NEWTON,
        max_iter        = MAX_ITER,
        precision_limit = PRECISION,
        max_insertions  = MAX_INSERTIONS,
    )

    print(f"Searching for X-point near (R={R_XPOINT_GUESS}, Z={Z_XPOINT_GUESS}) ...")
    if not mf.find_x_point(R_XPOINT_GUESS, Z_XPOINT_GUESS):
        raise RuntimeError("X-point Newton iteration did not converge.")
    xp = mf.x_point
    print(f"X-point: R = {xp[0]:.6f} m,  Z = {xp[1]:.6f} m")

    print(f"Computing primary segment ({N_INTERVALS + 1} points) ...")
    seg = mf.primary_segment(N_INTERVALS)
    print(f"  Primary: {seg.shape[0]} points")

    for i in range(1, N_SEGMENTS):
        print(f"Mapping segment {i + 1} / {N_SEGMENTS} ...")
        _, seg = mf.new_segment(seg, L_LIM, THETA_LIM)
        print(f"  Segment {i + 1}: {seg.shape[0]} points")

    total_pts = sum(s.shape[0] for s in mf.output_data)
    print(f"Done. {len(mf.output_data)} segments, {total_pts} total points.")
    return mf


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    source = build_source()

    # Each manifold needs its own Maglit instance: the Manifold constructor sets
    # inverse_map on the tracer, so stable and unstable cannot share one tracer.
    # stability=0 → forward map → unstable manifold
    # stability=1 → inverse map → stable manifold
    unstable = grow_manifold(build_tracer(source), stability=0)
    stable   = grow_manifold(build_tracer(source), stability=1)

    unstable.save("unstable.npz")
    unstable.save("unstable.dat")
    stable.save("stable.npz")
    stable.save("stable.dat")
    print("\nSaved: unstable.npz  unstable.dat  stable.npz  stable.dat")


if __name__ == "__main__":
    main()
