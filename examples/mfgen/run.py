"""
Invariant manifold computation — runs the integration and saves outputs.

Usage
-----

Computes both the stable (stability=0) and unstable (stability=1) manifolds
of the X-point and saves them to .npz and .dat files for use with plot.py.

Data paths default to tests/data/ relative to the repo root.
Override DATA_DIR at the top of the script to point to a different dataset.

The number of mapped segments (N_SEGMENTS) is set to 6 for faster computing
upon testing. Resolving the full topology of the manifolds may require resolving
more segments depending on the value set to EPSILON.
"""

from pathlib import Path

from maglib import M3DC1Source, Maglit, Manifold

# ── Data paths ────────────────────────────────────────────────────────────────

_REPO = Path(__file__).resolve().parents[2]
DATA_DIR = _REPO / "tests" / "data"

HDF5_PATH = str(DATA_DIR / "C1.h5")
WALL_PATH = str(DATA_DIR / "tcabr_first_wall.txt")

# ── Tracing parameters (mirror mfgen_input.txt defaults) ──────────────────────

TIMESLICE      = 1        # 0 = vacuum, 1 = plasma response
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
THETA_LIM      = 20.0     # turning-angle refinement threshold (deg)

# Initial X-point guess (from TCABR equilibrium metadata)
R_XPOINT_GUESS = 0.497999
Z_XPOINT_GUESS = -0.218603


# ── Helpers ───────────────────────────────────────────────────────────────────

def build_source():
    print(f"Loading M3DC1Source: {HDF5_PATH}  (timeslice {TIMESLICE})")
    source = M3DC1Source(HDF5_PATH, TIMESLICE)
    if not source.is_valid():
        raise RuntimeError(f"Failed to load M3DC1Source from {HDF5_PATH}")
    return source


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
    stable   = grow_manifold(build_tracer(source), stability=0)
    unstable = grow_manifold(build_tracer(source), stability=1)

    stable.save("stable.npz")
    stable.save("stable.dat")
    unstable.save("unstable.npz")
    unstable.save("unstable.dat")
    print("\nSaved: stable.npz  stable.dat  unstable.npz  unstable.dat")


if __name__ == "__main__":
    main()
