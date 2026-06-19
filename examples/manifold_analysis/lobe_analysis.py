"""
Lobe analysis — characterizes lobes formed by the intersection of the perturbed
stable and unstable invariant manifolds with the equilibrium separatrix.

In a tokamak with magnetic perturbations (MPs), the stable and unstable manifolds
of the X-point split away from the equilibrium separatrix. The lobes of the resulting
homoclinic tangle are the regions bounded between consecutive intersection points of
a perturbed manifold with the equilibrium separatrix. Their area, effective width,
and poloidal position characterise the radial transport across the separatrix and
the heat load footprint pattern on the divertor targets.

Usage
-----
    # Step 1: trace the equilibrium separatrix (timeslice=-1, closes the loop)
    PYTHONPATH=<maglib_root>/python python examples/manifold_analysis/trace_separatrix.py
    # → separatrix.npz, separatrix.dat

    # Step 2: grow the perturbed manifolds (timeslice=1)
    PYTHONPATH=<maglib_root>/python python examples/mfgen/run.py
    # → stable.npz, unstable.npz

    # Step 3: run this analysis
    PYTHONPATH=<maglib_root>/python python examples/manifold_analysis/lobe_analysis.py

Inputs
------
    separatrix.npz               — closed equilibrium separatrix from trace_separatrix.py
    stable.npz / unstable.npz   — segment arrays from mfgen/run.py (timeslice=1)
    tcabr_first_wall.txt        — first wall contour [R, Z] (for plotting only)

All path defaults are resolved relative to the repository root.
The .npz files are expected in the working directory (where the scripts wrote them).

Workflow
--------
    1. Load and stack segments from all three .npz files into (N, 2) curves.
    2. Simplify all three curves with Ramer-Douglas-Peucker (maglib.simplify) to
       remove clustered points in folded regions before intersection search.
    3. Read magnetic axis from HDF5 (scalars/xmag, scalars/zmag).
    4. Run maglib.lobe_map(separatrix, manifold, mag_axis, x_point) for each
       perturbed manifold. Consecutive intersection points along the manifold
       bound each lobe.
    5. Print a summary table and produce an R-Z plot.
"""

from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt
import maglib

# ── Paths ─────────────────────────────────────────────────────────────────────

_REPO    = Path(__file__).resolve().parents[2]
DATA_DIR = _REPO / "tests" / "data"

HDF5_PATH  = str(DATA_DIR / "C1.h5")
WALL_PATH  = str(DATA_DIR / "tcabr_first_wall.txt")

SEPARATRIX_NPZ = "separatrix.npz"
STABLE_NPZ     = "stable.npz"
UNSTABLE_NPZ   = "unstable.npz"

# ── Parameters ────────────────────────────────────────────────────────────────

# RDP tolerance for simplify(). Dense clusters in highly folded regions can
# cause spurious near-parallel intersection artefacts.
SIMPLIFY_TOL = 5e-4   # metres

# X-point position (matches mfgen/run.py default).
R_XPOINT = 0.497999
Z_XPOINT = -0.218603


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_manifold(path: str) -> np.ndarray:
    """Stack all segments from an .npz file into one (N, 2) curve, dropping NaN/inf rows."""
    data = np.load(path)
    segs = [data[k] for k in sorted(data.files)]
    arr = np.vstack(segs)
    valid = np.all(np.isfinite(arr), axis=1)
    n_dropped = (~valid).sum()
    if n_dropped:
        print(f"  (dropped {n_dropped} non-finite rows from {path})")
    return arr[valid]


def load_wall(path: str) -> np.ndarray:
    """Load wall contour and close it if the first and last points differ."""
    wall = np.loadtxt(path)
    if not np.allclose(wall[0], wall[-1]):
        wall = np.vstack([wall, wall[0]])
    return wall


def read_mag_axis(hdf5_path: str, timeslice: int = 1) -> np.ndarray:
    """Read magnetic axis [R, Z] from M3DC1 HDF5 scalars."""
    with h5py.File(hdf5_path, 'r') as f:
        R = float(f['scalars/xmag'][timeslice])
        Z = float(f['scalars/zmag'][timeslice])
    return np.array([R, Z])


def print_lobe_table(lobes: np.ndarray, label: str) -> None:
    print(f"\n{label} — {lobes.shape[0]} lobe(s) found")
    if lobes.shape[0] == 0:
        return
    header = f"  {'#':>3}  {'R_mid':>8}  {'Z_mid':>8}  {'angle':>8}  {'perim':>8}  {'area':>10}  {'h':>8}"
    units  = f"  {'':>3}  {'(m)':>8}  {'(m)':>8}  {'(rad)':>8}  {'(m)':>8}  {'(m^2)':>10}  {'(m)':>8}"
    print(header)
    print(units)
    print("  " + "-" * (len(header) - 2))
    for i, row in enumerate(lobes):
        print(f"  {i:>3}  {row[0]:8.4f}  {row[1]:8.4f}  {row[2]:8.4f}  "
              f"{row[3]:8.4f}  {row[4]:10.6f}  {row[5]:8.4f}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # ── 1. Load manifold data ─────────────────────────────────────────────────
    print(f"Loading stable manifold   from {STABLE_NPZ}")
    stable_full   = load_manifold(STABLE_NPZ)
    print(f"  {stable_full.shape[0]} points")

    print(f"Loading unstable manifold from {UNSTABLE_NPZ}")
    unstable_full = load_manifold(UNSTABLE_NPZ)
    print(f"  {unstable_full.shape[0]} points")

    # ── 2. Equilibrium separatrix ─────────────────────────────────────────────
    # Pre-computed by trace_separatrix.py (mfgen tracing, timeslice=-1).
    #
    # Alternative for datasets where the HDF5 flux surface extracts cleanly:
    #   separatrix_full = maglib.get_separatrix(HDF5_PATH, psi_n_level=0.999)
    print(f"Loading equilibrium separatrix from {SEPARATRIX_NPZ}")
    separatrix_full = load_manifold(SEPARATRIX_NPZ)
    print(f"  {separatrix_full.shape[0]} points")

    # ── 3. Wall (for plotting context only) ───────────────────────────────────
    wall = load_wall(WALL_PATH)

    # ── 4. Simplify ───────────────────────────────────────────────────────────
    print(f"\nSimplifying with RDP tolerance = {SIMPLIFY_TOL*1e3:.1f} mm ...")
    stable_simp, unstable_simp, sep_simp = maglib.simplify(
        [stable_full, unstable_full, separatrix_full], SIMPLIFY_TOL)
    print(f"  stable:      {stable_full.shape[0]} → {stable_simp.shape[0]} points")
    print(f"  unstable:    {unstable_full.shape[0]} → {unstable_simp.shape[0]} points")
    print(f"  separatrix:  {separatrix_full.shape[0]} → {sep_simp.shape[0]} points")

    # ── 5. Magnetic reference points ──────────────────────────────────────────
    mag_axis = read_mag_axis(HDF5_PATH, timeslice=1)
    x_point  = np.array([R_XPOINT, Z_XPOINT])
    print(f"\nMagnetic axis: R = {mag_axis[0]:.5f} m,  Z = {mag_axis[1]:.5f} m")
    print(f"X-point:       R = {x_point[0]:.5f} m,  Z = {x_point[1]:.5f} m")

    # ── 6. Lobe analysis ──────────────────────────────────────────────────────
    # Equilibrium separatrix is the reference boundary; each perturbed manifold
    # is the curve whose intersections with the separatrix bound the lobes.
    print("\nRunning lobe_map (stable manifold vs equilibrium separatrix) ...")
    stable_lobes = maglib.lobe_map(sep_simp, stable_simp, mag_axis, x_point)

    print("Running lobe_map (unstable manifold vs equilibrium separatrix) ...")
    unstable_lobes = maglib.lobe_map(sep_simp, unstable_simp, mag_axis, x_point)

    if stable_lobes.shape[0] == 0 and unstable_lobes.shape[0] == 0:
        print(
            "\nNo lobes found — the manifolds did not intersect the equilibrium separatrix.\n"
            "Increase N_SEGMENTS in examples/mfgen/run.py and re-run."
        )

    print_lobe_table(stable_lobes,   "Stable manifold lobes")
    print_lobe_table(unstable_lobes, "Unstable manifold lobes")

    # ── 7. Plot ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5, 8), dpi=90)
    ax.set_aspect('equal', adjustable='box')

    ax.plot(wall[:, 0],               wall[:, 1],               color='black',
            linewidth=1.0, label='wall')
    ax.plot(separatrix_full[:, 0],    separatrix_full[:, 1],    color='grey',
            linewidth=1.5, linestyle='--', label='separatrix (eq.)')
    ax.plot(stable_full[:, 0],   stable_full[:, 1],   color='steelblue',
            linewidth=1.2, label='stable')
    ax.plot(unstable_full[:, 0], unstable_full[:, 1], color='firebrick',
            linewidth=1.2, label='unstable')

    if stable_lobes.shape[0] > 0:
        ax.scatter(stable_lobes[:, 0], stable_lobes[:, 1],
                   color='steelblue', edgecolors='k', s=40, zorder=5,
                   label='stable lobe centroids')
        for i, row in enumerate(stable_lobes):
            ax.annotate(f"s{i}", xy=(row[0], row[1]),
                        xytext=(4, 4), textcoords='offset points', fontsize=7)
        ix = np.unique(stable_lobes[:, 6:10].reshape(-1, 2), axis=0)
        ax.scatter(ix[:, 0], ix[:, 1],
                   color='steelblue', marker='D', s=25, zorder=6,
                   label='stable intersections')

    if unstable_lobes.shape[0] > 0:
        ax.scatter(unstable_lobes[:, 0], unstable_lobes[:, 1],
                   color='firebrick', edgecolors='k', s=40, zorder=5,
                   label='unstable lobe centroids')
        for i, row in enumerate(unstable_lobes):
            ax.annotate(f"u{i}", xy=(row[0], row[1]),
                        xytext=(4, 4), textcoords='offset points', fontsize=7)
        ix = np.unique(unstable_lobes[:, 6:10].reshape(-1, 2), axis=0)
        ax.scatter(ix[:, 0], ix[:, 1],
                   color='firebrick', marker='D', s=25, zorder=6,
                   label='unstable intersections')

    ax.scatter(*x_point, marker='x', color='black', s=60, zorder=6, label='X-point')

    ax.set_xlabel("R (m)", fontsize=11)
    ax.set_ylabel("Z (m)", fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
