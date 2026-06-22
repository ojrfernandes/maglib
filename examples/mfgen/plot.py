"""
Manifold plotting — demonstrates all input modes of maglib.plot_manifold().

Usage
-----
    # First generate the output files:
    PYTHONPATH=<maglib_root>/python python examples/mfgen/run.py

    # Then plot:
    PYTHONPATH=<maglib_root>/python python examples/mfgen/plot.py

The primary goal is to visualise the full stable and unstable manifolds
together as two distinct structures, each in a single colour.

plot_manifold() accepts four source types:

  1. List of .dat paths  — one complete curve per file; primary workflow
  2. Single .dat path    — one complete curve; quick single-manifold look
  3. List of arrays      — loaded from .npz; same result as (1) after stacking
  4. Manifold object     — only available immediately after run.py; same as (1)

For types (1), (2), (3): the .dat / .npz files produced by run.py are enough.
"""

from pathlib import Path
import numpy as np
import maglib

WALL_PATH = str(Path(__file__).resolve().parents[2] / "tests" / "data" / "tcabr_first_wall.txt")


# ── Case 1: list of .dat paths (primary workflow) ─────────────────────────────
# Each file is the full manifold (all segments concatenated), plotted as one
# curve per file. This is the standard way to overlay stable and unstable.

def plot_both_from_files():
    print("Case 1: stable + unstable from .dat files")
    maglib.plot_manifold(
        ["stable.dat", "unstable.dat"],
        wall=WALL_PATH,
        labels=["stable", "unstable"],
        figsize=(6, 8),
    )


# ── Case 2: single .dat path ──────────────────────────────────────────────────
# Quick look at one manifold without loading anything into Python.

def plot_single_file():
    print("Case 2: single .dat path")
    maglib.plot_manifold("stable.dat", wall=WALL_PATH,
                         labels=["stable"],
                         figsize=(6, 8))


# ── Case 3: list of arrays loaded from .npz ───────────────────────────────────
# .npz stores one array per segment (keys seg_0, seg_1, ...). Stack all
# segments of each manifold into a single (N, 2) array to reproduce the same
# single-colour-per-manifold result as case 1.

def plot_both_from_npz():
    print("Case 3: stable + unstable from .npz files")

    def load_full(path):
        data = np.load(path)
        segs = [data[k] for k in sorted(data.files)]
        return np.vstack(segs)

    stable_arr   = load_full("stable.npz")
    unstable_arr = load_full("unstable.npz")

    maglib.plot_manifold(
        [stable_arr, unstable_arr],
        wall=WALL_PATH,
        labels=["stable", "unstable"],
        figsize=(6, 8),
    )


# ── Case 4: Manifold object (requires a live computation) ─────────────────────
# Only available immediately after growing the manifold in the same session.
# All segments share one colour and one legend label — same visual result as
# the cases above.
#
#   from run import build_source, build_tracer, grow_manifold
#   source   = build_source()
#   unstable = grow_manifold(build_tracer(source), stability=0)
#   stable   = grow_manifold(build_tracer(source), stability=1)
#   maglib.plot_manifold([stable, unstable],
#                        wall=WALL_PATH, labels=["stable", "unstable"])
#
# Note: plot_manifold accepts a list of Manifold objects, one entry per colour.


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    plot_both_from_files()   # primary workflow
    plot_single_file()
    plot_both_from_npz()


if __name__ == "__main__":
    main()
