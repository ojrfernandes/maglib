"""
crop_manifold.py — interactive tool to crop a manifold file to a chosen number of points.

Usage
-----
    python -m maglib.crop_manifold <file> [N] [--wall PATH]

    <file>    .npz, .dat, or .txt manifold file
    N         initial number of points to display (default: all)
    --wall    first-wall contour file (default: tests/data/tcabr_first_wall.txt)

Workflow
--------
    1. The file is loaded and stacked into a single (N_total, 2) curve.
    2. The first N points are plotted (blue) alongside the discarded tail (grey)
       and the first wall.
    3. Enter a new number to update the crop and replot.
    4. Type 'ok' to overwrite the file with the retained points, 'q' to quit
       without saving.

File formats
------------
    .npz    — numpy compressed archive with keys seg_0, seg_1, ...
              Saved back as single-key archive (key "seg_0").
    .dat / .txt — space-separated text, 3 columns: seg R Z
              Saved back with header and segment index 0 for all rows.
"""

import sys
import argparse
import shutil
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ── Default paths ─────────────────────────────────────────────────────────────

_REPO     = Path(__file__).resolve().parents[2]
WALL_PATH = str(_REPO / "tests" / "data" / "tcabr_first_wall.txt")


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_file(path: Path) -> np.ndarray:
    """Load a manifold file; return finite-only stacked (N, 2) [R, Z] array."""
    suffix = path.suffix.lower()
    if suffix == ".npz":
        data = np.load(path)
        arr = np.vstack([data[k] for k in sorted(data.files)])
    elif suffix in (".dat", ".txt"):
        raw = np.loadtxt(path, comments="#")
        if raw.ndim == 1:
            raw = raw[np.newaxis, :]
        arr = raw[:, 1:3]            # columns: seg R Z  →  take R, Z
    else:
        raise ValueError(f"Unsupported extension '{suffix}'. Use .npz, .dat, or .txt.")
    valid = np.all(np.isfinite(arr), axis=1)
    dropped = (~valid).sum()
    if dropped:
        print(f"  Dropped {dropped} non-finite rows.")
    return arr[valid]


def save_file(path: Path, arr: np.ndarray) -> None:
    """Write (N, 2) array back to file, preserving original format."""
    suffix = path.suffix.lower()
    if suffix == ".npz":
        np.savez_compressed(path, seg_0=arr)
    elif suffix in (".dat", ".txt"):
        with open(path, "w") as f:
            f.write("# seg                  R                   Z\n")
            for row in arr:
                f.write(f"0  {row[0]:.16f}  {row[1]:.16f}\n")
    else:
        raise ValueError(f"Unsupported extension '{suffix}'.")


def load_wall(path: str) -> np.ndarray | None:
    try:
        wall = np.loadtxt(path)
        if not np.allclose(wall[0], wall[-1]):
            wall = np.vstack([wall, wall[0]])
        return wall
    except Exception as e:
        print(f"Warning: could not load wall from {path} ({e}).")
        return None


# ── Interactive plot ──────────────────────────────────────────────────────────

def draw(ax, wall: np.ndarray | None, arr: np.ndarray, n: int) -> None:
    ax.clear()
    ax.set_aspect("equal", adjustable="box")

    if wall is not None:
        ax.plot(wall[:, 0], wall[:, 1], color="black", linewidth=1.0, label="wall")

    if n < len(arr):
        ax.plot(arr[n:, 0], arr[n:, 1],
                color="lightgrey", linewidth=0.8, zorder=1,
                label=f"discarded ({len(arr)-n} pts)")

    ax.plot(arr[:n, 0], arr[:n, 1],
            color="steelblue", linewidth=1.2, zorder=2, label=f"kept ({n} pts)")

    ax.scatter(arr[n - 1, 0], arr[n - 1, 1],
               color="crimson", s=30, zorder=5, label="cut point")

    pct = 100.0 * n / len(arr)
    ax.set_title(f"{n} / {len(arr)} points  ({pct:.1f}%)", fontsize=10)
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.legend(fontsize=8, loc="upper right")
    ax.figure.tight_layout()
    ax.figure.canvas.draw()
    ax.figure.canvas.flush_events()


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interactive manifold cropping tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("file", help="Manifold file (.npz, .dat, .txt)")
    parser.add_argument("N", nargs="?", type=int, default=None,
                        help="Initial number of points to display (default: all)")
    parser.add_argument("--wall", default=WALL_PATH,
                        help="First-wall contour file (default: tests/data/tcabr_first_wall.txt)")
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        sys.exit(f"File not found: {path}")

    print(f"Loading {path} ...")
    arr   = load_file(path)
    total = len(arr)
    print(f"  {total} finite points.")

    wall = load_wall(args.wall)

    n = args.N if args.N is not None else total
    n = max(1, min(n, total))

    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass

    plt.ion()
    fig, ax = plt.subplots(figsize=(6, 8), dpi=90)
    draw(ax, wall, arr, n)
    plt.show(block=False)
    plt.pause(0.05)

    print("\nCommands:")
    print("  <number>   update the crop and replot")
    print("  ok / save  overwrite the file with the current crop and exit")
    print("  q / quit   exit without saving")
    print()

    while True:
        try:
            raw = input(f"  Points [{n}/{total}]: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nAborted (no changes saved).")
            break

        if not raw:
            continue

        if raw.lower() in ("ok", "save", "yes", "y"):
            backup = path.with_suffix(path.suffix + ".bak")
            shutil.copy2(path, backup)
            save_file(path, arr[:n])
            print(f"Saved {n} points to {path}  (backup: {backup.name})")
            break

        if raw.lower() in ("q", "quit", "exit", "n", "no"):
            print("Exited without saving.")
            break

        try:
            new_n = int(raw)
        except ValueError:
            print("  Not understood. Enter a number, 'ok' to save, or 'q' to quit.")
            continue

        if not (1 <= new_n <= total):
            print(f"  Enter a value between 1 and {total}.")
            continue

        n = new_n
        draw(ax, wall, arr, n)
        plt.pause(0.05)

    plt.ioff()
    plt.close("all")


if __name__ == "__main__":
    main()
