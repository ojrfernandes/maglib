# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
# First-time setup: set Fusion-IO build directory
export FUSION_IO_DIR=/path/to/your/fusion-io/build

# Configure and build (from repo root)
mkdir -p build && cd build
cmake ..
make
```

Executables are placed in `build/bin/` (fpgen, mfgen, lbmap). The `maglit` library is built as a static library at `build/maglit/libmaglit.a`.

### Python bindings

The Python extension (`_maglib.so`) is built automatically if Python development headers are found. It is placed in `python/maglib/` (alongside `__init__.py`) so importing works by adding a single directory to `PYTHONPATH`:

```bash
export PYTHONPATH=/path/to/maglib/python:$PYTHONPATH
python -c "import maglib; print(maglib.Sode)"
```

To disable: `cmake .. -DBUILD_PYTHON_BINDINGS=OFF`

pybind11 is fetched automatically from GitHub the first time CMake runs (requires internet access). Subsequent runs use the cached download.

## Tests

Tests require sample M3D-C1 data in `tests/data/` (download from https://doi.org/10.6084/m9.figshare.30593324.v1).

```bash
cd build

# Run full C++ test suite (CTest)
make test

# Run individual C++ test binaries
./tests/test_sode
./tests/test_collider
./tests/test_maglit
./tests/test_fpgen
./tests/test_mfgen
```

C++ tests use GoogleTest (fetched from GitHub automatically if not found on the system).

### Python tests

```bash
# From repo root (after building the Python extension)
export PYTHONPATH=/path/to/maglib/python:$PYTHONPATH
pytest tests/python/
```

Python tests live in `tests/python/` and use pytest. `conftest.py` adds the `python/` directory to `sys.path` automatically so no manual `PYTHONPATH` setup is needed when running pytest from the repo root.

## Architecture

This is a plasma physics library for studying magnetic field line topology in tokamaks, built around M3D-C1 simulation data via the Fusion-IO framework.

### Core library: `maglit/`

The `maglit` class is the central abstraction used by all tools. It wraps Fusion-IO to read magnetic field data from HDF5 files and integrates field lines in cylindrical coordinates (R, φ, Z).

- **`sode`** — adaptive ODE solver with RK56/RK78 Butcher tables (Fehlberg, Cash-Karp, Dormand-Prince). Solves the field-line ODE `dR/dφ, dZ/dφ = B_R/B_φ, B_Z/B_φ`.
- **`collider`** — loads a vessel wall polygon from a `.txt` file and tests whether a point (R, Z) is inside/outside using winding number. Used to detect when field lines strike the wall.
- **`maglit`** — owns a `sode` solver and a `collider` boundary. Exposes `step()` for single-turn integration and `psin_eval()` / `psi_eval()` for flux surface queries.

### Executables

Each executable reads parameters from a `*_input.txt` file in the working directory (not the binary's directory). The input file path is hardcoded as a relative path in each `run.cpp`.

- **`fpgen/`** — Magnetic Footprint Generator. Traces a 2D grid of field lines from a target plate segment (defined by two (R,Z) points) through `nPhi` toroidal angles. Uses OpenMP to parallelize over grid points, one `maglit` instance per thread. Outputs connection length, minimum ψ_N, and turn number per point.

- **`mfgen/`** — Manifold Generator. Computes stable/unstable invariant manifolds of the X-point by iteratively refining segments via the Poincaré map. Uses Armadillo for linear algebra (Jacobian eigendecomposition to find the unstable direction). The refinement strategy inserts new points when segment arc length or turning angle exceeds configurable thresholds.

- **`lbmap/`** — Lobe and Boundary Mapping Tool. Takes pre-computed manifold curves (from mfgen output), finds their intersections, and computes lobe geometry (area, perimeter, H-parameter, poloidal angle). Does not depend on `maglit` directly — operates purely on geometry from input files.

### Python visualization: `magplot/`

A Python package (`pip install -e magplot/`) with two plotting functions:
- `plot_fp()` — plots fpgen output as 2D heat maps (connection length, ψ_N, and a combined normalized metric)
- `plot_mf()` — plots mfgen manifold curves in the (R, Z) poloidal plane

### Input file format

All executables use the same key-value input format parsed by their respective `input_read` classes:
```
key = value   # inline comments supported with #
```

### External dependency: Fusion-IO

`FUSION_IO_DIR` must point to a Fusion-IO build containing `lib/libfusionio.so` and `lib/libm3dc1.so` (shared libraries — CMake builds may produce static libs with different names; rename as needed). Headers must be at `$FUSION_IO_DIR/include`.
