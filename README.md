# maglib

C++ library and Python package for magnetic field-line integration and tokamak plasma analysis. Maglib provides:

- **Field-line tracing**
- **Magnetic footprints generation**
- **Invariant manifolds tracing**
- **Plotting and post-processing utilities**

Developed at the Plasma Physics Laboratory, Institute of Physics — University of São Paulo, Brazil.

> **Compatibility note:** maglib has been tested on M3D-C1 simulation data for the TCABR tokamak. Other Fusion-IO-compatible codes / experiments should work but have not been validated.

---

## Prerequisites

### System
- **CMake** ≥ 3.10
- **C++17** compiler (GCC or Clang)
- **HDF5** (C headers and libraries)
- **OpenMP**

### Fusion-IO
maglib requires [Fusion-IO](https://github.com/nferraro/fusion-io), developed by Nate Ferraro at PPPL.

> **Important:** maglib expects two shared libraries named `libfusionio.so` and `libm3dc1.so` in the Fusion-IO build directory. These are the default library names when building Fusion-IO with Makefiles. When building Fusion-IO with CMake, verify these names before proceeding.

### Python (optional, for bindings and post-processing)
- Python ≥ 3.8 with development headers
- **numpy**
- **matplotlib**
- **shapely**

pybind11 is fetched automatically by CMake — no manual installation needed.

---

## Installation

### 1. Set the Fusion-IO path

```bash
export FUSION_IO_DIR=/path/to/fusion-io/build
```

Add this to `~/.bashrc` or `~/.zshrc` to persist across sessions.

### 2. Build

```bash
git clone https://github.com/ojrfernandes/maglib.git
cd maglib
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The Python package is built automatically if Python development headers are found. To disable it:

```bash
cmake -B build -DBUILD_PYTHON_BINDINGS=OFF
```

### 3. Set up the Python package

```bash
export PYTHONPATH=/path/to/maglib/python:$PYTHONPATH

pip install numpy matplotlib shapely
```

### 4. (Optional) Add C++ binaries to PATH

```bash
export PATH=/path/to/maglib/build/bin:$PATH
```

---

## Python API

### Field source and tracer

```python
import maglib

# Load field data (timeslice: -1 = equilibrium, 0 = vacuum, 1 = plasma response)
source = maglib.M3DC1Source("C1.h5", 1)

# Create and configure a field-line integrator
tracer = maglib.Maglit(source)
tracer.configure(dphi_init=1e-2, dphi_min=1e-6, dphi_max=1e-2)
tracer.set_monitor("first_wall.txt")   # optional wall-collision monitor
```

### Magnetic footprint

```python
fp = maglib.Footprint(source, R1_min, Z1_min, R1_max, Z1_max,
                      nRZ=100, nPhi=8, max_turns=200)

# Serial run
fp.run([tracer])

# Parallel run — one independent source+tracer pair per OMP thread
# (Fusion-IO is not thread-safe; sources must not be shared across threads)
src2 = maglib.M3DC1Source("C1.h5", 1)
tracer2 = maglib.Maglit(src2)
tracer2.configure(1e-2, 1e-6, 1e-2)
fp.run([tracer, tracer2])

data = fp.output_data   # numpy (N, 6): R0, Z0, phi0, length, psiMin, turn
fp.save("footprint.npz")
maglib.plot_footprint(fp)
```

### Invariant manifolds

```python
# stability: 0 = stable manifold (forward map), 1 = unstable manifold (backward map)
mf = maglib.Manifold(tracer, phi=0.0, stability=0)
mf.configure(epsilon=1e-8, h=1e-8, tol=1e-14,
             max_iter=50, precision_limit=1e-14, max_insertions=50)
mf.find_x_point(r_guess=0.498, z_guess=-0.219)

# One-shot
mf.run(n_intervals=9, n_segments=10, method=1, l_lim=0.005, theta_lim=20.0)

# Step-by-step (interpolant method)
seg = mf.primary_segment(n_intervals=9)
for _ in range(9):
    _, seg = mf.new_segment(seg, l_lim=0.005, theta_lim=20.0)

data = mf.output_data   # list of (N_i, 2) arrays, columns [R, Z]
mf.save("manifold.npz")
maglib.plot_manifold(mf)
```

### Equilibrium separatrix

Two methods are available. Method A (HDF5 contour) is faster but may fail in some scenarios:

```python
# Method A — extract ψ_N contour from the equilibrium mesh (no field evaluation)
sep = maglib.get_separatrix("C1.h5", psi_n_level=0.999)  # shape (N, 2)
```

Method B (mfgen tracing) generally robust for all datasets:

```python
# Method B — trace the separatrix by iterating the Poincaré map at timeslice=-1
eq_source = maglib.M3DC1Source("C1.h5", -1)
eq_tracer = maglib.Maglit(eq_source)
eq_tracer.configure(1e-2, 1e-6, 1e-2)
eq_tracer.set_monitor("first_wall.txt")

sep, x_point = maglib.trace_separatrix(
    eq_tracer, phi=0.0,
    r_xpoint=0.498, z_xpoint=-0.219,
    n_segments=6,   # increase until the curve visually closes
    verbose=True,
)
# sep: shape (N, 2), x_point: shape (2,)
```

### Lobe analysis

Lobes are regions bounded between consecutive intersection points of a perturbed manifold with the equilibrium separatrix:

```python
import numpy as np

# Simplify curves before intersection search (Ramer-Douglas-Peucker)
sep_s, stable_s, unstable_s = maglib.simplify([sep, stable, unstable], tolerance=5e-4)

mag_axis = np.array([R_axis, Z_axis])
x_point  = np.array([R_xp,   Z_xp])

lobes = maglib.lobe_map(sep_s, stable_s, mag_axis, x_point)
# shape (N_lobes, 10):
#   Rmid, Zmid        — arc midpoint on the separatrix between the two intersections
#   angle [rad]       — poloidal angle relative to the X-point
#   perimeter [m]     — lobe boundary length
#   area [m²]         — lobe area
#   h_param [m]       — effective width (2·area / base)
#   R1, Z1, R2, Z2   — bounding intersection points
```

### Manifold post-processing utilities

```python
# Simplify segments (RDP, no interpolation)
simplified = maglib.simplify(source, tolerance)   # source: Manifold or list[ndarray]

# Convert to shapely LineString objects
lines = maglib.to_linestrings(source)
```

### Interactive crop tool

Trim a manifold file to a chosen number of points with live replot:

```bash
python -m maglib.crop_manifold stable.npz 800 [--wall first_wall.txt]
```

Enter a new number to update the cut, `ok` to overwrite the file (a `.bak` backup is created automatically), `q` to quit without saving. Supports `.npz`, `.dat`, and `.txt` formats.

---

## C++ binaries

`fpgen` and `mfgen` are thin CLI wrappers over the same library code. Each reads a plain-text input file from the working directory:

```bash
fpgen   # reads fpgen_input.txt
mfgen   # reads mfgen_input.txt
```

See `fpgen/fpgen_input.txt` and `mfgen/mfgen_input.txt` for documented input formats.

---

## Examples

End-to-end workflow scripts are in `examples/`:

| Script | Description |
|---|---|
| `sode/run.py` | ODE integrator usage |
| `fpgen/run.py` | Magnetic footprint computation |
| `fpgen/plot.py` | Footprint plotting |
| `mfgen/run.py` | Invariant manifold computation (stable + unstable) |
| `manifold_analysis/trace_separatrix.py` | Equilibrium separatrix tracing via `maglib.trace_separatrix()` |
| `manifold_analysis/lobe_analysis.py` | Lobe geometry from manifold–separatrix intersections |

Run any example with:

```bash
PYTHONPATH=/path/to/maglib/python python examples/<script>.py
```

---

## Testing

### Test data

C++ and Python tests require sample M3D-C1 data hosted on Figshare:

📦 https://doi.org/10.6084/m9.figshare.30593324.v1

Download and extract into `tests/data/`:

```
tests/
└── data/
    ├── C1.h5
    └── tcabr_first_wall.txt
```

### Run tests

```bash
# C++ tests (requires build)
ctest --test-dir build

# Python tests
pytest tests/python/
```

GoogleTest is fetched automatically by CMake. pytest must be installed: `pip install pytest`.

---

## Repository layout

```
maglib/
├── maglit/          Core integrator (sode, collider, maglit, M3DC1Source, FieldSource)
├── fpgen/           Footprint generator CLI (footprint.h/cpp, input_read, run.cpp)
├── mfgen/           Manifold generator CLI (manifold.h/cpp, input_read, run.cpp)
├── python/
│   ├── src/         pybind11 binding sources (one file per component)
│   └── maglib/      Python package
│       ├── __init__.py
│       ├── plot.py             — plot_footprint, plot_manifold
│       ├── manifold_tools.py  — simplify, to_linestrings, lobe_map,
│       │                        get_separatrix, trace_separatrix
│       └── crop_manifold.py   — interactive manifold crop CLI
├── examples/
│   ├── sode/
│   ├── fpgen/
│   ├── mfgen/
│   └── manifold_analysis/
└── tests/
    ├── data/        Test dataset
    ├── *.cpp        C++ tests (GoogleTest)
    └── python/      Python tests (pytest)
```

---

## Troubleshooting

- **CMake can't find HDF5** — install HDF5 development headers (`libhdf5-dev` on Debian/Ubuntu).
- **Fusion-IO linking errors** — verify `libfusionio.so` and `libm3dc1.so` exist in `$FUSION_IO_DIR/lib`.
- **`import maglib` fails** — check that `PYTHONPATH` includes the `python/` directory and that the build completed with Python bindings enabled.
- **`lobe_map` / `simplify` ImportError** — install shapely: `pip install shapely`.
