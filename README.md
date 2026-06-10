# maglib

C++ library and Python package for magnetic field-line integration and tokamak plasmas analysis. Built on top of [Fusion-IO](https://github.com/nferraro/fusion-io) for field data, maglib provides:

- **Field-line tracing**
- **Magnetic footprints** 
- **Invariant manifolds**
- **Magnetic Lobe analysis**
- **Plotting**

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

> **Important:** maglib expects two shared libraries named `libfusionio.so` and `libm3dc1.so` in the Fusion-IO build directory. These are the default library names when "manually" building Fusion-IO with Makefiles. When building Fusion-IO with CMake, verify these names before proceeding.

### Python (optional, for bindings)
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

The Python package is built automatically if Python development headers are found. To disable it explicitly:

```bash
cmake -B build -DBUILD_PYTHON_BINDINGS=OFF
```

### 3. Set up the Python package

```bash
# Make the package importable
export PYTHONPATH=/path/to/maglib/python:$PYTHONPATH

# Install Python dependencies
pip install numpy
pip install matplotlib
pip install shapely
```

### 4. (Optional) Add C++ binaries to PATH

```bash
export PATH=/path/to/maglib/build/bin:$PATH
```

---

## Usage

### Python API

```python
import maglib

# Load field data
src = maglib.M3DC1Source("C1.h5", 0.0)
src.open()

# Create a tracer
tracer = maglib.Maglit(src, r=1.12, z=0.0, phi=0.0)

# Magnetic footprint — serial
fp = maglib.Footprint(tracer, R1_min, Z1_min, R1_max, Z1_max,
                      nRZ=100, nPhi=8, max_turns=200)
fp.run([tracer])
data = fp.output_data   # numpy (N, 6): R0, Z0, phi0, length, psiMin, turn
fp.save("footprint.npz")
maglib.plot_footprint(fp)

# Magnetic footprint — parallel (one source+tracer per thread; Fusion-IO is not thread-safe)
src2 = maglib.M3DC1Source("C1.h5", 0.0); src2.open()
tracer2 = maglib.Maglit(src2, r=1.12, z=0.0, phi=0.0)
fp.run([tracer, tracer2])   # 2 OpenMP threads; list length sets the thread count

# Invariant manifold — single Poincaré section
mf = maglib.Manifold(tracer, phi=0.0, stability=0)   # stability: 0 = unstable, 1 = stable
mf.configure(epsilon=1e-4, h=1e-3, tol=1e-8,
             max_iter=1000, precision_limit=1e-12, max_insertions=500)
mf.find_x_point(r_guess=0.50, z_guess=-0.22)
mf.run(n_intervals=10, n_segments=5,
       method=0, l_lim=0.5, theta_lim=0.3)
data = mf.output_data   # list of (N_i, 2) arrays [R, Z] per segment
mf.save("manifold.npz")
maglib.plot_manifold(mf)

# Invariant manifold — multiple Poincaré sections (equivalent to mfgen nSections)
import numpy as np
for phi_deg in np.linspace(0, 90, 4):
    mf = maglib.Manifold(tracer, phi=np.deg2rad(phi_deg), stability=0)
    mf.configure(epsilon=1e-4, h=1e-3, tol=1e-8,
                 max_iter=1000, precision_limit=1e-12, max_insertions=500)
    mf.find_x_point(r_guess=0.50, z_guess=-0.22)
    mf.run(n_intervals=10, n_segments=5, method=0, l_lim=0.5, theta_lim=0.3)
    mf.save(f"manifold_{int(phi_deg)}.npz")
# For parallel execution across sections use threading.Thread with one source+tracer per thread.

# Lobe analysis (requires shapely)
lobes = maglib.lobe_map(equilibrium, perturbed, mag_axis)
# lobes: numpy (N_lobes, 6) — Rmid, Zmid, angle [rad], perimeter, area, h_param
```

### C++ binaries

`fpgen` and `mfgen` are thin CLI wrappers over the same library code.
Each reads a plain-text input file:

```bash
fpgen   # reads fpgen_input.txt  in the working directory
mfgen   # reads mfgen_input.txt  in the working directory
```

See `fpgen/fpgen_input.txt` and `mfgen/mfgen_input.txt` for documented input formats.

---

## Testing

### Test data

C++ and Python tests require sample M3D-C1 data hosted on Figshare:

📦 https://doi.org/10.6084/m9.figshare.30593324.v1

Download and extract into `tests/data/` so the layout is:

```
tests/
└── data/
    ├── C1.h5
    └── ...
```

### Test dependencies

GoogleTest (C++) and pytest (Python) are the only additional dependencies needed to run the test suite.

- **GoogleTest** is fetched automatically by CMake — no manual installation needed.
- **pytest** must be installed manually: `pip install pytest`

### Run tests

```bash
# C++ tests
ctest --test-dir build

# Python tests
pytest tests/python/
```

---

## Repository layout

```
maglib/
├── maglit/          Core integrator (sode, collider, maglit, M3DC1Source, FieldSource)
├── fpgen/           Footprint generator CLI (footprint.h/cpp, input_read, run.cpp)
├── mfgen/           Manifold generator CLI (manifold.h/cpp, input_read, run.cpp)
├── python/
│   ├── src/         pybind11 binding sources (one file per component)
│   └── maglib/      Python package (__init__.py, plot.py, manifold_tools.py)
└── tests/
    ├── data/        Test dataset
    ├── *.cpp        C++ tests (GoogleTest)
    └── python/      Python tests (pytest)
```

---

## Troubleshooting

- **CMake can't find HDF5** — ensure HDF5 development headers are installed (`libhdf5-dev` on Debian/Ubuntu).
- **Fusion-IO linking errors** — verify `libfusionio.so` and `libm3dc1.so` exist in `$FUSION_IO_DIR/lib`.
- **`import maglib` fails** — check that `PYTHONPATH` includes the `python/` directory and that the build completed with Python bindings enabled.
- **`lobe_map` / `simplify` ImportError** — install shapely: `pip install shapely`.
