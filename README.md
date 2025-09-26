# Maglib

Maglib is a library designed for computational tasks in plasma physics, leveraging the Fusion-IO framework for data extraction. It provides tools for generating high-resolution figures of magnetic footprints and invariant manifolds. 

**Current Status**: Maglib is a work in progress and has currently only been tested on data from M3D-C1 simulations for the TCABR Tokamak. 

Developed at the Plasma Physics Laboratory, Institute of Physics - Universidade de São Paulo, Brazil.

## Prerequisites

Before installing Maglib, ensure you have the following dependencies:

### System Requirements
- **CMake** (≥ 3.10): Build system generator  
  📄 [CMake Documentation](https://cmake.org/)

- **HDF5**: High-performance data management library (with C headers)  
  📄 [HDF5 Documentation](https://www.hdfgroup.org/)

- **OpenMP**: API for parallel programming  
  📄 [OpenMP Documentation](https://www.openmp.org/)

- **Armadillo**: C++ library for linear algebra & scientific computing  
  📄 [Armadillo Documentation](https://arma.sourceforge.net/)

### Fusion-IO Dependency
Maglib requires **Fusion-IO**, an interface for extracting data from plasma simulation codes developed by Nate Ferraro at Princeton Plasma Physics Laboratory.

📥 [Fusion-IO GitHub Repository](https://github.com/nferraro/fusion-io)

⚠️ **Important Compatibility Note**: Maglib expects Fusion-IO libraries to be compiled as two shared libraries named `libfusionio.so` and `libm3dc1.so`. When installing Fusion-IO using CMake, these may be compiled as static libraries or with different names. Please ensure the correct shared library names before proceeding with installation.

## Installation

### Step 1: Set Environment Variable
Add your Fusion-IO build directory to your environment variables:

```bash
export FUSION_IO_DIR=/path/to/your/fusion-io/build
```
or add this line to your `~/.bashrc` or `~/.zshrc` for persistence across sessions.

### Step 2: Clone and Build
```bash
# Clone the repository (if not already done)
git clone https://github.com/ojrfernandes/maglib.git maglib
cd maglib

# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make
```

### Step 3: (optional) Set the executables to your system PATH in your .bashrc
```bash
# set maglib directory
export MAGLIB_ROOT=/path/to/your/maglib/root
 
# add executables to path
export PATH=$MAGLIB_ROOT/build/bin:$PATH
```

## Troubleshooting

**Common Issues:**

- **CMake can't find dependencies**: Ensure all prerequisites are properly installed and in your system PATH
- **Fusion-IO linking errors**: Verify that `libfusionio.so` and `libm3dc1.so` exist in your Fusion-IO build/lib directory
- **Build errors**: Check that all dependencies have development headers installed (especially HDF5)