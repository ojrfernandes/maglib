# Maglib
Maglib is a library designed for computational tasks in plasma physics, leveraging the Fusion-IO framework for data extraction. It provides tools for generating high resolution figures of magnetic footprints and manifolds. Maglib is a work in progess and currently has only been tested on data from M3D-C1 simulations for the TCABR Tokamak. Maglib is being developped at the Plasma Physics Laboratory from Universidade de São Paulo - Brazil.

## Prerequisites
### Fusion-IO Dependency
maglib requires an installation of Fusion-IO, an interface for extracting data from various plasma simulation codes. Fusion-IO was developed by Nate Ferraro at the Princeton Plasma Physics Laboratory and can be found at the [Fusion-IO GitHub Repository](https://github.com/nferraro/fusion-io).

⚠️ Compatibility Warning:
This program **is not** compatible with the latest versions of Fusion-IO. For the best results, we recommend using a Fusion-IO version from around July 2023. Commit 022a77f is our recommentation.

To clone the compatible version, run:
```bash 
git clone https://github.com/nferraro/fusion-io.git
cd fusion-io
git checkout 022a77f
```
Keep in mind that Fusion-IO CMake install has been implemented prior to this version, hence, manual instalation is necessary. 

## Sode 
The Sode module (Solver for Ordinary Differential Equations) provides a general purpose integrator for systems of ordinary differential equations, including various adaptative Runge-Kutta methods.

## Maglit
The Maglit module (Magnetic Line Integrator) provides an interface between Sode and the Fusion-IO modules, making it possible to handle systems which require information from the fields at every integration step e.g. magnetic field lines.
