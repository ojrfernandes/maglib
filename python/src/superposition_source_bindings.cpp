#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "field_source.h"
#include "superposition_source.h"

namespace py = pybind11;
using namespace py::literals;

void bind_superposition_source(py::module_ &m) {
    py::class_<SuperpositionSource, FieldSource>(m, "SuperpositionSource", R"doc(
FieldSource that linearly superposes N M3DC1 perturbation components onto a
shared axisymmetric equilibrium:

    B(R, φ, Z) = B_eq(R, φ, Z) + Σ_i A_i · [B_i(R, φ − δ_i, Z) − B_eq(R, φ, Z)]

The phase shift δ_i is in **radians** and is applied as φ → φ − δ_i when
evaluating component i. The equilibrium is loaded automatically (timeslice = -1)
from the first component's file.

Example::

    import maglib

    source = maglib.SuperpositionSource()
    source.add_component("/path/IM_C1.h5", timeslice=1, phase_shift=  0.0, amplitude=1.0)
    source.add_component("/path/IL_C1.h5", timeslice=1, phase_shift=-100.0, amplitude=1.0)
    source.add_component("/path/IU_C1.h5", timeslice=1, phase_shift= +80.0, amplitude=1.0)

    if not source.is_valid():
        raise RuntimeError("failed to open one or more field sources")

    tracer = maglib.Maglit(source)
)doc")
        .def(py::init<>())

        .def("add_component",
             [](SuperpositionSource &self, const std::string &path, int timeslice,
                double phase_shift_deg, double amplitude) {
                 self.add_component(path, timeslice, phase_shift_deg * M_PI / 180.0, amplitude);
             },
             "path"_a, "timeslice"_a, "phase_shift"_a, "amplitude"_a,
             R"doc(Add a perturbation component.

Parameters
----------
path : str
    Path to the M3DC1 HDF5 file for this component.
timeslice : int
    0 = vacuum, 1 = full single-fluid response. The equilibrium (timeslice = -1)
    is opened automatically from the first component's file.
phase_shift : float
    Toroidal phase shift δ_i in **degrees**. Applied as φ → φ − δ_i.
amplitude : float
    Linear scale factor A_i. May be negative for anti-phase contributions.
)doc")

        .def("is_valid", &SuperpositionSource::is_valid,
             "Return True if all component sources were opened successfully.")

        .def("num_components", &SuperpositionSource::num_components,
             "Return the number of perturbation components.");
}
