#include <memory>
#include <pybind11/pybind11.h>
#include "equilibrium_source.h"
#include "field_source.h"
#include "grid_perturbation_source.h"
#include "perturbed_source.h"

namespace py = pybind11;
using namespace py::literals;

void bind_equilibrium_sources(py::module_ &m) {
    py::class_<EquilibriumSource, FieldSource>(m, "EquilibriumSource", R"doc(
Axisymmetric equilibrium from a 'magnetics' XML file or a G-EQDSK file
(format detected from content).

B_R = -sigma dpsi/dZ / (2 pi R), B_Z = sigma dpsi/dR / (2 pi R), B_phi = F(psi_N) / R
inside the separatrix and F_vac / R elsewhere. eval_psi returns the total flux in Wb.

Example::

    src = maglib.EquilibriumSource("equilibrium.xml")
    if not src.is_valid():
        raise RuntimeError("failed to load equilibrium")
)doc")
        .def(py::init<const std::string &>(), "path"_a)
        .def("is_valid", &EquilibriumSource::is_valid,
             "Return True if the file was read and is consistent.");

    py::class_<GridPerturbationSource, FieldSource>(m, "GridPerturbationSource", R"doc(
Perturbation field on a uniform (phi, R, Z) grid from a restore_field() stream binary.

The stored toroidal range may be one field period 2 pi / N; evaluation wraps phi into it.
eval_B fails outside the R-Z box; eval_psi / eval_psin always fail (no flux in the file).
The file is memory-mapped, so several instances (one per thread) share one copy.
)doc")
        .def(py::init<const std::string &>(), "path"_a)
        .def("is_valid", &GridPerturbationSource::is_valid,
             "Return True if the file was read and its header is consistent.")
        .def("period", &GridPerturbationSource::period,
             "Toroidal period of the stored data in radians.");

    py::class_<PerturbedSource, FieldSource>(m, "PerturbedSource", R"doc(
Equilibrium plus a scaled, toroidally shifted gridded perturbation:

    B(R, phi, Z) = B_eq(R, phi, Z) + amplitude * B_pert(R, phi - phase_shift, Z)

phase_shift is in **radians**. psi and psi_N come from the equilibrium.

Example::

    src = maglib.PerturbedSource("equilibrium.xml", "perturbation.dat")
    tracer = maglib.Maglit(src)
)doc")
        .def(py::init([](const std::string &equilibrium_path, const std::string &perturbation_path,
                         double amplitude, double phase_shift) {
                 auto eq = std::make_unique<EquilibriumSource>(equilibrium_path);
                 if (!eq->is_valid())
                     throw std::runtime_error("failed to load equilibrium: " + equilibrium_path);
                 auto pert = std::make_unique<GridPerturbationSource>(perturbation_path);
                 if (!pert->is_valid())
                     throw std::runtime_error("failed to load perturbation: " + perturbation_path);
                 return std::make_unique<PerturbedSource>(std::move(eq), std::move(pert),
                                                          amplitude, phase_shift);
             }),
             "equilibrium_path"_a, "perturbation_path"_a, "amplitude"_a = 1.0, "phase_shift"_a = 0.0);
}
