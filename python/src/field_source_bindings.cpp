#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "field_source.h"
#include "m3dc1_source.h"

namespace py = pybind11;
using namespace py::literals;

void bind_field_source(py::module_ &m) {
    py::class_<FieldSource>(m, "FieldSource", R"doc(
Abstract base class for magnetic field data sources.

Concrete implementations (e.g. M3DC1Source) provide eval_B, eval_psin,
and eval_psi. One instance per thread is required for multi-threaded use.
)doc")
        .def("eval_B", [](FieldSource &self, double R, double phi, double Z) {
            double B[3];
            bool ok = self.eval_B(R, phi, Z, B);
            auto arr = py::array_t<double>(3);
            std::copy(B, B + 3, arr.mutable_data());
            return py::make_tuple(ok, arr);
        }, "R"_a, "phi"_a, "Z"_a,
        "Evaluate magnetic field. Returns (success, B) where B is a (3,) array [B_R, B_phi, B_Z].")

        .def("eval_psin", [](FieldSource &self, double R, double phi, double Z) {
            double psin;
            bool ok = self.eval_psin(R, phi, Z, psin);
            return py::make_tuple(ok, psin);
        }, "R"_a, "phi"_a, "Z"_a,
        "Evaluate normalised poloidal flux. Returns (success, psin).")

        .def("eval_psi", [](FieldSource &self, double R, double phi, double Z) {
            double psi;
            bool ok = self.eval_psi(R, phi, Z, psi);
            return py::make_tuple(ok, psi);
        }, "R"_a, "phi"_a, "Z"_a,
        "Evaluate poloidal flux. Returns (success, psi).");

    py::class_<M3DC1Source, FieldSource>(m, "M3DC1Source", R"doc(
M3D-C1 magnetic field source backed by a Fusion-IO HDF5 file.

Example::

    import maglib

    source = maglib.M3DC1Source("C1.h5", timeslice=-1)
    if not source.is_valid():
        raise RuntimeError("failed to open field source")

    ok, B = source.eval_B(0.62, 0.0, 0.0)
)doc")
        .def(py::init<const char *, int>(), "path"_a, "timeslice"_a,
             "Open an M3D-C1 HDF5 file. timeslice=-1 selects the last available slice.")

        .def("is_valid", &M3DC1Source::is_valid,
             "Return True if the file was opened and the magnetic field loaded successfully.");
}
