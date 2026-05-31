#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "field_source.h"
#include "maglit.h"
#include "sode.h"

namespace py = pybind11;
using namespace py::literals;

void bind_maglit(py::module_ &m) {
    py::class_<maglit>(m, "Maglit", R"doc(
Magnetic field-line integrator.

Integrates the field-line ODE dR/dphi = B_R/B_phi, dZ/dphi = B_Z/B_phi in
cylindrical coordinates using an adaptive Runge-Kutta solver.

The FieldSource must outlive this object.

Example::

    import maglib

    source = maglib.M3DC1Source("C1.h5", timeslice=-1)
    tracer = maglib.Maglit(source)
    tracer.configure(0.01, 1e-6, 0.1)

    R, Z, phi = 0.62, 0.0, 0.0
    phi_max = 2 * 3.14159265
    R, Z, phi, status = tracer.step(R, Z, phi, phi_max)
)doc")
        .def(py::init<FieldSource &>(), "source"_a, py::keep_alive<1, 2>(),
             "Create integrator. source must outlive this Maglit instance.")

        .def("configure", &maglit::configure,
             "dphi_init"_a, "dphi_min"_a, "dphi_max"_a,
             "Set initial, minimum, and maximum toroidal step size.")

        .def("inverse_map", &maglit::inverse_map, "inverse"_a,
             "Reverse the integration direction (True = backward map).")

        .def("reset", &maglit::reset,
             "Reset integrator state before starting a new field line.")

        .def("set_monitor", &maglit::set_monitor, "path"_a,
             "Load a vessel wall polygon from a text file and activate "
             "boundary detection in step().")

        .def("step", [](maglit &self,
                        double R, double Z, double phi,
                        double phi_max, int dir) {
            int status = self.step(R, Z, phi, phi_max, dir);
            return py::make_tuple(R, Z, phi, static_cast<sode_status>(status));
        },
        "R"_a, "Z"_a, "phi"_a, "phi_max"_a, "dir"_a = 0,
        R"doc(
Advance the field line by one adaptive step.

Parameters
----------
R, Z : float
    Current position in the poloidal plane (metres).
phi : float
    Current toroidal angle (radians).
phi_max : float
    Integration will not advance past this angle.
dir : int
    Boundary monitor direction: 0 = ignore, 1 = stop on inside→outside,
    -1 = stop on outside→inside.

Returns
-------
(R, Z, phi, status) : (float, float, float, SodeStatus)
)doc")

        .def("calc_mag_field", [](maglit &self, double R, double phi, double Z) {
            double x[3] = {R, phi, Z};
            double B[3];
            bool ok = self.calc_mag_field(x, B);
            auto arr = py::array_t<double>(3);
            std::copy(B, B + 3, arr.mutable_data());
            return py::make_tuple(ok, arr);
        },
        "R"_a, "phi"_a, "Z"_a,
        "Evaluate the magnetic field at (R, phi, Z). "
        "Returns (success, B) where B is a (3,) array [B_R, B_phi, B_Z].")

        .def("psin_eval", [](maglit &self, double R, double phi, double Z) {
            double psin;
            self.psin_eval(R, phi, Z, &psin);
            return psin;
        },
        "R"_a, "phi"_a, "Z"_a,
        "Return the normalised poloidal flux psi_N at (R, phi, Z).")

        .def("psi_eval", [](maglit &self, double R, double phi, double Z) {
            double psi;
            self.psi_eval(R, phi, Z, &psi);
            return psi;
        },
        "R"_a, "phi"_a, "Z"_a,
        "Return the poloidal flux psi at (R, phi, Z).")

        .def_readonly("boundary", &maglit::boundary,
            "Vessel wall boundary (Collider). Populated after set_monitor() is called.");
}
