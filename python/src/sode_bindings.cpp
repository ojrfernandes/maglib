#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "sode_wrapper.h"

namespace py = pybind11;
using namespace py::literals;

static py::array_t<double> to_numpy(const std::vector<double> &v) {
    auto arr = py::array_t<double>(static_cast<py::ssize_t>(v.size()));
    std::copy(v.begin(), v.end(), arr.mutable_data());
    return arr;
}

void bind_sode(py::module_ &m) {
    py::enum_<sode_type>(m, "SodeMethod",
        "Runge-Kutta method used by the adaptive ODE solver.")
        .value("RK56_FB", SODE_RK56_FB, "Fehlberg RK5(6)")
        .value("RK56_CK", SODE_RK56_CK, "Cash-Karp RK5(6)")
        .value("RK78_DP", SODE_RK78_DP, "Dormand-Prince RK7(8)")
        .export_values();

    py::enum_<sode_status>(m, "SodeStatus",
        "Return status from Sode.step() or Sode.integrate().")
        .value("OK",              SODE_OK)
        .value("SUCCESS_TIME",    SODE_SUCCESS_TIME,    "Reached t_end")
        .value("SUCCESS_MONITOR", SODE_SUCCESS_MONITOR, "Event detected")
        .value("CONTINUE_GOOD",   SODE_CONTINUE_GOOD_STEP)
        .value("CONTINUE_BAD",    SODE_CONTINUE_BAD_STEP)
        .value("FAILED",          SODE_FAILED,          "Step size fell below h_min")
        .value("BAD_FUNC",        SODE_BAD_FUNC,        "System function returned an error")
        .export_values();

    py::class_<SodeSolver>(m, "Sode", R"doc(
Adaptive-step ODE solver with event detection.

Wraps the sode C++ engine (Fehlberg / Cash-Karp / Dormand-Prince embedded
Runge-Kutta pairs) with a Python-friendly callback interface.

Example::

    import numpy as np
    import maglib

    solver = maglib.Sode(maglib.SodeMethod.RK56_CK, dim=2)
    solver.configure_steps(1e-3, 1e-6, 0.1)
    solver.configure_tol(1e-10, 1e-12, 1e-12, 0.9)

    def sho(x, t):
        return np.array([x[1], -x[0]])

    solver.set_system(sho)
    x, t, status = solver.integrate(np.array([1.0, 0.0]), t0=0.0, t_end=2*np.pi)
)doc")
        .def(py::init<sode_type, int>(), "method"_a, "dim"_a,
             "Create solver. method: SodeMethod enum. dim: number of state variables.")

        .def("configure_steps", &SodeSolver::configure_steps,
             "h_init"_a, "h_min"_a, "h_max"_a,
             "Set initial, minimum, and maximum step size.")

        .def("configure_tol", &SodeSolver::configure_tol,
             "tol_sol"_a, "tol_mon"_a, "tol_end"_a, "damp"_a,
             "Set solution tolerance, monitor bisection tolerance, "
             "end-time tolerance, and step-size damping factor.")

        .def("set_system", [](SodeSolver &self, py::function fn) {
            self.set_system([fn](const std::vector<double> &x, double t) -> std::vector<double> {
                py::array_t<double> x_arr(static_cast<py::ssize_t>(x.size()), x.data());
                py::object          result = fn(x_arr, t);
                return result.cast<std::vector<double>>();
            });
        }, "fn"_a,
        "Set ODE right-hand side: fn(x: np.ndarray, t: float) -> np.ndarray. "
        "The returned array must have the same length as x.")

        .def("set_monitor", [](SodeSolver &self, py::function fn) {
            self.set_monitor([fn](const std::vector<double> &x, double t) -> bool {
                py::array_t<double> x_arr(static_cast<py::ssize_t>(x.size()), x.data());
                return fn(x_arr, t).cast<bool>();
            });
        }, "fn"_a,
        "Set event monitor: fn(x: np.ndarray, t: float) -> bool. "
        "Integration stops when the boolean value changes in the direction "
        "specified by monitor_dir in step() / integrate().")

        .def("reset", &SodeSolver::reset,
             "Reset integrator state (call before integrating a new trajectory).")

        .def("set_verbose", &SodeSolver::set_verbose,
             "Enable verbose diagnostic output.")

        .def("step", [](SodeSolver &self,
                        py::array_t<double, py::array::c_style | py::array::forcecast> x_arr,
                        double t, double t_end, int monitor_dir) {
            auto                buf = x_arr.request();
            std::vector<double> x(static_cast<double *>(buf.ptr),
                                   static_cast<double *>(buf.ptr) + buf.size);
            auto [new_t, status] = self.step(x, t, t_end, monitor_dir);
            return py::make_tuple(to_numpy(x), new_t, static_cast<sode_status>(status));
        },
        "x"_a, "t"_a, "t_end"_a, "monitor_dir"_a = 0,
        R"doc(
Take a single adaptive step.

Parameters
----------
x : array_like
    Current state vector (length dim).
t : float
    Current time.
t_end : float
    Target end time (step will not exceed this).
monitor_dir : int
    0 = no event detection; 1 = stop on False→True; -1 = stop on True→False.

Returns
-------
(x_new, t_new, status) : (np.ndarray, float, SodeStatus)
)doc")

        .def("integrate", [](SodeSolver &self,
                              py::array_t<double, py::array::c_style | py::array::forcecast> x0_arr,
                              double t0, double t_end, int monitor_dir) {
            auto                buf = x0_arr.request();
            std::vector<double> x0(static_cast<double *>(buf.ptr),
                                    static_cast<double *>(buf.ptr) + buf.size);
            auto [x, t, status] = self.integrate(x0, t0, t_end, monitor_dir);
            return py::make_tuple(to_numpy(x), t, static_cast<sode_status>(status));
        },
        "x0"_a, "t0"_a, "t_end"_a, "monitor_dir"_a = 0,
        R"doc(
Integrate from t0 to t_end (or until an event is detected).

Parameters
----------
x0 : array_like
    Initial state vector (length dim).
t0 : float
    Initial time.
t_end : float
    End time.
monitor_dir : int
    0 = no event detection; 1 = stop on False→True; -1 = stop on True→False.

Returns
-------
(x, t, status) : (np.ndarray, float, SodeStatus)
    Final state, time, and terminal status (SUCCESS_TIME or SUCCESS_MONITOR).
)doc");
}
