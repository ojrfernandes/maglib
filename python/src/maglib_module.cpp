#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "collider.h"
#include "field_source.h"
#include "footprint.h"
#include "m3dc1_source.h"
#include "maglit.h"
#include "sode_wrapper.h"

namespace py = pybind11;
using namespace py::literals;

// Helper: copy a std::vector<double> into a new 1-D numpy array (owns its memory).
static py::array_t<double> to_numpy(const std::vector<double> &v) {
    auto arr = py::array_t<double>(static_cast<py::ssize_t>(v.size()));
    std::copy(v.begin(), v.end(), arr.mutable_data());
    return arr;
}

PYBIND11_MODULE(_maglib, m) {
    m.doc() = "maglib — magnetic field line integration and ODE utilities";

    // ── Integration method ────────────────────────────────────────────────────
    py::enum_<sode_type>(m, "SodeMethod",
        "Runge-Kutta method used by the adaptive ODE solver.")
        .value("RK56_FB", SODE_RK56_FB, "Fehlberg RK5(6)")
        .value("RK56_CK", SODE_RK56_CK, "Cash-Karp RK5(6)")
        .value("RK78_DP", SODE_RK78_DP, "Dormand-Prince RK7(8)")
        .export_values();

    // ── Step status ───────────────────────────────────────────────────────────
    py::enum_<sode_status>(m, "SodeStatus",
        "Return status from Sode.step() or Sode.integrate().")
        .value("OK",             SODE_OK)
        .value("SUCCESS_TIME",   SODE_SUCCESS_TIME,    "Reached t_end")
        .value("SUCCESS_MONITOR",SODE_SUCCESS_MONITOR, "Event detected")
        .value("CONTINUE_GOOD",  SODE_CONTINUE_GOOD_STEP)
        .value("CONTINUE_BAD",   SODE_CONTINUE_BAD_STEP)
        .value("FAILED",         SODE_FAILED,          "Step size fell below h_min")
        .value("BAD_FUNC",       SODE_BAD_FUNC,        "System function returned an error")
        .export_values();

    // ── Sode class ────────────────────────────────────────────────────────────
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
                // GIL is held: we are always called from Python context.
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

    // ── Collider class ────────────────────────────────────────────────────────
    py::class_<collider>(m, "Collider", R"doc(
Vessel wall boundary tester.

Loads a tokamak vessel wall polygon from a two-column text file (R Z per line)
and tests whether points are inside the enclosed region using the winding-number
algorithm.

Example::

    import maglib

    wall = maglib.Collider()
    wall.load_shape("tcabr_first_wall.txt")
    print(wall.inside(0.62, 0.0))   # True — inside TCABR
    print(wall.inside(2.0,  0.0))   # False — outside
)doc")
        .def(py::init<>(), "Create an unloaded Collider.")

        .def("load_shape", &collider::load_shape, "path"_a,
             "Load vessel wall vertices from a whitespace-separated R Z text file. "
             "Returns True on success.")

        .def("inside", &collider::inside, "R"_a, "Z"_a,
             "Return True if the point (R, Z) is inside (or on the boundary of) the vessel.")

        .def("is_loaded", &collider::is_loaded,
             "Return True if a vessel shape has been successfully loaded.")

        .def("get_vertices", [](const collider &self) {
            const auto &verts = self.get_vertices();
            auto arr = py::array_t<double>({static_cast<py::ssize_t>(verts.size()),
                                            static_cast<py::ssize_t>(2)});
            auto buf = arr.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(verts.size()); ++i) {
                buf(i, 0) = verts[static_cast<size_t>(i)].first;
                buf(i, 1) = verts[static_cast<size_t>(i)].second;
            }
            return arr;
        }, "Return vessel wall vertices as a (N, 2) numpy array of (R, Z) coordinates.");

    // ── FieldSource abstract base ─────────────────────────────────────────────
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

    // ── M3DC1Source ───────────────────────────────────────────────────────────
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

    // ── Maglit ────────────────────────────────────────────────────────────────
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

    // ── Footprint ─────────────────────────────────────────────────────────────

    // Helper: copy outputData into a (nRZ*nPhi, 6) numpy array.
    auto fp_to_numpy = [](const footprint &fp) {
        const auto   &data  = fp.outputData;
        py::ssize_t   nrows = static_cast<py::ssize_t>(data.size());
        auto          arr   = py::array_t<double>({nrows, static_cast<py::ssize_t>(6)});
        auto          buf   = arr.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < nrows; ++i)
            for (py::ssize_t j = 0; j < 6; ++j)
                buf(i, j) = data[static_cast<size_t>(i)][static_cast<size_t>(j)];
        return arr;
    };

    py::class_<footprint>(m, "Footprint", R"doc(
Magnetic footprint generator.

Traces a 2-D grid of field lines from a target plate segment through nPhi
toroidal angles and records connection length, minimum psi_N, and turn count
per point.

Pass one Maglit tracer per thread to run() for parallel execution; a
single-element list gives serial execution.

Example::

    import maglib

    sources = [maglib.M3DC1Source("C1.h5", 1) for _ in range(4)]
    tracers = [maglib.Maglit(s) for s in sources]
    for t in tracers:
        t.configure(0.01, 1e-6, 0.1)
        t.set_monitor("wall.txt")

    fp = maglib.Footprint(
        manifold=1,
        grid_R1=0.435, grid_Z1=-0.239,
        grid_R2=0.435, grid_Z2=-0.232,
        nRZ=100, nPhi=20, max_turns=1000,
    )
    fp.run(tracers)
    data = fp.output_data   # numpy array (nRZ*nPhi, 6)
    fp.save("output.dat")
)doc")
        .def(py::init<int, double, double, double, double, int, int, int>(),
             "manifold"_a, "grid_R1"_a, "grid_Z1"_a,
             "grid_R2"_a, "grid_Z2"_a,
             "nRZ"_a, "nPhi"_a, "max_turns"_a,
             R"doc(
Create a Footprint grid.

Parameters
----------
manifold : int
    0 = stable (backward map), 1 = unstable (forward map).
grid_R1, grid_Z1 : float
    First endpoint of the target plate segment (metres).
grid_R2, grid_Z2 : float
    Second endpoint of the target plate segment (metres).
nRZ : int
    Number of grid points along the plate segment.
nPhi : int
    Number of toroidal starting angles.
max_turns : int
    Maximum toroidal turns before a field line is considered lost.
)doc")

        .def("run", [](footprint &self, py::list tracers) {
            std::vector<maglit*> ptrs;
            ptrs.reserve(py::len(tracers));
            for (auto &item : tracers)
                ptrs.push_back(item.cast<maglit*>());
            py::gil_scoped_release release;
            self.run(ptrs);
        }, "tracers"_a,
        "Run the grid. Pass one Maglit per thread for parallel execution, "
        "or a single-element list for serial execution. The GIL is released "
        "during integration.")

        .def_property_readonly("output_data", fp_to_numpy,
            R"doc(
Output as a (nRZ*nPhi, 6) numpy array.

Columns: R0, Z0, phi0, length, psiMin, turn
)doc")

        .def("save", [fp_to_numpy](const footprint &self, const std::string &path) {
            // Infer format from extension
            auto ends_with = [&](const std::string &ext) {
                return path.size() >= ext.size() &&
                       path.substr(path.size() - ext.size()) == ext;
            };

            if (ends_with(".npy")) {
                py::module_::import("numpy").attr("save")(path, fp_to_numpy(self));
            } else if (ends_with(".npz")) {
                // numpy adds .npz automatically if not present; pass path directly
                py::module_::import("numpy").attr("savez_compressed")(
                    path, "data"_a = fp_to_numpy(self));
            } else {
                // delegate .dat / .txt / .csv to C++
                if (!self.save(path))
                    throw std::runtime_error("Failed to save: " + path);
            }
        }, "path"_a,
        R"doc(
Save output_data to file. Format is inferred from the extension:

  .dat / .txt  — space-separated text with header
  .csv         — comma-separated text with header
  .npy         — numpy binary (numpy.save)
  .npz         — numpy compressed archive (numpy.savez_compressed),
                 array stored under key 'data'
)doc");
}
