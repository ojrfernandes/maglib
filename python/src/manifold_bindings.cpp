#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "maglit.h"
#include "manifold.h"

namespace py = pybind11;
using namespace py::literals;

void bind_manifold(py::module_ &m) {

    // Helper: std::vector<point> → numpy (N, 2) float64 array, columns [R, Z]
    auto seg_to_numpy = [](const std::vector<point> &seg) {
        auto arr = py::array_t<double>({static_cast<py::ssize_t>(seg.size()),
                                        static_cast<py::ssize_t>(2)});
        auto buf = arr.mutable_unchecked<2>();
        for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(seg.size()); ++i) {
            buf(i, 0) = seg[static_cast<size_t>(i)].R;
            buf(i, 1) = seg[static_cast<size_t>(i)].Z;
        }
        return arr;
    };

    // Helper: numpy (N, 2) array → std::vector<point>
    auto numpy_to_seg = [](py::array_t<double, py::array::c_style | py::array::forcecast> arr)
                           -> std::vector<point> {
        if (arr.ndim() != 2 || arr.shape(1) != 2)
            throw std::invalid_argument("segment must be a (N, 2) array");
        auto buf = arr.unchecked<2>();
        std::vector<point> seg(static_cast<size_t>(buf.shape(0)));
        for (py::ssize_t i = 0; i < buf.shape(0); ++i)
            seg[static_cast<size_t>(i)] = {buf(i, 0), buf(i, 1)};
        return seg;
    };

    py::class_<manifold>(m, "Manifold", R"doc(
Invariant manifold calculator for the X-point of the Poincaré map.

Computes stable and unstable invariant manifolds by iteratively applying
the Poincaré map and refining segments.

The Maglit tracer must outlive this object.

Example::

    import numpy as np
    import maglib

    source = maglib.M3DC1Source("C1.h5", 1)
    tracer = maglib.Maglit(source)
    tracer.configure(1e-2, 1e-6, 1e-2)

    mf = maglib.Manifold(tracer, phi=0.0, stability=0)
    mf.configure(epsilon=1e-6, h=1e-8, tol=1e-14,
                 max_iter=50, precision_limit=1e-14, max_insertions=100)

    mf.find_x_point(0.498, -0.219)
    seg0 = mf.primary_segment(10)
    prev, seg1 = mf.new_segment(seg0, phi=0.0, l_lim=0.005, theta_lim=20.0)
)doc")

        .def(py::init<maglit &, double, int>(),
             "tracer"_a, "phi"_a, "stability"_a,
             py::keep_alive<1, 2>(),
             R"doc(
Create a Manifold object.

Parameters
----------
tracer : Maglit
    Field-line integrator. Must outlive this object.
phi : float
    Toroidal angle of the Poincaré section (radians).
stability : int
    0 = stable manifold (forward map), 1 = unstable manifold (backward map).
)doc")

        .def("configure", &manifold::configure,
             "epsilon"_a, "h"_a, "tol"_a, "max_iter"_a,
             "precision_limit"_a, "max_insertions"_a,
             R"doc(
Set numerical parameters.

Parameters
----------
epsilon : float
    Distance from the X-point to the pivot point.
h : float
    Step size for numerical Jacobian derivatives.
tol : float
    Convergence tolerance for Newton's method.
max_iter : int
    Maximum Newton iterations.
precision_limit : float
    Minimum arc length below which point insertion is skipped.
max_insertions : int
    Maximum insertions per refinement pass before giving up.
)doc")

        .def("set_verbose", &manifold::setVerbose,
             "Enable verbose diagnostic output to stdout.")

        .def("find_x_point",
            [](manifold &self, double r_guess, double z_guess) {
                py::gil_scoped_release release;
                return self.find_xPoint(r_guess, z_guess);
            },
            "r_guess"_a, "z_guess"_a,
            R"doc(
Find the nearest 1-period fixed point (X-point) by Newton's method.

Parameters
----------
r_guess, z_guess : float
    Initial guess for the X-point position (metres).

Returns
-------
bool
    True if convergence was achieved within max_iter iterations.
)doc")

        .def_property("x_point",
            [](const manifold &self) {
                auto arr = py::array_t<double>(2);
                arr.mutable_data()[0] = self.xPoint.R;
                arr.mutable_data()[1] = self.xPoint.Z;
                return arr;
            },
            [](manifold &self,
               py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
                if (arr.size() != 2)
                    throw std::invalid_argument("x_point must be a (2,) array [R, Z]");
                self.xPoint = {arr.data()[0], arr.data()[1]};
            },
            "X-point coordinates as a (2,) numpy array [R, Z]. "
            "Set by find_x_point(); can also be assigned directly.")

        .def("primary_segment",
            [seg_to_numpy](manifold &self, size_t n_intervals) {
                std::vector<point> seg;
                {
                    py::gil_scoped_release release;
                    seg = self.primarySegment(n_intervals);
                }
                return seg_to_numpy(seg);
            },
            "n_intervals"_a,
            R"doc(
Compute the primary segment of the manifold.

Requires x_point to be set (call find_x_point first).

Parameters
----------
n_intervals : int
    Number of intervals; produces n_intervals+1 points.

Returns
-------
np.ndarray, shape (n_intervals+1, 2)
    Segment points. Columns: [R, Z].
)doc")

        .def("new_segment",
            [seg_to_numpy, numpy_to_seg](manifold &self,
                py::array_t<double, py::array::c_style | py::array::forcecast> prev_arr,
                double l_lim, double theta_lim) {
                std::vector<point> prev = numpy_to_seg(prev_arr);
                std::vector<point> new_seg;
                {
                    py::gil_scoped_release release;
                    new_seg = self.newSegment(prev, l_lim, theta_lim);
                }
                return py::make_tuple(seg_to_numpy(prev), seg_to_numpy(new_seg));
            },
            "prev_seg"_a, "l_lim"_a, "theta_lim"_a,
            R"doc(
Compute a refined segment from a previous segment (interpolant method).

The Poincaré section angle is the one set at construction time.

Parameters
----------
prev_seg : np.ndarray, shape (N, 2)
    Input segment; may be refined with additional points.
l_lim : float
    Arc-length threshold for segment refinement (metres).
theta_lim : float
    Turning-angle threshold in degrees.

Returns
-------
(updated_prev_seg, new_seg) : tuple of np.ndarray, shapes (M, 2) and (K, 2)
    updated_prev_seg is prev_seg after any inserted refinement points.
    new_seg is the image of prev_seg under one application of the Poincaré map.
)doc")

        .def("new_segment",
            [seg_to_numpy, numpy_to_seg](manifold &self,
                py::array_t<double, py::array::c_style | py::array::forcecast> prev_arr,
                int n_seg, double l_lim, double theta_lim) {
                std::vector<point> prev = numpy_to_seg(prev_arr);
                std::vector<point> new_seg;
                {
                    py::gil_scoped_release release;
                    new_seg = self.newSegment(prev, n_seg, l_lim, theta_lim);
                }
                return py::make_tuple(seg_to_numpy(prev), seg_to_numpy(new_seg));
            },
            "prev_seg"_a, "n_seg"_a, "l_lim"_a, "theta_lim"_a,
            R"doc(
Compute a refined segment by applying the map n_seg times (exact-map method).

The Poincaré section angle is the one set at construction time.

Parameters
----------
prev_seg : np.ndarray, shape (N, 2)
    Primary segment (used as source; not modified).
n_seg : int
    Number of Poincaré map applications.
l_lim : float
    Arc-length threshold (metres).
theta_lim : float
    Turning-angle threshold in degrees.

Returns
-------
(prev_seg, new_seg) : tuple of np.ndarray, shapes (N, 2) and (K, 2)
)doc")

        .def("run",
            [](manifold &self, size_t n_intervals, int n_segments, int method,
               double l_lim, double theta_lim) {
                py::gil_scoped_release release;
                self.run(n_intervals, n_segments, method, l_lim, theta_lim);
            },
            "n_intervals"_a, "n_segments"_a, "method"_a, "l_lim"_a, "theta_lim"_a,
            R"doc(
Compute all manifold segments in one call.

Requires x_point to be set (call find_x_point first). Results are
accumulated in output_data and can be saved with save().

Parameters
----------
n_intervals : int
    Intervals for the primary segment; produces n_intervals+1 points.
n_segments : int
    Total number of segments including the primary (>= 1).
method : int
    0 = exact-map (each segment mapped directly from the primary),
    1 = interpolant (each segment mapped from the previous one).
l_lim : float
    Arc-length threshold for refinement (metres).
theta_lim : float
    Turning-angle threshold in degrees.
)doc")

        .def("progress_bar", &manifold::progressBar, "j"_a, "n_seg"_a,
             "Print a progress bar for segment j of n_seg to stdout.")

        .def_property_readonly("output_data",
            [seg_to_numpy](const manifold &self) {
                py::list result;
                for (const auto &seg : self.get_output_data())
                    result.append(seg_to_numpy(seg));
                return result;
            },
            R"doc(
Accumulated segment data as a list of (N_i, 2) float64 arrays.

Each call to primary_segment() appends its result here; each call to
new_segment() appends the new segment (second element of the returned tuple).
Columns: [R, Z].
)doc")

        .def("save", [seg_to_numpy](const manifold &self, const std::string &path) {
            auto ends_with = [&](const std::string &ext) {
                return path.size() >= ext.size() &&
                       path.substr(path.size() - ext.size()) == ext;
            };
            if (ends_with(".npz")) {
                py::object savez = py::module_::import("numpy").attr("savez_compressed");
                py::tuple  args  = py::make_tuple(py::str(path));
                py::dict   kw;
                for (size_t i = 0; i < self.get_output_data().size(); ++i)
                    kw[("seg_" + std::to_string(i)).c_str()] = seg_to_numpy(self.get_output_data()[i]);
                py::object r = py::reinterpret_steal<py::object>(
                    PyObject_Call(savez.ptr(), args.ptr(), kw.ptr()));
                if (!r) throw py::error_already_set();
            } else {
                if (!self.save(path))
                    throw std::runtime_error("Failed to save: " + path);
            }
        }, "path"_a,
        R"doc(
Save accumulated segment data to file. Format is inferred from the extension:

  .dat / .txt  — space-separated text with header (columns: seg, R, Z)
  .csv         — comma-separated text with header (columns: seg, R, Z)
  .npz         — numpy compressed archive; arrays stored as seg_0, seg_1, ...

Raises RuntimeError for unsupported extensions.
)doc");
}
