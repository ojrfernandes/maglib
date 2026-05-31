#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include "collider.h"

namespace py = pybind11;
using namespace py::literals;

void bind_collider(py::module_ &m) {
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
}
