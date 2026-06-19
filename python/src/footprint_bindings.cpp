#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "footprint.h"
#include "maglit.h"

namespace py = pybind11;
using namespace py::literals;

void bind_footprint(py::module_ &m) {
    auto fp_to_numpy = [](const footprint &fp) {
        const auto   &data  = fp.get_output_data();
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
        stability=1,
        grid_R1=0.435, grid_Z1=-0.239,
        grid_R2=0.435, grid_Z2=-0.232,
        nRZ=100, nPhi=20, max_turns=1000,
    )
    fp.run(tracers)
    data = fp.output_data   # numpy array (nRZ*nPhi, 6)
    fp.save("output.dat")
)doc")
        .def(py::init<int, double, double, double, double, int, int, int>(),
             "stability"_a, "grid_R1"_a, "grid_Z1"_a,
             "grid_R2"_a, "grid_Z2"_a,
             "nRZ"_a, "nPhi"_a, "max_turns"_a,
             R"doc(
Create a Footprint grid.

Parameters
----------
stability : int
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
            auto ends_with = [&](const std::string &ext) {
                return path.size() >= ext.size() &&
                       path.substr(path.size() - ext.size()) == ext;
            };

            if (ends_with(".npy")) {
                py::module_::import("numpy").attr("save")(path, fp_to_numpy(self));
            } else if (ends_with(".npz")) {
                py::module_::import("numpy").attr("savez_compressed")(
                    path, "data"_a = fp_to_numpy(self));
            } else {
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
