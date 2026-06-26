#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_sode(py::module_ &m);
void bind_collider(py::module_ &m);
void bind_field_source(py::module_ &m);
void bind_superposition_source(py::module_ &m);
void bind_maglit(py::module_ &m);
void bind_footprint(py::module_ &m);
void bind_manifold(py::module_ &m);
