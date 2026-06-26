#include <pybind11/pybind11.h>
#include "bindings.h"

namespace py = pybind11;

PYBIND11_MODULE(_maglib, m) {
    m.doc() = "maglib — magnetic field line integration and ODE utilities";
    bind_sode(m);
    bind_collider(m);
    bind_field_source(m);
    bind_superposition_source(m);
    bind_maglit(m);
    bind_footprint(m);
    bind_manifold(m);
}
