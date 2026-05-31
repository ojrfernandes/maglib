# python/CMakeLists.txt
# Builds the _maglib Python extension module.
#
# sode.cpp is compiled directly here — this avoids linking against the full
# maglit static library, which requires Fusion-IO and HDF5. As more C++ classes
# are exposed to Python those will be added here incrementally.

pybind11_add_module(_maglib
    src/sode_wrapper.cpp
    src/maglib_module.cpp
    ${CMAKE_SOURCE_DIR}/maglit/sode.cpp
    ${CMAKE_SOURCE_DIR}/maglit/collider.cpp
)

target_include_directories(_maglib PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_features(_maglib PRIVATE cxx_std_17)

# Place the compiled extension alongside __init__.py so that adding
# python/ to PYTHONPATH is all that is needed to import maglib.
set_target_properties(_maglib PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/maglib
)
