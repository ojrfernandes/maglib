# python/CMakeLists.txt
# Builds the _maglib Python extension module.
#
# C++ sources are compiled directly into the extension rather than linking
# against the maglit static library, so each source's dependencies are
# listed explicitly here. Fusion-IO and HDF5 are required from this point
# onward (needed by maglit.cpp / m3dc1_source.cpp).

pybind11_add_module(_maglib
    src/sode_wrapper.cpp
    src/maglib_module.cpp
    src/sode_bindings.cpp
    src/collider_bindings.cpp
    src/field_source_bindings.cpp
    src/superposition_source_bindings.cpp
    src/maglit_bindings.cpp
    src/footprint_bindings.cpp
    src/manifold_bindings.cpp
    ${CMAKE_SOURCE_DIR}/maglit/sode.cpp
    ${CMAKE_SOURCE_DIR}/maglit/collider.cpp
    ${CMAKE_SOURCE_DIR}/maglit/maglit.cpp
    ${CMAKE_SOURCE_DIR}/maglit/m3dc1_source.cpp
    ${CMAKE_SOURCE_DIR}/maglit/superposition_source.cpp
    ${CMAKE_SOURCE_DIR}/fpgen/src/footprint.cpp
    ${CMAKE_SOURCE_DIR}/mfgen/src/manifold.cpp
)

target_link_libraries(_maglib PRIVATE
    ${HDF5_LIBRARIES}
    ${FUSION_IO_LIB_DIR}/libfusionio.so
    ${FUSION_IO_LIB_DIR}/libm3dc1.so
    OpenMP::OpenMP_CXX
)

target_include_directories(_maglib PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
    ${CMAKE_SOURCE_DIR}/fpgen/src
    ${CMAKE_SOURCE_DIR}/mfgen/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
    ${FUSION_IO_INCLUDE_DIR}
)

target_compile_features(_maglib PRIVATE cxx_std_17)

# Place the compiled extension alongside __init__.py so that adding
# python/ to PYTHONPATH is all that is needed to import maglib.
set_target_properties(_maglib PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/maglib
)

# Install: binaries go to ${prefix}/bin, Python package to ${prefix}/python/maglib/.
# Add ${prefix}/python to PYTHONPATH to import maglib after installation.
install(TARGETS _maglib DESTINATION python/maglib)
install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/maglib/__init__.py
    ${CMAKE_CURRENT_SOURCE_DIR}/maglib/plot.py
    ${CMAKE_CURRENT_SOURCE_DIR}/maglib/manifold_tools.py
    DESTINATION python/maglib
)
