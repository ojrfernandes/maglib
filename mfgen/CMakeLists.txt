# mfgen/CMakeLists.txt

# Create executable
add_executable(mfgen
    src/input_read.cpp
    src/manifold.cpp
    src/run.cpp
)

# Link libraries using modern CMake
target_link_libraries(mfgen PRIVATE
    maglit
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# Include directories
target_include_directories(mfgen PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)

install(TARGETS mfgen DESTINATION bin)

# Copy input file to build directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/mfgen_input.txt
    ${CMAKE_BINARY_DIR}/bin/mfgen_input.txt
    COPYONLY
)