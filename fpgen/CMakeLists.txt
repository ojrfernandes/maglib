# fpgen/CMakeLists.txt

# Create executable
add_executable(fpgen
    src/input_read.cpp
    src/footprint.cpp
    src/run.cpp
)

# Link libraries using modern CMake
target_link_libraries(fpgen PRIVATE
    maglit                     # Our library (includes sode automatically)
    ${HDF5_LIBRARIES}          # HDF5 from find_package
    OpenMP::OpenMP_CXX         # OpenMP from find_package
)

# Include directories
target_include_directories(fpgen PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)

install(TARGETS fpgen DESTINATION bin)

# Copy input file to build directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/fpgen_input.txt
    ${CMAKE_BINARY_DIR}/bin/fpgen_input.txt
    COPYONLY
)