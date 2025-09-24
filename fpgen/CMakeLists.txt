# fpgen/CMakeLists.txt

# Collect source files
file(GLOB FPGEN_SOURCES "src/*.cpp")

# Create executable
add_executable(fpgen ${FPGEN_SOURCES})

# Link libraries using modern CMake
target_link_libraries(fpgen PRIVATE
    maglit                     # Our library (includes sode automatically)
    ${HDF5_LIBRARIES}          # HDF5 from find_package
    OpenMP::OpenMP_CXX         # OpenMP from find_package
    ${ARMADILLO_LIBRARIES}     # Armadillo from find_package
)

# Include directories
target_include_directories(fpgen PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
)

# Copy input file to build directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/fp_input.txt
    ${CMAKE_CURRENT_BINARY_DIR}/fp_input.txt
    COPYONLY
)