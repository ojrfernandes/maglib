# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# Define the project
project(maglib)

# Enable C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#=============================================================================

# Find standard dependencies
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(Armadillo REQUIRED)
find_package(OpenMP REQUIRED)

# Custom dependency - fusion-io
set(FUSION_IO_DIR "$ENV{FUSION_IO_DIR}" CACHE PATH "Path to fusion-io build directory")
if(NOT FUSION_IO_DIR)
    message(FATAL_ERROR "FUSION_IO_DIR must be set either as environment variable or cmake option")
endif()

# Verify fusion-io structure
if(NOT EXISTS "${FUSION_IO_DIR}/lib")
    message(FATAL_ERROR "Cannot find fusion-io libraries at ${FUSION_IO_DIR}/lib")
endif()

# Set fusion-io paths for subdirectories to use
set(FUSION_IO_INCLUDE_DIR "${FUSION_IO_DIR}/include")
set(FUSION_IO_LIB_DIR "${FUSION_IO_DIR}/lib")

# Status messages
message(STATUS "FUSION_IO_DIR: ${FUSION_IO_DIR}")
message(STATUS "FUSION_IO_INCLUDE_DIR: ${FUSION_IO_INCLUDE_DIR}")
message(STATUS "FUSION_IO_LIB_DIR: ${FUSION_IO_LIB_DIR}")

#=============================================================================

# Include subdirectories
add_subdirectory(maglit)
add_subdirectory(sode)
add_subdirectory(fpgen)
add_subdirectory(mfgen)
add_subdirectory(lbmap)