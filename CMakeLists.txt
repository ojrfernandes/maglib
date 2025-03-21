# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# Define the project
project(maglib)

# Enable C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#=============================================================================
# Edit the following paths to match your system

set(FIO_DIR "/home/jfernandes/fusion-io/build" CACHE STRING "Fusion IO directory")
set(HDF5_DIR "/home/jfernandes/CMake-hdf5-1.14.0/hdf5-1.14.0/hdf5" CACHE STRING "HDF5 directory")

#=============================================================================

# Include subdirectories
add_subdirectory(maglit)
add_subdirectory(sode)
add_subdirectory(fpgen)
add_subdirectory(mfgen)
add_subdirectory(lbmap)
