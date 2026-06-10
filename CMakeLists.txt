# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# Define the project
project(maglib)

# Enable C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#=============================================================================

# Find standard dependencies
find_package(HDF5 REQUIRED COMPONENTS C)
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

# Testing setup
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    enable_testing()
    
    # Try to find system GoogleTest first
    find_package(GTest QUIET)
    
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.15.2
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
        set(BUILD_GTEST ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)

        set(GTEST_LIBRARIES gtest)
        set(GTEST_MAIN_LIBRARIES gtest_main)
    else()
        message(STATUS "Using system GoogleTest")
        set(GTEST_LIBRARIES GTest::gtest)
        set(GTEST_MAIN_LIBRARIES GTest::gtest_main)
    endif()
    
    # Add tests subdirectory
    add_subdirectory(tests)
endif()

#=============================================================================

# Set output directories for all executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)

#=============================================================================

# Python bindings (requires Python development headers; skipped if not found)
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" ON)
if(BUILD_PYTHON_BINDINGS)
    find_package(Python COMPONENTS Interpreter Development QUIET)
    if(Python_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG        v2.13.6
        )
        FetchContent_MakeAvailable(pybind11)
        add_subdirectory(python)
        message(STATUS "Python bindings enabled (Python ${Python_VERSION})")
    else()
        message(STATUS "Python development headers not found — skipping Python bindings")
    endif()
endif()

#=============================================================================

# Include subdirectories
add_subdirectory(maglit)
add_subdirectory(fpgen)
add_subdirectory(mfgen)
add_subdirectory(lbmap)