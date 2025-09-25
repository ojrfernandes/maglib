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

# Testing setup
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    enable_testing()
    
    # Try to find system GoogleTest first
    find_package(GTest QUIET)
    
    if(NOT GTest_FOUND)
        # Check if we already have googletest cloned
        set(GTEST_ROOT ${CMAKE_SOURCE_DIR}/external/googletest)
        
        if(NOT EXISTS ${GTEST_ROOT})
            message(STATUS "GoogleTest not found, cloning from GitHub...")
            
            # Create external directory
            file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/external)
            
            # Clone googletest
            execute_process(
                COMMAND git clone https://github.com/google/googletest.git ${GTEST_ROOT}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external
                RESULT_VARIABLE GIT_RESULT
            )
            
            if(NOT GIT_RESULT EQUAL "0")
                message(FATAL_ERROR "Failed to clone GoogleTest from GitHub")
            endif()
            
            message(STATUS "GoogleTest cloned successfully")
        else()
            message(STATUS "Using existing GoogleTest at ${GTEST_ROOT}")
        endif()
        
        # Build googletest as part of our project
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
        set(BUILD_GTEST ON CACHE BOOL "" FORCE)
        add_subdirectory(${GTEST_ROOT} ${CMAKE_BINARY_DIR}/googletest EXCLUDE_FROM_ALL)
        
        # Set up variables for our tests
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

# Include subdirectories
add_subdirectory(maglit)
add_subdirectory(sode)
add_subdirectory(fpgen)
add_subdirectory(mfgen)
add_subdirectory(lbmap)