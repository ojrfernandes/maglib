# mfgen/CMakeLists.txt

# Collect source files
file(GLOB MFGEN_SOURCES "src/*.cpp")

# Create executable
add_executable(mfgen ${MFGEN_SOURCES})

# Link libraries using modern CMake
target_link_libraries(mfgen PRIVATE
    maglit                     
    ${HDF5_LIBRARIES}          
    OpenMP::OpenMP_CXX         
    ${ARMADILLO_LIBRARIES}     
)

# Include directories
target_include_directories(mfgen PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
)

# Copy input file to build directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/mf_input.txt
    ${CMAKE_CURRENT_BINARY_DIR}/mf_input.txt
    COPYONLY
)