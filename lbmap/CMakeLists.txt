# lbmap/CMakeLists.txt

# Collect source files
file(GLOB LBMAP_SOURCES "src/*.cpp")

# Create executable
add_executable(lbmap ${LBMAP_SOURCES})

# Link libraries using modern CMake
target_link_libraries(lbmap PRIVATE
    ${HDF5_LIBRARIES}          
)

# Include directories
target_include_directories(lbmap PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
)

# Copy input file to build directory
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/lb_input.txt
    ${CMAKE_CURRENT_BINARY_DIR}/lb_input.txt
    COPYONLY
)