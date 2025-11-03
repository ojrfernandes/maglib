# tests/CMakeLists.txt

# Test for sode library
add_executable(test_sode test_sode.cpp)
target_link_libraries(test_sode 
    PRIVATE 
    sode
    ${GTEST_MAIN_LIBRARIES}
    ${GTEST_LIBRARIES}
)
target_include_directories(test_sode PRIVATE
    ${CMAKE_SOURCE_DIR}/sode
)

# Test for collider class in maglit library
add_executable(test_collider test_collider.cpp)
target_link_libraries(test_collider 
    PRIVATE
    maglit
    ${GTEST_MAIN_LIBRARIES}
    ${GTEST_LIBRARIES}
)
target_include_directories(test_collider PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
)

# Test for maglit library
add_executable(test_maglit test_maglit.cpp)
target_link_libraries(test_maglit 
    PRIVATE 
    maglit
    sode
    ${GTEST_MAIN_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${HDF5_LIBRARIES}
)
target_include_directories(test_maglit PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
    ${CMAKE_SOURCE_DIR}/sode
    ${HDF5_INCLUDE_DIRS}
)
target_compile_definitions(test_maglit PRIVATE
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data"
)

# Integration tests (for programs)
add_executable(test_fpgen
    test_fpgen.cpp
    ${CMAKE_SOURCE_DIR}/fpgen/src/input_read.cpp
    ${CMAKE_SOURCE_DIR}/fpgen/src/tcabr_collider.cpp
    ${CMAKE_SOURCE_DIR}/fpgen/src/footprint.cpp)
    
target_link_libraries(test_fpgen 
    PRIVATE 
    maglit
    sode
    ${GTEST_MAIN_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
)
target_include_directories(test_fpgen PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
    ${CMAKE_SOURCE_DIR}/sode
    ${HDF5_INCLUDE_DIRS}
)

# Register tests with CTest
add_test(NAME SodeTests COMMAND test_sode)
add_test(NAME ColliderTests COMMAND test_collider)
add_test(NAME MaglitTests COMMAND test_maglit)
add_test(NAME FpgenTests COMMAND test_fpgen)