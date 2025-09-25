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

# Integration tests (for programs)
add_executable(test_integration test_integration.cpp)
target_link_libraries(test_integration 
    PRIVATE 
    maglit
    sode
    ${GTEST_MAIN_LIBRARIES}
    ${GTEST_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    OpenMP::OpenMP_CXX
)
target_include_directories(test_integration PRIVATE
    ${CMAKE_SOURCE_DIR}/maglit
    ${CMAKE_SOURCE_DIR}/sode
    ${HDF5_INCLUDE_DIRS}
    ${ARMADILLO_INCLUDE_DIRS}
)

# Register tests with CTest
add_test(NAME SodeTests COMMAND test_sode)
add_test(NAME MaglitTests COMMAND test_maglit)
add_test(NAME IntegrationTests COMMAND test_integration)