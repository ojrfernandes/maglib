# sode/CMakeLists.txt

add_library(sode STATIC
    sode.cpp
    sode.h
)

target_include_directories(sode PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)