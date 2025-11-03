# maglit/CMakeLists.txt
add_library(maglit STATIC 
    maglit.cpp
    maglit.h
    collider.cpp
    collider.h
)

# Link dependencies
target_link_libraries(maglit 
    PUBLIC sode                          
    PRIVATE ${HDF5_LIBRARIES}           
            ${FUSION_IO_LIB_DIR}/libfusionio.so  
            ${FUSION_IO_LIB_DIR}/libm3dc1.so    
)

# Include directories
target_include_directories(maglit 
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}    
           ${FUSION_IO_INCLUDE_DIR}       
    PRIVATE ${HDF5_INCLUDE_DIRS}         
)