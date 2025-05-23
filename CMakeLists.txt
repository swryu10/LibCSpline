cmake_minimum_required(VERSION 3.31)

project(CSpline)
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_library(CSpline STATIC)
target_sources(CSpline PRIVATE ${PROJECT_SOURCE_DIR}/src/InterCSpline1D.cpp
                               ${PROJECT_SOURCE_DIR}/src/InterCSpline2D.cpp
                               ${PROJECT_SOURCE_DIR}/src/InterCSpline3D.cpp)
target_include_directories(CSpline PUBLIC "${PROJECT_SOURCE_DIR}/include")

project(test_CSpline1D)
add_executable(test_CSpline1D ${PROJECT_SOURCE_DIR}/src/test_CSpline1D.cpp)
set_target_properties(test_CSpline1D PROPERTIES OUTPUT_NAME "test_CSpline1D.exec")
target_link_libraries(test_CSpline1D PUBLIC CSpline)

project(test_CSpline2D)
add_executable(test_CSpline2D ${PROJECT_SOURCE_DIR}/src/test_CSpline2D.cpp)
set_target_properties(test_CSpline2D PROPERTIES OUTPUT_NAME "test_CSpline2D.exec")
target_link_libraries(test_CSpline2D PUBLIC CSpline)

project(test_CSpline3D)
add_executable(test_CSpline3D ${PROJECT_SOURCE_DIR}/src/test_CSpline3D.cpp)
set_target_properties(test_CSpline3D PROPERTIES OUTPUT_NAME "test_CSpline3D.exec")
target_link_libraries(test_CSpline3D PUBLIC CSpline)
