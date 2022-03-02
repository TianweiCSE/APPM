cmake_minimum_required(VERSION 3.0)

project(appm LANGUAGES CXX)

add_executable(appm)

find_package (HDF5 REQUIRED COMPONENTS C CXX HL)

target_include_directories(appm PRIVATE "E:/Packages/Eigen/eigen-3.3.7")
target_include_directories(appm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(appm PRIVATE ${HDF5_INCLUDE_DIR})

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC)
target_sources(appm PRIVATE ${SRC})

target_link_libraries(appm PRIVATE  ${HDF5_LIBRARIES})

set(CMAKE_CXX_FLAGS "-fopenmp")

# set(PROJECT_BINARY_DIR "D:/ytw/Thesis/build")
message(STATUS "Project binary dir: ${PROJECT_BINARY_DIR}")
message(STATUS "HDF5 root: ${HDF5_ROOT}")
message(STATUS "HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIR}")
message(STATUS "HDF5 CXX lib: ${HDF5_LIBRARIES}")
message(STATUS "${CMAKE_CXX_FLAGS}")

