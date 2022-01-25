cmake_minimum_required(VERSION 3.0)

project(appm LANGUAGES CXX)

add_executable(appm)

find_package (HDF5 REQUIRED COMPONENTS CXX)
find_package (Eigen3 REQUIRED CONFIG) 

target_include_directories(appm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(appm PRIVATE ${HDF5_INCLUDE_DIR})

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC)
target_sources(appm PRIVATE ${SRC})

target_link_libraries(appm PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(appm PRIVATE Eigen3::Eigen) 

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -fopenmp -march=native -std=c++11")

