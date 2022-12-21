cmake_minimum_required(VERSION 3.0)

project(appm LANGUAGES CXX)

add_executable(appm)
set_property(TARGET appm PROPERTY CXX_STANDARD 17)

find_package (HDF5 REQUIRED COMPONENTS C CXX HL)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# Please change the path to eigen accourding to your setup 
target_include_directories(appm PRIVATE "E:/Packages/Eigen/eigen-3.3.7")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
find_package (Eigen3 REQUIRED CONFIG) 
target_link_libraries(appm PRIVATE Eigen3::Eigen)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

target_include_directories(appm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(appm PRIVATE ${HDF5_INCLUDE_DIR})

target_link_libraries(appm PRIVATE ${HDF5_LIBRARIES})

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC)
target_sources(appm PRIVATE ${SRC})

set(CMAKE_CXX_FLAGS "-fopenmp -O3 -DNDEBUG -march=native")

