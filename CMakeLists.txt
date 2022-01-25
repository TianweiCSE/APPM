cmake_minimum_required(VERSION 3.0)

project(appm LANGUAGES CXX)

add_compile_definitions(EIGEN_USE_MKL_ALL)

add_executable(appm)

find_package (HDF5 REQUIRED COMPONENTS CXX)
find_package (Eigen3 REQUIRED CONFIG) 
# find_package (MKL REQUIRED CONFIG) # This command cannot find the path on Euler

target_include_directories(appm PRIVATE $ENV{MKLROOT}/include)
target_include_directories(appm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(appm PRIVATE ${HDF5_INCLUDE_DIR})

link_directories($ENV{MKLROOT}/lib/intel64)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC)
target_sources(appm PRIVATE ${SRC})

target_link_libraries(appm PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(appm PRIVATE Eigen3::Eigen)
target_link_libraries(appm PRIVATE mkl_intel_ilp64 mkl_intel_thread mkl_core iomp5)


set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -fopenmp -march=native -std=c++11")

