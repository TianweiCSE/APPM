cmake_minimum_required(VERSION 3.0)

project(appm LANGUAGES CXX)

add_compile_definitions(EIGEN_USE_MKL_ALL)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -march=native -fopenmp")

add_executable(appm)
set_property(TARGET appm PROPERTY CXX_STANDARD 17)

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
target_link_libraries(appm PRIVATE mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
#target_link_libraries(appm PUBLIC -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)

