project(Groups)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(SlicerExecutionModel REQUIRED)
find_package(MeshLib REQUIRED)

include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories(Mesh GroupwiseRegistration)
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )


include_directories(src)

include_directories(wrapper)
add_subdirectory(wrapper)
