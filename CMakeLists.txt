project(Groups)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(SlicerExecutionModel REQUIRED)

include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories(Mesh GroupwiseRegistration)
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

#Convert my XML into a C/C++ file->Creates GRccmTestProjectCLP.h in the binary directory
generateCLP(Groups_SRC wrapper/GroupwiseRegistration.xml)

# MESH library
set(MESH_LIB_DIR /NIRAL/work/ilwoolyu/MeshLib)
include_directories(${MESH_LIB_DIR}/include ${MESH_LIB_DIR}/include/Util)
link_directories(${MESH_LIB_DIR}/lib)

set (Registration_SOURCES src/GroupwiseRegistration.cpp)
include_directories(src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable(Groups ${Registration_SOURCES} wrapper/main.cpp ${Groups_SRC})
target_link_libraries(Groups ${LAPACK_LIBRARIES} Mesh)
