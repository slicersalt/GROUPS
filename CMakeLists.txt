project(Groups)
cmake_minimum_required(VERSION 2.8)
enable_language(Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(SlicerExecutionModel REQUIRED)

include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories(/System/Library/Frameworks/Accelerate.framework/Headers)
include_directories(Mesh GroupwiseRegistration)
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

#Convert my XML into a C/C++ file->Creates GRccmTestProjectCLP.h in the binary directory
generateCLP(GRTestProject_SRC GRTestProject.xml )

set (Registration_SOURCES GroupwiseRegistration/GroupwiseRegistration.cpp
	GroupwiseRegistration/GroupwiseRegistration.h
	GroupwiseRegistration/SphericalHarmonics.h GroupwiseRegistration/newuoa.h
	GroupwiseRegistration/itkTimer.cpp GroupwiseRegistration/itkTimer.h)

set (Mesh_SOURCES Mesh/Mesh.cpp Mesh/MeshIO.cpp Mesh/MeshIO.h Mesh/ObjIO.cpp Mesh/ObjIO.h Mesh/VtkIO.cpp Mesh/VtkIO.h)

add_executable(Groups ${Registration_SOURCES} ${Mesh_SOURCES} GRTestProject.cpp ${GRTestProject_SRC})
target_link_libraries(Groups ${LAPACK_LIBRARIES})
