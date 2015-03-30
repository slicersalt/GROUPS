project(GRTestProject)
cmake_minimum_required(VERSION 2.8)

#Find GenerateCLP
#eg: /devel/linux/SlicerExecutionModel/SlicerExecutionModel-build_ITKv4.3.0_stat_v3Compatible_linux64
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})

#Where is my source code 
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

#Convert my XML into a C/C++ file->Creates GRTestProjectCLP.h in the binary directory
generateCLP( GRTestProject_SRC GRTestProject.xml )

#Create my executable
add_executable( GRTestProject GRTestProject.cpp ${GRTestProject_SRC} )
