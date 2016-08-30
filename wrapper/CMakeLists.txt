
#Convert my XML into a C/C++ file->Creates GRccmTestProjectCLP.h in the binary directory
generateCLP(GRTestProject_SRC GRTestProject.xml)

set (Registration_SOURCES ../src/GroupwiseRegistration.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
add_executable(Groups ${Registration_SOURCES} GRTestProject.cpp ${GRTestProject_SRC})
target_link_libraries(Groups ${LAPACK_LIBRARIES} Mesh)
