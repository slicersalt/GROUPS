#-----------------------------------------------------------------------------
set(MODULE_NAME Groups)

#-----------------------------------------------------------------------------

# SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

# ITK
#set(${PROJECT_NAME}_ITK_COMPONENTS
#  ITKIOImageBase
#  ITKSmoothing
#  )
#find_package(ITK 4.6 COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
#set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # See Libs/ITKFactoryRegistration/CMakeLists.txt
#include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  )

set(MODULE_SRCS
  GroupwiseRegistration.cpp
  )

if (WIN32)
list(APPEND MODULE_SRCS dirent.c)
endif()

set(MODULE_TARGET_LIBRARIES
#  ${ITK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  Mesh
  ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
#  add_subdirectory(Testing)
endif()
