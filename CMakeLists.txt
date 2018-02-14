cmake_minimum_required(VERSION 2.8.3)
cmake_policy(VERSION 2.8.3)

project(GROUPS)

set(GROUPS_VERSION_MAJOR 1)
set(GROUPS_VERSION_MINOR 0)
set(GROUPS_VERSION_PATCH 0)

set(LOCAL_PROJECT_NAME GROUPS)

#-----------------------------------------------------------------------------
# Update CMake module path
# We need to update the CMake Module path in this main CMakeLists.txt file
# so that we can include SlicerExtensionsConfigureMacros which is in the current
# ${Project}/CMake folder
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/SuperBuild
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_NAME GROUPS)
set(EXTENSION_HOMEPAGE "")
set(EXTENSION_CATEGORY "GROUPS")
set(EXTENSION_CONTRIBUTORS "Mahmoud Mostapha (UNC)")
set(EXTENSION_DESCRIPTION "Cortical correspondence method employing group-wise registration in a spherical parametrization space for the use in neuroimaging studies")

set(EXTENSION_ICONURL "")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "") 

set(EXTENSION_BUILD_SUBDIRECTORY ${LOCAL_PROJECT_NAME}-build)
set(EXTENSION_STATUS Beta)

set(SUPERBUILD_TOPLEVEL_PROJECT ${LOCAL_PROJECT_NAME})
#-----------------------------------------------------------------------------
# Slicer
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
include(SlicerExtensionsConfigureMacros)
mark_as_superbuild(Slicer_DIR)

# Slicer Execution Model
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})

# Project version number.
set(${LOCAL_PROJECT_NAME}_VERSION_MAJOR "0")
set(${LOCAL_PROJECT_NAME}_VERSION_MINOR "1")
set(${LOCAL_PROJECT_NAME}_VERSION_PATCH "0")

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  ${LOCAL_PROJECT_NAME}_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with ${LOCAL_PROJECT_NAME}_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the ${LOCAL_PROJECT_NAME}, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------

option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
endif()

add_subdirectory(Modules)

#-----------------------------------------------------------------------------
# CPACK
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
include(${Slicer_EXTENSION_CPACK})