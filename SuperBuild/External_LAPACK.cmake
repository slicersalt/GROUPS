
set(proj LAPACK)

# Set dependency list
set(${proj}_DEPENDS
  ""
  )

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(LAPACK_DIR CACHE)
  find_package(LAPACK REQUIRED)
  unset(LAPACKE_DIR CACHE)
  find_package(LAPACKE REQUIRED)
endif()

# Sanity checks
if(DEFINED LAPACK_DIR AND NOT EXISTS ${LAPACK_DIR})
  message(FATAL_ERROR "LAPACK_DIR [${LAPACK_DIR}] variable is defined but corresponds to nonexistent directory")
endif()
if(DEFINED LAPACKE_DIR AND NOT EXISTS ${LAPACKE_DIR})
  message(FATAL_ERROR "LAPACKE_DIR [${LAPACKE_DIR}] variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED LAPACK_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_LAPACK
   AND NOT DEFINED LAPACKE_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_LAPACKE)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/Reference-LAPACK/lapack.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "c5471e8be2754345ec650c8d3b08c6eb680ebd0c"
    QUIET
    )

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      # Output directories
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_BIN_DIR}
      -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_LIB_DIR}
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
      # Install directories
      -DLAPACK_INSTALL_RUNTIME_DIR:STRING=${Slicer_INSTALL_THIRDPARTY_LIB_DIR}
      -DLAPACK_INSTALL_LIBRARY_DIR:STRING=${Slicer_INSTALL_THIRDPARTY_LIB_DIR}
      # Options
      -DBUILD_TESTING:BOOL=OFF
      -DCBLAS:BOOL=OFF
      -DLAPACKE:BOOL=ON
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDS}
    )
  set(LAPACK_DIR ${EP_BINARY_DIR})
  set(LAPACKE_DIR ${EP_BINARY_DIR})

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(LAPACK_DIR:PATH)
mark_as_superbuild(LAPACKE_DIR:PATH)
