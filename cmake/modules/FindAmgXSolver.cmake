# - Try to find AmgX
# Once done this will define
#
#  AMGXSOLVER_FOUND        - system has PETSc
#  AMGXSOLVER_INCLUDES     - the PETSc include directories
#  AMGXSOLVER_LIBRARIES    - Link these to use PETSc
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_policy(VERSION 3.3)

include(DuneMPI)

set(AMGX_VALID_COMPONENTS
    C
    CXX)

  # search AMGX header
find_path(AMGX_INCLUDE_DIR amgx_c.h amgx_capi.h
  PATHS ${AMGX_ROOT} ${AMGX_DIR}
  PATH_SUFFIXES include include/amgx
  NO_DEFAULT_PATH
  DOC "Include directory of AMGX")

find_library(AMGX_LIBRARY
    NAMES "amgx" "amgxsh"
    PATHS ${AMGX_ROOT}
    PATH_SUFFIXES "lib" "lib" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    NO_DEFAULT_PATH)

find_path(AMGX_WRAPPER_INCLUDE_DIR AmgXSolver.hpp
  PATHS ${AMGX_WRAPPER_ROOT} ${AMGX_WRAPPER_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Include directory of AMGX_WRAPPER")

find_library(AMGX_WRAPPER_LIBRARY
    NAMES "AmgXWrapper"
    PATHS ${AMGX_WRAPPER_ROOT}
    PATH_SUFFIXES "lib" "lib" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    NO_DEFAULT_PATH)

message("amg = ${AMGX_LIBRARY}")
message("amg wrapper = ${AMGX_WRAPPER_LIBRARY}")

if( AMGX_WRAPPER_LIBRARY )
  set(AMGXSOLVER_INCLUDE_DIRS ${AMGX_INCLUDE_DIR} ${AMGX_WRAPPER_INCLUDE_DIR})
  set(AMGXSOLVER_LIBRARIES ${AMGX_LIBRARY} ${AMGX_WRAPPER_LIBRARY})

  set(AMGXSOLVER_FOUND TRUE)
endif()

#set HAVE_AMGXSOLVER for config.h
set(HAVE_AMGXSOLVER ${AMGXSOLVER_FOUND})

# register all AMGX related flags
if(AMGXSOLVER_FOUND)
  dune_register_package_flags(LIBRARIES "${AMGXSOLVER_LIBRARIES}"
                              INCLUDE_DIRS "${AMGXSOLVER_INCLUDE_DIRS}")
endif()
