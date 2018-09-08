# - Try to find AMGX
# Once done this will define
#
#  AMGX_FOUND        - system has PETSc
#  AMGX_INCLUDES     - the PETSc include directories
#  AMGX_LIBRARIES    - Link these to use PETSc
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_policy(VERSION 3.3)

include(DuneMPI)

set(AMGX_VALID_COMPONENTS
    C
    CXX)

# search metis header
find_path(AMGX_INCLUDE_DIR amgx_c.h
  PATHS ${AMGX_ROOT} ${AMGX_DIR}
  PATH_SUFFIXES include include/amgx
  NO_DEFAULT_PATH
  DOC "Include directory of AMGX")

find_library(AMGX_LIBRARY
    NAMES "amgx" "amgxsh"
    PATHS ${AMGX_ROOT}
    PATH_SUFFIXES "lib" "lib" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    NO_DEFAULT_PATH)

if( AMGX_LIBRARY )
  set(AMGX_INCLUDE_DIRS ${AMGX_INCLUDE_DIR})
  set(AMGX_LIBRARIES ${AMGX_LIBRARY})
  set(AMGX_FOUND TRUE)
endif()

#set HAVE_AMGX for config.h
set(HAVE_AMGX ${AMGX_FOUND})

# register all METIS related flags
if(AMGX_FOUND)
  dune_register_package_flags(LIBRARIES "${AMGX_LIBRARIES}"
                              INCLUDE_DIRS "${AMGX_INCLUDE_DIRS}")
endif()
