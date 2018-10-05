# - Try to find AmgX
# Once done this will define
#
#  VIENNACL_FOUND        - system has ViennaCL
#  VIENNACL_INCLUDES     - the ViennaCL include directories
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

cmake_policy(VERSION 3.3)

include(DuneMPI)

# search VIENNACL header
find_path(VIENNACL_INCLUDE_DIR version.hpp
  PATHS ${VIENNACL_ROOT} ${VIENNACL_DIR}
  PATH_SUFFIXES include include/viennacl
  NO_DEFAULT_PATH
  DOC "Include directory of VIENNACL")

if( VIENNACL_INCLUDE_DIR )
  set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR})

  set(VIENNACL_FOUND TRUE)
endif()

#set HAVE_VIENNACL for config.h
set(HAVE_VIENNACL ${VIENNACL_FOUND})

# register all VIENNACL related flags
if(VIENNACL_FOUND)
  dune_register_package_flags(INCLUDE_DIRS "${VIENNACL_INCLUDE_DIRS}")
endif()
