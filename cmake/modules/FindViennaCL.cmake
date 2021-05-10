# - Try to find AmgX
# Once done this will define
#
#  VIENNACL_FOUND        - system has ViennaCL
#  VIENNACL_INCLUDES     - the ViennaCL include directories

cmake_policy(VERSION 3.3)

include(AddMPIFlags)

find_package(OpenCL)

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

#set HAVE_OPENCL for config.h
set(HAVE_OPENCL ${OpenCL_FOUND})

# register all VIENNACL related flags
if(VIENNACL_FOUND)
  dune_register_package_flags(INCLUDE_DIRS "${VIENNACL_INCLUDE_DIRS}")
endif()

# register all OpenCL related flags
if(OpenCL_FOUND)
  message("opencl = ${OpenCL_LIBRARIES}")
  dune_register_package_flags(
    LIBRARIES "${OpenCL_LIBRARIES}"
    INCLUDE_DIRS "${OpenCL_INCLUDE_DIRS}")
endif()
