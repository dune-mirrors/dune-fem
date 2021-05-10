# - Try to find AmgX
# Once done this will define
#
#  AMGXSOLVER_FOUND        - system has PETSc
#  AMGXSOLVER_INCLUDES     - the PETSc include directories
#  AMGXSOLVER_LIBRARIES    - Link these to use PETSc

cmake_policy(VERSION 3.3)

include(AddMPIFlags)

set(AMGX_VALID_COMPONENTS
    C
    CXX)

SET(CUDA_DIR "$ENV{CUDA_DIR}" CACHE PATH "The path to CUDA.")
find_package(CUDA)
find_package(PETSc)

if(${CUDA_FOUND} AND ${PETSC_FOUND})
    SET(CUDA_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cusolver_LIBRARY})
    #MESSAGE("-- Finding CUDA - Success")
    message("Cuda = ${CUDA_LIBRARY_DIRS} ${CUDA_LIBRARIES}")

  # search AMGX header
  find_path(AMGX_INCLUDE_DIR amgx_c.h amgx_capi.h
    PATHS ${AMGX_ROOT} ${AMGX_DIR}
    PATH_SUFFIXES include include/amgx
    NO_DEFAULT_PATH
    DOC "Include directory of AMGX")

  find_library(AMGX_LIBRARY
      NAMES "amgxsh" "amgx"
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
    set(AMGXSOLVER_INCLUDE_DIRS ${AMGX_WRAPPER_INCLUDE_DIR} ${AMGX_INCLUDE_DIR})
    set(AMGXSOLVER_LIBRARIES ${AMGXSOLVER_LIBRARIES} ${AMGX_WRAPPER_LIBRARY})
    set(AMGXSOLVER_LIBRARIES ${AMGXSOLVER_LIBRARIES} ${AMGX_LIBRARY})
    set(AMGXSOLVER_LIBRARIES ${AMGXSOLVER_LIBRARIES} ${CUDA_LIBRARIES})
    set(AMGXSOLVER_FOUND TRUE)
    message("AMG Libraries = ${AMGXSOLVER_LIBRARIES}")
  endif()

  #set HAVE_AMGXSOLVER for config.h
  set(HAVE_AMGXSOLVER ${AMGXSOLVER_FOUND})

  # register all AMGX related flags
  if(AMGXSOLVER_FOUND)
    dune_register_package_flags(LIBRARIES "${AMGXSOLVER_LIBRARIES}"
                                INCLUDE_DIRS "${AMGXSOLVER_INCLUDE_DIRS}")
  endif()

# endif CUDA_FOUND AND PETSC_FOUND
endif()
