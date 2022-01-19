# Try to find PSCToolkit
# https://psctoolkit.github.io/
#
#  PSCTOOLKIT_FOUND         - system has PSCToolkit
#  PSCTOOLKIT_INCLUDE_DIRS  - the PSCToolkit include directories
#  PSCTOOLKIT_LIBRARIES     - Link these to use PSCToolkit

set(PSCTOOLKIT_VALID_COMPONENTS
    C
    CXX)

include(AddMPIFlags)

# required for PSCToolkit
find_package(METIS)

# search PSCTOOLKIT header
find_path(PSCTOOLKIT_INCLUDE_DIRS psb_krylov_cbind.h psb_prec_cbind.h psb_c_base.h psb_util_cbind.h
  PATHS ${PSCTOOLKIT_ROOT} ${PSCTOOLKIT_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Include directory of PSCTOOLKIT")

foreach (psclib psb_base psb_cbind psb_krylov psb_prec psb_util)
  find_library(PSCTOOLKIT_LIB_${psclib}
    NAMES ${psclib}
    PATHS ${PSCTOOLKIT_ROOT} ${PSCTOOLKIT_DIR}
    PATH_SUFFIXES "lib" "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    NO_DEFAULT_PATH)
  if(PSCTOOLKIT_LIB_${psclib})
    set(PSCTOOLKIT_LIBRARIES ${PSCTOOLKIT_LIBRARIES} ${PSCTOOLKIT_LIB_${psclib}})
  endif()
endforeach()

find_package_handle_standard_args(PSCToolkit
  REQUIRED_VARS PSCTOOLKIT_INCLUDE_DIRS PSCTOOLKIT_LIBRARIES)
