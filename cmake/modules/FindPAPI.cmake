# Module that checks whether PAPI is available.
#
# Variables used by this module which you may want to set:
# PAPI_ROOT        Path list to search for PAPI
#
# Sets the following variables
#
# PAPI_FOUND             True if PAPI was found and usable
# HAVE_PAPI              True if PAPI was found and usable
# PAPI_INCLUDE_DIRS      Path to the PAPI include dirs
# PAPI_LIBRARIES         Name of the PAPI libraries
#

set(PAPI_ROOT "" CACHE PATH "Path list to search for PAPI")

mark_as_advanced(PAPI_ROOT)

#look for header files at positions given by the user
find_path(PAPI_INCLUDE_DIR
  NAMES "papi.h"
  PATHS ${PAPI_ROOT}
  PATH_SUFFIXES "include" 
  NO_DEFAULT_PATH
)
#now also look for default paths
find_path(PAPI_INCLUDE_DIR
  NAMES "papi.h"
  PATH_SUFFIXES "include" 
)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -DENABLE_PAPI")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${PAPI_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
check_include_files(papi.h PAPI_HEADER_USABLE)

#look for library at positions given by the user
find_library(PAPI_LIBRARY
  NAMES "papi"
  PATHS ${PAPI_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)
#now  also include the default paths
find_library(PAPI_LIBRARY
  NAMES "papi"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)

# check if library papi works
include(CheckSymbolExists)
if(PAPI_LIBRARY)
  get_filename_component(PAPI_LIB_PATH ${PAPI_LIBRARY} PATH)
  check_library_exists(papi PAPI_flops ${PAPI_LIBRARY} PAPI_LIB_WORKS)
endif(PAPI_LIBRARY)

cmake_pop_check_state()



# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "PAPI"
  DEFAULT_MSG
  PAPI_INCLUDE_DIR
  PAPI_LIBRARY
  PAPI_HEADER_USABLE
  PAPI_LIB_WORKS
)

mark_as_advanced(PAPI_INCLUDE_DIR PAPI_LIBRARY)

# if both headers and library are found, store results
if(PAPI_FOUND)
  set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
  set(PAPI_LIBRARIES ${PAPI_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of PAPI succeded:\n"
    "Include directory: ${PAPI_INCLUDE_DIRS}\n"
    "Library directory: ${PAPI_LIBRARIES}\n\n")
  set(PAPI_DUNE_COMPILE_FLAGS "-I${PAPI_INCLUDE_DIRS} -DENABLE_PAPI=1"
    CACHE STRING "Compile Flags used by DUNE when compiling with PAPI programs")
  set(PAPI_DUNE_LIBRARIES ${PAPI_LIBRARIES} 
    CACHE STRING "Libraries used by DUNE when linking PAPI programs")
else(PAPI_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of PAPI failed:\n"
    "Include directory: ${PAPI_INCLUDE_DIRS}\n"
    "Library directory: ${PAPI_LIBRARIES}\n\n")
endif(PAPI_FOUND)

#set HAVE_PAPI for config.h
set(HAVE_PAPI ${PAPI_FOUND})

#add all papi related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_papi_flags
if(PAPI_FOUND)
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "${PAPI_DUNE_COMPILE_FLAGS}")
  foreach(dir "${PAPI_INCLUDE_DIRS}")
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()


