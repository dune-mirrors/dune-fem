# Module that checks whether SIONLIB is available.
#
# Variables used by this module which you may want to set:
# SIONLIB_ROOT        Path list to search for SIONLIB
# SIONLIB_SUFFIX      suffix to the library name , e.g. gcc or something
# SIONLIB_INCLUDEDIR  directory with SIONLIB headers inside
# SIONLIB_LIBDIR      directory with SIONLIB libraries inside
#
# Sets the following variables
#
# SIONLIB_FOUND          True if SIONLIB was found and usable
# HAVE_SIONLIB           True if SIONLIB was found and usable
# SIONLIB_INCLUDE_DIRS   Path to the SIONLIB include dirs
# SIONLIB_LIBRARIES      Name of the SIONLIB libraries
#

set(SIONLIB_ROOT "" CACHE PATH "Path list to search for SIONLIB")
set(SIONLIB_SUFFIX "_lib64" CACHE STRING "suffix to the library name , e.g. gcc or something")
set(SIONLIB_INCLUDEDIR "" CACHE PATH "directory with SIONLIB headers inside")
set(SIONLIB_LIBDIR "" CACHE PATH "directory with SIONLIB libraries inside")

mark_as_advanced(SIONLIB_ROOT SIONLIB_SUFFIX SIONLIB_INCLUDEDIR SIONLIB_LIBDIR)


#look for header files at positions given by the user
find_path(SIONLIB_INCLUDE_DIR
  NAMES "sion.h"
  PATHS ${SIONLIB_ROOT} ${SIONLIB_INCLUDEDIR}
  PATH_SUFFIXES "include"
  NO_DEFAULT_PATH
)
#now also look for default paths
find_path(SIONLIB_INCLUDE_DIR
  NAMES "sion.h"
  PATH_SUFFIXES "include"
)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()
set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} ${MPI_DUNE_COMPILE_FLAGS} -DENABLE_SIONLIB")
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_DUNE_INCLUDE_PATH} ${SIONLIB_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${MPI_DUNE_LIBRARIES})
include(CheckIncludeFiles)
check_include_files(sion.h SIONLIB_HEADER_USABLE)


#look for library at positions given by the user
find_library(SIONLIB_LIBRARY
  NAMES "sion"
  PATHS ${SIONLIB_ROOT} ${SIONLIB_LIBDIR}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)
#now  also include the default paths
find_library(SIONLIB_LIBRARY
  NAMES "sion"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)

# check if library sion/sionser works
include(CheckSymbolExists)
if(SIONLIB_LIBRARY)
  get_filename_component(SIONLIB_LIB_PATH ${SIONLIB_LIBRARY} PATH)
  check_library_exists("sion${SIONLIB_SUFFIX}" main ${SIONLIB_LIBRARY} SIONLIB_LIB_WORKS)
  check_library_exists("sionser${SIONLIB_SUFFIX}" main ${SIONLIB_LIBRARY} SIONLIB_LIB_SIONSER_WORKS)
endif(SIONLIB_LIBRARY)

cmake_pop_check_state()



# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "SIONlib"
  DEFAULT_MSG
  SIONLIB_INCLUDE_DIR
  SIONLIB_LIBRARY
  SIONLIB_HEADER_USABLE
  SIONLIB_LIB_WORKS
  SIONLIB_LIB_SIONSER_WORKS
)

mark_as_advanced(SIONLIB_INCLUDE_DIR SIONLIB_LIBRARY)

# if both headers and library are found, store results
if(SIONLIB_FOUND)
  set(SIONLIB_INCLUDE_DIRS ${SIONLIB_INCLUDE_DIR})
  set(SIONLIB_LIBRARIES ${SIONLIB_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of SIONLIB succeded:\n"
    "Include directory: ${SIONLIB_INCLUDE_DIRS}\n"
    "Library directory: ${SIONLIB_LIBRARIES}\n\n")
  set(SIONLIB_DUNE_COMPILE_FLAGS "-I${SIONLIB_INCLUDE_DIRS} -DENABLE_SIONLIB=1"
    CACHE STRING "Compile Flags used by DUNE when compiling with SIONLIB programs")
  set(SIONLIB_DUNE_LIBRARIES ${SIONLIB_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking SIONLIB programs")
else(SIONLIB_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of SIONLIB failed:\n"
    "Include directory: ${SIONLIB_INCLUDE_DIRS}\n"
    "Library directory: ${SIONLIB_LIBRARIES}\n\n")
endif(SIONLIB_FOUND)

#set HAVE_SIONLIB for config.h
set(HAVE_SIONLIB ${SIONLIB_FOUND})

#add all sionlib related flags to ALL_PKG_FLAGS, this must happen regardless of a target using add_dune_sionlib_flags
if(SIONLIB_FOUND)
  set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "${SIONLIB_DUNE_COMPILE_FLAGS}")
  foreach(dir "${SIONLIB_INCLUDE_DIRS}")
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
endif()


