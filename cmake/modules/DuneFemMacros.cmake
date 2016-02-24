
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
  CACHE FILEPATH "path to local libs in dune-fem" )
  mark_as_advanced(LOCAL_LIBS)

set(USE_OPENMP OFF CACHE BOOL "whether we are using OpenMP.")
# if open mp should be used perform cmake check
if(USE_OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    # add flags to compiler flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)

#find endian headers
set(ENDIAN_HEADER_ROOT "" CACHE STRING "path of endian header")
check_include_file_cxx(endian.h HAVE_ENDIAN_HEADER_HH)
if( HAVE_ENDIAN_HEADER_HH )
  find_path(SYSTEM_ENDIAN_HEADER_PATH
    NAMES endian.h
    PATHS ${ENDIAN_HEADER_ROOT}
    DOC "Path where endian.h was found"
  )
  if(EXISTS SYSTEM_ENDIAN_HEADER_PATH)
    set(SYSTEM_ENDIAN_HEADER "<endian.h>")
  endif(EXISTS SYSTEM_ENDIAN_HEADER_PATH)
else(HAVE_ENDIAN_HEADER_HH)
  check_include_file_cxx(machine/endian.h HAVE_ENDIAN_MACHINE_HEADER_HH)
  if( HAVE_ENDIAN_MACHINE_HEADER_HH )
    find_path(SYSTEM_ENDIAN_HEADER_PATH
      NAMES endian.h
      PATHS "${ENDIAN_HEADER_ROOT}/machine"
      PATH_SUFFIXES "machine"
      DOC "Path where machine/endian.h was found"
    )
    if(EXISTS SYSTEM_ENDIAN_HEADER_PATH)
      set(SYSTEM_ENDIAN_HEADER "<machine/endian.h>")
    endif(EXISTS SYSTEM_ENDIAN_HEADER_PATH)
  endif(HAVE_ENDIAN_MACHINE_HEADER_HH)
endif(HAVE_ENDIAN_HEADER_HH)
mark_as_advanced(ENDIAN_HEADER_ROOT SYSTEM_ENDIAN_HEADER_PATH)

include(CheckCXXSourceCompiles)

# check for phtreads
include(FindPThreads)

###############################
# end pthreads
################################

find_package(SIONlib)
include(AddSIONlibFlags)
find_package(PAPI)
include(AddPAPIFlags)

# download PETSc cmake files
set(PETSC_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/petsc/")
if(NOT EXISTS "${PETSC_CMAKE_MODULES}")
  message (STATUS "Downloading cmake-modules from Jed Brown into ${PETSC_CMAKE_MODULES}")
  execute_process(COMMAND git clone https://github.com/jedbrown/cmake-modules.git ${PETSC_CMAKE_MODULES}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
else()
  message (STATUS "Updating cmake-modules from Jed Brown into ${PETSC_CMAKE_MODULES}")
  execute_process(COMMAND git pull
                  WORKING_DIRECTORY ${PETSC_CMAKE_MODULES})
endif ()

# check for PETSc if Jed Browns cmake-modules have been successfully downloaded
# needs to set environment variables PETSC_DIR and optional PETSC_ARCH
if(EXISTS "${PETSC_CMAKE_MODULES}")
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PETSC_CMAKE_MODULES}")
  set(PETSC_DIR  $ENV{PETSC_DIR})
  set(PETSC_ARCH $ENV{PETSC_ARCH})
  find_package(PETSc)
  if(PETSC_FOUND)
    # set HAVE_PETSC for config.h
    set(HAVE_PETSC ${PETSC_FOUND})
    # log result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining location of PETSC succeded:\n"
      "Include directory: ${PETSC_INCLUDES}\n"
      "Library directory: ${PETSC_LIBRARIES}\n\n")
    set(PETSC_DUNE_COMPILE_FLAGS "-DENABLE_PETSC=1 -I${PETSC_INCLUDES}"
      CACHE STRING "Compile Flags used by DUNE when compiling with PETSc programs")
    set(PETSC_DUNE_LIBRARIES ${PETSC_LIBRARIES}
      CACHE STRING "Libraries used by DUNE when linking PETSc programs")
    foreach(dir ${PETSC_INCLUDES})
      set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
    endforeach()
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_LIBS "${PETSC_LIBRARIES}")
    # register includes and libs to global flags
    dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_PETSC=1"
                                INCLUDE_DIRS ${PETSC_INCLUDES}
                                LIBRARIES ${PETSC_LIBRARIES})

  else()
    # log errornous result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determing location of PETSc failed:\n"
      "Include directory: ${PETSC_INCLUDES}\n"
      "Library directory: ${PETSC_LIBRARIES}\n\n")
  endif()
endif()

# check eigen
pkg_check_modules(EIGEN eigen3)
if(EIGEN_FOUND)
  foreach( ARG ${EIGEN_INCLUDE_DIRS} )
    if(IS_DIRECTORY ${ARG} )
      add_definitions("-isystem${ARG}")
      include_directories(${ARG})
    else()
      message( STATUS "Include directory ${EIGEN_INCLUDE_DIRS} (needed for eigen) does not exist" )
    endif()
  endforeach()
  set( HAVE_EIGEN 1 )
endif(EIGEN_FOUND)

# check for XDR (deprecated)
find_package(XDR)

####### abbreviations
include(FemShort)

####### hacks
include(CommandLineHacks)

# check SuiteSparse support
find_package(SuiteSparse OPTIONAL_COMPONENTS UMFPACK SPQR LDL)
include(AddSuiteSparseFlags)
