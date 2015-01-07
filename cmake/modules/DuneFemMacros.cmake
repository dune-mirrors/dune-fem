
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
  CACHE FILEPATH "path to local libs in dune-fem" )
  mark_as_advanced(LOCAL_LIBS)

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

########################
# pthreads....
########################
# we are using the cmake default implementation for threads at the moment
include(FindThreads)
set(HAVE_PTHREAD 0)
set(USE_PTHREADS 0 CACHE BOOL "whether we are using pthreads.")
if(CMAKE_USE_PTHREADS_INIT AND NOT USE_PTHREADS)
  set(HAVE_PTHREAD 1)

  set(CMAKE_THREAD_PREFER_PTHREAD 1)

  include(CMakePushCheckState)
  cmake_push_check_state()

  #new settings
  set(CMAKE_REQUIRED_LIBS "${CMAKE_REQUIRED_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}" )
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${PTHREAD_CFLAGS}" )

  check_cxx_source_compiles(
                "include <pthread.h>
                static const int maxThreads = 2 ;
                static int& threadNumber () {
                    static __thread int tn;
                    return tn;
                }
                static int threadnumber[ maxThreads ];
                pthread_barrier_t barrier;
                static void* setThreadNumber(void *number) {
                  const int num = *((int*) number);
                  threadNumber () = num ;
                  //sleep (1);
                  threadnumber[ num ] = threadNumber ();
                  pthread_barrier_wait( &barrier );
                  return 0;
                }
                int main () {
                  pthread_t thread[ maxThreads ];
                  int num[ maxThreads ];
                  for (int i=0; i<maxThreads; ++i )
                  {
                    num[ i ] = i;
                    threadnumber[ i ] = -1;
                  }
                  pthread_barrier_init( &barrier, 0, maxThreads );
                  for( int i=0; i<maxThreads; ++ i)
                    pthread_create (&thread[ i ], 0, &setThreadNumber, (void *) &num[ i ]);

                  for( int i=0; i<maxThreads; ++ i)
                    pthread_join (thread[ i ], 0);
                  int result = 0;
                  for( int i=0; i<maxThreads; ++i )
                    if( threadnumber[ i ] != i ) result = 1;
                  return result;
                }"
                DUNE_CV_PTHREAD_TLS_COMPILES)

  #message( "DUNE_CV_PTHREAD_TLS_COMPILES is set to ${DUNE_CV_PTHREAD_TLS_COMPILES}" )

  if(DUNE_CV_PTHREAD_TLS_COMPILES)
    set(HAVE_PTHREAD_TLS 1)
  endif(DUNE_CV_PTHREAD_TLS_COMPILES)

  cmake_pop_check_state()

endif(CMAKE_USE_PTHREADS_INIT AND NOT USE_PTHREADS)

###############################
# end pthreads
################################

find_package(SIONlib)
include(AddSIONlibFlags)
find_package(PAPI)
include(AddPAPIFlags)

# download PETSc cmake files
set(PETSC_CMAKE_MODULES "${PROJECT_SOURCE_DIR}/cmake/modules/petsc/")
if(NOT EXISTS "${PETSC_CMAKE_MODULES}")
  message (STATUS "Downloading cmake-modules from Jed Brown into ${PETSC_CMAKE_MODULES}")
  execute_process(COMMAND git clone https://github.com/jedbrown/cmake-modules.git ${PETSC_CMAKE_MODULES}
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake/modules/)
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
    set(PETSC_DUNE_COMPILE_FLAGS "-I${PETSC_INCLUDES}"
      CACHE STRING "Compile Flags used by DUNE when compiling with PETSc programs")
    set(PETSC_DUNE_LIBRARIES ${PETSC_LIBRARIES}
      CACHE STRING "Libraries used by DUNE when linking PETSc programs")
    foreach(dir ${PETSC_INCLUDES})
      set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
    endforeach()
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_LIBS "${PETSC_LIBRARIES}")
  else()
    # log errornous result
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Determing location of PETSc failed:\n"
      "Include directory: ${PETSC_INCLUDES}\n"
      "Library directory: ${PETSC_LIBRARIES}\n\n")
  endif()
endif()

# check for XDR (deprecated)
find_package(XDR)

####### abbreviations
include(FemShort)

####### hacks
include(CommandLineHacks)
