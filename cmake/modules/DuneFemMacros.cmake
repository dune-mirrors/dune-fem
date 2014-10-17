
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
		CACHE STRING "path to local libs in dune-fem" )
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
  #todo: think about this compile definition
  # < and > may break cmake tools on some systems...
  set(SYSTEM_ENDIAN_HEADER "<endian.h>")
  add_definitions("-DSYSTEM_ENDIAN_HEADER=${SYSTEM_ENDIAN_HEADER}")
else(HAVE_ENDIAN_HEADER_HH)
  check_include_file_cxx(maschine/endian.h HAVE_ENDIAN_MACHINE_HEADER_HH)
  if( HAVE_ENDIAN_MACHINE_HEADER_HH )
    find_path(SYSTEM_ENDIAN_HEADER_PATH
      NAMES endian.h
      PATHS "${ENDIAN_HEADER_ROOT}/machine"
      DOC "Path where machine/endian.h was found"
    )
    #todo: think about this compile definition
    # < and > may break cmake tools on some systems...
    set(SYSTEM_ENDIAN_HEADER "<machine/endian.h>")
    add_definitions("-DSYSTEM_ENDIAN_HEADER=${SYSTEM_ENDIAN_HEADER}")
  endif(HAVE_ENDIAN_MACHINE_HEADER_HH)
endif(HAVE_ENDIAN_HEADER_HH)
mark_as_advanced(ENDIAN_HEADER_ROOT)

#todo: there seems to be no cmake equivalent for
#xdr issues in dune-common. Thus, the following lines are only for completeness
#Please check whether it works and/or delete this stuff.
set(XDR_ROOT "" CACHE STRING "root of header rpc.h")
find_path(DUNE_PATH_XDR
  NAMES rpc.h
  PATHS ${XDR_ROOT}
  PATH_SUFFIXES "rpc"
  DOC "path to header rpc.h"
)
if(EXISTS DUNE_PATH_XDR)
  include(CheckFunctionExists)
  check_function_exists(xdr_uint64_t XDR_UINT64_FUNC_EXISTS)
  if( XDR_UINT64_FUNC_EXISTS )
    add_definitions("-DXDR_UINT64_FUNC=xdr_uint64_t")
  else(XDR_UINT64_FUNC_EXISTS )
    check_function_exists(xdr_u_int64_t XDR_U_INT64_FUNC_EXISTS)
    if( XDR_U_INT64_FUNC_EXISTS )
      add_definitions("-DXDR_UINT64_FUNC=xdr_u_int64_t")
    endif(XDR_U_INT64_FUNC_EXISTS )
  endif(XDR_UINT64_FUNC_EXISTS )
else()
  message("can not find rpc.h")
endif(EXISTS DUNE_PATH_XDR)
mark_as_advanced(XDR_ROOT)

include(CheckCXXSourceCompiles)


########################
# pthreads....
########################
set(USE_PTHREADS 0 CACHE BOOL "whether we are using pthreads.")
if(USE_PTHREADS)
  message(AUTHOR_WARNING "TODO. Please check pthread issues. Not all systems are supported, yet")

  # we are using the cmake default implemenation
  set(CMAKE_THREAD_PREFER_PTHREAD 1)
  include(FindThreads)
  
  if(Threads_FOUND)
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
      add_definitions("-DHAVE_PTHREAD_TLS=1")
    endif(DUNE_CV_PTHREAD_TLS_COMPILES)
    
    cmake_pop_check_state()

  endif(Threads_FOUND)

endif(USE_PTHREADS)

###############################
# end pthreads
################################


find_package( sionlib )


message(AUTHOR_WARNING "TODO. Improve module test.")

