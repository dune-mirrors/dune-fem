
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
		CACHE STRING "path to local libs in dune-fem" )

#find endian headers
check_include_file_cxx(endian.h HAVE_ENDIAN_HEADER_HH)
if( HAVE_ENDIAN_HEADER_HH )
  find_path(SYSTEM_ENDIAN_HEADER_PATH
    NAMES endian.h )
  set(SYSTEM_ENDIAN_HEADER "<endian.h>"
    CACHE STRING "Path where endian.h was found")
else(HAVE_ENDIAN_HEADER_HH)
  check_include_file_cxx(maschine/endian.h HAVE_ENDIAN_MACHINE_HEADER_HH)
  if( HAVE_ENDIAN_MACHINE_HEADER_HH )
    find_path(SYSTEM_ENDIAN_HEADER_PATH
      NAMES endian.h )
    set(SYSTEM_ENDIAN_HEADER "<machine/endian.h>"
      CACHE STRING "Path where endian.h was found")
  endif(HAVE_ENDIAN_MACHINE_HEADER_HH)
endif(HAVE_ENDIAN_HEADER_HH)
add_definitions("-DSYSTEM_ENDIAN_HEADER=${SYSTEM_ENDIAN_HEADER}")

#todo: there seems to be no cmake equivalent for
#xdr issues in dune-common. Thus, the following lines are only for completeness
#Please check whether it works and/or delete this stuff.
find_path(DUNE_PATH_XDR
  NAMES rpc/rpc.h )
if(DUNE_PATH_XDR)
  include(CheckFunctionExists)
  check_function_exists(xdr_uint64_t XDR_UINT64_FUNC_EXISTS)
  if( XDR_UINT64_FUNC_EXISTS )
    set(XDR_UINT64_FUNC "xdr_uint64_t"
      CACHE STRING "xdr_unit64_t routine")
  else(XDR_UINT64_FUNC_EXISTS )
    check_function_exists(xdr_u_int64_t XDR_U_INT64_FUNC_EXISTS)
    if( XDR_U_INT64_FUNC_EXISTS )
      set(XDR_UINT64_FUNC "xdr_u_int64_t"
        CACHE STRING "xdr_unit64_t routine")
    endif(XDR_U_INT64_FUNC_EXISTS )
  endif(XDR_UINT64_FUNC_EXISTS )
endif(DUNE_PATH_XDR)


include(CheckCXXSourceCompiles)


########################
# pthreads....
########################
if(USE_PTHREADS)
  message(AUTHOR_WARNING "TODO. Please check pthread issues. Not all systems are supported, yet")

  # we are using the cmake default implemenation
  set(CMAKE_THREAD_PREFER_PTHREAD 1)
  include(FindThreads)
  
  if(Threads_FOUND)
    #save old settings
    set(CMAKE_REQUIRED_LIBS_OLD CMAKE_REQUIRED_LIBRARIES )
    set(CMAKE_REQUIRED_FLAGS_OLD CMAKE_REQUIRED_FLAGS )
    
    #new settings
    set(CMAKE_REQUIRED_LIBS "${CMAKE_REQUIRED_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}" )
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${PTHREAD_CFLAGS}" )
    
    #debug information
    #message( "CMAKE_REQUIRED_LIBS is temporary set to ${CMAKE_REQUIRED_LIBS}" ) 
    #message( "CMAKE_REQUIRED_FLAGS is temporary set to ${CMAKE_REQUIRED_FLAGS}" ) 
    
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
    
    #restore old settings
    set(CMAKE_REQUIRED_LIBS CMAKE_REQUIRED_LIBRARIES_OLD )
    set(CMAKE_REQUIRED_FLAGS CMAKE_REQUIRED_FLAGS_OLD )
    
    #message( "DUNE_CV_PTHREAD_TLS_COMPILES is set to ${DUNE_CV_PTHREAD_TLS_COMPILES}" )
    
    if(DUNE_CV_PTHREAD_TLS_COMPILES)
      set(HAVE_PTHREAD_TLS 1
        CACHE BOOL "This is true if the keyword __thread can be used")
      message( "We found thread local storage...")
    endif(DUNE_CV_PTHREAD_TLS_COMPILES)

  endif(Threads_FOUND)

endif(USE_PTHREADS)

###############################
# end pthreads
################################


find_package( sionlib )


message(AUTHOR_WARNING "TODO. Improve module test.")

