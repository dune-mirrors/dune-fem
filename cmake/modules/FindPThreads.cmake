########################
# pthreads....
########################
# we are using the cmake default implementation for threads at the moment
include(FindThreads)
set(HAVE_PTHREAD 0)
set(USE_PTHREADS OFF CACHE BOOL "whether we are using pthreads.")
if(CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREAD)
  set(HAVE_PTHREAD 1)

  set(CMAKE_THREAD_PREFER_PTHREAD 1)

  include(CMakePushCheckState)
  cmake_push_check_state()

  #new settings
  set(CMAKE_REQUIRED_LIBS "${CMAKE_REQUIRED_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}" )
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${PTHREAD_CFLAGS}" )

  if(USE_PTHREADS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PTHREADS=1")
  endif(USE_PTHREADS)

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

endif(CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREAD)
