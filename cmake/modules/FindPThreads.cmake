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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PTHREADS")
    #string(FIND ${CMAKE_CXX_FLAGS} "-DUSE_PTHREADS=1" PTHREADS_FLAGS_FOUND)
    #if(${PTHREADS_FLAGS_FOUND} STREQUAL "-1")
    #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PTHREADS=1")
    #endif()
  endif(USE_PTHREADS)

  cmake_pop_check_state()

endif(CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREAD)
