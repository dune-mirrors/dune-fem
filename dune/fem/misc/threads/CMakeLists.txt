set(threadsdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/misc/threads)
set(SUBDIRS  )
set(threads_HEADERS  domainthreaditerator.hh  threaditerator.hh
                                                                        threadmanager.hh  threadpartitioner.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${threads_HEADERS} DESTINATION ${threadsdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
