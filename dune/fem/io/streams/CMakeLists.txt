# don't follow the full GNU-standard; we definitely need automake >= 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
set(SUBDIRS  test)
set(streamsdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/streams)
set(streams_HEADERS  asciistreams.hh
                  binarystreams.hh
                  sionlibstreams.hh
                  standardstreams.hh
                  streams.hh
                  streams_inline.hh
                  tuples.hh
                  virtualstreams.hh
                  xdrstreams.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${streams_HEADERS} DESTINATION ${streamsdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
