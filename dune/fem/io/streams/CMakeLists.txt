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
dune_add_subdirs(${SUBDIRS})
