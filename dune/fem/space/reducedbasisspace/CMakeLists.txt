# don't follow the full GNU-standard; we definitely need automake >= 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
set(SUBDIRS  test)
set(reducedbasisspacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/reducedbasisspace)
set(reducedbasisspace_HEADERS  basefunctionset.hh
                            commdatahandle.hh
                            mapper.hh
                            reducedbasisspace.hh
                            reducedbasisspace_inline.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${reducedbasisspace_HEADERS} DESTINATION ${reducedbasisspacedir})

