set(basisfunctionsetdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/basisfunctionset)
set(basisfunctionset_HEADERS  basisfunctionset.hh
                           default.hh
                           functor.hh
                           proxy.hh
                           simple.hh
                           transformation.hh
                           vectorial.hh)
set(SUBDIRS  test)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${basisfunctionset_HEADERS} DESTINATION ${basisfunctionsetdir})

dune_add_subdirs(${SUBDIRS})
