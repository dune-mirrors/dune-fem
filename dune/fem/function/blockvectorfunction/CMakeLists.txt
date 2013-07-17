set(blockvectorfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/blockvectorfunction)
set(blockvectorfunction_HEADERS  blockvectorfunction.hh
blockvectorfunction_inline.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${blockvectorfunction_HEADERS} DESTINATION ${blockvectorfunctiondir})