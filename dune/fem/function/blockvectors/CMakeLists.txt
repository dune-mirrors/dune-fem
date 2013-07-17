set(blockvectorsdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/blockvectors)
set(blockvectors_HEADERS  referenceblockvector.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${blockvectors_HEADERS} DESTINATION ${blockvectorsdir})