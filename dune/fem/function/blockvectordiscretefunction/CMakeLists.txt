set(blockvectordiscretefunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/blockvectordiscretefunction)
set(blockvectordiscretefunction_HEADERS  blockvectordiscretefunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${blockvectordiscretefunction_HEADERS} DESTINATION ${blockvectordiscretefunctiondir})