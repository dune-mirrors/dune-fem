set(combinedfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/combinedfunction)
set(combinedfunction_HEADERS  combinedfunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${combinedfunction_HEADERS} DESTINATION ${combinedfunctiondir})