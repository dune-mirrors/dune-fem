set(vectorfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/vectorfunction)
set(vectorfunction_HEADERS  managedvectorfunction.hh
                         vectorfunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${vectorfunction_HEADERS} DESTINATION ${vectorfunctiondir})