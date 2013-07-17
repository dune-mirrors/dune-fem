set(adaptivefunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/adaptivefunction)
set(adaptivefunction_HEADERS  adaptivefunction.hh
                           adaptiveimp.cc
                           adaptiveimp.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${adaptivefunction_HEADERS} DESTINATION ${adaptivefunctiondir})