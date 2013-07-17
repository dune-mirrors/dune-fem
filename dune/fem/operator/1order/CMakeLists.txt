set(1orderdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/1order)
set(1order_HEADERS  dgmasspass.hh localmassmatrix.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${1order_HEADERS} DESTINATION ${1orderdir})