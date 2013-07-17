set(0orderdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/0order)
set(0order_HEADERS  )
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${0order_HEADERS} DESTINATION ${0orderdir})