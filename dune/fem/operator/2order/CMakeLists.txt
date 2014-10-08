set(2orderdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/2order)
set(2order_HEADERS  dgprimaloperator.hh
                    dgmatrixsetup.hh 
                    lagrangematrixsetup.hh) 
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${2order_HEADERS} DESTINATION ${2orderdir})
