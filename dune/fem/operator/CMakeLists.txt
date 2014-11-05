set(operatordir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator)
set(operator_HEADERS  discreteoperatorimp.hh
                   integrationoperator.hh
                   lagrangeinterpolation.hh)
set(SUBDIRS  common 0order 1order 2order matrix model projection)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${operator_HEADERS} DESTINATION ${operatordir})
dune_add_subdirs(${SUBDIRS})
