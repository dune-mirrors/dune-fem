set(operatordir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator)
set(operator_HEADERS  discreteoperatorimp.hh
                   integrationoperator.hh
                   lagrangeinterpolation.hh)
set(SUBDIRS  common 0order 1order 2order matrix model projection)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${operator_HEADERS} DESTINATION ${operatordir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
