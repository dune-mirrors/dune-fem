set(SUBDIRS  common
             localfunction
             adaptivefunction
             attachedfunction
             combinedfunction
             blockvectorfunction
             productfunction
             vectorfunction
             blockvectordiscretefunction
             blockvectors
             petscdiscretefunction)
set(functiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function)
set(function_HEADERS  adaptivefunction.hh
                      attachedfunction.hh
                      blockvectorfunction.hh
                      combinedfunction.hh
                      productfunction.hh
                      vectorfunction.hh
                      subfunction.hh
                      blockvectordiscretefunction.hh
                      petscdiscretefunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${function_HEADERS} DESTINATION ${functiondir})
dune_add_subdirs(${SUBDIRS})
