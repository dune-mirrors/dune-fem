set(petscdiscretefunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/petscdiscretefunction)
set(petscdiscretefunction_HEADERS   petscdiscretefunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${petscdiscretefunction_HEADERS} DESTINATION ${petscdiscretefunctiondir})