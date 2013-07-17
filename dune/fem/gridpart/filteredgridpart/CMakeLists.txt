set(filteredgridpartdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/filteredgridpart)
set(filteredgridpart_HEADERS  capabilities.hh  datahandle.hh
                           intersectioniterator.hh  iterator.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${filteredgridpart_HEADERS} DESTINATION ${filteredgridpartdir})