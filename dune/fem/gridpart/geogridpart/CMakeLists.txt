set(geogridpartdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/geogridpart)
set(geogridpart_HEADERS  capabilities.hh  cornerstorage.hh  datahandle.hh
                      entity.hh  entitypointer.hh  geometry.hh
                      intersection.hh  intersectioniterator.hh  iterator.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${geogridpart_HEADERS} DESTINATION ${geogridpartdir})